# main.py — FastAPI backend for Chemdem  v3.0
# Run: uvicorn main:app --reload --port 8000

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional

from model import predict as _legacy_predict, predict_from_smiles
from lookup import validate_compound
from reference_search import search_all_references
from side_products import predict_side_products
from rdkit import Chem
from rdkit.Chem import Descriptors

app = FastAPI(title="Chemdem API", version="3.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# ── Request / Response Models ─────────────────────────────────────────────────

class SmilesRequest(BaseModel):
    amine_smiles:    str
    squarate_smiles: str  = "CCOC1=C(OCC)C(=O)C1=O"
    temperature:     str  = "r.t."
    catalyst:        bool = False


class FullRequest(BaseModel):
    """Full prediction + reference search + side products."""
    amine_smiles:        str
    squarate_smiles:     str  = "CCOC1=C(OCC)C(=O)C1=O"
    temperature:         str  = "r.t."
    catalyst:            bool = False
    solvent:             str  = "EtOH"
    amine_identifier:    str  = ""   # optional: raw user input for compound lookup
    squarate_identifier: str  = ""


class LegacyRequest(BaseModel):
    amine_type:  str
    substituent: str
    position:    str
    temperature: str
    catalyst:    bool


class ValidateRequest(BaseModel):
    smiles: str


class LookupRequest(BaseModel):
    identifier: str   # SMILES, CAS, IUPAC name, or common name


# ── Utility ───────────────────────────────────────────────────────────────────

def _confidence_summary(prediction: dict, references: dict) -> dict:
    """
    Build a human-readable confidence summary with source attribution.
    """
    match_type = prediction.get("match_type", "inferred")
    exact      = references["local"]["exact"]
    similar    = references["local"]["similar"]

    if exact:
        level   = "high"
        label   = "High — exact match found in internal dataset"
        explain = (
            f"This exact reaction (same amine, same squarate, same solvent) has been "
            f"reported experimentally. Yield data comes from the internal dataset "
            f"({references['local']['sources'][0] if references['local']['sources'] else 'internal'})."
        )
    elif similar and any(r["similarity"] >= 0.70 for r in similar):
        level   = "medium"
        label   = "Medium — similar reactions found in internal dataset"
        explain = (
            "No exact match was found, but structurally similar reactions are available "
            "in the internal dataset. The predicted yield is based on k-NN extrapolation "
            "from these analogues."
        )
    elif similar:
        level   = "medium"
        label   = "Medium — analogous reactions found (lower similarity)"
        explain = (
            "Reactions with moderate structural similarity were found. The predicted yield "
            "is an extrapolation with greater uncertainty than a high-similarity match."
        )
    else:
        level   = "low"
        label   = "Low — no reference reactions found"
        explain = (
            "No structurally similar reactions were found in the internal dataset. "
            "The predicted yield is based on chemical class rules only and should be "
            "treated as a rough estimate."
        )

    sources_used = []
    sources_used.append("RDKit (structure validation and Morgan fingerprint matching)")
    if references["local"]["sources"]:
        sources_used += references["local"]["sources"]
    sources_used.append("PubChem (compound validation and CAS lookup)")
    sources_used.append("OPSIN (IUPAC name conversion, if applicable)")
    sources_used.append("ORD — not yet connected (architecture ready)")
    sources_used.append("USPTO patent dataset — not yet connected (architecture ready)")

    return {
        "level":        level,
        "label":        label,
        "explanation":  explain,
        "sources_used": sources_used,
    }


def _status_label(prediction: dict, refs: dict) -> str:
    y = prediction.get("yield_percent", 0)
    exact = refs["local"]["exact"]
    if exact:
        avg = sum(r["yield_percent"] for r in exact) / len(exact)
        if avg >= 50: return "Likely to work"
        if avg >= 20: return "May work — moderate yield"
        return "Unlikely to work"
    if y >= 50:  return "Likely to work"
    if y >= 20:  return "May work — moderate yield expected"
    return "Unlikely to work"


# ── Endpoints ─────────────────────────────────────────────────────────────────

@app.post("/predict/full", tags=["v3"])
async def predict_full(req: FullRequest):
    """
    Comprehensive prediction endpoint (v3).

    Returns:
      amine_validation    — PubChem + OPSIN + RDKit compound check
      squarate_validation — same for the squarate
      references          — local dataset + external stubs
      prediction          — k-NN yield + product SMILES
      side_products       — rule-based side product list
      confidence          — scored + explained
      summary             — status label + yield breakdown
    """
    # ── 1. Compound validation (async, runs in parallel) ───────────────────
    import asyncio
    amine_id    = req.amine_identifier or req.amine_smiles
    squarate_id = req.squarate_identifier or req.squarate_smiles

    amine_val, sq_val = await asyncio.gather(
        validate_compound(amine_id),
        validate_compound(squarate_id),
    )

    # ── 2. Core prediction ────────────────────────────────────────────────
    pred = predict_from_smiles(
        amine_smiles    = req.amine_smiles,
        squarate_smiles = req.squarate_smiles,
        temperature     = req.temperature,
        catalyst        = req.catalyst,
    )

    # ── 3. Reference search ───────────────────────────────────────────────
    refs = search_all_references(req.amine_smiles)

    # ── 4. Side products ──────────────────────────────────────────────────
    at  = pred.get("amine_matched", "").split(" / ")[0] if pred.get("amine_matched") else "aniline"
    sp  = predict_side_products(
        amine_smiles     = req.amine_smiles,
        amine_type       = at,
        temperature      = req.temperature,
        catalyst         = req.catalyst,
        solvent          = req.solvent,
        monoamide_smiles = pred.get("product_smiles", ""),
    )

    # ── 5. Confidence + summary ───────────────────────────────────────────
    conf    = _confidence_summary(pred, refs)
    status  = _status_label(pred, refs)

    # Yield breakdown
    exact_yields = [r["yield_percent"] for r in refs["local"]["exact"] if r["yield_percent"] > 0]
    sim_yields   = [r["yield_percent"] for r in refs["local"]["similar"] if r["yield_percent"] > 0]

    yield_display = {
        "reported":        round(sum(exact_yields)/len(exact_yields), 1) if exact_yields else None,
        "similar_min":     min(sim_yields) if sim_yields else None,
        "similar_max":     max(sim_yields) if sim_yields else None,
        "predicted":       pred.get("yield_percent"),
        "reported_source": refs["local"]["sources"][0] if exact_yields and refs["local"]["sources"] else None,
    }

    return {
        "amine_validation":    amine_val,
        "squarate_validation": sq_val,
        "references":          refs,
        "prediction":          pred,
        "side_products":       sp,
        "confidence":          conf,
        "summary": {
            "status":        status,
            "yield_display": yield_display,
        },
    }


@app.post("/lookup", tags=["v3"])
async def lookup_compound(req: LookupRequest):
    """
    Resolve any compound identifier to validated structure data.
    Accepts SMILES, CAS number, IUPAC name, or common name.
    """
    return await validate_compound(req.identifier)


@app.post("/predict/smiles", tags=["v2"])
def predict_smiles(req: SmilesRequest):
    """SMILES-based prediction (v2) — fast, no reference search."""
    result = predict_from_smiles(
        amine_smiles    = req.amine_smiles,
        squarate_smiles = req.squarate_smiles,
        temperature     = req.temperature,
        catalyst        = req.catalyst,
    )
    if result.get("error") == "invalid_smiles":
        raise HTTPException(status_code=422, detail="Invalid amine SMILES string.")
    return result


@app.post("/predict", tags=["v1-legacy"])
def predict_legacy(req: LegacyRequest):
    """Legacy dropdown-encoded prediction (v1) — React Native mobile app."""
    return _legacy_predict(
        amine_type  = req.amine_type,
        substituent = req.substituent,
        position    = req.position,
        temperature = req.temperature,
        catalyst    = req.catalyst,
    )


@app.post("/validate", tags=["utility"])
def validate_smiles_endpoint(req: ValidateRequest):
    """Validate a SMILES string, return canonical form + formula + MW."""
    mol = Chem.MolFromSmiles(req.smiles.strip())
    if mol is None:
        return {"valid": False, "canonical_smiles": None, "formula": None,
                "mol_weight": None, "error": "Could not parse SMILES."}
    return {
        "valid":           True,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "formula":         Chem.rdMolDescriptors.CalcMolFormula(mol),
        "mol_weight":      round(Descriptors.ExactMolWt(mol), 3),
        "error":           None,
    }


@app.get("/options", tags=["utility"])
def get_options():
    return {
        "amine_types": ["aniline", "benzylamine", "heterocyclic"],
        "substituents": {
            "aniline":      ["none", "F", "Cl", "Br", "OH", "OMe", "CF3", "tolyl"],
            "benzylamine":  ["none", "Cl", "Br", "fused_ring"],
            "heterocyclic": [
                "pyrrolidine", "azetidine", "piperidine", "azepane",
                "thiazolidine", "morpholine", "dimethylmorpholine",
                "thiomorpholine", "thiomorpholine_dioxide",
                "boc_piperazine", "boc_dimethylpiperazine",
                "bromo_tetrahydroquinoline", "piperazine"
            ],
        },
        "positions":    ["ortho", "meta", "para", "none"],
        "temperatures": ["r.t.", "reflux"],
    }


@app.get("/health", tags=["utility"])
def health():
    return {
        "status":       "ok",
        "version":      "3.0",
        "dataset_size": 42,
        "endpoints": [
            "/predict/full", "/predict/smiles", "/predict",
            "/lookup", "/validate", "/options",
        ],
    }
