# lookup.py
# Compound validation pipeline for Chemdem
# Sources (in order of query):
#   1. RDKit   — SMILES parsing & canonicalisation
#   2. OPSIN   — IUPAC name → SMILES  (opsin.ch.cam.ac.uk)
#   3. PubChem — name / CAS / SMILES lookup + CAS synonym extraction
#
# All data is labelled with its source.  Mismatches between sources trigger a warning.

import re
import httpx
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
OPSIN_BASE   = "https://www.ebi.ac.uk/opsin/ws"   # opsin.ch.cam.ac.uk → EBI (permanent redirect)
_TIMEOUT     = 8.0   # seconds per HTTP call
_PC_PROPS    = "MolecularFormula,MolecularWeight,ConnectivitySMILES,IUPACName,CID"

_CAS_RE    = re.compile(r"^\d{2,7}-\d{2}-\d$")
_SMILES_RE = re.compile(r"[=#@\[\]()\\\/]|(?<![A-Z])[cnosp](?![a-zA-Z])")


# ── Heuristics ────────────────────────────────────────────────────────────────

def _looks_like_smiles(s: str) -> bool:
    return bool(_SMILES_RE.search(s)) or (len(s) > 3 and any(c.isdigit() for c in s) and s[0].isupper())

def _looks_like_cas(s: str) -> bool:
    return bool(_CAS_RE.match(s.strip()))

def _canonical_or_none(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None

def _rdkit_props(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        "canonical_smiles": Chem.MolToSmiles(mol),
        "formula":          rdMolDescriptors.CalcMolFormula(mol),
        "mol_weight":       round(Descriptors.ExactMolWt(mol), 3),
    }


# ── HTTP helpers (all async, use httpx) ──────────────────────────────────────

async def _opsin_convert(name: str) -> str | None:
    """IUPAC name → SMILES via OPSIN web API (EBI mirror)."""
    try:
        async with httpx.AsyncClient(timeout=_TIMEOUT, follow_redirects=True) as c:
            r = await c.get(
                f"{OPSIN_BASE}/{name}",
                headers={"Accept": "chemical/x-daylight-smiles"},
            )
            if r.status_code == 200 and r.text.strip():
                return _canonical_or_none(r.text.strip())
    except Exception:
        pass
    return None


async def _pubchem_fetch(url: str, **kwargs) -> dict | None:
    """GET/POST PubChem REST, return first Properties item or None."""
    try:
        async with httpx.AsyncClient(timeout=_TIMEOUT) as c:
            r = await (c.post(url, **kwargs) if kwargs else c.get(url))
            if r.status_code != 200:
                return None
            body = r.json()
            props = (
                body.get("PropertyTable", {}).get("Properties") or
                body.get("InformationList", {}).get("Information") or
                []
            )
            return props[0] if props else None
    except Exception:
        return None


async def _pubchem_by_name(name: str) -> dict | None:
    url = f"{PUBCHEM_BASE}/compound/name/{name}/property/{_PC_PROPS}/JSON"
    return await _pubchem_fetch(url)


async def _pubchem_by_smiles(smiles: str) -> dict | None:
    url = f"{PUBCHEM_BASE}/compound/smiles/property/{_PC_PROPS}/JSON"
    return await _pubchem_fetch(url, data={"smiles": smiles})


async def _pubchem_get_cas(cid: int) -> str | None:
    """Fetch CAS number from PubChem synonym list."""
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON"
    data = await _pubchem_fetch(url)
    if not data:
        return None
    for syn in data.get("Synonym", []):
        if _CAS_RE.match(str(syn)):
            return syn
    return None


# ── Main validation function ──────────────────────────────────────────────────

async def validate_compound(identifier: str) -> dict:
    """
    Resolve any identifier (SMILES, CAS, IUPAC name, common name) to a
    validated compound record.

    Returns:
        valid            bool
        input_type       'smiles' | 'cas' | 'name'
        canonical_smiles str | None
        iupac_name       str | None
        cas              str | None
        formula          str | None
        mol_weight       float | None
        sources          list[{name, type, detail}]
        warnings         list[str]
        cross_checks     list[{check, status, note}]
    """
    s = identifier.strip()

    result: dict = {
        "input":           s,
        "valid":           False,
        "input_type":      None,
        "canonical_smiles": None,
        "iupac_name":      None,
        "cas":             None,
        "formula":         None,
        "mol_weight":      None,
        "sources":         [],
        "warnings":        [],
        "cross_checks":    [],
    }

    pc_data       = None
    opsin_smiles  = None
    rdkit_smiles  = None

    # ── Branch 1: SMILES ─────────────────────────────────────────────────────
    if _looks_like_smiles(s):
        mol = Chem.MolFromSmiles(s)
        if mol:
            result["input_type"] = "smiles"
            props = _rdkit_props(s)
            result.update(props)
            rdkit_smiles = result["canonical_smiles"]
            result["sources"].append({
                "name":   "RDKit",
                "type":   "smiles_validation",
                "detail": "SMILES parsed and canonicalised",
            })
            pc_data = await _pubchem_by_smiles(rdkit_smiles)

    # ── Branch 2: CAS number ─────────────────────────────────────────────────
    elif _looks_like_cas(s):
        result["input_type"] = "cas"
        result["cas"]        = s
        pc_data              = await _pubchem_by_name(s)

    # ── Branch 3: Name (IUPAC or common) ─────────────────────────────────────
    else:
        result["input_type"] = "name"
        opsin_smiles = await _opsin_convert(s)
        if opsin_smiles:
            props = _rdkit_props(opsin_smiles)
            result.update(props)
            result["iupac_name"] = s
            result["sources"].append({
                "name":   "OPSIN",
                "type":   "iupac_conversion",
                "detail": f"IUPAC → SMILES: {opsin_smiles}",
            })
        pc_data = await _pubchem_by_name(s)

    # ── Merge PubChem data ────────────────────────────────────────────────────
    if pc_data:
        pc_smi    = pc_data.get("ConnectivitySMILES") or pc_data.get("CanonicalSMILES", "")
        pc_iupac  = pc_data.get("IUPACName", "")
        pc_cid    = pc_data.get("CID")
        pc_formula = pc_data.get("MolecularFormula", "")
        pc_mw     = pc_data.get("MolecularWeight")

        result["sources"].append({
            "name":   "PubChem",
            "type":   "database_lookup",
            "detail": f"CID {pc_cid}" if pc_cid else "name/SMILES search",
        })

        # Fill missing fields
        if not result.get("formula"):
            result["formula"] = pc_formula
        if not result.get("mol_weight") and pc_mw:
            result["mol_weight"] = float(pc_mw)
        if not result.get("iupac_name") and pc_iupac:
            result["iupac_name"] = pc_iupac

        # If we had no SMILES yet, use PubChem's
        if not result.get("canonical_smiles") and pc_smi:
            props = _rdkit_props(pc_smi)
            result.update(props)

        # Cross-check SMILES consistency
        if rdkit_smiles and pc_smi:
            a = _canonical_or_none(rdkit_smiles)
            b = _canonical_or_none(pc_smi)
            if a and b:
                if a == b:
                    result["cross_checks"].append({
                        "check":  "SMILES consistency (RDKit vs PubChem)",
                        "status": "pass",
                        "note":   "Both sources agree on structure",
                    })
                else:
                    result["warnings"].append(
                        f"Structure mismatch: your input resolves to {a}, "
                        f"but PubChem ({pc_iupac or 'name search'}) has {b}. "
                        "Verify the identifier."
                    )
                    result["cross_checks"].append({
                        "check":  "SMILES consistency (RDKit vs PubChem)",
                        "status": "warning",
                        "note":   f"Input: {a} | PubChem: {b}",
                    })

        # OPSIN vs PubChem cross-check
        if opsin_smiles and pc_smi:
            a = _canonical_or_none(opsin_smiles)
            b = _canonical_or_none(pc_smi)
            if a and b:
                status = "pass" if a == b else "warning"
                result["cross_checks"].append({
                    "check":  "SMILES consistency (OPSIN vs PubChem)",
                    "status": status,
                    "note":   "Agree" if a == b else f"OPSIN: {a} | PubChem: {b}",
                })
                if a != b:
                    result["warnings"].append(
                        f"OPSIN parsed '{s}' as {a}, but PubChem has {b}. "
                        "The name may be ambiguous."
                    )

        # Fetch CAS if not already known
        if not result.get("cas") and pc_cid:
            cas = await _pubchem_get_cas(int(pc_cid))
            if cas:
                result["cas"] = cas
                result["sources"].append({
                    "name":   "PubChem Synonyms",
                    "type":   "cas_extraction",
                    "detail": f"CAS {cas} from synonym list (CID {pc_cid})",
                })

        result["valid"] = True

    elif result.get("canonical_smiles"):
        result["valid"] = True
    else:
        result["warnings"].append(
            f"Could not resolve '{s}' from RDKit, OPSIN, or PubChem. "
            "Check the spelling or try a CAS number / SMILES."
        )

    return result
