# reference_search.py
# Reaction reference search for Chemdem.
#
# Search tiers (in order):
#   Tier 1 — Local dataset (Papers 1 & 2 — exact experimental data)
#   Tier 2 — PubChem compound cross-reference (CID-based reaction hints)
#   Tier 3 — ORD / USPTO stubs (architecture ready, datasets not yet downloaded)
#
# Each result is clearly labelled with its source and confidence tier.

import math
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from model import DATASET, COMPOUND_DB

# ── Build a search index from COMPOUND_DB at import time ─────────────────────

def _canonical(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None

def _morgan_fp(smiles, r=2, n=2048):
    mol = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(mol, r, nBits=n) if mol else None

# Pre-compute canonical SMILES + Morgan FPs for every training compound
_COMPOUND_INDEX = []   # (at, sub, pos, canonical_smiles, fp)
for _at, _sub, _pos, _smi in COMPOUND_DB:
    _can = _canonical(_smi)
    _fp  = _morgan_fp(_can) if _can else None
    if _can and _fp:
        _COMPOUND_INDEX.append((_at, _sub, _pos, _can, _fp))

# Map (amine_type, substituent, position) → DATASET rows
_DATASET_MAP = {}
for _row in DATASET:
    _key = (_row[0], _row[1], _row[2])
    _DATASET_MAP.setdefault(_key, []).append(_row)

# Source metadata
_LOCAL_SOURCES = {
    "paper_1": {
        "name": "Long et al., SynOpen 2023",
        "type": "Internal dataset — peer-reviewed journal",
        "doi":  "10.1055/a-XXXX-XXXX",
        "url":  "https://www.thieme-connect.com/products/ejournals/journal/10.1055/s-00000121",
    },
    "paper_2": {
        "name": "Long et al., Bioorg. Med. Chem. 2024",
        "type": "Internal dataset — peer-reviewed journal",
        "doi":  None,
        "url":  None,
    },
}

# Paper 2 heterocyclic entries (amine_type == 'heterocyclic')
_PAPER2_SUBS = {
    "pyrrolidine", "azetidine", "piperidine", "azepane", "thiazolidine",
    "morpholine", "dimethylmorpholine", "thiomorpholine", "thiomorpholine_dioxide",
    "boc_piperazine", "boc_dimethylpiperazine", "bromo_tetrahydroquinoline", "piperazine",
}


def _paper_for_row(row) -> str:
    return "paper_2" if row[1] in _PAPER2_SUBS else "paper_1"


def _row_to_ref(row, sim: float, match_tier: str) -> dict:
    """Convert a DATASET row to a structured reference dict."""
    at, sub, pos, temp_i, cat_i, yld = row
    paper_key = _paper_for_row(row)
    src = _LOCAL_SOURCES[paper_key]
    return {
        "amine_type":   at,
        "substituent":  sub,
        "position":     pos,
        "temperature":  "reflux" if temp_i else "r.t.",
        "catalyst":     "Zn(OTf)₂ (10–13 mol%)" if cat_i else "none",
        "yield_percent": yld,
        "solvent":      "EtOH",   # all training data used EtOH
        "match_tier":   match_tier,   # 'exact' | 'high' | 'medium' | 'low'
        "similarity":   round(sim, 3),
        "source": {
            "label":    src["name"],
            "type":     src["type"],
            "doi":      src["doi"],
        },
    }


# ── Tier 1: Local dataset search ─────────────────────────────────────────────

def search_local(amine_smiles: str, top_k: int = 5) -> dict:
    """
    Search the local training dataset for reactions involving amines
    similar to the query.

    Returns:
        exact     list[dict]   — same canonical SMILES in training set
        similar   list[dict]   — Tanimoto-sorted similar amines (>0.30)
        sources   list[str]    — source labels used
    """
    can = _canonical(amine_smiles)
    q_fp = _morgan_fp(can) if can else None

    exact_rows   = []
    similar_hits = []   # (similarity, row)

    for at, sub, pos, db_can, db_fp in _COMPOUND_INDEX:
        # Exact match
        if can and can == db_can:
            key = (at, sub, pos)
            for row in _DATASET_MAP.get(key, []):
                exact_rows.append(_row_to_ref(row, 1.0, "exact"))
            continue

        # Similarity search
        if q_fp is not None and db_fp is not None:
            sim = DataStructs.TanimotoSimilarity(q_fp, db_fp)
            if sim >= 0.30:
                key = (at, sub, pos)
                for row in _DATASET_MAP.get(key, []):
                    similar_hits.append((sim, row))

    # Sort similar hits by decreasing similarity, deduplicate
    similar_hits.sort(key=lambda x: -x[0])
    seen = set()
    similar_refs = []
    for sim, row in similar_hits:
        k = (row[0], row[1], row[2], row[3], row[4])
        if k not in seen:
            seen.add(k)
            tier = "high" if sim >= 0.70 else "medium" if sim >= 0.50 else "low"
            similar_refs.append(_row_to_ref(row, sim, tier))
        if len(similar_refs) >= top_k:
            break

    sources = []
    if exact_rows or similar_refs:
        sources = [
            "Long et al., SynOpen 2023 (internal dataset)",
            "Long et al., Bioorg. Med. Chem. 2024 (internal dataset)",
        ]

    return {
        "exact":   exact_rows,
        "similar": similar_refs,
        "sources": sources,
    }


def compute_yield_range(refs: list[dict]) -> dict | None:
    """Compute min / max / mean yield from a list of reference dicts."""
    yields = [r["yield_percent"] for r in refs if r["yield_percent"] > 0]
    if not yields:
        return None
    return {
        "min":  min(yields),
        "max":  max(yields),
        "mean": round(sum(yields) / len(yields), 1),
        "n":    len(yields),
    }


# ── Tier 2: PubChem reaction hints ───────────────────────────────────────────
# PubChem does not expose a general reaction search API.
# We use it only for compound-level cross-referencing (done in lookup.py).
# This stub is here for architectural completeness.

def search_pubchem_reactions(amine_cid: int | None) -> dict:
    """
    Placeholder — PubChem BioAssay / Patent cross-reference.
    Returns empty results with a clear status note.
    """
    return {
        "status":  "not_implemented",
        "message": (
            "PubChem does not provide a general reaction search API. "
            "Compound-level validation (CAS, SMILES, IUPAC) is performed via PubChem; "
            "reaction references come from the internal dataset."
        ),
        "results": [],
        "source": {
            "label": "PubChem",
            "type":  "compound_database_only",
            "url":   "https://pubchem.ncbi.nlm.nih.gov",
        },
    }


# ── Tier 3: ORD stub ──────────────────────────────────────────────────────────
# The Open Reaction Database is available at https://open-reaction-database.org
# and on GitHub: https://github.com/open-reaction-database/ord-data
# The full dataset is ~50 GB in protobuf format and requires the ord-schema package.
#
# To enable this tier:
#   1.  pip install ord-schema
#   2.  Download the squaric acid reaction subset from ORD (filter by reaction SMARTS)
#   3.  Convert to JSON and load here
#   4.  Implement the search function below

def search_ord(amine_smiles: str) -> dict:
    """
    Placeholder — Open Reaction Database search.
    Requires the ORD dataset to be downloaded and indexed locally.
    """
    return {
        "status":  "not_connected",
        "message": (
            "ORD search is not yet connected. "
            "To enable: download the ORD dataset from https://github.com/open-reaction-database/ord-data, "
            "filter for squaric acid reactions, and load the JSON index here."
        ),
        "results": [],
        "source": {
            "label": "Open Reaction Database (ORD)",
            "type":  "open_reaction_database",
            "url":   "https://open-reaction-database.org",
        },
    }


# ── Tier 4: USPTO stub ────────────────────────────────────────────────────────
# The USPTO patent reaction dataset (Lowe, 2012/2017) contains ~1M reactions
# extracted from US patents.  Accessing it requires downloading the CML/RXN files.
#
# To enable this tier:
#   1.  Download from https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents/5104873
#   2.  Extract and index reactions containing squaric acid cores
#   3.  Implement a SMARTS-based similarity search here

def search_uspto(amine_smiles: str) -> dict:
    """
    Placeholder — USPTO patent reaction dataset search.
    Requires the dataset to be downloaded and indexed locally.
    """
    return {
        "status":  "not_connected",
        "message": (
            "USPTO patent reaction search is not yet connected. "
            "To enable: download the Lowe USPTO reaction dataset from Figshare "
            "(https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents/5104873) "
            "and index the squaric acid subset."
        ),
        "results": [],
        "source": {
            "label": "USPTO Patent Reaction Dataset (Lowe 2012/2017)",
            "type":  "patent_database",
            "url":   "https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents/5104873",
        },
    }


# ── Main search function ──────────────────────────────────────────────────────

def search_all_references(amine_smiles: str) -> dict:
    """
    Run all search tiers and return a unified reference report.

    Keys:
        local          dict   — Tier 1 results
        pubchem        dict   — Tier 2 (compound-level only)
        ord            dict   — Tier 3 (stub)
        uspto          dict   — Tier 4 (stub)
        reported_yield dict | None  — statistics over exact + high-similarity matches
        best_source    str    — highest-quality source that returned data
    """
    local   = search_local(amine_smiles)
    pc      = search_pubchem_reactions(None)
    ord_res = search_ord(amine_smiles)
    uspto   = search_uspto(amine_smiles)

    # Aggregate for yield statistics — exact matches first, then high-similarity
    primary = local["exact"] or local["similar"][:3]
    yield_stats = compute_yield_range(primary)

    # Best source label
    if local["exact"]:
        best_source = "Internal dataset — exact match"
    elif local["similar"]:
        best_source = "Internal dataset — similar reactions"
    else:
        best_source = "No reference reactions found"

    return {
        "local":          local,
        "pubchem":        pc,
        "ord":            ord_res,
        "uspto":          uspto,
        "reported_yield": yield_stats,
        "best_source":    best_source,
    }
