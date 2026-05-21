# model.py  v2.0
# SMILES-based prediction engine for squaric acid monoamide synthesis
# Accepts SMILES inputs → RDKit validation → Morgan FP matching → k-NN → product SMILES
#
# Sources:
#   Paper 1: Long et al., SynOpen 2023      (aniline/benzylamine series, 29 reactions)
#   Paper 2: Long et al., Bioorg Med Chem 2024  (heterocyclic series, 13 reactions)

import math
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning, module="rdkit")

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


# ── Training Dataset ──────────────────────────────────────────────────────────
# (amine_type, substituent, position, temperature, catalyst, yield_pct)
# temperature: 0 = r.t.   1 = reflux
# catalyst:    0 = none   1 = Zn(OTf)2

DATASET = [
    # ── Paper 1: Aniline / Benzylamine series ────────────────────────────────
    ("benzylamine",  "none",                    "none",  0, 0, 33),
    ("aniline",      "none",                    "none",  0, 1, 37),
    ("aniline",      "F",                       "para",  0, 1, 75),
    ("aniline",      "Cl",                      "meta",  0, 1, 26),
    ("aniline",      "Cl",                      "para",  0, 1, 60),
    ("aniline",      "Br",                      "meta",  0, 1, 64),
    ("aniline",      "Br",                      "para",  0, 1, 66),
    ("aniline",      "OH",                      "para",  0, 0, 44),
    ("aniline",      "OH",                      "meta",  0, 0, 54),
    ("aniline",      "OH",                      "para",  0, 0, 45),
    ("aniline",      "OMe",                     "ortho", 0, 1, 18),
    ("aniline",      "OMe",                     "para",  0, 1, 94),
    ("aniline",      "OMe",                     "meta",  0, 1, 49),
    ("aniline",      "tolyl",                   "para",  0, 0, 29),
    ("aniline",      "tolyl",                   "meta",  0, 0, 89),
    ("aniline",      "tolyl",                   "ortho", 0, 0, 75),
    ("aniline",      "CF3",                     "para",  0, 0, 54),
    ("aniline",      "CF3",                     "meta",  0, 1, 48),
    ("aniline",      "Cl",                      "ortho", 0, 1, 38),
    ("benzylamine",  "Cl",                      "ortho", 0, 0, 25),
    ("benzylamine",  "Cl",                      "para",  1, 0, 99),
    ("benzylamine",  "Br",                      "ortho", 0, 0, 34),
    ("benzylamine",  "Br",                      "para",  0, 0, 57),
    ("benzylamine",  "Br",                      "meta",  0, 0, 41),
    ("benzylamine",  "Br",                      "para",  0, 0, 91),
    ("aniline",      "OH",                      "ortho", 1, 1, 31),
    ("benzylamine",  "fused_ring",              "none",  1, 0, 48),
    # Paper 1 — confirmed FAILS
    ("aniline",      "Cl",                      "ortho", 1, 1,  0),
    ("aniline",      "Br",                      "ortho", 1, 1,  0),

    # ── Paper 2: Heterocyclic amine series ───────────────────────────────────
    ("heterocyclic", "pyrrolidine",             "none",  0, 0, 10),
    ("heterocyclic", "azetidine",               "none",  0, 0, 40),
    ("heterocyclic", "piperidine",              "none",  0, 1, 54),
    ("heterocyclic", "azepane",                 "none",  1, 0, 21),
    ("heterocyclic", "thiazolidine",            "none",  0, 0, 21),
    ("heterocyclic", "morpholine",              "none",  0, 0, 72),
    ("heterocyclic", "dimethylmorpholine",      "none",  0, 0, 76),
    ("heterocyclic", "thiomorpholine",          "none",  1, 0, 95),
    ("heterocyclic", "thiomorpholine_dioxide",  "none",  1, 0, 27),
    ("heterocyclic", "boc_piperazine",          "none",  0, 0, 73),
    ("heterocyclic", "boc_dimethylpiperazine",  "none",  0, 0, 59),
    ("heterocyclic", "bromo_tetrahydroquinoline","none", 0, 0, 21),
    # Paper 2 — confirmed FAIL
    ("heterocyclic", "piperazine",              "none",  0, 0,  0),
]


# ── Compound Database (label → SMILES) ───────────────────────────────────────
# Maps (amine_type, substituent, position) to a canonical SMILES for Morgan FP matching.
# SMILES are PubChem-verified (ConnectivitySMILES, validated in enrich_dataset.py).

COMPOUND_DB = [
    # (amine_type, substituent, position, smiles)
    # ── Paper 1 ──────────────────────────────────────────────────────────────
    ("benzylamine",  "none",                    "none",  "NCC1=CC=CC=C1"),
    ("aniline",      "none",                    "none",  "Nc1ccccc1"),
    ("aniline",      "F",                       "para",  "Nc1ccc(F)cc1"),
    ("aniline",      "Cl",                      "meta",  "Nc1cccc(Cl)c1"),
    ("aniline",      "Cl",                      "para",  "Nc1ccc(Cl)cc1"),
    ("aniline",      "Cl",                      "ortho", "Nc1ccccc1Cl"),
    ("aniline",      "Br",                      "meta",  "Nc1cccc(Br)c1"),
    ("aniline",      "Br",                      "para",  "Nc1ccc(Br)cc1"),
    ("aniline",      "Br",                      "ortho", "Nc1ccccc1Br"),
    ("aniline",      "OH",                      "para",  "Nc1ccc(O)cc1"),
    ("aniline",      "OH",                      "meta",  "Nc1cccc(O)c1"),
    ("aniline",      "OH",                      "ortho", "Nc1ccccc1O"),
    ("aniline",      "OMe",                     "ortho", "Nc1ccccc1OC"),
    ("aniline",      "OMe",                     "para",  "Nc1ccc(OC)cc1"),
    ("aniline",      "OMe",                     "meta",  "Nc1cccc(OC)c1"),
    ("aniline",      "tolyl",                   "para",  "Nc1ccc(C)cc1"),
    ("aniline",      "tolyl",                   "meta",  "Nc1cccc(C)c1"),
    ("aniline",      "tolyl",                   "ortho", "Nc1ccccc1C"),
    ("aniline",      "CF3",                     "para",  "Nc1ccc(C(F)(F)F)cc1"),
    ("aniline",      "CF3",                     "meta",  "Nc1cccc(C(F)(F)F)c1"),
    ("benzylamine",  "Cl",                      "ortho", "NCC1=CC=CC=C1Cl"),
    ("benzylamine",  "Cl",                      "para",  "NCC1=CC=C(Cl)C=C1"),
    ("benzylamine",  "Br",                      "ortho", "NCC1=CC=CC=C1Br"),
    ("benzylamine",  "Br",                      "para",  "NCC1=CC=C(Br)C=C1"),
    ("benzylamine",  "Br",                      "meta",  "NCC1=CC=CC(Br)=C1"),
    # ── Paper 2 ──────────────────────────────────────────────────────────────
    ("heterocyclic", "pyrrolidine",             "none",  "C1CCNC1"),
    ("heterocyclic", "azetidine",               "none",  "C1CNC1"),
    ("heterocyclic", "piperidine",              "none",  "C1CCNCC1"),
    ("heterocyclic", "azepane",                 "none",  "C1CCNCCC1"),
    ("heterocyclic", "thiazolidine",            "none",  "C1CNCS1"),
    ("heterocyclic", "morpholine",              "none",  "C1CNCCO1"),
    ("heterocyclic", "dimethylmorpholine",      "none",  "CC1CNCC(C)O1"),
    ("heterocyclic", "thiomorpholine",          "none",  "C1CNCCS1"),
    ("heterocyclic", "thiomorpholine_dioxide",  "none",  "O=S1(=O)CCNCC1"),
    ("heterocyclic", "boc_piperazine",          "none",  "CC(C)(C)OC(=O)N1CCNCC1"),
    ("heterocyclic", "boc_dimethylpiperazine",  "none",  "CC(C)(C)OC(=O)N1CC(C)NCC1"),
    ("heterocyclic", "bromo_tetrahydroquinoline","none", "Brc1ccc2c(c1)CCNC2"),
    ("heterocyclic", "piperazine",              "none",  "C1CNCCN1"),
]


# ── Feature Encoding ──────────────────────────────────────────────────────────

AMINE_ENC      = {"aniline": 0, "benzylamine": 1, "heterocyclic": 2}
POSITION_ENC   = {"ortho": 0, "meta": 1, "para": 2, "none": 2}
SUBSTITUENT_ENC = {
    "none": 0, "F": 1, "Cl": 2, "Br": 3, "OH": 4,
    "OMe": 5, "CF3": 6, "tolyl": 7, "fused_ring": 8,
    "pyrrolidine": 9, "azetidine": 10, "piperidine": 11, "azepane": 12,
    "thiazolidine": 13, "morpholine": 14, "dimethylmorpholine": 15,
    "thiomorpholine": 16, "thiomorpholine_dioxide": 17,
    "boc_piperazine": 18, "boc_dimethylpiperazine": 19,
    "bromo_tetrahydroquinoline": 20, "piperazine": 21,
}


def _encode_categorical(amine_type, substituent, position, temperature, catalyst):
    return [
        AMINE_ENC.get(amine_type, 0),
        SUBSTITUENT_ENC.get(substituent, 0),
        POSITION_ENC.get(position, 2),
        int(temperature == "reflux"),
        int(bool(catalyst)),
    ]


# ── Pre-compute Morgan fingerprints for COMPOUND_DB at module load ────────────

def _canonical(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def _morgan_fp(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)


# Index: (amine_type, substituent, position, canonical_smiles, morgan_fp, enc[])
_DB_INDEX = []
for _at, _sub, _pos, _smi in COMPOUND_DB:
    _can = _canonical(_smi)
    _fp  = _morgan_fp(_can) if _can else None
    if _can and _fp:
        _enc = _encode_categorical(_at, _sub, _pos, "r.t.", False)  # temp/cat slot filled at query time
        _DB_INDEX.append((_at, _sub, _pos, _can, _fp, _enc))


# ── Hard-Fail Detection (SMARTS) ──────────────────────────────────────────────

_HARD_FAIL_PATTERNS = [
    (
        Chem.MolFromSmarts("Nc1ccccc1Cl"),       # 2-chloroaniline
        "Ortho-Cl aniline fails due to steric crowding + poor nucleophilicity. "
        "No reaction even at reflux + 48 h + Zn(OTf)₂ catalyst.",
        "Switch to meta or para position for significantly higher yield.",
    ),
    (
        Chem.MolFromSmarts("Nc1ccccc1Br"),       # 2-bromoaniline
        "Ortho-Br aniline fails due to steric crowding + poor nucleophilicity. "
        "No reaction even at reflux + 48 h + Zn(OTf)₂ catalyst.",
        "Switch to meta or para position for significantly higher yield.",
    ),
    (
        Chem.MolFromSmarts("[NH]1CC[NH]CC1"),    # unprotected piperazine
        "Unprotected piperazine fails — both N atoms compete, "
        "giving unwanted double substitution rather than the monoamide.",
        "Use Boc-protected piperazine instead (N-Boc-piperazine gives 73% yield).",
    ),
]


def _check_hard_fail(amine_smiles: str):
    """Returns (is_fail: bool, message: str, recommendation: str)."""
    mol = Chem.MolFromSmiles(amine_smiles)
    if mol is None:
        return False, "", ""
    for pat, msg, rec in _HARD_FAIL_PATTERNS:
        if pat and mol.HasSubstructMatch(pat):
            return True, msg, rec
    return False, "", ""


# ── Amine Class Inference (SMARTS fallback) ───────────────────────────────────

_PAT_ANILINE     = Chem.MolFromSmarts("Nc1ccccc1")
_PAT_BENZYLAMINE = Chem.MolFromSmarts("NCc1ccccc1")
_PAT_NH2         = Chem.MolFromSmarts("[NH2]")
_PAT_NH          = Chem.MolFromSmarts("[NH1;R]")  # ring NH


def _infer_amine_class(smiles: str):
    """
    Rough SMARTS classification for compounds not in COMPOUND_DB.
    Returns (amine_type, substituent, position).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "aniline", "none", "none"
    if mol.HasSubstructMatch(_PAT_BENZYLAMINE):
        return "benzylamine", "none", "none"
    if mol.HasSubstructMatch(_PAT_ANILINE):
        return "aniline", "none", "none"
    if mol.HasSubstructMatch(_PAT_NH):
        return "heterocyclic", "morpholine", "none"   # closest generic cyclic amine
    return "aniline", "none", "none"


# ── Morgan-FP Similarity Matching ────────────────────────────────────────────

def _match_to_training(amine_smiles: str):
    """
    Find the closest training compound by Morgan fingerprint Tanimoto similarity.
    Returns (amine_type, substituent, position, match_type, similarity, matched_name).
    match_type: 'exact' | 'similar' | 'inferred'
    """
    can = _canonical(amine_smiles)
    if can is None:
        return None, None, None, "invalid", 0.0, "—"

    # 1. Exact canonical SMILES match
    for at, sub, pos, db_can, db_fp, _ in _DB_INDEX:
        if can == db_can:
            label = f"{at} / {sub}" + (f" / {pos}" if pos != "none" else "")
            return at, sub, pos, "exact", 1.0, label

    # 2. Morgan FP Tanimoto similarity
    q_fp = _morgan_fp(can)
    if q_fp is None:
        at, sub, pos = _infer_amine_class(can)
        return at, sub, pos, "inferred", 0.0, f"{at} (inferred)"

    best_sim, best_at, best_sub, best_pos, best_label = -1.0, None, None, None, "—"
    for at, sub, pos, db_can, db_fp, _ in _DB_INDEX:
        sim = DataStructs.TanimotoSimilarity(q_fp, db_fp)
        if sim > best_sim:
            best_sim = sim
            best_at, best_sub, best_pos = at, sub, pos
            best_label = f"{at} / {sub}" + (f" / {pos}" if pos != "none" else "")

    if best_sim >= 0.30:
        return best_at, best_sub, best_pos, "similar", round(best_sim, 3), best_label

    # 3. SMARTS fallback
    at, sub, pos = _infer_amine_class(can)
    return at, sub, pos, "inferred", round(best_sim, 3), f"{at} (inferred)"


# ── k-NN Yield Estimator ──────────────────────────────────────────────────────

def _euclidean(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _knn_predict(query_features, k=3):
    """
    Inverse-distance weighted k-NN regressor over DATASET.
    Returns (predicted_yield: float, neighbours: list).
    """
    encoded = []
    for row in DATASET:
        at, sub, pos, temp_i, cat_i, yld = row
        feats = _encode_categorical(
            at, sub, pos,
            "reflux" if temp_i else "r.t.",
            bool(cat_i),
        )
        encoded.append((feats, yld, row))

    distances = [(
        _euclidean(query_features, feats), yld, row
    ) for feats, yld, row in encoded]
    distances.sort(key=lambda x: x[0])
    neighbours = distances[:k]

    total_w, weighted_sum = 0.0, 0.0
    for d, yld, _ in neighbours:
        w = 1.0 / (d + 0.001)
        weighted_sum += w * yld
        total_w += w

    return round(weighted_sum / total_w, 1), neighbours


# ── Product SMILES via Reaction SMARTS ───────────────────────────────────────

# Squarate ring SMARTS — RDKit perceives the squarate ring as aromatic (bond 1.5),
# so we use lowercase [c] for ring atoms.
# C:2 carries OCC that is displaced; C:3 retains its OCC; C:4/C:5 are ketones.
_RXN_PRIMARY = AllChem.ReactionFromSmarts(
    "[NH2:1].[c:2]1([O:6]CC)[c:3](OCC)[c:4](=O)[c:5]1=O"
    ">>"
    "[NH1:1][c:2]1[c:3](OCC)[c:4](=O)[c:5]1=O"
)
_RXN_SECONDARY = AllChem.ReactionFromSmarts(
    "[NH1:1].[c:2]1([O:6]CC)[c:3](OCC)[c:4](=O)[c:5]1=O"
    ">>"
    "[N:1][c:2]1[c:3](OCC)[c:4](=O)[c:5]1=O"
)


def _compute_product_smiles(amine_smiles: str, squarate_smiles: str) -> str:
    """
    Run a 2-component reaction SMARTS to produce the squaric acid monoamide.
    Returns canonical product SMILES, or '' on failure.
    """
    try:
        amine_mol = Chem.MolFromSmiles(amine_smiles)
        sq_mol    = Chem.MolFromSmiles(squarate_smiles)
        if amine_mol is None or sq_mol is None:
            return ""

        is_primary   = amine_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]"))
        is_secondary = amine_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH1]"))

        rxn = _RXN_PRIMARY if is_primary else (_RXN_SECONDARY if is_secondary else None)
        if rxn is None:
            return ""

        products = rxn.RunReactants((amine_mol, sq_mol))
        if not products:
            return ""

        prod_mol = products[0][0]
        Chem.SanitizeMol(prod_mol)
        return Chem.MolToSmiles(prod_mol)
    except Exception:
        return ""


# ── Soft-Rule Warnings ────────────────────────────────────────────────────────

def _soft_warnings(amine_type, position, catalyst):
    warnings = []
    if amine_type == "aniline" and not catalyst:
        warnings.append(
            "⚠️ Anilines typically require Zn(OTf)₂ catalyst "
            "(10–13 mol%) for good yield."
        )
    if position == "ortho":
        warnings.append(
            "⚠️ Ortho-substituted amines often give lower yields "
            "due to steric hindrance."
        )
    return " ".join(warnings) if warnings else None


# ── Recommendation Builder ────────────────────────────────────────────────────

def _build_recommendation(amine_type, position, temperature, catalyst, predicted_yield):
    tips = []
    if predicted_yield < 40:
        if amine_type == "aniline" and not catalyst:
            tips.append("Add Zn(OTf)₂ (10–13 mol%) to initiate reaction.")
        if temperature != "reflux":
            tips.append("Try heating to reflux to drive the reaction forward.")
        if position == "ortho":
            tips.append("Consider switching to para position for higher yield.")
    elif predicted_yield >= 80:
        tips.append("Conditions look optimal based on similar reactions in the training set.")
    return " ".join(tips) if tips else "Conditions look reasonable."


# ── Main Entry Point (SMILES-based) ───────────────────────────────────────────

def predict_from_smiles(
    amine_smiles:    str,
    squarate_smiles: str  = "CCOC1=C(OCC)C(=O)C1=O",
    temperature:     str  = "r.t.",
    catalyst:        bool = False,
) -> dict:
    """
    Predict reaction outcome given SMILES inputs.

    Args:
        amine_smiles    : SMILES of the amine coupling partner
        squarate_smiles : SMILES of the squarate ester (default = diethyl squarate)
        temperature     : 'r.t.' or 'reflux'
        catalyst        : True if Zn(OTf)2 is added

    Returns dict with keys:
        success, yield_percent, confidence, match_type, similarity,
        amine_matched, product_smiles, message, warning, recommendation,
        similar_reactions, error (optional)
    """

    # ── 1. Validate amine SMILES ─────────────────────────────────────────────
    amine_can = _canonical(amine_smiles)
    if amine_can is None:
        return {
            "success": False, "yield_percent": 0.0, "confidence": "none",
            "match_type": "invalid", "similarity": 0.0, "amine_matched": "—",
            "product_smiles": "", "message": "❌ Invalid amine SMILES.",
            "warning": None, "recommendation": "Check SMILES for typos.",
            "similar_reactions": [], "error": "invalid_smiles",
        }

    # ── 2. Hard-fail check (SMARTS) ──────────────────────────────────────────
    is_fail, fail_msg, fail_rec = _check_hard_fail(amine_can)
    if is_fail:
        return {
            "success": False, "yield_percent": 0.0, "confidence": "high",
            "match_type": "hard_fail", "similarity": 1.0, "amine_matched": amine_can,
            "product_smiles": "",
            "message": f"❌ {fail_msg}",
            "warning": None, "recommendation": fail_rec,
            "similar_reactions": [],
        }

    # ── 3. Match amine to training compound ──────────────────────────────────
    at, sub, pos, match_type, similarity, matched_label = _match_to_training(amine_can)

    # ── 4. k-NN yield prediction ─────────────────────────────────────────────
    query_feats = _encode_categorical(at, sub, pos, temperature, catalyst)
    predicted_yield, neighbours = _knn_predict(query_feats, k=3)

    success = predicted_yield >= 20.0

    # ── 5. Confidence (nearest-neighbour distance + match quality) ───────────
    nearest_dist = neighbours[0][0]
    if match_type == "exact" and nearest_dist < 0.5:
        confidence = "high"
    elif nearest_dist < 1.5 and match_type in ("exact", "similar"):
        confidence = "medium"
    else:
        confidence = "low (limited data for this combination)"

    # ── 6. Product SMILES ────────────────────────────────────────────────────
    product_smiles = ""
    if success:
        product_smiles = _compute_product_smiles(amine_can, squarate_smiles)

    # ── 7. Assemble response ─────────────────────────────────────────────────
    warning = _soft_warnings(at, pos, catalyst)
    recommendation = _build_recommendation(at, pos, temperature, catalyst, predicted_yield)

    icon = "✅" if success else "❌"
    message = (
        f"{icon} Predicted yield: {predicted_yield}% "
        f"(confidence: {confidence})"
    )

    similar = [
        {
            "amine":       r[0],
            "substituent": r[1],
            "position":    r[2],
            "temperature": "reflux" if r[3] else "r.t.",
            "catalyst":    bool(r[4]),
            "yield":       r[5],
        }
        for _, _, r in neighbours
    ]

    return {
        "success":          success,
        "yield_percent":    predicted_yield,
        "confidence":       confidence,
        "match_type":       match_type,       # 'exact' | 'similar' | 'inferred' | 'hard_fail'
        "similarity":       similarity,
        "amine_matched":    matched_label,
        "product_smiles":   product_smiles,
        "message":          message,
        "warning":          warning,
        "recommendation":   recommendation,
        "similar_reactions": similar,
    }


# ── Legacy Entry Point (dropdown-encoded, kept for mobile app compat) ─────────

_LEGACY_AMINE_MAP  = {"aniline": 0, "benzylamine": 1}
_LEGACY_POS_MAP    = {"ortho": 0, "meta": 1, "para": 2, "none": 2}
_LEGACY_SUB_MAP    = {
    "none": 0, "F": 1, "Cl": 2, "Br": 3, "OH": 4,
    "OMe": 5, "CF3": 6, "tolyl": 7, "fused_ring": 8,
    "pyrrolidine": 9, "azetidine": 10, "piperidine": 11, "azepane": 12,
    "thiazolidine": 13, "morpholine": 14, "dimethylmorpholine": 15,
    "thiomorpholine": 16, "thiomorpholine_dioxide": 17,
    "boc_piperazine": 18, "boc_dimethylpiperazine": 19,
    "bromo_tetrahydroquinoline": 20, "piperazine": 21,
}

_LEGACY_HARD_FAILS = {
    ("aniline",      "Cl",         "ortho"): (
        "Ortho-Cl aniline fails due to steric crowding. No reaction even at reflux.",
        "Change position to meta or para.",
    ),
    ("aniline",      "Br",         "ortho"): (
        "Ortho-Br aniline fails due to steric crowding. No reaction even at reflux.",
        "Change position to meta or para.",
    ),
    ("heterocyclic", "piperazine", "none"): (
        "Unprotected piperazine fails — both N atoms compete.",
        "Use Boc-protected piperazine instead.",
    ),
}


def predict(amine_type: str, substituent: str, position: str,
            temperature: str, catalyst: bool) -> dict:
    """
    Legacy predict() — accepts dropdown-encoded inputs.
    Kept for backward compatibility with the React Native mobile app.
    """
    key = (amine_type, substituent, position)
    if key in _LEGACY_HARD_FAILS:
        msg, rec = _LEGACY_HARD_FAILS[key]
        return {
            "success": False, "yield_percent": 0.0, "confidence": "high",
            "message": f"❌ {msg}", "warning": None,
            "recommendation": rec, "similar_reactions": [],
        }

    warning = None
    if amine_type == "aniline" and not catalyst:
        warning = "⚠️ Anilines typically require Zn(OTf)₂ catalyst."

    query = [
        _LEGACY_AMINE_MAP.get(amine_type, 0),
        _LEGACY_SUB_MAP.get(substituent, 0),
        _LEGACY_POS_MAP.get(position, 2),
        int(temperature == "reflux"),
        int(bool(catalyst)),
    ]
    predicted_yield, neighbours = _knn_predict(query, k=3)
    success = predicted_yield >= 20.0

    nearest_dist = neighbours[0][0]
    if nearest_dist < 0.5:
        confidence = "high"
    elif nearest_dist < 1.5:
        confidence = "medium"
    else:
        confidence = "low (limited data for this combination)"

    icon = "✅" if success else "❌"
    similar = [
        {"amine": r[0], "substituent": r[1], "position": r[2], "yield": r[5]}
        for _, _, r in neighbours
    ]

    return {
        "success":          success,
        "yield_percent":    predicted_yield,
        "confidence":       confidence,
        "message":          f"{icon} Predicted yield: {predicted_yield}% (confidence: {confidence})",
        "warning":          warning,
        "recommendation":   _build_recommendation(amine_type, position, temperature, catalyst, predicted_yield),
        "similar_reactions": similar,
    }
