# side_products.py
# Rule-based side-product prediction for squaric acid monoamide synthesis.
#
# Chemistry background:
#   DES (diethyl squarate) + R-NH2  →  monoamide + EtOH   (target)
#   Competing pathways:
#     1. Bis-squaramide  (second substitution with excess / reactive amine)
#     2. Squaric acid    (hydrolysis by trace water)
#     3. Squaric acid monoamide oligomers  (very rare, concentrated solutions)
#
# Each entry returns:
#   name           str
#   smiles         str | None
#   probability    'high' | 'medium' | 'low'
#   source         'known chemistry' | 'predicted'
#   explanation    str
#   risk_factors   list[str]   — what INCREASES risk
#   mitigation     list[str]   — how to AVOID it

from rdkit import Chem
from rdkit.Chem import AllChem

# ── SMARTS patterns ───────────────────────────────────────────────────────────
_PAT_NH2      = Chem.MolFromSmarts("[NH2]")
_PAT_NH1      = Chem.MolFromSmarts("[NH1]")
_PAT_ANILINE  = Chem.MolFromSmarts("Nc1ccccc1")
_PAT_BENZYLAM = Chem.MolFromSmarts("NCc1ccccc1")
_PAT_PIPERAZ  = Chem.MolFromSmarts("[NH]1CC[NH]CC1")      # unprotected piperazine
_PAT_EDG      = Chem.MolFromSmarts("Nc1ccc([OH,OC,C])cc1")   # para ED group on aniline
# EWG patterns split by atom type — comma-OR doesn't parse in aromatic SMARTS context
_PAT_EWG_HAL  = Chem.MolFromSmarts("Nc1ccc([F,Cl,Br])cc1")
_PAT_EWG_CF3  = Chem.MolFromSmarts("Nc1ccccc1C(F)(F)F")

# ── Reaction SMARTS for bis-squaramide ───────────────────────────────────────
# Replace the remaining OEt on the monoamide with the amine's N fragment.
# Uses a single generic SMARTS that matches [NH2] or [NH1] amines.
_RXN_BIS = AllChem.ReactionFromSmarts(
    "[N:1].[c:2]1([O:6]CC)[c:3]([N:8])[c:4](=O)[c:5]1=O"
    ">>"
    "[N:1][c:2]1[c:3]([N:8])[c:4](=O)[c:5]1=O"
)


def _compute_bis_smiles(amine_smiles: str, monoamide_smiles: str) -> str:
    """
    Compute the bis-squaramide SMILES by substituting the remaining OEt
    on the monoamide with a second equivalent of the amine.
    Returns canonical SMILES or '' on failure.
    """
    try:
        amine_mol = Chem.MolFromSmiles(amine_smiles)
        mono_mol  = Chem.MolFromSmiles(monoamide_smiles)
        if amine_mol is None or mono_mol is None:
            return ""
        prods = _RXN_BIS.RunReactants((amine_mol, mono_mol))
        if prods:
            p = prods[0][0]
            Chem.SanitizeMol(p)
            return Chem.MolToSmiles(p)
    except Exception:
        pass
    return ""


def _assess_bis_risk(amine_mol, amine_type: str, temperature: str, catalyst: bool) -> str:
    """Return 'high' | 'medium' | 'low' risk for bis-squaramide formation."""
    # Unprotected piperazine always gives bis-squaramide → hard fail catches this first,
    # but include here for completeness
    if amine_mol.HasSubstructMatch(_PAT_PIPERAZ):
        return "high"

    # Cyclic secondary amines (morpholine, piperidine, etc.)
    # — only one N-H, can't undergo second substitution on same N
    # — second substitution would need another molecule's N, unlikely in dilute conditions
    is_nh1 = amine_mol.HasSubstructMatch(_PAT_NH1) and not amine_mol.HasSubstructMatch(_PAT_NH2)
    if is_nh1:
        return "low" if temperature != "reflux" else "medium"

    # Primary amines at reflux with electron-donating groups: higher risk
    if amine_mol.HasSubstructMatch(_PAT_NH2):
        if amine_mol.HasSubstructMatch(_PAT_EDG) and temperature == "reflux":
            return "medium"
        if amine_type == "benzylamine":
            return "medium"
        if temperature == "reflux":
            return "medium"
        return "low"

    return "low"


def predict_side_products(
    amine_smiles:    str,
    amine_type:      str,   # 'aniline' | 'benzylamine' | 'heterocyclic'
    temperature:     str,   # 'r.t.' | 'reflux'
    catalyst:        bool,
    solvent:         str,   # 'EtOH' | 'MeOH' | 'DMF' | 'DMSO' | 'THF'
    monoamide_smiles: str = "",
) -> list[dict]:
    """
    Predict likely side products for squaric acid monoamide synthesis.

    Returns a list of side-product dicts, each with:
      name, smiles, probability, source, explanation, risk_factors, mitigation
    """
    amine_mol = Chem.MolFromSmiles(amine_smiles)
    if amine_mol is None:
        return []

    side_products = []

    # ── 1. Bis-squaramide (double substitution) ───────────────────────────────
    bis_risk = _assess_bis_risk(amine_mol, amine_type, temperature, catalyst)
    bis_smiles = ""
    if monoamide_smiles:
        bis_smiles = _compute_bis_smiles(amine_smiles, monoamide_smiles)

    risk_factors_bis = []
    if temperature == "reflux":
        risk_factors_bis.append("Reflux conditions accelerate the second substitution")
    if amine_mol.HasSubstructMatch(_PAT_EDG):
        risk_factors_bis.append("Electron-donating substituent increases amine nucleophilicity")
    if amine_type == "benzylamine":
        risk_factors_bis.append("Benzylamines are more nucleophilic than anilines")
    if not catalyst and amine_type == "aniline":
        risk_factors_bis.append("Absence of Zn(OTf)₂ may slow first substitution, allowing competition")

    side_products.append({
        "name":        "Bis-squaramide (double substitution product)",
        "smiles":      bis_smiles,
        "probability": bis_risk,
        "source":      "known chemistry",
        "explanation": (
            "The squaric acid monoamide still contains one ethyl ester group that "
            "can undergo a second substitution with another equivalent of the amine. "
            "This is the most common side reaction in squaric acid monoamide synthesis. "
            "Using exactly 1.0 equiv of amine and monitoring reaction progress by TLC "
            "are key to minimising this product."
        ),
        "risk_factors": risk_factors_bis or ["No elevated risk factors identified under these conditions"],
        "mitigation": [
            "Use exactly 1.0 equiv of amine (avoid excess)",
            "Monitor by TLC and stop reaction when monoamide spot appears",
            "Lower temperature to room temperature if possible",
            "Use slow dropwise addition of amine to squarate",
        ],
    })

    # ── 2. Squaric acid (hydrolysis) ──────────────────────────────────────────
    protic_solvents = {"MeOH", "EtOH"}
    hydrolysis_risk = "medium" if solvent in protic_solvents else "low"

    side_products.append({
        "name":        "Squaric acid (hydrolysis product)",
        "smiles":      "OC1=C(O)C(=O)C1=O",
        "probability": hydrolysis_risk,
        "source":      "known chemistry",
        "explanation": (
            "Diethyl squarate is susceptible to hydrolysis in the presence of water or "
            "protic solvents. Squaric acid forms when both ethyl ester groups are "
            "replaced by hydroxyl groups. This is particularly relevant in MeOH/EtOH "
            "solvents that may carry trace moisture."
        ),
        "risk_factors": (
            [f"{solvent} is a protic solvent — trace water accelerates hydrolysis"]
            if solvent in protic_solvents
            else ["Risk is low under standard anhydrous conditions"]
        ),
        "mitigation": [
            "Use anhydrous (dry) solvent",
            "Perform reaction under inert atmosphere (N₂ or Ar)",
            "Avoid prolonged reaction times after product formation",
            "Use fresh, dry DES starting material",
        ],
    })

    # ── 3. Recovered / partially reacted DES ─────────────────────────────────
    # This is a risk when amine nucleophilicity is low or catalyst is missing
    sm_risk = "low"
    sm_factors = []
    if amine_type == "aniline" and not catalyst:
        sm_risk = "medium"
        sm_factors.append("Anilines require Zn(OTf)₂ to activate the squarate electrophile")
    if temperature == "r.t." and amine_type == "aniline":
        sm_factors.append("Room temperature is borderline for anilines without catalyst")

    side_products.append({
        "name":        "Unreacted diethyl squarate (recovered starting material)",
        "smiles":      "CCOC1=C(OCC)C(=O)C1=O",
        "probability": sm_risk,
        "source":      "known chemistry",
        "explanation": (
            "If the amine is insufficiently nucleophilic or the reaction conditions are "
            "too mild, the squarate starting material may be recovered unreacted. "
            "This is particularly relevant for electron-poor anilines without a Lewis "
            "acid catalyst. This is not a true side product but indicates incomplete conversion."
        ),
        "risk_factors": sm_factors or ["No elevated risk factors under these conditions"],
        "mitigation": [
            "Add Zn(OTf)₂ (10–13 mol%) for anilines",
            "Extend reaction time and check conversion by TLC",
            "Consider heating to reflux if room temperature gives incomplete conversion",
        ],
    })

    # ── 4. N-oxide or oxidation artefact (rare) ───────────────────────────────
    # Only relevant for morpholine-type amines where N can be oxidised
    has_o_adjacent_n = Chem.MolFromSmarts("[N]1CCOCC1")  # morpholine-type
    if amine_mol.HasSubstructMatch(has_o_adjacent_n):
        side_products.append({
            "name":        "N-oxide (morpholine N-oxide if oxidation occurs)",
            "smiles":      "O=N1CCOCC1",
            "probability": "low",
            "source":      "predicted",
            "explanation": (
                "Morpholine and similar oxygen-containing cyclic amines can form "
                "N-oxides under oxidative workup conditions. This is uncommon under "
                "standard reaction conditions but may appear as a minor impurity."
            ),
            "risk_factors": [
                "Morpholine-type amines are susceptible to N-oxidation",
                "Oxidative workup or air exposure",
            ],
            "mitigation": [
                "Perform workup under inert atmosphere",
                "Avoid oxidising reagents or conditions during workup",
            ],
        })

    # ── 5. Squaric acid monoethyl ester (partially hydrolysed intermediate) ──
    side_products.append({
        "name":        "Squaric acid monoethyl ester (reaction intermediate)",
        "smiles":      "CCOC1=C(O)C(=O)C1=O",
        "probability": "low",
        "source":      "predicted",
        "explanation": (
            "Partial hydrolysis of one ester group can give the squaric acid "
            "monoethyl ester, which can then co-elute with the product. "
            "This intermediate is typically only observed under aqueous conditions "
            "or prolonged reaction times."
        ),
        "risk_factors": [
            "Aqueous workup",
            "Trace moisture in protic solvents",
        ],
        "mitigation": [
            "Use anhydrous conditions",
            "Minimise aqueous workup steps",
        ],
    })

    return side_products
