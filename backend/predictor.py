# predictor.py
# Public interface wrappers — keeps main.py clean.

from model import predict as _predict_legacy, predict_from_smiles as _predict_smiles


def predict_reaction(amine_type: str, substituent: str, position: str,
                     temperature: str, catalyst: bool) -> dict:
    """Legacy wrapper — dropdown-encoded inputs."""
    return _predict_legacy(
        amine_type=amine_type,
        substituent=substituent,
        position=position,
        temperature=temperature,
        catalyst=catalyst,
    )


def predict_reaction_smiles(amine_smiles: str, squarate_smiles: str,
                            temperature: str, catalyst: bool) -> dict:
    """v2 wrapper — SMILES inputs."""
    return _predict_smiles(
        amine_smiles=amine_smiles,
        squarate_smiles=squarate_smiles,
        temperature=temperature,
        catalyst=catalyst,
    )
