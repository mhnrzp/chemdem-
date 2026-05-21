# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

Chemdem predicts the outcome (yield, product SMILES, side products, confidence) of **squaric acid monoamide synthesis** — specifically `Amine + Diethyl Squarate (DES) → Monosquarate-amide + EtOH`. All training data is from two papers by the Wren group (SynOpen 2023, Bioorg Med Chem 2024), totalling 42 reactions.

The project has three layers:
- **`backend/`** — FastAPI prediction server (Python, RDKit)
- **`landing/predict.html`** — Single-file React 18 web UI (no build step, CDN deps)
- **`app/`** — Legacy React Native / Expo mobile skeleton (uses the old `/predict` endpoint only)

---

## Running the Backend

```bash
cd backend
pip install -r requirements.txt      # fastapi uvicorn[standard] pydantic python-dotenv rdkit httpx
python -m uvicorn main:app --reload --port 8000
```

Health check: `GET http://localhost:8000/health`

The mobile app's `api.ts` points to `http://localhost:8000` — change to the machine's LAN IP when testing on a real device.

---

## Backend Architecture

### Endpoint hierarchy (`main.py`)

| Endpoint | Version | Description |
|---|---|---|
| `POST /predict/full` | v3 (current) | Full pipeline: validation + k-NN + references + side products + confidence |
| `POST /predict/smiles` | v2 | k-NN only, fast, no reference search |
| `POST /predict` | v1 legacy | Dropdown-encoded inputs for the React Native app |
| `POST /lookup` | v3 | Resolve any identifier to a validated structure |
| `POST /validate` | utility | SMILES → canonical form + formula + MW |

### Prediction pipeline (`model.py`)

1. **Hard-fail SMARTS check** — ortho-Cl aniline, ortho-Br aniline, unprotected piperazine → immediate fail, no k-NN run
2. **Morgan FP matching** — `GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)` + Tanimoto against `COMPOUND_DB` (38 entries). Threshold 0.30 → `similar`; exact canonical SMILES → `exact`; below threshold → SMARTS class inference → `inferred`
3. **k-NN yield** — inverse-distance weighted, k=3, 5-feature encoding `[amine_type, substituent, position, temperature, catalyst]` over `DATASET` (42 rows)
4. **Product SMILES** — reaction SMARTS. **Critical**: RDKit perceives the squarate ring as aromatic (bond type 1.5), so ring atoms must be `[c:n]` not `[C:n]` in SMARTS. `_RXN_PRIMARY` handles `[NH2]`, `_RXN_SECONDARY` handles `[NH1]`.

### Reference search (`reference_search.py`)

Four tiers, always all run and clearly labelled:
- **Tier 1** — Local dataset (Morgan FP similarity, same logic as model.py)
- **Tier 2** — PubChem (compound-level cross-reference only, no reaction search API exists)
- **Tier 3** — ORD stub (`status: "not_connected"`)
- **Tier 4** — USPTO stub (`status: "not_connected"`)

### Compound lookup (`lookup.py`)

Resolves any identifier (SMILES, CAS, IUPAC name, common name) via:
1. RDKit (SMILES input direct)
2. OPSIN at `https://www.ebi.ac.uk/opsin/ws/{name}` with `Accept: chemical/x-daylight-smiles` — **must use `follow_redirects=True`** in httpx (old `opsin.ch.cam.ac.uk` URL permanently redirected)
3. PubChem REST: `ConnectivitySMILES` key (not `CanonicalSMILES`)

Cross-checks OPSIN vs PubChem SMILES and emits `warnings[]` on mismatch.

### Side products (`side_products.py`)

Rule-based, always returns 5 entries: bis-squaramide, squaric acid, recovered DES, N-oxide (morpholine only), monoethyl ester. Each has `probability` (high/medium/low), `risk_factors[]`, and `mitigation[]`. The EWG SMARTS pattern `Nc1ccc([F,Cl,Br,C(F)(F)F])cc1` emits a non-fatal parse warning at startup — known issue.

---

## Frontend Architecture (`landing/predict.html`)

Single HTML file, no build step. React 18 + Babel standalone + SmilesDrawer 2.0.1, all from CDN.

**State lives in `App`:**
- `sm1` / `sm2` — compound panel state `{ iupac, cas, compound }` where `compound` is a full compound object
- `result` — raw API response (either `/predict/full` shape or offline fallback shape)
- `tab` — `'input'` | `'prediction'`

**Key normalisation pattern** — `/predict/full` nests prediction under `result.prediction`; the offline k-NN fallback puts it at the top level. Every component that reads yield/status/product_smiles must handle both shapes:
```js
const pred = result.prediction || result;
const productSmiles = result?.prediction?.product_smiles || result?.product_smiles || '';
```

The `_reactionResult` block in `App` derives a single `_isNoReaction` / `_productLabel` used by the product box. **Never read `result.success` directly** — it is `undefined` for API results (the field is at `result.prediction.success`).

**`AutocompleteInput`** does a two-tier search:
1. Instant local search over the 38 hard-coded `COMPOUNDS` array
2. Debounced (600 ms) `POST /lookup` call for any input ≥ 3 chars — returns a compound with `fromPubchem: true`, no `enc` field

External (PubChem-resolved) compounds have no `enc` vector, so the offline k-NN fallback (`runPrediction`) returns `_no_enc: true` for them — a dedicated error UI handles this case.

**`ResultPanel`** derives status from `summ?.status` → `pred.yield_percent` thresholds → `pred.success` in that priority order. It never assumes a field exists.

---

## Dataset & Chemistry Facts

- Default squarate: `CCOC1=C(OCC)C(=O)C1=O` (diethyl squarate, DES), CAS 5231-87-8
- All training reactions used EtOH solvent at 1:1 DES:amine ratio
- Catalyst is Zn(OTf)₂ at 10–13 mol%; required for anilines, not for benzylamines
- Hard fails (experimentally confirmed): ortho-Cl aniline, ortho-Br aniline, unprotected piperazine
- Best yields: para-substituted amines > meta > ortho; benzylamines > anilines > heterocyclics
- Most active HDAC8 compound: para-F aniline monoamide (IC₅₀ = 1.20 μM)

---

## Source Labels in Responses

| `paper_key` | Citation |
|---|---|
| `paper_1` | Long et al., SynOpen 2023 (aniline/benzylamine series, 29 reactions) |
| `paper_2` | Long et al., Bioorg. Med. Chem. 2024 (heterocyclic series, 13 reactions) |

The DOI for Paper 1 in `reference_search.py` is currently a placeholder (`10.1055/a-XXXX-XXXX`) — needs the real DOI.
