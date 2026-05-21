"""
extract_library.py
Reads ChemInventory-Export.xlsx, deduplicates by CAS, prints JS-ready data
for the compounds relevant to the Chemdem predict.html interface.
"""
import pandas as pd
import json, re

LIBRARY_PATH = r'C:\Users\mahan\Downloads\ChemInventory-Export.xlsx'

df = pd.read_excel(LIBRARY_PATH)
df.columns = [c.strip() for c in df.columns]
print("Columns:", list(df.columns))
print(f"Total rows: {len(df)}")

# Normalise column names
col_map = {}
for c in df.columns:
    lc = c.lower().replace(' ', '_').replace('-', '_')
    col_map[c] = lc
df.rename(columns=col_map, inplace=True)

print("\nNormalised columns:", list(df.columns))

# Print first few rows
print("\nFirst 5 rows:")
print(df.head())

# Find SMILES and formula columns
smiles_col = [c for c in df.columns if 'smiles' in c.lower()]
formula_col = [c for c in df.columns if 'formula' in c.lower()]
cas_col = [c for c in df.columns if 'cas' in c.lower()]
name_col = [c for c in df.columns if 'name' in c.lower() or 'substance' in c.lower()]
mw_col = [c for c in df.columns if 'weight' in c.lower() or 'mw' in c.lower()]

print(f"\nSMILES col: {smiles_col}")
print(f"Formula col: {formula_col}")
print(f"CAS col: {cas_col}")
print(f"Name col: {name_col}")
print(f"MW col: {mw_col}")

# Use first match for each
SC = smiles_col[0] if smiles_col else None
FC = formula_col[0] if formula_col else None
CASC = cas_col[0] if cas_col else None
NC = name_col[0] if name_col else None

print(f"\nUsing: name={NC}, cas={CASC}, smiles={SC}, formula={FC}")

# Build deduplicated compound dict keyed by CAS
compounds = {}
for _, row in df.iterrows():
    cas  = str(row[CASC]).strip() if CASC else ''
    name = str(row[NC]).strip()   if NC   else ''
    smi  = str(row[SC]).strip()   if SC   else ''
    frm  = str(row[FC]).strip()   if FC   else ''

    # Skip if no CAS
    if not cas or cas in ('nan', '', '-'):
        continue
    # Skip unknown SMILES
    if smi.lower() in ('nan', '', 'unknown', '-', 'none'):
        smi = ''

    # Deduplicate: keep first occurrence (or prefer the one with SMILES)
    if cas not in compounds:
        compounds[cas] = {'name': name, 'cas': cas, 'smiles': smi, 'formula': frm}
    else:
        # Update only if current has no SMILES but new one does
        if not compounds[cas]['smiles'] and smi:
            compounds[cas]['smiles'] = smi
            compounds[cas]['formula'] = frm

print(f"\nUnique CAS entries: {len(compounds)}")

# Print all compounds (for review)
print("\n=== ALL LIBRARY COMPOUNDS ===")
for cas, c in sorted(compounds.items()):
    smiles_preview = c['smiles'][:40] + '...' if len(c['smiles']) > 40 else c['smiles']
    print(f"  CAS {cas:15s} | {c['formula']:12s} | {c['name'][:40]:40s} | {smiles_preview}")

# Output as JSON for integration
out = list(compounds.values())
with open(r'C:\Users\mahan\Downloads\Chemdem\data\library_compounds.json', 'w') as f:
    json.dump(out, f, indent=2)
print(f"\n✓ Saved {len(out)} compounds to library_compounds.json")
