# Dataset Notes — Table 1 (Wren Group, SynOpen 2023)

## Summary
- 27 successful reactions (compounds 16–43)
- 2 confirmed FAILED reactions (ortho-Cl and ortho-Br anilines)
- All reactions: Amine + Diethyl Squarate (DES) → Monosquarate-amide
- Solvent: EtOH (except compound 16: MeOH)

## Key Patterns Learned From Data

### Substituent Position Effect on Yield
| Position | Avg Yield | Notes |
|---|---|---|
| Para | ~70% | Best position consistently |
| Meta | ~48% | Moderate |
| Ortho | ~18–38% | Low — steric clash (Cl/Br = FAIL) |

### Amine Type Effect
| Amine | Catalyst needed? | Avg Yield |
|---|---|---|
| Benzylamine | No | ~57% (up to 99% at reflux) |
| Aniline | Usually yes | ~47% |

### Catalyst Effect (Zn(OTf)₂)
- Anilines without catalyst: often no reaction
- With 10–20 mol%: reaction initiates
- Benzylamines: react without catalyst (more nucleophilic)

### Temperature Effect
- r.t. is standard
- Reflux can boost yield dramatically (compound 37: 99%)
- But reflux not always transferable between substrates

### Complete FAIL Conditions
- Ortho-Cl aniline → NO REACTION even at reflux + 48h + catalyst
- Ortho-Br aniline → NO REACTION even at reflux + 48h + catalyst
- Reason: poor nucleophilicity + steric crowding combined

## Paper 3 — Biological Activity Data (HDAC8 IC₅₀, ChemMedChem 2026)
Active compounds (single-digit micromolar):
- 9b (para-F aniline monosquarate) → HDAC8 IC₅₀ = 1.20 μM ⭐ BEST
- 14b (N-methyl bis-squaramide benzylamine) → 1.3 μM
- 13b (bis-squaramide benzylamine) → 1.8 μM
- 13c (bis-squaramide para-Br aniline) → 2.0 μM
- All aliphatic/heterocyclic compounds → >35 μM (inactive)

Key SAR: para > meta > ortho for both yield AND HDAC8 activity
Benzylamine > aniline in both synthesis and biology

## Missing Data (Need From Lecturer)
- Exact substituent identities for compounds 30, 31, 32 (tolyl position unclear)
- Any unpublished reaction attempts
- Reactions with different DES:amine ratios (all 1:1 in this paper)
- Squaric acid (not DES) as starting material results
