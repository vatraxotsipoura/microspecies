# pKa — Microspecies Ionization Curve Calculator

Generate Chemicalize-style microspecies distribution plots from a SMILES string and user-supplied pKa values.

**Input:** SMILES + pKa values (+ optional logP) → **Output:** interactive HTML with chart, molecule structures, annotated pKa diagram, and (with logP) a logD vs pH chart.

---

## Quick Start

```bash
conda activate pka

# Pseudoephedrine (2 ionizable sites)
python ionization.py "C[C@@H](NC)[C@H](O)c1ccccc1" \
    --pka "9.52:base:NH,13.89:acid:OH" \
    --name "Pseudoephedrine" \
    -o pseudoephedrine.html

# Salbutamol (4 ionizable sites)
python ionization.py "CC(C)(C)NCC(C1=CC(=C(C=C1)O)CO)O" \
    --pka "9.5:base:NH,9.8:acid:phenol-OH,12.4:acid:chain-OH,13.1:acid:benzyl-OH" \
    --name "Salbutamol" \
    -o salbutamol.html

# Butorphanol with logD vs pH (needs logP of the neutral species)
python ionization.py "C1CC[C@]2([C@H]3CC4=C([C@]2(C1)CCN3CC5CCC5)C=C(C=C4)O)O" \
    --pka "8.19:base:N:12,10.0:acid:phenol-OH:22" \
    --logp 3.77 --name "Butorphanol" -o butorphanol.html

# Glycine (amino acid, zwitterion)
python ionization.py "NCC(=O)O" \
    --pka "2.34:acid:COOH,9.60:base:NH2" \
    --name "Glycine"
```

---

## Installation

Requires a conda environment with Python 3.11:

```bash
conda create -n pka python=3.11 -y
conda activate pka
pip install rdkit numpy plotly Pillow
```

No GPU, no PyTorch, no heavy ML frameworks needed.

---

## Usage

```
python ionization.py "SMILES" --pka "PKA_SPEC" [--name NAME] [-o FILE] [--logp LOGP] [--logp-type neutral|zwitterion]
                    [--cation-offset X] [--anion-offset Y] [--zwitterion-offset Z] [--ref POINTS]
```

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `SMILES` | Yes | SMILES string of the molecule (positional) |
| `--pka`  | Yes | Comma-separated pKa specifications (see format below) |
| `--name` | No | Molecule name (used in the HTML title) |
| `-o` / `--output` | No | Output HTML file path (default: auto-generated temp file) |
| `--logp` | No | logP of the neutral species — adds a logD vs pH chart |
| `--logp-type` | No | `neutral` (default) or `zwitterion` — what the `--logp` value refers to |
| `--cation-offset` | No | logP penalty per protonated base group (default `3.0`) |
| `--anion-offset` | No | logP penalty per deprotonated acid group (default `4.0`) |
| `--zwitterion-offset` | No | logP penalty per (+, −) pair on one microspecies (default `3.0`) |
| `--ref` | No | Literature logD points `pH:logD,...` to compare against (needs `--logp`) — prints residuals + RMSE and overlays the points on the logD chart |

### pKa Specification Format

Each ionizable site is specified as:

```
pKa:type[:label[:atom_idx]]
```

| Field | Required | Values | Description |
|-------|----------|--------|-------------|
| `pKa` | Yes | float | The pKa value |
| `type` | Yes | `acid` or `base` | Whether the site is acidic (HA → A⁻ + H⁺) or basic (BH⁺ → B + H⁺) |
| `label` | No | string | Human-readable label (e.g. "NH", "OH", "COOH") |
| `atom_idx` | No | int | RDKit atom index — override for auto-detection |

Multiple sites are separated by commas:

```
"9.52:base:NH,13.89:acid:OH"
"2.34:acid:COOH:4,9.60:base:NH2:0"    # glycine NCC(=O)O with explicit atom indices
```

### Atom Index Assignment

The program needs to know which atom in the molecule each pKa value belongs to, in order to render correct protonation states.

**Option A — Auto-detection (default):**
The program finds candidate atoms and ranks them by group type:
- **Base sites**: aliphatic amines / amidines / imines, then aniline-type N, then
  pyridine-type aromatic N. Amide, sulfonamide, nitro, nitrile, N-oxide and
  quaternary N are never used.
- **Acid sites**: carboxylic/sulfonic/phosphonic O-H, then phenols, then thiols
  and sulfonamide/imide N-H, then alcohols.

Acid sites are matched in order of increasing pKa to the most acidic candidates,
base sites in order of decreasing pKa to the most basic ones (ties follow SMILES
order). Atoms given explicitly are never reused. This works well for simple
molecules; for two groups of the same kind (e.g. two COOH) check the result or
give indices.

All sites are checked before anything is computed: indices must exist and be
distinct, acid sites must carry an O-H/S-H/N-H, and base sites must be neutral
nitrogens. Problems stop the run with an error message.

**Option B — Manual override (fallback):**
If auto-detection assigns the wrong atom (e.g., multiple similar groups), specify the RDKit atom index explicitly as the 4th field:

```
"9.8:acid:phenol-OH:12,13.1:acid:chain-OH:5"
```

To find atom indices, you can inspect the molecule in RDKit:
```python
from rdkit import Chem
mol = Chem.MolFromSmiles("your_smiles")
for atom in mol.GetAtoms():
    print(atom.GetIdx(), atom.GetSymbol())
```

---

## Output

The program generates a **standalone HTML file** (no server needed) containing:

1. **Interactive Plotly chart** — pH (0–14) vs microspecies distribution (%).
   Hover over any point to see exact pH and percentage. Each curve is color-coded.

2. **Molecule structures** — one 2D structure per significant species (those
   reaching ≥1% at any pH), rendered by RDKit. Each structure shows the correct
   protonation state and has a colored border matching its curve.

3. **pKa table** — summary of all ionizable sites with type and value.

4. **Annotated structure** — the base molecule with pKa values displayed next
   to each ionizable group.

5. **logD vs pH chart** (only with `--logp`) — the logD curve and logP as a
   reference line (logD₇.₄ is given above the chart). With `--ref`, literature
   points are overlaid.

---

## How It Works

### Henderson-Hasselbalch Microspecies Model

For a molecule with N ionizable sites, there are 2^N possible protonation
states (microspecies). At each pH, the fractional population of each
microspecies is computed as:

For each site *i* with pKa *Kᵢ*:

$$
\text{fraction protonated} = \frac{1}{1 + 10^{(\text{pH} - K_i)}}
$$

$$
\text{fraction deprotonated} = \frac{10^{(\text{pH} - K_i)}}{1 + 10^{(\text{pH} - K_i)}}
$$

This formula is identical for both acid and base sites — the pKa always refers
to the conjugate acid form (BH⁺ → B + H⁺ or HA → A⁻ + H⁺). The difference
is only in labeling: for an acid, the protonated form is neutral; for a base,
the protonated form is charged.

Each microspecies fraction is the product of the per-site fractions (assuming
independent equilibria), and all 2^N fractions sum to 100% at every pH.

### logD (lipophilicity vs pH)

With `--logp`, logD is the fraction-weighted partition coefficient over all
microspecies:

$$
\log D(\text{pH}) = \log_{10}\Big(\sum_i f_i(\text{pH}) \cdot P_i\Big)
$$

- The fully uncharged microspecies has logPᵢ = logP.
- Each ionized group lowers logPᵢ, modelling ion-pair partitioning
  (Avdeef's "diff 3–4" rule, octanol / 0.15 M KCl): a protonated base by
  `--cation-offset` (default 3.0), a deprotonated acid by `--anion-offset`
  (default 4.0). Each (+, −) pair on the same microspecies counts once as
  `--zwitterion-offset` (default 3.0) instead. This makes logD level off at
  extreme pH, as measured curves do.
- The defaults were fitted to measured profiles (see the validation below);
  if a measured ion-pair logP is known, pass logP − logP(ion) as the offset.
- `--logp-type zwitterion` says the `--logp` value is the zwitterion's, for
  compounds that exist almost only as zwitterions (cetirizine, amino acids).
  The tool warns when the zwitterion outnumbers the uncharged form >10× and
  `--logp-type` was left at `neutral`.

Ampholytes (an acid and a base site, e.g. phenol + amine):

- The model treats sites as independent. For close pKa values, prefer the
  **microscopic** pKa values that start from the uncharged form (uncharged ⇌
  cation for the base, uncharged ⇌ anion for the acid) over the macroscopic
  ones. For the morphinans below this cut the logD₇.₄ error from up to 0.6 to
  ≤ 0.05.
- With macroscopic values the error grows with the zwitterion/uncharged ratio
  (≈ 0.1 for morphine, 0.6 for naltrexone and oxymorphone).

Limitations:

- Permanently charged groups (quaternary ammonium) have no pKa and can't be
  modelled.
- Penalties add up per charged group, so species with many charges reach
  unrealistically low logD (e.g. −10 and below); measured logD rarely goes
  below about −3.

Accuracy notes:

- Errors in logP shift the whole curve; errors in pKa shift it sideways, which
  becomes a 1:1 logD error where the molecule is mostly ionized.
- The plateau height far from the pKa depends on counterion and ionic strength
  (±1 log unit); near pH 7.4 and the pKa values the offset barely matters.
- Literature "logP" values are often actually logD₇.₄ or predicted clogP —
  check before using them.

### Validating logD against literature

Pass measured logD values with `--ref`; the tool prints calc − ref at each pH
and the RMSE, and plots the points on the logD chart.

Reference compounds used to choose the default offsets:

| Compound | Type | logP | pKa | Reference data | Best-fit offset |
|---|---|---|---|---|---|
| Propranolol | base | 3.48 | 9.53 | Potentiometric profile, ion-pair logP 0.78 (Avdeef/Scherrer, 0.15 M KCl) | 2.70 (measured) |
| Lidocaine | base | 2.45 | 7.92 | pH-metric logD points (digitized from ChemAxon logP/logD docs) | 3.0 |
| Ibuprofen | acid | 3.97 | 4.42 | Observed logD points (digitized from ChemAxon logP/logD docs) | 4.1 |

RMSE (calc − ref) for a single offset vs. the type-specific defaults:

| Compound | 3.5 for all ions (old) | cation 3.0 / anion 4.0 (default) |
|---|---|---|
| Propranolol | 0.37 | 0.14 |
| Lidocaine | 0.21 | 0.12 |
| Ibuprofen | 0.23 | 0.11 |

```bash
python ionization.py "CC(C)NCC(COc1cccc2ccccc12)O" --pka "9.53:base:NH" --logp 3.48 \
    --ref "2:0.78,4:0.78,6:0.84,7:1.17,7.4:1.45,8:1.97,8.5:2.42,9:2.84,9.53:3.18,10:3.35,11:3.47,12:3.48" \
    --name "Propranolol" -o propranolol.html

python ionization.py "CC(C)Cc1ccc(cc1)C(C)C(=O)O" --pka "4.42:acid:COOH:14" --logp 3.97 \
    --ref "2:3.98,4:3.71,5:3.10,6:2.26,7:1.30,8:0.54,9:0.00" \
    --name "Ibuprofen" -o ibuprofen.html

python ionization.py "CCN(CC)CC(=O)Nc1c(C)cccc1C" --pka "7.92:base:NEt2:2" --logp 2.45 \
    --ref "3.89:-0.56,5.38:0.21,5.87:0.57,6.65:1.13,7.0:1.42,7.94:1.98" \
    --name "Lidocaine" -o lidocaine.html
```

Notes:

- Digitized points are read off small published plots (±0.05 log units).
- Shake-flask logD for lidocaine from the same source sits ~1 unit below the
  pH-metric points across pH 5–7; no offset choice fits it, so the
  measurement method matters as much as the model.
- Shake-flask logD₇.₄ for propranolol (≈1.2–1.3, Wenlock et al. 2011) is within
  0.2 of the calculated 1.40.
- Ion-pair partitioning varies by compound (Scherrer 2009: > 3 log units
  across amines), so the plateau far from the pKa stays the least certain part
  of the curve.

Zwitterion-capable compounds (Mazák & Noszál, ChemistryOpen 2019, species-specific
logP of cation, zwitterion, uncharged form and anion; reference curve built from
their microconstants):

| Compound | logP (uncharged) | zwitterion penalty | logD₇.₄ error, macro / micro pKa | pH 2 plateau error |
|---|---|---|---|---|
| Morphine | 0.93 | 3.03 | +0.13 / +0.02 | +0.04 |
| Naloxone | 2.18 | 3.03 | +0.24 / +0.03 | +1.13 |
| Naltrexone | 2.24 | 3.03 | +0.60 / +0.05 | +1.00 |
| Oxymorphone | 1.22 | 3.03 | +0.58 / +0.04 | +0.46 |

- All four zwitterions sit 3.03 below the uncharged form, which sets the
  `--zwitterion-offset` default. For these compounds the uncharged form
  dominates partitioning, so the zwitterion penalty barely moves logD.
- Their cations sit 3.0–4.1 below the uncharged form, hence the 0.5–1.1 error
  on the low-pH plateau for naloxone, naltrexone and oxymorphone.
- Cetirizine (Pagliara et al. 1998) is the opposite case: almost purely
  zwitterionic from pH 3.5–7.5 with logD = logP(zwitterion) = 1.5. Passing 1.5
  as an uncharged logP gives logD ≈ −1.5 (3 units off, and the tool warns);
  with `--logp-type zwitterion` the plateau is reproduced.

### Accuracy

The math is exact — it's the same equation Chemicalize and other tools use.
The accuracy of the plot depends entirely on the pKa values you provide.
Good sources for pKa values:

- [MolGpKa](https://xundrug.cn/molgpka/) — free online micro-pKa prediction from SMILES
- [Chemicalize](https://chemicalize.com/) (free tier, uses ChemAxon)
- [DrugBank](https://go.drugbank.com/)
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- Literature / experimental data

Any source that gives you per-site pKa values works. You just need to know
which group each value belongs to (acid vs base) and its pKa.

---

## File Structure

```
pka/
├── README.md          # this file
├── ionization.py      # CLI entry point
├── microspecies.py    # Henderson-Hasselbalch math engine
├── plot.py            # Plotly chart + RDKit structure rendering
├── test_micro.py      # quick test for the math engine
└── test_plot.py       # quick test for plot generation
```

### `microspecies.py`

Core math module. Key functions:

- `compute_microspecies(sites, ph_min, ph_max, ph_step)` — compute all 2^N
  species fractions across pH range.
- `auto_detect_atom_indices(smiles, sites)` — find RDKit atom indices for
  ionizable groups.
- `microspecies_smiles(smiles, sites, state)` — generate SMILES for a specific
  protonation state (adding/removing H, setting formal charges).
- `net_charge(sites, state)` — calculate net charge for a microspecies.

### `plot.py`

Visualization module. Key functions:

- `build_chart(ph_values, species)` — create the Plotly figure.
- `build_html(ph_values, species, smiles, sites, title)` — full HTML page
  with chart, structures, table, and annotated diagram.
- `save_and_open(html, filepath)` — write HTML and open in browser.

### `ionization.py`

CLI wrapper. Parses `--pka` strings, calls the engine and plotter.

---

## Examples

### Pseudoephedrine (2 sites: base + acid)

```bash
python ionization.py "C[C@@H](NC)[C@H](O)c1ccccc1" \
    --pka "9.52:base:NH,13.89:acid:OH" \
    --name "Pseudoephedrine" -o pseudoephedrine.html
```

| pH range | Dominant species | Charge |
|----------|-----------------|--------|
| 0–9 | NH⁺ / OH | +1 |
| 10–13 | NH / OH | 0 |
| 13.5+ | NH / O⁻ | −1 |

### Salbutamol (4 sites: 1 base + 3 acids)

```bash
python ionization.py "CC(C)(C)NCC(C1=CC(=C(C=C1)O)CO)O" \
    --pka "9.5:base:NH,9.8:acid:phenol-OH,12.4:acid:chain-OH,13.1:acid:benzyl-OH" \
    --name "Salbutamol" -o salbutamol.html
```

### Aspartic Acid (3 sites: 2 acids + 1 base)

```bash
python ionization.py "OC(=O)CC(N)C(=O)O" \
    --pka "1.88:acid:COOH1,3.65:acid:COOH2,9.60:base:NH2" \
    --name "Aspartic acid" -o aspartic.html
```
