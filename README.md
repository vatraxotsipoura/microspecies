# pKa — Microspecies Ionization Curve Calculator

Generate Chemicalize-style microspecies distribution plots from a SMILES string and user-supplied pKa values.

**Input:** SMILES + pKa values → **Output:** interactive HTML with chart, molecule structures, and annotated pKa diagram.

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
python ionization.py "SMILES" --pka "PKA_SPEC" [--name NAME] [-o FILE]
```

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `SMILES` | Yes | SMILES string of the molecule (positional) |
| `--pka`  | Yes | Comma-separated pKa specifications (see format below) |
| `--name` | No | Molecule name (used in the HTML title) |
| `-o` / `--output` | No | Output HTML file path (default: auto-generated temp file) |

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
"2.34:acid:COOH:8,9.60:base:NH2:0"    # with explicit atom indices
```

### Atom Index Assignment

The program needs to know which atom in the molecule each pKa value belongs to, in order to render correct protonation states.

**Option A — Auto-detection (default):**
The program automatically finds candidate atoms:
- **Base sites**: non-aromatic nitrogen atoms (N)
- **Acid sites**: oxygen atoms bonded to hydrogen (O-H), sulfur atoms bonded to hydrogen (S-H)

Sites are matched in the order they appear in the `--pka` argument to the order candidates are found in the molecule. This works well for simple molecules with unambiguous ionizable groups.

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
