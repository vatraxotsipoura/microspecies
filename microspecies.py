"""
Microspecies distribution engine using Henderson-Hasselbalch.

Given N ionizable sites (each with a pKa and type acid/base),
computes the fractional population of all 2^N protonation states
across a pH range.
"""

import itertools
import numpy as np
from dataclasses import dataclass, field

from rdkit import Chem


@dataclass
class IonizableSite:
    pka: float
    site_type: str  # "acid" or "base"
    label: str = ""  # e.g. "NH", "OH", "COOH"
    atom_idx: int | None = None  # RDKit atom index (auto-detected if None)


def compute_microspecies(sites: list[IonizableSite],
                         ph_min: float = 0.0,
                         ph_max: float = 14.0,
                         ph_step: float = 0.05):
    """
    Compute microspecies fractions across pH range.

    For each site at a given pH:
      - Acid (HA ⇌ A⁻ + H⁺):
          fraction protonated (HA) = 1 / (1 + 10^(pH - pKa))
          fraction deprotonated (A⁻) = 10^(pH - pKa) / (1 + 10^(pH - pKa))
      - Base (BH⁺ ⇌ B + H⁺):
          fraction protonated (BH⁺) = 1 / (1 + 10^(pH - pKa))
          fraction deprotonated (B) = 10^(pH - pKa) / (1 + 10^(pH - pKa))

    Note: for both acid and base, the math is identical — it's always
    pKa of the conjugate acid. The difference is only in labeling
    (which form is "default" / neutral).

    Each microspecies is a tuple of (protonated?, protonated?, ...)
    for each site. The fraction is the product of per-site terms.

    Returns:
        ph_values: 1D array of pH values
        species: list of dicts with keys:
            - "protonated": tuple of bools (one per site)
            - "label": human-readable label
            - "fractions": 1D array of fractional populations (0-100%)
    """
    n = len(sites)
    ph_values = np.arange(ph_min, ph_max + ph_step / 2, ph_step)

    # All 2^N combinations: True = protonated, False = deprotonated
    states = list(itertools.product([True, False], repeat=n))

    species = []
    for state in states:
        # Compute fraction for this microspecies at each pH
        fraction = np.ones_like(ph_values)
        for i, (site, is_protonated) in enumerate(zip(sites, state)):
            ratio = 10.0 ** (ph_values - site.pka)
            if is_protonated:
                fraction *= 1.0 / (1.0 + ratio)
            else:
                fraction *= ratio / (1.0 + ratio)

        label = _make_label(sites, state)
        species.append({
            "protonated": state,
            "label": label,
            "fractions": fraction * 100.0,  # percent
        })

    return ph_values, species


def _make_label(sites: list[IonizableSite], state: tuple[bool, ...]) -> str:
    """
    Generate a human-readable label for a microspecies.

    Examples for pseudoephedrine (NH base pKa 9.52, OH acid pKa 13.89):
      (True, True)   → "NH₂⁺ / OH"   (both protonated)
      (True, False)   → "NH₂⁺ / O⁻"
      (False, True)  → "NH / OH"      (neutral form)
      (False, False)  → "NH / O⁻"
    """
    parts = []
    for site, is_protonated in zip(sites, state):
        if site.label:
            base_label = site.label
        else:
            base_label = site.site_type.upper()

        if site.site_type == "acid":
            # Acid: protonated = neutral form (HA), deprotonated = charged (A⁻)
            if is_protonated:
                parts.append(base_label)
            else:
                parts.append(f"{base_label}(−)")
        else:
            # Base: protonated = charged (BH⁺), deprotonated = neutral (B)
            if is_protonated:
                parts.append(f"{base_label}(+)")
            else:
                parts.append(base_label)

    return " / ".join(parts)


def net_charge(sites: list[IonizableSite], state: tuple[bool, ...]) -> int:
    """Calculate net charge for a microspecies state."""
    charge = 0
    for site, is_protonated in zip(sites, state):
        if site.site_type == "acid" and not is_protonated:
            charge -= 1  # A⁻
        elif site.site_type == "base" and is_protonated:
            charge += 1  # BH⁺
    return charge


def _has_double_bond_to_o(atom) -> bool:
    """True if `atom` carries a C=O, S=O or P=O (or C=S) double bond."""
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if bond.GetBondType() == Chem.BondType.DOUBLE and other.GetSymbol() in ("O", "S"):
            return True
    return False


def _acid_rank(atom) -> int | None:
    """
    Rank an atom as an acid candidate (lower = more acidic), or None.

    0: carboxylic / sulfonic / phosphonic O-H   1: phenol O-H
    2: thiol S-H, sulfonamide / imide N-H        3: alcohol O-H
    """
    if atom.GetTotalNumHs() == 0 or atom.GetFormalCharge() != 0:
        return None
    sym = atom.GetSymbol()
    heavy = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
    if sym == "O":
        if len(heavy) != 1:
            return None
        nb = heavy[0]
        if _has_double_bond_to_o(nb):
            return 0
        if nb.GetIsAromatic():
            return 1
        return 3
    if sym == "S":
        return 2
    if sym == "N" and not atom.GetIsAromatic():
        sulfonyl = any(nb.GetSymbol() == "S" and _has_double_bond_to_o(nb) for nb in heavy)
        n_acyl = sum(_has_double_bond_to_o(nb) for nb in heavy)
        if sulfonyl or n_acyl >= 2:
            return 2
    return None


def _base_rank(atom) -> int | None:
    """
    Rank an atom as a base candidate (lower = more basic), or None.

    0: aliphatic amine, amidine/guanidine/imine N   1: aniline-type N
    2: pyridine-type aromatic N
    Amide, sulfonamide, nitro, nitrile, N-oxide and quaternary N are excluded.
    """
    if atom.GetSymbol() != "N" or atom.GetFormalCharge() != 0:
        return None
    heavy = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
    if any(nb.GetSymbol() == "O" for nb in heavy):
        return None
    if any(b.GetBondType() == Chem.BondType.TRIPLE for b in atom.GetBonds()):
        return None
    if atom.GetIsAromatic():
        if atom.GetTotalNumHs() == 0 and atom.GetDegree() == 2:
            return 2
        return None
    if any(_has_double_bond_to_o(nb) for nb in heavy):
        return None  # amide / sulfonamide / urea
    if any(nb.GetIsAromatic() for nb in heavy):
        return 1
    return 0


def auto_detect_atom_indices(smiles: str, sites: list[IonizableSite]) -> list[IonizableSite]:
    """
    Auto-detect atom indices for ionizable sites that don't have one set.

    Candidates are ranked by group type (see _acid_rank / _base_rank, ties in
    SMILES order). Acid sites are matched in order of increasing pKa to the
    most acidic candidates, base sites in order of decreasing pKa to the most
    basic ones. Atoms already assigned explicitly are never reused.

    Returns a new list with atom_idx populated.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    taken = {s.atom_idx for s in sites if s.atom_idx is not None}
    candidates = {"acid": [], "base": []}
    for atom in mol.GetAtoms():
        for site_type, rank_fn in (("acid", _acid_rank), ("base", _base_rank)):
            rank = rank_fn(atom)
            if rank is not None:
                candidates[site_type].append((rank, atom.GetIdx()))
    for c in candidates.values():
        c.sort()

    updated = [IonizableSite(pka=s.pka, site_type=s.site_type,
                             label=s.label, atom_idx=s.atom_idx) for s in sites]
    for site_type, strongest_first in (("acid", False), ("base", True)):
        pending = sorted((s for s in updated
                          if s.site_type == site_type and s.atom_idx is None),
                         key=lambda s: s.pka, reverse=strongest_first)
        for site in pending:
            pool = [idx for _, idx in candidates[site_type] if idx not in taken]
            if not pool:
                group = "O-H/S-H/N-H" if site_type == "acid" else "basic N"
                raise ValueError(
                    f"Not enough {group} groups found for {site_type} site "
                    f"'{site.label or site.pka}'. Specify atom_idx manually."
                )
            site.atom_idx = pool[0]
            taken.add(site.atom_idx)
    return updated


def validate_sites(smiles: str, sites: list[IonizableSite]) -> None:
    """
    Check that resolved sites point at sensible, distinct atoms.
    Raises ValueError with a readable message otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    seen = {}
    for site in sites:
        name = site.label or f"{site.site_type} pKa {site.pka:g}"
        idx = site.atom_idx
        if idx is None or not 0 <= idx < mol.GetNumAtoms():
            raise ValueError(f"Site '{name}': atom index {idx} is out of range "
                             f"(molecule has atoms 0-{mol.GetNumAtoms() - 1}).")
        if idx in seen:
            raise ValueError(f"Sites '{seen[idx]}' and '{name}' both use atom {idx}.")
        seen[idx] = name
        atom = mol.GetAtomWithIdx(idx)
        desc = f"atom {idx} ({atom.GetSymbol()}, {atom.GetTotalNumHs()} H)"
        if site.site_type == "acid" and (atom.GetSymbol() not in ("O", "S", "N")
                                         or atom.GetTotalNumHs() == 0):
            raise ValueError(f"Acid site '{name}': {desc} has no ionizable O-H/S-H/N-H.")
        if site.site_type == "base" and (atom.GetSymbol() != "N"
                                         or atom.GetFormalCharge() != 0):
            raise ValueError(f"Base site '{name}': {desc} is not a neutral nitrogen.")


def microspecies_smiles(smiles: str,
                        sites: list[IonizableSite],
                        state: tuple[bool, ...]) -> str:
    """
    Generate SMILES for a specific protonation microspecies.

    Modifies the base molecule by adding/removing protons and
    adjusting formal charges at the relevant atoms.
    """
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles))
    if mol is None:
        return smiles
    mol = Chem.RWMol(Chem.AddHs(mol))

    # Collect removals first (removing atoms shifts indices)
    h_atoms_to_remove = []

    for site, is_protonated in zip(sites, state):
        if site.atom_idx is None:
            continue
        atom = mol.GetAtomWithIdx(site.atom_idx)

        if site.site_type == "base" and is_protonated:
            # BH⁺: add an explicit H, set charge +1
            h_idx = mol.AddAtom(Chem.Atom(1))
            mol.AddBond(site.atom_idx, h_idx, Chem.BondType.SINGLE)
            atom.SetFormalCharge(1)

        elif site.site_type == "acid" and not is_protonated:
            # A⁻: remove one H neighbor, set charge -1
            for nb in atom.GetNeighbors():
                if nb.GetSymbol() == "H":
                    h_atoms_to_remove.append(nb.GetIdx())
                    break
            atom.SetFormalCharge(-1)

    # Remove H atoms in reverse order to preserve indices
    for h_idx in sorted(h_atoms_to_remove, reverse=True):
        mol.RemoveAtom(h_idx)

    try:
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        return smiles


def compute_logd(sites: list[IonizableSite],
                 species: list[dict],
                 logp: float,
                 cation_offset: float = 3.0,
                 anion_offset: float = 4.0,
                 zwitterion_offset: float = 3.0) -> np.ndarray:
    """
    Compute logD across the pH grid of `species`.

        logD(pH) = log10( Σ_i f_i(pH) · P_i )

    The neutral (fully uncharged) microspecies partitions with P = 10^logp.
    Each ionized group lowers a microspecies' logP (ion-pair partitioning,
    Avdeef's "diff 3-4" rule for octanol/0.15 M KCl): a protonated base by
    `cation_offset`, a deprotonated acid by `anion_offset`. Each (+, -) pair
    on the same microspecies counts once as `zwitterion_offset` instead, since
    the internal charges partly neutralize each other.

    Defaults are fitted to measured profiles: lidocaine 3.0 and propranolol
    2.7 (bases), ibuprofen 4.1 (acid), and 3.03 for the zwitterions of
    morphine, naloxone, naltrexone and oxymorphone (Mazak & Noszal 2019).
    """
    p_total = np.zeros_like(species[0]["fractions"])
    for s in species:
        n_cation = sum(site.site_type == "base" and is_prot
                       for site, is_prot in zip(sites, s["protonated"]))
        n_anion = sum(site.site_type == "acid" and not is_prot
                      for site, is_prot in zip(sites, s["protonated"]))
        n_pair = min(n_cation, n_anion)
        logp_i = (logp - zwitterion_offset * n_pair
                  - cation_offset * (n_cation - n_pair)
                  - anion_offset * (n_anion - n_pair))
        p_total += (s["fractions"] / 100.0) * 10.0 ** logp_i
    return np.log10(p_total)
