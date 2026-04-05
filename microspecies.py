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


def auto_detect_atom_indices(smiles: str, sites: list[IonizableSite]) -> list[IonizableSite]:
    """
    Auto-detect atom indices for ionizable sites that don't have one set.

    Strategy:
      - base sites: find N atoms (not aromatic, not amide) — assign in order
      - acid sites: find O-H groups (alcohols, carboxylic acids) — assign in order

    Returns a new list with atom_idx populated.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)

    # Collect candidate atoms
    base_candidates = []  # N atoms
    acid_candidates = []  # O atoms bonded to H

    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        idx = atom.GetIdx()
        if sym == "N" and not atom.GetIsAromatic():
            base_candidates.append(idx)
        elif sym == "O":
            # Check if bonded to at least one H
            has_h = any(nb.GetSymbol() == "H" for nb in atom.GetNeighbors())
            if has_h:
                acid_candidates.append(idx)
        elif sym == "S":
            has_h = any(nb.GetSymbol() == "H" for nb in atom.GetNeighbors())
            if has_h:
                acid_candidates.append(idx)

    base_iter = iter(base_candidates)
    acid_iter = iter(acid_candidates)

    updated = []
    for site in sites:
        if site.atom_idx is not None:
            updated.append(site)
            continue
        new_site = IonizableSite(
            pka=site.pka, site_type=site.site_type,
            label=site.label, atom_idx=site.atom_idx,
        )
        if site.site_type == "base":
            try:
                new_site.atom_idx = next(base_iter)
            except StopIteration:
                raise ValueError(
                    f"Not enough N atoms found for base site '{site.label}'. "
                    "Specify atom_idx manually."
                )
        else:  # acid
            try:
                new_site.atom_idx = next(acid_iter)
            except StopIteration:
                raise ValueError(
                    f"Not enough O-H groups found for acid site '{site.label}'. "
                    "Specify atom_idx manually."
                )
        updated.append(new_site)
    return updated


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
