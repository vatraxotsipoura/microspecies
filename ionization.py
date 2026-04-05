"""
CLI entry point for microspecies ionization curve plotting.

Usage:
    python ionization.py "SMILES" --pka "9.52:base:NH,13.89:acid:OH"
    python ionization.py "SMILES" --pka "9.52:base:NH,13.89:acid:OH" --name "Pseudoephedrine"
    python ionization.py "SMILES" --pka "9.52:base:NH:4,13.89:acid:OH:12"  # with atom indices
    python ionization.py "SMILES" --pka "2.34:acid,9.60:base" -o glycine.html
"""

import argparse
import sys

from microspecies import IonizableSite, compute_microspecies, net_charge
from plot import build_html, save_and_open


def parse_pka_string(pka_str: str) -> list[IonizableSite]:
    """
    Parse pKa specification string.

    Format: "pka:type[:label[:atom_idx]], ..."
    Examples:
        "9.52:base:NH,13.89:acid:OH"
        "9.52:base,13.89:acid"
        "9.52:base:NH:4,13.89:acid:OH:12"
    """
    sites = []
    for part in pka_str.split(","):
        fields = part.strip().split(":")
        if len(fields) < 2:
            print(f"Error: invalid pKa spec '{part}'. Expected 'pka:type[:label[:atom_idx]]'")
            sys.exit(1)

        pka = float(fields[0])
        site_type = fields[1].lower()
        if site_type not in ("acid", "base"):
            print(f"Error: site_type must be 'acid' or 'base', got '{site_type}'")
            sys.exit(1)

        label = fields[2] if len(fields) > 2 else ""
        atom_idx = int(fields[3]) if len(fields) > 3 else None

        sites.append(IonizableSite(pka=pka, site_type=site_type,
                                   label=label, atom_idx=atom_idx))
    return sites


def main():
    parser = argparse.ArgumentParser(
        description="Plot microspecies ionization curves from SMILES + pKa values.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python ionization.py "C[C@@H](NC)[C@H](O)c1ccccc1" --pka "9.52:base:NH,13.89:acid:OH" --name "Pseudoephedrine"
  python ionization.py "NCC(=O)O" --pka "2.34:acid:COOH,9.60:base:NH2" --name "Glycine"
  python ionization.py "OC(=O)CC(N)C(=O)O" --pka "1.88:acid:COOH1,3.65:acid:COOH2,9.60:base:NH2" --name "Aspartic acid"
        """,
    )
    parser.add_argument("smiles", help="SMILES string of the molecule")
    parser.add_argument("--pka", required=True,
                        help="Comma-separated pKa specs: 'pka:type[:label[:atom_idx]],...'")
    parser.add_argument("--name", default=None,
                        help="Molecule name (used in title)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output HTML file path (default: auto-generated temp file)")

    args = parser.parse_args()

    # Parse pKa sites
    sites = parse_pka_string(args.pka)

    # Build title
    if args.name:
        title = f"{args.name} — Microspecies Distribution"
    else:
        title = "Microspecies Distribution"

    # Compute
    print(f"SMILES: {args.smiles}")
    print(f"Sites:  {len(sites)}")
    for s in sites:
        idx_str = f" (atom {s.atom_idx})" if s.atom_idx is not None else ""
        print(f"  pKa {s.pka:6.2f}  {s.site_type:5s}  {s.label}{idx_str}")
    print()

    ph, species = compute_microspecies(sites)

    # Print summary
    print(f"{'Species':<25s} {'Charge':>6s}  {'Peak %':>6s}  {'at pH':>5s}")
    print("-" * 50)
    for s in species:
        charge = net_charge(sites, s["protonated"])
        peak = s["fractions"].max()
        peak_ph = ph[s["fractions"].argmax()]
        if peak >= 0.1:
            charge_str = f"{charge:+d}"
            print(f"{s['label']:<25s} {charge_str:>6s}  {peak:>5.1f}%  {peak_ph:>5.1f}")
    print()

    # Build and save HTML
    html = build_html(ph, species, args.smiles, sites, title)
    path = save_and_open(html, args.output)
    print(f"Saved to: {path}")


if __name__ == "__main__":
    main()
