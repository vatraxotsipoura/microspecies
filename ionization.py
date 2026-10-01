"""
CLI entry point for microspecies ionization curve plotting.

Usage:
    python ionization.py "SMILES" --pka "9.52:base:NH,13.89:acid:OH"
    python ionization.py "SMILES" --pka "9.52:base:NH,13.89:acid:OH" --name "Pseudoephedrine"
    python ionization.py "SMILES" --pka "9.52:base:NH:4,13.89:acid:OH:12"  # with atom indices
    python ionization.py "SMILES" --pka "2.34:acid,9.60:base" -o glycine.html
    python ionization.py "SMILES" --pka "8.19:base" --logp 3.77            # adds logD vs pH chart
    python ionization.py "SMILES" --pka "9.53:base" --logp 3.48 --ref "2:0.78,7.4:1.28"  # compare to literature
"""

import argparse
import sys

import numpy as np

from microspecies import (IonizableSite, auto_detect_atom_indices, compute_logd,
                          compute_microspecies, net_charge, validate_sites)
from plot import build_html, logd_at, save_and_open


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

        try:
            pka = float(fields[0])
        except ValueError:
            print(f"Error: invalid pKa value '{fields[0]}' in '{part}'")
            sys.exit(1)
        site_type = fields[1].lower()
        if site_type not in ("acid", "base"):
            print(f"Error: site_type must be 'acid' or 'base', got '{site_type}'")
            sys.exit(1)

        label = fields[2] if len(fields) > 2 else ""
        try:
            atom_idx = int(fields[3]) if len(fields) > 3 else None
        except ValueError:
            print(f"Error: invalid atom index '{fields[3]}' in '{part}'")
            sys.exit(1)

        sites.append(IonizableSite(pka=pka, site_type=site_type,
                                   label=label, atom_idx=atom_idx))
    return sites


def parse_ref_string(ref_str: str) -> list[tuple[float, float]]:
    """Parse literature logD points: "pH:logD,pH:logD,..." """
    points = []
    for part in ref_str.split(","):
        fields = part.strip().split(":")
        if len(fields) != 2:
            print(f"Error: invalid --ref point '{part}'. Expected 'pH:logD'")
            sys.exit(1)
        try:
            ph, logd = float(fields[0]), float(fields[1])
        except ValueError:
            print(f"Error: invalid --ref point '{part}'. Expected numbers 'pH:logD'")
            sys.exit(1)
        if not 0.0 <= ph <= 14.0:
            print(f"Error: --ref pH {ph:g} is outside the computed range 0-14")
            sys.exit(1)
        points.append((ph, logd))
    return points


def _warn_if_zwitterion_dominates(sites: list[IonizableSite], species: list[dict]) -> None:
    """
    Warn when the zwitterion outnumbers the uncharged form by >10x: literature
    "logP" for such compounds is usually the zwitterion's, not the neutral's.
    """
    neutral = zwitterion = 0.0
    for s in species:
        cation = any(site.site_type == "base" and prot
                     for site, prot in zip(sites, s["protonated"]))
        anion = any(site.site_type == "acid" and not prot
                    for site, prot in zip(sites, s["protonated"]))
        if net_charge(sites, s["protonated"]) != 0:
            continue
        if cation and anion:
            zwitterion = max(zwitterion, s["fractions"].max())
        elif not cation and not anion:
            neutral = max(neutral, s["fractions"].max())
    if zwitterion > 10 * neutral:
        print(f"Warning: the zwitterion (peak {zwitterion:.1f}%) far outnumbers the "
              f"uncharged form (peak {neutral:.2g}%).")
        print("         Literature logP for such compounds is usually the zwitterion's;")
        print("         if so, rerun with --logp-type zwitterion.")
        print()


def main():
    # Species labels contain '−'; Windows consoles may not use UTF-8
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

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
    parser.add_argument("--logp", type=float, default=None,
                        help="logP of the neutral species (see --logp-type); adds a logD vs pH chart")
    parser.add_argument("--logp-type", choices=("neutral", "zwitterion"), default="neutral",
                        help="What --logp refers to: the uncharged form (default) or the "
                             "zwitterion, for compounds that exist almost only as zwitterions")
    parser.add_argument("--cation-offset", type=float, default=3.0,
                        help="logP penalty per protonated base group (default 3.0)")
    parser.add_argument("--anion-offset", type=float, default=4.0,
                        help="logP penalty per deprotonated acid group (default 4.0)")
    parser.add_argument("--zwitterion-offset", type=float, default=3.0,
                        help="logP penalty per (+,-) pair on one microspecies (default 3.0)")
    parser.add_argument("--ref", default=None,
                        help="Literature logD points to compare against: 'pH:logD,...' (needs --logp)")

    args = parser.parse_args()
    if args.ref and args.logp is None:
        print("Error: --ref needs --logp")
        sys.exit(1)
    ref_points = parse_ref_string(args.ref) if args.ref else None

    # Parse pKa sites and resolve/validate their atoms before computing anything
    sites = parse_pka_string(args.pka)
    try:
        sites = auto_detect_atom_indices(args.smiles, sites)
        validate_sites(args.smiles, sites)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    if args.logp_type == "zwitterion" and not (
            any(s.site_type == "acid" for s in sites) and any(s.site_type == "base" for s in sites)):
        print("Error: --logp-type zwitterion needs at least one acid and one base site")
        sys.exit(1)

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
    width = max(25, *(len(s["label"]) for s in species))
    print(f"{'Species':<{width}s} {'Charge':>6s}  {'Peak %':>6s}  {'at pH':>5s}")
    print("-" * (width + 25))
    for s in species:
        charge = net_charge(sites, s["protonated"])
        peak = s["fractions"].max()
        peak_ph = ph[s["fractions"].argmax()]
        if peak >= 0.1:
            charge_str = f"{charge:+d}" if charge else "0"
            print(f"{s['label']:<{width}s} {charge_str:>6s}  {peak:>5.1f}%  {peak_ph:>5.1f}")
    print()

    # logD
    logd = None
    if args.logp is not None:
        # The model is anchored on the uncharged form; a zwitterion logP is
        # converted back with the zwitterion penalty.
        logp_neutral = args.logp
        if args.logp_type == "zwitterion":
            logp_neutral = args.logp + args.zwitterion_offset
        else:
            _warn_if_zwitterion_dominates(sites, species)
        logd = {
            "logp": args.logp,
            "logp_type": args.logp_type,
            "cation_offset": args.cation_offset,
            "anion_offset": args.anion_offset,
            "zwitterion_offset": args.zwitterion_offset,
            "curve": compute_logd(sites, species, logp_neutral, args.cation_offset,
                                  args.anion_offset, args.zwitterion_offset),
        }
        print(f"logP ({args.logp_type}) = {args.logp:.2f}")
        print(f"{'pH':>5s}  {'logD':>6s}")
        for p in (1.0, 2.0, 5.5, 7.4, 9.0, 12.0):
            print(f"{p:>5.1f}  {logd_at(ph, logd['curve'], p):>6.2f}")
        print()

    # Comparison with literature logD points
    if ref_points:
        logd["ref"] = ref_points
        print(f"{'pH':>5s}  {'ref':>6s}  {'calc':>6s}  {'diff':>6s}")
        diffs = []
        for p, ref in logd["ref"]:
            calc = logd_at(ph, logd["curve"], p)
            diffs.append(calc - ref)
            print(f"{p:>5.2f}  {ref:>6.2f}  {calc:>6.2f}  {calc - ref:>+6.2f}")
        rmse = float(np.sqrt(np.mean(np.square(diffs))))
        print(f"RMSE = {rmse:.2f}, max |diff| = {max(abs(d) for d in diffs):.2f}")
        print()

    # Build and save HTML
    html = build_html(ph, species, args.smiles, sites, title, logd)
    path = save_and_open(html, args.output)
    print(f"Saved to: {path}")


if __name__ == "__main__":
    main()
