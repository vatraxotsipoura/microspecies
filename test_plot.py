"""Quick test of the full pipeline: microspecies + plot."""
from microspecies import IonizableSite, compute_microspecies
from plot import build_html, save_and_open

smiles = "C[C@@H](NC)[C@H](O)c1ccccc1"  # pseudoephedrine
sites = [
    IonizableSite(pka=9.52, site_type="base", label="NH"),
    IonizableSite(pka=13.89, site_type="acid", label="OH"),
]

ph, species = compute_microspecies(sites)
html = build_html(ph, species, smiles, sites, title="Pseudoephedrine — Microspecies Distribution")
path = save_and_open(html, "pseudoephedrine.html")
print(f"Saved to: {path}")
