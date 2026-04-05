"""Quick test of the microspecies engine."""
from microspecies import IonizableSite, compute_microspecies, net_charge

sites = [
    IonizableSite(pka=9.52, site_type="base", label="NH"),
    IonizableSite(pka=13.89, site_type="acid", label="OH"),
]

ph, species = compute_microspecies(sites)

print(f"pH points: {len(ph)} ({ph[0]:.1f} to {ph[-1]:.1f})")
print(f"Microspecies: {len(species)}\n")

for s in species:
    charge = net_charge(sites, s["protonated"])
    peak = s["fractions"].max()
    peak_ph = ph[s["fractions"].argmax()]
    print(f"  {s['label']:20s}  charge={charge:+d}  peak={peak:.1f}% at pH {peak_ph:.1f}")

# Verify: at pH 7 the dominant species should be NH(+)/OH (cation)
idx_7 = abs(ph - 7.0).argmin()
print(f"\nAt pH 7.0:")
for s in species:
    print(f"  {s['label']:20s}  {s['fractions'][idx_7]:.1f}%")
