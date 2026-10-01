"""
Plotly chart + RDKit structure rendering for microspecies distributions.

Generates a standalone HTML file with:
  - Interactive Plotly line chart (pH vs microspecies %)
  - Row of 2D molecule structures below, one per significant species,
    each with a colored border matching its curve.
"""

import base64
import io
import webbrowser
import tempfile
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw

from microspecies import IonizableSite, net_charge, auto_detect_atom_indices, microspecies_smiles

# Chemicalize-style color palette
COLORS = [
    "#1f77b4",  # blue
    "#d62728",  # red
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
]


def _mol_to_base64_png(smiles: str, size: tuple[int, int] = (200, 170)) -> str:
    """Render a molecule from SMILES to a base64-encoded PNG."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    img = Draw.MolToImage(mol, size=size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _add_colored_border(b64_png: str, color: str, border_width: int = 4) -> str:
    """Add a colored border around a base64 PNG image."""
    img_data = base64.b64decode(b64_png)
    img = Image.open(io.BytesIO(img_data))
    bordered = Image.new(
        "RGB",
        (img.width + 2 * border_width, img.height + 2 * border_width),
        color,
    )
    bordered.paste(img, (border_width, border_width))
    buf = io.BytesIO()
    bordered.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _filter_significant_species(species: list[dict],
                                 threshold: float = 1.0) -> list[dict]:
    """Keep only species that reach at least `threshold` % at some pH."""
    return [s for s in species if s["fractions"].max() >= threshold]


def _annotated_mol_png(smiles: str, sites: list[IonizableSite],
                       size: tuple[int, int] = (400, 350)) -> str:
    """
    Render the molecule with pKa values annotated next to each ionizable group.
    Returns base64-encoded PNG.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""

    # Set atom notes with pKa values at the resolved atom indices
    resolved = auto_detect_atom_indices(smiles, sites)
    for site in resolved:
        if site.atom_idx is not None:
            atom = mol.GetAtomWithIdx(site.atom_idx)
            atom.SetProp("atomNote", f"pKa {site.pka:.2f}")

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.annotationFontScale = 0.7
    opts.additionalAtomLabelPadding = 0.1
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png_bytes = drawer.GetDrawingText()
    return base64.b64encode(png_bytes).decode("ascii")


def build_chart(ph_values: np.ndarray,
                species: list[dict],
                title: str = "Microspecies Distribution") -> go.Figure:
    """Build the Plotly line chart."""
    fig = go.Figure()
    significant = _filter_significant_species(species)

    for i, s in enumerate(significant):
        color = COLORS[i % len(COLORS)]
        charge = net_charge.__wrapped__(s) if hasattr(net_charge, '__wrapped__') else None
        fig.add_trace(go.Scatter(
            x=ph_values,
            y=s["fractions"],
            mode="lines",
            name=s["label"],
            line=dict(color=color, width=2.5),
            hovertemplate="pH %{x:.1f}<br>%{y:.1f}%<extra>" + s["label"] + "</extra>",
        ))

    fig.update_layout(
        xaxis=dict(
            title="pH",
            range=[ph_values[0], ph_values[-1]],
            dtick=2,
            gridcolor="#eee",
        ),
        yaxis=dict(
            title="Microspecies distribution (%)",
            range=[0, 105],
            dtick=10,
            gridcolor="#eee",
        ),
        plot_bgcolor="white",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5,
            font=dict(size=11),
        ),
        margin=dict(l=60, r=30, t=80, b=60),
        width=900,
        height=500,
    )
    return fig


def logd_at(ph_values: np.ndarray, logd: np.ndarray, ph: float) -> float:
    """Interpolate logD at a given pH."""
    return float(np.interp(ph, ph_values, logd))


def build_logd_chart(ph_values: np.ndarray, logd: dict) -> go.Figure:
    """
    Build the logD vs pH chart.

    `logd` holds "logp", "logp_type", the offsets and "curve" (logD at each pH).
    """
    fig = go.Figure()
    curve = logd["curve"]
    name = "logD"

    fig.add_trace(go.Scatter(
        x=ph_values,
        y=curve,
        mode="lines",
        name=name,
        line=dict(color="#1f77b4", width=2.5),
        hovertemplate="pH %{x:.1f}<br>logD %{y:.2f}<extra></extra>",
    ))

    if logd.get("ref"):
        ref_ph, ref_logd = zip(*logd["ref"])
        fig.add_trace(go.Scatter(
            x=ref_ph,
            y=ref_logd,
            mode="markers",
            name="Literature",
            marker=dict(color="#ff7f0e", size=9, symbol="diamond",
                        line=dict(color="white", width=1)),
            hovertemplate="pH %{x:.2f}<br>logD %{y:.2f}<extra>literature</extra>",
        ))
        curve_min = min(curve.min(), min(ref_logd))
    else:
        curve_min = curve.min()

    logp_type = logd.get("logp_type", "neutral")
    fig.add_hline(y=logd["logp"], line=dict(color="#2ca02c", width=1, dash="dot"),
                  annotation_text=f"logP ({logp_type}) = {logd['logp']:.2f}",
                  annotation_position="top left")

    y_min = np.floor(curve_min) - 0.5
    y_max = np.ceil(max(logd["logp"], curve.max())) + 1

    fig.update_layout(
        xaxis=dict(title="pH", range=[ph_values[0], ph_values[-1]],
                   dtick=2, gridcolor="#eee"),
        yaxis=dict(title="logD", range=[y_min, y_max], dtick=1, gridcolor="#eee"),
        plot_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="center", x=0.5, font=dict(size=11)),
        margin=dict(l=60, r=30, t=80, b=60),
        width=900,
        height=450,
    )
    return fig


def build_html(ph_values: np.ndarray,
               species: list[dict],
               smiles: str,
               sites: list[IonizableSite],
               title: str = "Microspecies Distribution",
               logd: dict | None = None) -> str:
    """
    Build a complete standalone HTML page with chart + structure images.
    If `logd` is given (see build_logd_chart), a logD vs pH chart is added.
    """
    fig = build_chart(ph_values, species, title)
    chart_html = fig.to_html(full_html=False, include_plotlyjs="cdn")

    logd_section = ""
    if logd is not None:
        logd_fig = build_logd_chart(ph_values, logd)
        logd_chart = logd_fig.to_html(full_html=False, include_plotlyjs=False)
        model = (f"Ion-pair partitioning: logP − {logd['cation_offset']:g} per cationic "
                 f"group, − {logd['anion_offset']:g} per anionic group, "
                 f"− {logd.get('zwitterion_offset', 3.0):g} per (+,−) pair.")
        logd74 = logd_at(ph_values, logd["curve"], 7.4)
        logp_type = logd.get("logp_type", "neutral")
        logd_section = f"""
    <div class="chart-container">
        <h3 style="margin-top:0; color: #333;">Lipophilicity — logD vs pH</h3>
        <div style="color:#666; font-size:13px;">
            logP ({logp_type}) = {logd['logp']:.2f} · logD<sub>7.4</sub> = {logd74:.2f} · {model}
        </div>
        {logd_chart}
    </div>"""

    # Auto-detect atom indices if needed
    resolved_sites = auto_detect_atom_indices(smiles, sites)

    # Render molecule structure for each significant species
    significant = _filter_significant_species(species)

    structure_cards = []
    for i, s in enumerate(significant):
        color = COLORS[i % len(COLORS)]
        species_smi = microspecies_smiles(smiles, resolved_sites, s["protonated"])
        species_img = _mol_to_base64_png(species_smi)
        bordered_img = _add_colored_border(species_img, color) if species_img else ""

        charge = 0
        for site, is_prot in zip(sites, s["protonated"]):
            if site.site_type == "acid" and not is_prot:
                charge -= 1
            elif site.site_type == "base" and is_prot:
                charge += 1

        charge_str = f"+{charge}" if charge > 0 else str(charge) if charge < 0 else "0"

        card = f"""
        <div style="text-align:center; margin: 0 15px;">
            <img src="data:image/png;base64,{bordered_img}"
                 style="border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.12);" />
            <div style="margin-top: 8px; font-weight: 600; color: {color}; font-size: 14px;">
                {s['label']}
            </div>
            <div style="margin-top: 2px; color: #666; font-size: 12px;">
                charge: {charge_str}
            </div>
        </div>
        """
        structure_cards.append(card)

    structures_row = "\n".join(structure_cards)

    # pKa summary table
    pka_rows = ""
    for site in sites:
        pka_rows += f"<tr><td>{site.label or site.site_type}</td><td>{site.site_type}</td><td>{site.pka:.2f}</td></tr>\n"

    # Annotated molecule image with pKa values
    annotated_img = _annotated_mol_png(smiles, sites)
    annotated_section = ""
    if annotated_img:
        annotated_section = f"""
    <div class="annotated-mol">
        <h3 style="margin-top:0; color: #333;">pKa Sites</h3>
        <img src="data:image/png;base64,{annotated_img}"
             style="border-radius: 8px;" />
    </div>"""

    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            max-width: 1000px;
            margin: 0 auto;
            padding: 20px;
            background: #fafafa;
        }}
        .chart-container {{
            background: white;
            border-radius: 12px;
            padding: 20px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.08);
            margin-bottom: 24px;
        }}
        .structures-row {{
            display: flex;
            justify-content: center;
            flex-wrap: wrap;
            gap: 10px;
            padding: 20px;
            background: white;
            border-radius: 12px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.08);
            margin-bottom: 24px;
        }}
        .pka-table {{
            background: white;
            border-radius: 12px;
            padding: 20px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.08);
            margin-bottom: 24px;
        }}
        .annotated-mol {{
            background: white;
            border-radius: 12px;
            padding: 20px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.08);
            text-align: center;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
        }}
        th, td {{
            padding: 8px 16px;
            text-align: left;
            border-bottom: 1px solid #eee;
        }}
        th {{
            color: #666;
            font-size: 13px;
            text-transform: uppercase;
        }}
        .smiles {{
            font-family: monospace;
            color: #888;
            font-size: 13px;
            margin-bottom: 16px;
        }}
    </style>
</head>
<body>
    <h2 style="margin: 0 0 4px 0; color: #222;">{title}</h2>
    <div class="smiles">SMILES: {smiles}</div>

    <div class="chart-container">
        {chart_html}
    </div>

    <div class="structures-row">
        {structures_row}
    </div>
{logd_section}

    <div class="pka-table">
        <h3 style="margin-top:0; color: #333;">pKa Values</h3>
        <table>
            <thead>
                <tr><th>Site</th><th>Type</th><th>pKa</th></tr>
            </thead>
            <tbody>
                {pka_rows}
            </tbody>
        </table>
    </div>

    {annotated_section}
</body>
</html>"""
    return html


def save_and_open(html: str, filepath: str | None = None) -> str:
    """Save HTML to file and open in browser. Returns the file path."""
    if filepath is None:
        fd, filepath = tempfile.mkstemp(suffix=".html", prefix="microspecies_")
        import os
        os.close(fd)

    Path(filepath).write_text(html, encoding="utf-8")
    webbrowser.open(f"file:///{Path(filepath).resolve()}")
    return filepath
