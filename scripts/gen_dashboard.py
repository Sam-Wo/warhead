"""Build reports/screens_dashboard.html - the interactive, one-tab-per-screen
dashboard (EC90/IC50 potency, HCC/CRC selectivity, clinical/ADC context).

    PYTHONPATH=src py scripts/gen_dashboard.py
"""
from pathlib import Path
import tempfile

from warhead.analysis.nci60 import _annotated
from warhead.analysis.pdxe import crc_response_ranking, load_metrics
from warhead.analysis.screen_potency import rank_potency, selectivity
from warhead.reporting.dashboard import render_dashboard
from warhead.reporting.screen_data import load_dr_screens
from warhead.reporting.screen_overlap import (DR_SOURCES, render_overlap_venn,
                                              tested_sets, three_set_counts)
import pandas as pd

cans, counts, meta = load_dr_screens()
tested = tested_sets(cans)
meta = meta.set_index("source")

# overlap Venn (rendered to a PNG we inline into the Overlap tab) + region counts
venn_path = Path(tempfile.gettempdir()) / "warhead_venn.png"
render_overlap_venn(cans, out_path=venn_path)
venn_png = venn_path.read_bytes()
A, B, C = (tested[s] for s in DR_SOURCES)
overlap_counts = three_set_counts(A, B, C)
overlap_totals = {s: len(tested[s]) for s in DR_SOURCES}

screens = [{"label": "Overlap", "type": "overlap", "meta": {}}]  # placeholder meta (skipped below)
for src, short in [("GDSC2", "GDSC2"),
                   ("PRISM Repurposing (secondary)", "PRISM"),
                   ("CTRP v2", "CTRP v2")]:
    c = cans[src]
    screens.append({"label": short, "type": "dr", "meta": meta.loc[src].to_dict(),
                    "sel": {i: selectivity(c, i) for i in ("CRC", "HCC")},
                    "rank": {i: rank_potency(c, i, emax_max=0.5) for i in ("CRC", "HCC")}})

screens.append({"label": "PDXE", "type": "pdxe", "meta": meta.loc["PDXE (Novartis)"].to_dict(),
                "ranking": crc_response_ranking(load_metrics())})

sel = pd.read_pickle("data/interim/nci60_crc_selectivity.pkl").sort_values("delta_z", ascending=False)
screens.append({"label": "NCI-60", "type": "nci60", "meta": meta.loc["NCI-60"].to_dict(),
                "selectivity": _annotated(sel)})

out = render_dashboard(screens, out_path="reports/screens_dashboard.html", tested=tested,
                       venn_png=venn_png, overlap_counts=overlap_counts, overlap_totals=overlap_totals)
print("wrote", out)
