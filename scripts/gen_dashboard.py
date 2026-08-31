"""Build reports/screens_dashboard.html - the interactive, one-tab-per-screen
dashboard (EC90/IC50 potency, HCC/CRC selectivity, clinical/ADC context).

    PYTHONPATH=src py scripts/gen_dashboard.py
"""
from pathlib import Path
import tempfile

from warhead.analysis.clinical_tox import clinical_tox_table
from warhead.analysis.conjugation import (add_chemistry, add_efflux, add_window,
                                          delivery_scorecard, growth_lookup)
from warhead.analysis.nci60 import _annotated
from warhead.analysis.pdxe import crc_response_ranking, load_metrics
from warhead.analysis.screen_potency import rank_potency, selectivity
from warhead.io.pubchem import fetch_smiles
from warhead.reporting.dashboard import render_dashboard
from warhead.reporting.screen_curves import load_ctrp_curve_data, load_prism_curve_data
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

# dose-response curves for the two wide-window screens (shared frames with gen_curves).
# Both are drawn as median + IQR across cell lines: CTRP from measured wells (markers),
# PRISM from the cross-line spread of the per-line 4PL fits (smooth line).
prism_band, prism_curve_summary = load_prism_curve_data(top=20, with_band=True)
prism_curves = {"note": "PRISM 4PL, median + IQR across lines", "markers": False,
                "pooled": prism_band, "summary": prism_curve_summary}
ctrp_pooled, ctrp_curve_summary = load_ctrp_curve_data(top=20)
ctrp_curves = {"note": "CTRP measured, median + IQR across lines", "markers": True,
               "pooled": ctrp_pooled, "summary": ctrp_curve_summary}
curves_by_src = {"PRISM Repurposing (secondary)": prism_curves, "CTRP v2": ctrp_curves}

# conjugation-suitability scorecard (CTRP CRC top-20) for the synthesis tab
ctrp = cans["CTRP v2"]
_sc = delivery_scorecard(ctrp, rank_potency(ctrp, "CRC", emax_max=0.5), clinical_tox_table(),
                         indication="CRC", growth_fn=growth_lookup(), top=20)
_sc = add_window(add_efflux(add_chemistry(_sc, fetch_smiles(_sc["compound"].tolist())), ctrp))

screens = [{"label": "Overlap", "type": "overlap", "meta": {}},        # placeholder meta (skipped below)
           {"label": "Conjugation", "type": "conjugation", "meta": {}, "scorecard": _sc}]
for src, short in [("GDSC2", "GDSC2"),
                   ("PRISM Repurposing (secondary)", "PRISM"),
                   ("CTRP v2", "CTRP v2")]:
    c = cans[src]
    screens.append({"label": short, "type": "dr", "meta": meta.loc[src].to_dict(),
                    "sel": {i: selectivity(c, i) for i in ("CRC", "HCC")},
                    "rank": {i: rank_potency(c, i, emax_max=0.5) for i in ("CRC", "HCC")},
                    "curves": curves_by_src.get(src)})

screens.append({"label": "PDXE", "type": "pdxe", "meta": meta.loc["PDXE (Novartis)"].to_dict(),
                "ranking": crc_response_ranking(load_metrics())})

sel = pd.read_pickle("data/interim/nci60_crc_selectivity.pkl").sort_values("delta_z", ascending=False)
screens.append({"label": "NCI-60", "type": "nci60", "meta": meta.loc["NCI-60"].to_dict(),
                "selectivity": _annotated(sel)})

out = render_dashboard(screens, out_path="reports/screens_dashboard.html", tested=tested,
                       venn_png=venn_png, overlap_counts=overlap_counts, overlap_totals=overlap_totals)
print("wrote", out)
