"""Build reports/screens_dashboard.html - the interactive, one-tab-per-screen
dashboard (EC90/IC50 potency, HCC/CRC selectivity, clinical/ADC context).

    PYTHONPATH=src py scripts/gen_dashboard.py
"""
from pathlib import Path
import tempfile

import numpy as np

from warhead.analysis.clinical_tox import clinical_tox_table
from warhead.analysis.conjugation import (add_chemistry, add_efflux, add_window,
                                          delivery_scorecard, growth_lookup)
from warhead.analysis.indications import GASTRIC_CODES, label_indication, line_set
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

# conjugation-suitability scorecard (PRISM CRC top-20) for the synthesis tab.
# PRISM (not CTRP) for the G1 potency bar: its IC50 is calibrated so the reference
# payloads read true (exatecan 0.1 nM); CTRP reads ~4.7x weaker.
_pr = cans["PRISM Repurposing (secondary)"]
_sc = delivery_scorecard(_pr, rank_potency(_pr, "CRC", emax_max=0.5), clinical_tox_table(),
                         indication="CRC", growth_fn=growth_lookup(), top=20)
_sc = add_window(add_efflux(add_chemistry(_sc, fetch_smiles(_sc["compound"].tolist())), _pr))

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

# Extra-indication selectivity tabs (same analysis as CRC/HCC). CTRP/PRISM lines are
# relabelled via DepMap Oncotree; GDSC is rebuilt from gdsc2_ec90 with its TCGA code.
_FULL = {"GDSC2": "GDSC2", "PRISM": "PRISM Repurposing (secondary)", "CTRP v2": "CTRP v2"}
_g = pd.read_pickle("data/interim/gdsc2_ec90.pkl")


def _indication_tab(ind, screen_labels, *, disease=None, codes=None, gdsc_tcga=None, caveat=""):
    frames = {}
    for s in screen_labels:
        if s == "GDSC2":
            gmap = {"COREAD": "CRC", "LIHC": "HCC", gdsc_tcga: ind}
            frames[s] = pd.DataFrame({
                "source": "GDSC2", "compound": _g.drug_name, "target": _g.target, "moa": _g.pathway,
                "model_id": _g.cell_line, "indication": _g.tcga_desc.map(gmap).fillna("other"),
                "ic50_nM": _g.ic50_uM * 1e3, "ec90_nM": _g.ec90_uM * 1e3, "emax": np.nan,
                "ec90_extrapolated": _g.ec90_range.eq("extrapolated"), "clinical_phase": pd.NA})
        else:
            base = cans[_FULL[s]]
            ls = line_set(base["model_id"].unique(), disease_substr=disease, codes=codes)
            frames[s] = label_indication(base, ls, ind)
    return {"label": ind, "type": "indication", "meta": {}, "indication": ind, "screens": screen_labels,
            "caveat": caveat,
            "sel": {s: selectivity(frames[s], ind) for s in screen_labels},
            "rank": {s: rank_potency(frames[s], ind, emax_max=0.5) for s in screen_labels}}


screens.append(_indication_tab(
    "AML", ["CTRP v2", "GDSC2"], disease="Acute Myeloid Leukemia", gdsc_tcga="LAML",
    caveat=('Same selectivity analysis as the CRC/HCC tabs, on the two screens with AML lines (CTRP 30, '
            'GDSC 26). PRISM is absent - its secondary subset has no blood-cancer lines. '
            '<b style="color:#B87C22">Caveat:</b> leukaemia lines are globally hypersensitive in vitro '
            '(the cloud shifts up, median &Delta; +0.2 to +0.5 logs), so plain significance flags ~90% of '
            'compounds. Read the <b>top-right</b> of each panel (largest &Delta;) - those recover the real '
            'AML dependencies: Aurora&nbsp;B (barasertib), MCL1 (AZD5991), PLK1, BRD4, KIF11, cytarabine.')))
screens.append(_indication_tab(
    "Gastric", ["CTRP v2", "GDSC2", "PRISM"], codes=GASTRIC_CODES, gdsc_tcga="STAD",
    caveat=('Same selectivity analysis as CRC/HCC, on the three screens with gastric (stomach-'
            'adenocarcinoma) lines (CTRP 30, GDSC 24, PRISM 17). Gastric is a solid tumour, so there is no '
            'global-sensitivity confound as in AML - read the <b>top-right</b> (selective &amp; potent) for '
            'gastric-preferential hits. PRISM has the fewest lines, so its &Delta; estimates are the noisiest.')))

# gate-cascade reference (defines the G1-G6 the Conjugation tab applies)
screens.append({"label": "Cascade", "type": "cascade", "meta": {}})

out = render_dashboard(screens, out_path="reports/screens_dashboard.html", tested=tested,
                       venn_png=venn_png, overlap_counts=overlap_counts, overlap_totals=overlap_totals)
print("wrote", out)
