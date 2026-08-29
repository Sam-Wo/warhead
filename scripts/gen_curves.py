"""Top-20 dose-response curves for the two wide-window screens:

  reports/ctrp_top20_curves.pdf   measured median+IQR, pooled from the 16-pt raw
                                  wells (data/interim/ctrp_curves.csv) with
                                  pool_measured_curves() to remove grid zig-zag;
  reports/prism_top20_curves.pdf  fitted 4PL with a free lower asymptote (Emax),
                                  params median-aggregated per compound so the
                                  drawn curve and its IC50/EC90 markers agree.

Both rank compounds by CRC potency (rank_potency, emax_max=0.5). The same curve
frames feed the interactive dashboard (scripts/gen_dashboard.py). Run:

    PYTHONPATH=src py scripts/gen_curves.py
"""
from warhead.reporting.screen_curves import (load_ctrp_curve_data, load_prism_curve_data,
                                             render_fitted_curves, render_measured_curves)

pooled, ctrp_summary = load_ctrp_curve_data(top=20)
print("CTRP curves:", len(ctrp_summary), "compounds")
print("wrote", render_measured_curves(pooled, ctrp_summary, source="CTRP v2",
                                       out_path="reports/ctrp_top20_curves.pdf"))

prism_summary = load_prism_curve_data(top=20)
print("PRISM curves:", len(prism_summary), "compounds")
print("wrote", render_fitted_curves(prism_summary, source="PRISM Repurposing (secondary)",
                                     out_path="reports/prism_top20_curves.pdf"))
