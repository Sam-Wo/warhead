"""WARHEAD command line.

    warhead demo            run the G1->G2b slice on the synthetic fixture and
                            write reports/proliferation_independence.{pdf,xlsx}
    warhead exatecan        run the G3b exatecan-partner search on the fixture and
                            write reports/exatecan_partner.{pdf,xlsx}
    warhead collateral      run the G2c collateral-lethality scan (CRC + HCC) and
                            write reports/collateral_lethality_*.{pdf,xlsx}
    warhead gdsc            GDSC EC90 potency + HCC/CRC selectivity (real data) ->
                            reports/gdsc_ec90_selectivity.pdf + ranking xlsx
    warhead g2b  --raw-dir  run the G2b slice on real PRISM + DepMap data
    warhead info            print the resolved gate thresholds and repo paths
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from . import __version__
from .config import REPO_ROOT, REPORTS, ensure_dirs, load_gates


def _write_outputs(state, out_dir: Path, config: dict) -> dict[str, Path]:
    from .reporting import render_proliferation_report

    out_dir.mkdir(parents=True, exist_ok=True)
    pdf = render_proliferation_report(
        state.g2b_stats, state.sensitivity, _model_meta_cache["meta"],
        out_path=out_dir / "proliferation_independence.pdf", config=config,
    )
    xlsx = out_dir / "proliferation_independence.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        state.g2b_stats.to_excel(xw, sheet_name="g2b_stats", index=False)
        state.g2b.passed.to_excel(xw, sheet_name="g2b_pass", index=False)
        state.g2b.failed.to_excel(xw, sheet_name="g2b_fail", index=False)
        pd.DataFrame(state.provenance).to_excel(xw, sheet_name="provenance", index=False)
    return {"pdf": pdf, "xlsx": xlsx}


# module-level stash so _write_outputs can see the model_meta used
_model_meta_cache: dict = {}


def _summarise(state) -> None:
    st = state.g2b_stats
    print(f"\nrefit/QC: {state.provenance[0]}")
    print(f"G1: {state.g1.n_pass}/{state.g1.n_in} compounds pass potency gate")
    print(f"G2b: {state.g2b.n_pass}/{state.g2b.n_in} compounds proliferation-independent")
    if st is not None and len(st):
        show = st.sort_values("std_slope")[["compound_id", "std_slope", "q", "n_lines"]]
        if "prolif_class" in st.columns:
            show = st.sort_values("std_slope")[["compound_id", "std_slope", "q", "n_lines", "prolif_class"]]
        with pd.option_context("display.width", 120, "display.max_columns", None):
            print("\nper-compound (sorted by proliferation dependence):")
            print(show.to_string(index=False))


def cmd_demo(args: argparse.Namespace) -> int:
    from .cascade import run_g2b_slice
    from .fixtures import make

    ensure_dirs()
    cfg = load_gates()
    data = make(seed=args.seed, n_lines=args.n_lines)
    _model_meta_cache["meta"] = data.model_meta
    state = run_g2b_slice(data.dose_response, data.model_meta, config=cfg,
                          apply_g1_filter=not args.no_g1_filter)
    _summarise(state)
    out = _write_outputs(state, Path(args.out), cfg)
    print(f"\nwrote {out['pdf']}")
    print(f"wrote {out['xlsx']}")

    # Report how well the refit recovered the planted proliferation betas.
    truth = data.compound_truth.set_index("compound_id")["prolif_beta"]
    st = state.g2b_stats.set_index("compound_id")
    common = st.index.intersection(truth.index)
    if len(common) > 2:
        corr = pd.Series({c: st.loc[c, "slope"] for c in common}).corr(truth.loc[common])
        print(f"\nrecovery check: corr(recovered slope, planted beta) = {corr:.3f} "
              f"over {len(common)} compounds")
    return 0


def cmd_exatecan(args: argparse.Namespace) -> int:
    from .cascade import run_g3b_slice
    from .fixtures import make
    from .reporting import render_exatecan_report

    ensure_dirs()
    cfg = load_gates()
    data = make(seed=args.seed, n_lines=args.n_lines)
    state = run_g3b_slice(data.dose_response, data.expression, config=cfg,
                          apply_g1_filter=not args.no_g1_filter)

    g3b = state.g3b
    print(f"\nG1: {state.g1.n_pass}/{state.g1.n_in} compounds pass potency gate")
    print(f"G3b: {int(g3b['is_partner_candidate'].sum())}/{len(g3b)} orthogonal partner candidates")
    cols = ["rank", "compound_id", "slfn11_slope", "orthogonality", "resistant_potency",
            "is_partner_candidate"]
    with pd.option_context("display.width", 130, "display.max_columns", None):
        print("\nexatecan-partner ranking (top orthogonal potency first):")
        print(g3b[cols].to_string(index=False))
    print(f"\n=> empirical top partner: {g3b.iloc[0]['compound_id']}")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf = render_exatecan_report(g3b, state.sensitivity, data.expression,
                                 out_path=out_dir / "exatecan_partner.pdf", config=cfg)
    xlsx = out_dir / "exatecan_partner.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        g3b.to_excel(xw, sheet_name="partner_ranking", index=False)
        pd.DataFrame(state.provenance).to_excel(xw, sheet_name="provenance", index=False)
    print(f"\nwrote {pdf}\nwrote {xlsx}")
    return 0


def cmd_collateral(args: argparse.Namespace) -> int:
    from .cascade import run_g2c_slice
    from .fixtures import make
    from .reporting import render_collateral_report

    ensure_dirs()
    cfg = load_gates()
    data = make(seed=args.seed, n_lines=args.n_lines)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    frames = {}
    for indication in ("CRC", "HCC"):
        state = run_g2c_slice(data.chronos, data.copy_number, data.tcga_recurrence,
                              data.model_meta, indication=indication, config=cfg)
        targets = state.g2c
        frames[indication] = targets
        prov = state.provenance[-1]
        print(f"\n=== G2c {indication} ===  candidate targets = {prov['n_candidate_targets']}"
              f"  |  {prov['positive_control']} recovered = {prov['positive_control_recovered']}")
        cols = ["gene", "delta", "q", "loss_frequency", "common_essential", "candidate", "collateral_score"]
        with pd.option_context("display.width", 130, "display.max_columns", None):
            print(targets[cols].to_string(index=False))
        pdf = render_collateral_report(
            targets, data.chronos, data.copy_number, indication=indication,
            out_path=out_dir / f"collateral_lethality_{indication}.pdf", config=cfg,
        )
        print(f"wrote {pdf}")

    xlsx = out_dir / "collateral_lethality_crc_hcc.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        for indication, targets in frames.items():
            targets.to_excel(xw, sheet_name=indication, index=False)
    print(f"\nwrote {xlsx}")
    return 0


def cmd_gdsc(args: argparse.Namespace) -> int:
    from .analysis.gdsc_ec90 import indication_ranking, selectivity
    from .io.gdsc import load_with_ec90
    from .reporting import render_gdsc_report

    ensure_dirs()
    df = load_with_ec90(Path(args.raw_dir) / "gdsc", dataset=args.dataset)
    print(f"{args.dataset}: {len(df)} curves | EC90 within tested range overall: "
          f"{(df['ec90_range'] == 'within').mean():.1%}")

    rankings, sels = {}, {}
    for ind in ("CRC", "HCC"):
        rankings[ind] = indication_ranking(df, ind)
        sels[ind] = selectivity(df, ind)
        top = rankings[ind].head(5)["drug_name"].tolist()
        nsel = int(sels[ind]["selective_potent"].sum()) if len(sels[ind]) else 0
        print(f"  {ind}: lowest-EC90 -> {', '.join(top)}  |  selective&potent = {nsel}")

    counts = {
        "total_compounds": int(df["drug_name"].nunique()),
        "total_lines": int(df["cell_line"].nunique()),
        "CRC_lines": int(df.loc[df["tcga_desc"] == "COREAD", "cell_line"].nunique()),
        "HCC_lines": int(df.loc[df["tcga_desc"] == "LIHC", "cell_line"].nunique()),
    }
    print(f"dataset: {counts['total_compounds']} compounds x {counts['total_lines']} lines "
          f"(CRC {counts['CRC_lines']}, HCC {counts['HCC_lines']})")
    out_dir = Path(args.out)
    pdf = render_gdsc_report(rankings, sels, out_path=out_dir / "gdsc_ec90_selectivity.pdf", counts=counts)
    xlsx = out_dir / "gdsc_ec90_ranking.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        for ind in ("CRC", "HCC"):
            rankings[ind].to_excel(xw, sheet_name=f"{ind}_ec90_rank", index=False)
            sels[ind].to_excel(xw, sheet_name=f"{ind}_selectivity", index=False)
    print(f"\nwrote {pdf}\nwrote {xlsx}")
    return 0


def cmd_gdsc_curves(args: argparse.Namespace) -> int:
    from .analysis.gdsc_curves import extract_raw_curves, pool_by_conc
    from .io.gdsc import load_with_ec90
    from .reporting.gdsc_curves import render_top15_curves, render_top15_fitted_curves

    ensure_dirs()
    gdir = Path(args.raw_dir) / "gdsc"
    df = load_with_ec90(gdir, dataset=args.dataset)

    # top-N most potent by median IC50 (all lines), summary per drug_id
    summ = (df.groupby(["drug_id", "drug_name"]).agg(
        target=("target", "first"),
        median_ic50_uM=("ic50_uM", "median"),
        median_ec90_uM=("ec90_uM", "median"),
        median_scal=("scal", "median"),
        min_conc_uM=("min_conc_uM", "median"),
        max_conc_uM=("max_conc_uM", "median"),
    ).reset_index().sort_values("median_ic50_uM").head(args.top))
    print(f"top {args.top} by median IC50: {', '.join(summ['drug_name'])}")

    raw_csv = gdir / "GDSC2_public_raw_data.csv"
    out_pdf = Path(args.out) / "gdsc_top15_curves.pdf"
    if raw_csv.exists() and not args.fitted:
        # measured curves from the raw well-level data
        drug_ids = [int(x) for x in summ["drug_id"]]
        print(f"streaming {raw_csv} ({raw_csv.stat().st_size/1e9:.2f} GB) for {len(drug_ids)} drugs...")
        pooled = pool_by_conc(extract_raw_curves(raw_csv, drug_ids))
        out = render_top15_curves(pooled, summ, out_path=out_pdf)
    else:
        # fitted-model curves (raw well data not present); shows tested window vs EC90
        if not raw_csv.exists():
            print("raw well-level file not present -> plotting GDSC FITTED curves "
                  "(download GDSC2_public_raw_data for measured points)")
        out = render_top15_fitted_curves(summ, out_path=out_pdf)
    print(f"wrote {out}")
    return 0


def _screen_reports(source_label, can, out_dir, stem):
    from .analysis.screen_potency import rank_potency, selectivity
    from .reporting.screen import render_screen_report
    R = {i: rank_potency(can, i, emax_max=0.5) for i in ("CRC", "HCC")}
    S = {i: selectivity(can, i) for i in ("CRC", "HCC")}
    counts = {"total_compounds": int(can["compound"].nunique()), "total_lines": int(can["model_id"].nunique()),
              "CRC_lines": int(can[can.indication == "CRC"].model_id.nunique()),
              "HCC_lines": int(can[can.indication == "HCC"].model_id.nunique())}
    pdf = render_screen_report(source_label, R, S, out_path=out_dir / f"{stem}.pdf", counts=counts)
    with pd.ExcelWriter(out_dir / f"{stem}.xlsx", engine="openpyxl") as xw:
        for i in ("CRC", "HCC"):
            R[i].to_excel(xw, sheet_name=f"{i}_ec90_rank", index=False)
            S[i].to_excel(xw, sheet_name=f"{i}_selectivity", index=False)
    print(f"{source_label}: {counts}  ->  wrote {pdf}")
    return R


def cmd_prism(args: argparse.Namespace) -> int:
    from .io.prism import load_canonical
    from .config import DATA_INTERIM
    ensure_dirs()
    can = load_canonical(Path(args.raw_dir) / "prism", Path(args.raw_dir) / "depmap" / "Model.csv")
    can.to_pickle(DATA_INTERIM / "prism_canonical.pkl")   # for warhead selectivity-html
    _screen_reports("PRISM Repurposing (secondary)", can, Path(args.out), "prism_ec90_selectivity")
    return 0


def cmd_ctrp(args: argparse.Namespace) -> int:
    from .io.ctrp import load_canonical, CTRP_EXPORT
    from .config import DATA_INTERIM
    ensure_dirs()
    can = load_canonical(args.export_csv or CTRP_EXPORT)
    can.to_pickle(DATA_INTERIM / "ctrp_canonical.pkl")    # for warhead selectivity-html
    _screen_reports("CTRP v2", can, Path(args.out), "ctrp_ec90_selectivity")
    return 0


def cmd_selectivity_html(args: argparse.Namespace) -> int:
    """Interactive (Plotly) HCC/CRC selectivity across whichever canonical screen
    pickles are available in data/interim/."""
    from .analysis.screen_potency import selectivity
    from .config import DATA_INTERIM
    from .reporting.interactive import render_selectivity_html
    ensure_dirs()
    sources = {"GDSC2": "gdsc_canonical.pkl", "PRISM Repurposing (secondary)": "prism_canonical.pkl",
               "CTRP v2": "ctrp_canonical.pkl"}
    sel = {}
    for src, fn in sources.items():
        p = DATA_INTERIM / fn
        if p.exists():
            can = pd.read_pickle(p)
            sel[src] = {i: selectivity(can, i) for i in ("CRC", "HCC")}
    if not sel:
        print("no canonical screen pickles in data/interim/ (run warhead prism / ctrp first)")
        return 1
    out = render_selectivity_html(sel, out_path=Path(args.out) / "selectivity_interactive.html")
    print(f"wrote {out}  (sources: {', '.join(sel)})")
    return 0


def cmd_pdxe(args: argparse.Namespace) -> int:
    from .analysis.pdxe import load_metrics, crc_response_ranking, crc_response_selectivity
    from .reporting.pdxe import render_pdxe_report
    ensure_dirs()
    m = load_metrics()
    rank, sel = crc_response_ranking(m), crc_response_selectivity(m)
    ncrc = int(m[m["Tumor Type"] == "CRC"].Model.nunique())
    pdf = render_pdxe_report(rank, sel, out_path=Path(args.out) / "pdxe_crc_response.pdf", n_crc_models=ncrc)
    with pd.ExcelWriter(Path(args.out) / "pdxe_crc_response.xlsx", engine="openpyxl") as xw:
        rank.to_excel(xw, sheet_name="crc_response_rank", index=False)
        sel.to_excel(xw, sheet_name="crc_response_selectivity", index=False)
    print(f"PDXE CRC ({ncrc} models) -> wrote {pdf}")
    return 0


def cmd_clinical_tox(args: argparse.Namespace) -> int:
    from .analysis.clinical_tox import clinical_tox_table
    from .reporting.clinical_tox_fig import render_clinical_tox
    ensure_dirs()
    t = clinical_tox_table()
    render_clinical_tox(t, out_path=Path(args.out) / "clinical_toxicity.pdf")
    t.to_excel(Path(args.out) / "clinical_toxicity.xlsx", index=False)
    print(f"clinical/toxicity table ({len(t)} compounds) -> wrote reports/clinical_toxicity.{{pdf,xlsx}}")
    return 0


def cmd_g2b_real(args: argparse.Namespace) -> int:
    from .analysis.gdsc_proliferation import run_real_g2b
    from .reporting import render_real_g2b_report

    ensure_dirs()
    cfg = load_gates()
    res = run_real_g2b(Path(args.raw_dir) / "gdsc", Path(args.raw_dir) / "depmap", config=cfg)
    st = res.stats
    print(f"real G2b (GDSC2 x DepMap growth): {res.n_lines} lines, {len(st)} compounds")
    print(f"proliferation-independent (pass G2b): {res.gate.n_pass}/{res.gate.n_in}")
    cols = ["compound_id", "target", "pathway", "std_slope", "q", "n_lines"]
    with pd.option_context("display.width", 150, "display.max_columns", None):
        print("\nMOST proliferation-DEPENDENT (top 10 positive slope):")
        print(st.sort_values("std_slope", ascending=False).head(10)[cols].to_string(index=False))
        print("\nMOST proliferation-INDEPENDENT among the potent set (flattest, not significant):")
        flat = st[st["q"] > cfg["g2"]["proliferation"]["fdr_alpha"]]
        print(flat.reindex(flat["std_slope"].abs().sort_values().index).head(10)[cols].to_string(index=False))

    out_dir = Path(args.out)
    pdf = render_real_g2b_report(
        res.stats, res.sensitivity, res.model_meta,
        out_path=out_dir / "proliferation_independence_real.pdf", config=cfg,
    )
    xlsx = out_dir / "proliferation_independence_real.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        res.stats.sort_values("std_slope", ascending=False).to_excel(xw, sheet_name="g2b_real", index=False)
    print(f"\nwrote {pdf}\nwrote {xlsx}")
    return 0


def cmd_g2b(args: argparse.Namespace) -> int:
    from .cascade import run_g2b_slice
    from .io.depmap import load_model_metadata
    from .io.prism import load_secondary_doseresponse

    ensure_dirs()
    cfg = load_gates()
    raw = Path(args.raw_dir)
    dose_response = load_secondary_doseresponse(raw / "prism")
    model_meta = load_model_metadata(raw / "depmap")
    _model_meta_cache["meta"] = model_meta
    state = run_g2b_slice(dose_response, model_meta, config=cfg,
                          apply_g1_filter=not args.no_g1_filter)
    _summarise(state)
    out = _write_outputs(state, Path(args.out), cfg)
    print(f"\nwrote {out['pdf']}\nwrote {out['xlsx']}")
    return 0


def cmd_info(args: argparse.Namespace) -> int:
    cfg = load_gates()
    print(f"WARHEAD {__version__}")
    print(f"repo root : {REPO_ROOT}")
    print(f"reports   : {REPORTS}")
    print("\nGate thresholds (config/gates.yaml):")
    import yaml
    print(yaml.safe_dump(cfg, sort_keys=False))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="warhead", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"warhead {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("demo", help="run the G1->G2b slice on the synthetic fixture")
    d.add_argument("--seed", type=int, default=7)
    d.add_argument("--n-lines", type=int, default=80)
    d.add_argument("--out", default=str(REPORTS))
    d.add_argument("--no-g1-filter", action="store_true",
                   help="run G2b on all compounds, not just G1 passers")
    d.set_defaults(func=cmd_demo)

    e = sub.add_parser("exatecan", help="G3b orthogonal-resistance (exatecan partner) on the fixture")
    e.add_argument("--seed", type=int, default=7)
    e.add_argument("--n-lines", type=int, default=80)
    e.add_argument("--out", default=str(REPORTS))
    e.add_argument("--no-g1-filter", action="store_true")
    e.set_defaults(func=cmd_exatecan)

    c = sub.add_parser("collateral", help="G2c collateral-lethality scan (CRC + HCC) on the fixture")
    c.add_argument("--seed", type=int, default=7)
    c.add_argument("--n-lines", type=int, default=80)
    c.add_argument("--out", default=str(REPORTS))
    c.set_defaults(func=cmd_collateral)

    gd = sub.add_parser("gdsc", help="GDSC EC90 ranking + HCC/CRC selectivity (real data)")
    gd.add_argument("--raw-dir", default="data/raw")
    gd.add_argument("--dataset", default="GDSC2", choices=["GDSC1", "GDSC2"])
    gd.add_argument("--out", default=str(REPORTS))
    gd.set_defaults(func=cmd_gdsc)

    gr = sub.add_parser("g2b-real", help="G2b proliferation independence on real GDSC2 x DepMap growth")
    gr.add_argument("--raw-dir", default="data/raw")
    gr.add_argument("--out", default=str(REPORTS))
    gr.set_defaults(func=cmd_g2b_real)

    gc = sub.add_parser("gdsc-curves", help="measured dose-response curves for the most potent GDSC2 compounds")
    gc.add_argument("--raw-dir", default="data/raw")
    gc.add_argument("--dataset", default="GDSC2", choices=["GDSC1", "GDSC2"])
    gc.add_argument("--top", type=int, default=15)
    gc.add_argument("--fitted", action="store_true", help="force fitted-model curves even if raw data is present")
    gc.add_argument("--out", default=str(REPORTS))
    gc.set_defaults(func=cmd_gdsc_curves)

    pm = sub.add_parser("prism", help="PRISM Repurposing EC90 ranking + HCC/CRC selectivity (real data)")
    pm.add_argument("--raw-dir", default="data/raw")
    pm.add_argument("--out", default=str(REPORTS))
    pm.set_defaults(func=cmd_prism)

    cp = sub.add_parser("ctrp", help="CTRP v2 EC90 ranking + HCC/CRC selectivity (via base-R export)")
    cp.add_argument("--export-csv", default=None, help="ctrp_export.csv from scripts/ctrp_export.R")
    cp.add_argument("--out", default=str(REPORTS))
    cp.set_defaults(func=cmd_ctrp)

    sh = sub.add_parser("selectivity-html", help="interactive Plotly HCC/CRC selectivity across screens")
    sh.add_argument("--out", default=str(REPORTS))
    sh.set_defaults(func=cmd_selectivity_html)

    px = sub.add_parser("pdxe", help="Novartis PDXE in-vivo CRC response ranking + selectivity")
    px.add_argument("--out", default=str(REPORTS))
    px.set_defaults(func=cmd_pdxe)

    ct = sub.add_parser("clinical-tox", help="curated clinical-validation + patient-toxicity table")
    ct.add_argument("--out", default=str(REPORTS))
    ct.set_defaults(func=cmd_clinical_tox)

    g = sub.add_parser("g2b", help="run the slice on real PRISM + DepMap data")
    g.add_argument("--raw-dir", default="data/raw")
    g.add_argument("--out", default=str(REPORTS))
    g.add_argument("--no-g1-filter", action="store_true")
    g.set_defaults(func=cmd_g2b)

    i = sub.add_parser("info", help="print resolved gate thresholds and paths")
    i.set_defaults(func=cmd_info)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv if argv is not None else sys.argv[1:])
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
