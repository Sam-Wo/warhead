"""WARHEAD command line.

    warhead demo            run the G1->G2b slice on the synthetic fixture and
                            write reports/proliferation_independence.{pdf,xlsx}
    warhead g2b  --raw-dir  run the same slice on real PRISM + DepMap data
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
    d.add_argument("--seed", type=int, default=20260827)
    d.add_argument("--n-lines", type=int, default=60)
    d.add_argument("--out", default=str(REPORTS))
    d.add_argument("--no-g1-filter", action="store_true",
                   help="run G2b on all compounds, not just G1 passers")
    d.set_defaults(func=cmd_demo)

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
