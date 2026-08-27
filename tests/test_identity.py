import pandas as pd

from warhead.identity import normalise_cellline_name, resolve_model_ids, resolution_report


def test_normalise():
    assert normalise_cellline_name("NCI-H23") == "NCIH23"
    assert normalise_cellline_name("hct 116") == "HCT116"


def test_resolve_model_ids_alias_and_passthrough():
    meta = pd.DataFrame(
        {
            "ModelID": ["ACH-000001", "ACH-000002"],
            "CellLineName": ["HCT-116", "SW620"],
            "StrippedCellLineName": ["HCT116", "SW620"],
        }
    )
    df = pd.DataFrame({"name": ["hct 116", "ACH-000002", "unknown-line"]})
    out = resolve_model_ids(df, "name", meta)
    assert out.loc[0, "ModelID"] == "ACH-000001"   # alias resolved
    assert out.loc[1, "ModelID"] == "ACH-000002"   # ACH passthrough
    assert pd.isna(out.loc[2, "ModelID"])           # unresolved, not fabricated

    rep = resolution_report(out["ModelID"])
    assert rep["resolved"] == 2 and rep["unresolved"] == 1
