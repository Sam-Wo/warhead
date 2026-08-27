"""G3b exatecan-partner search: the planted orthogonal partner (ATRi) must win,
the Top1i compounds must sink (they share the resistance mechanism)."""


def test_atri_is_top_partner(g3b_state):
    g = g3b_state.g3b.set_index("compound_id")
    # ATRi is the single best orthogonal partner (rank 1 by orthogonality).
    assert g["rank"].idxmin() == "ATRi_like"
    assert bool(g.loc["ATRi_like", "is_partner_candidate"])
    # It is potent specifically where SLFN11 is low (negative dependence).
    assert g.loc["ATRi_like", "slfn11_slope"] < -0.1


def test_top1i_compounds_are_worst_partners(g3b_state):
    g = g3b_state.g3b.set_index("compound_id")
    # Top1i compounds depend on SLFN11 (positive slope) -> not partners...
    for top1i in ("exatecan_like", "SN38_like"):
        assert g.loc[top1i, "slfn11_slope"] > 0.15
        assert not bool(g.loc[top1i, "is_partner_candidate"])
    # ...and ATRi clearly out-scores them on orthogonal potency.
    assert g.loc["ATRi_like", "orthogonality"] > g.loc["exatecan_like", "orthogonality"]
    assert g.loc["ATRi_like", "orthogonality"] > g.loc["SN38_like", "orthogonality"]


def test_efflux_controlled_partner_not_just_mdr(g3b_state):
    # The efflux substrate should not top the partner list once ABCB1 is removed.
    g = g3b_state.g3b.set_index("compound_id")
    assert g.loc["ATRi_like", "rank"] < g.loc["MDR_substrate_like", "rank"]
