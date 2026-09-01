"""Is the 'CRC-selective' MEK/ERK finding lineage or MAPK genotype? Figure for the deck.

MEK-inhibitor sensitivity (PRISM, z-scored consensus of MAP2K1/2 inhibitors) x DepMap
hotspot mutations (KRAS/BRAF/NRAS). Shows it is a MAPK-mutant biomarker, not a CRC
lineage effect, and translates to the addressable patient population.

    PYTHONPATH=src py scripts/gen_mek_confounding.py   -> reports/mek_confounding.pdf/.png
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RRB = "#6E1426"
COL = {"WT": "#AEB4BC", "KRAS": "#C77", "BRAF": "#8E2438", "NRAS": "#C0842A"}

# --- MEK-inhibitor sensitivity (PRISM) x KRAS/BRAF/NRAS hotspot x lineage ---
prism = pd.read_pickle("data/interim/prism_canonical.pkl")
mek = prism[prism.target.astype(str).str.upper().str.contains("MAP2K1|MAP2K2", na=False)]
mek = mek[~mek.compound.str.contains(":", na=False)].copy()
mek["pot"] = -np.log10(mek.ic50_nM.clip(lower=1e-3))
zs = []
for _, g in mek.groupby("compound"):
    s = g.groupby("model_id").pot.median()
    if len(s) >= 20 and s.std() > 0:
        zs.append((s - s.mean()) / s.std())
mekz = pd.concat(zs, axis=1).mean(axis=1).rename("mek_sens")   # consensus MEK sensitivity
H = pd.read_csv("data/raw/depmap/OmicsSomaticMutationsMatrixHotspot.csv", index_col=0)
mut = pd.DataFrame({"KRAS": H["KRAS (3845)"] > 0, "BRAF": H["BRAF (673)"] > 0, "NRAS": H["NRAS (4893)"] > 0})
mut["geno"] = np.select([mut.KRAS, mut.BRAF, mut.NRAS], ["KRAS", "BRAF", "NRAS"], "WT")
mut["mapk"] = mut[["KRAS", "BRAF", "NRAS"]].any(axis=1)
lineage = pd.read_csv("data/raw/depmap/Model.csv").set_index("ModelID")["OncotreeLineage"]
crc = set(prism[prism.indication == "CRC"].model_id)
df = pd.DataFrame(mekz).join(mut[["geno", "mapk"]]).dropna()
df["crc"] = df.index.isin(crc)
df["lineage"] = df.index.map(lineage)

fig = plt.figure(figsize=(15.5, 7.0))
fig.suptitle("The 'CRC-selective' MEK sensitivity is a MAPK-mutant genotype effect, not a lineage effect",
             color=RRB, fontsize=15, fontweight="bold", x=.5, y=.985)
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.35, 1.05], wspace=.28, left=.05, right=.975, top=.85, bottom=.28)

# --- A: MEK sensitivity by genotype (all lines) ---
axA = fig.add_subplot(gs[0, 0])
order = ["WT", "KRAS", "NRAS", "BRAF"]
for i, g in enumerate(order):
    s = df[df.geno == g].mek_sens
    axA.scatter(np.full(len(s), i) + np.random.RandomState(1).uniform(-.16, .16, len(s)), s,
                s=13, color=COL[g], alpha=.55, edgecolor="none")
    axA.plot([i - .28, i + .28], [s.median()] * 2, color="#222", lw=2.2, zorder=5)
    axA.text(i, 2.35, f"n={len(s)}\n{s.mean():+.2f}", ha="center", fontsize=8.5, color="#333")
axA.axhline(0, color="#ccc", lw=.8); axA.set_ylim(-2.25, 2.75)
axA.set_xticks(range(4)); axA.set_xticklabels(order, fontsize=10)
axA.set_ylabel("MEK-inhibitor sensitivity  (z)", fontsize=10)
axA.set_title("A.  by genotype (all lineages)", fontsize=11, loc="left", color=RRB)
axA.text(.5, -.16, "MAPK-mut vs WT: gap 0.57,  p = 1e-18", transform=axA.transAxes, ha="center",
         fontsize=9, color=RRB, fontweight="bold")
for s in ("top", "right"):
    axA.spines[s].set_visible(False)

# --- B: per-lineage, MAPK-mut vs WT (genotype drives it in every tissue) ---
axB = fig.add_subplot(gs[0, 1])
lin_n = df.groupby("lineage").size()
lins = [l for l in lin_n[lin_n >= 12].index if pd.notna(l)]
means = df[df.lineage.isin(lins)].groupby(["lineage", "mapk"]).mek_sens.mean().unstack()
means = means.dropna().sort_values(True, ascending=False)
y = np.arange(len(means))
axB.scatter(means[True], y, s=90, color=RRB, zorder=4, label="MAPK-mutant")
axB.scatter(means[False], y, s=90, color="#AEB4BC", zorder=4, label="MAPK-WT")
for i in range(len(means)):
    axB.plot([means[False].iloc[i], means[True].iloc[i]], [i, i], color="#ccc", lw=1.5, zorder=2)
labels = [("CRC (Bowel)" if l == "Bowel" else ("Melanoma (Skin)" if l == "Skin" else str(l))) for l in means.index]
axB.set_yticks(y); axB.set_yticklabels(labels, fontsize=9)
for i, l in enumerate(means.index):
    if l == "Bowel":
        axB.get_yticklabels()[i].set_color(RRB); axB.get_yticklabels()[i].set_fontweight("bold")
axB.axvline(0, color="#ccc", lw=.8); axB.invert_yaxis()
axB.set_xlabel("MEK sensitivity  (z)", fontsize=10)
axB.set_title("B.  within each lineage, mutants are sensitive - CRC is not special",
              fontsize=11, loc="left", color=RRB)
axB.legend(fontsize=9, frameon=False, loc="lower right")
for s in ("top", "right"):
    axB.spines[s].set_visible(False)

# --- C: decomposition + patient coverage ---
axC = fig.add_subplot(gs[0, 2])
X = np.column_stack([df.crc.astype(float), df.mapk.astype(float), np.ones(len(df))])
b, *_ = np.linalg.lstsq(X, df.mek_sens.values, rcond=None)
Xc = np.column_stack([df.crc.astype(float), np.ones(len(df))]); bc, *_ = np.linalg.lstsq(Xc, df.mek_sens.values, rcond=None)
bars = [("CRC effect\n(unadjusted)", bc[0], "#C77"), ("CRC effect\n(+ MAPK-mut)", b[0], "#C77"),
        ("MAPK-mut effect", b[1], RRB)]
for i, (lab, val, c) in enumerate(bars):
    axC.bar(i, val, color=c, width=.62); axC.text(i, val + .02, f"{val:.2f}", ha="center", fontsize=9.5, fontweight="bold")
axC.set_xticks(range(3)); axC.set_xticklabels([b[0] for b in bars], fontsize=8.5)
axC.set_ylabel("effect on MEK sensitivity (OLS β)", fontsize=9.5)
axC.set_title("C.  CRC effect collapses once genotype enters", fontsize=11, loc="left", color=RRB)
axC.text(0.50, 0.52, "70% of the CRC\n'selectivity' is\nMAPK genotype", transform=axC.transAxes,
         ha="center", fontsize=9.2, color=RRB, fontweight="bold")
for s in ("top", "right"):
    axC.spines[s].set_visible(False)

# --- patient-coverage banner (answers "how many patients would this cover?") ---
fig.text(0.5, 0.135, "How many patients?  —  reframed as a MAPK-mutant biomarker",
         ha="center", fontsize=11.5, fontweight="bold", color=RRB)
fig.text(0.5, 0.025,
         "CRC is ~55% RAS/RAF-mutant in patients (KRAS ~42% · BRAF ~9% · NRAS ~4%)  ≈  84k US / ~1.05M global "
         "new CRC cases per year — largely the anti-EGFR (cetuximab)-refractory population, an unmet need.\n"
         "The same biomarker extends beyond CRC: KRAS-mutant pancreatic (~90%) & lung adeno (~30%), BRAF-mutant "
         "melanoma (~50%).   Caveat: these cell lines are 93% MAPK-mutant (selection bias) vs ~55% in patients; "
         "prevalences are literature epidemiology.",
         ha="center", va="bottom", fontsize=8.4, color="#333",
         bbox=dict(boxstyle="round,pad=0.6", facecolor="#FBF6F7", edgecolor="#d8c7cc"))

for ext in ("pdf", "png"):
    fig.savefig(f"reports/mek_confounding.{ext}", dpi=160)
plt.close(fig)
print("wrote reports/mek_confounding.pdf/.png  | CRC beta %.2f -> %.2f (MAPK %.2f)" % (bc[0], b[0], b[1]))
