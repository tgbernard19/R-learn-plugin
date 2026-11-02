#!/usr/bin/env python3
"""
How different are the N-responsive root gene sets between time points after N
resupply? Varala et al. 2018 PNAS (GEO GSE97500), Arabidopsis roots,
N-starved seedlings resupplied with N (KNO3/NH4NO3) or a KCl mock.

Run from anywhere (python -I keeps the data directory off the import path):

    bash ../scripts/fetch_data.sh                 # downloads into ../data
    python3 -I varala_timepoint_redundancy.py     # writes figures/CSVs next to this file

Options: --data DIR, --out DIR, --cache DIR, --recompute (ignore cached DE fits).

Sections
  1. Load counts, read the design from the column names, cross-check every
     column against the per-sample GEO count files, QC (PCA).
  2. Spline sanity check (time x treatment model, moderated F-test), compared
     with the paper's ~2,174 dynamic N-responsive root genes.
  3. Per-time-point N vs KCl DE with pydeseq2 (one model, ~group).
  4. Overlap analysis: counts, cumulative union, Jaccard, greedy redundancy;
     all 9 points and the 4 points matching the planned Avena 15/30/60/120 min.
  5. GO biological-process enrichment for time-specific gene sets.
"""
import argparse
import glob
import json
import os
import re
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats, special
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*keyword-only.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--data", default=os.path.join(HERE, "..", "data"))
ap.add_argument("--out", default=HERE)
ap.add_argument("--cache", default=os.path.join(HERE, "..", "cache"))
ap.add_argument("--recompute", action="store_true")
ap.add_argument("--cpus", type=int, default=4)
args = ap.parse_args()
DATA, OUT, CACHE = (os.path.abspath(p) for p in (args.data, args.out, args.cache))
os.makedirs(OUT, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

ALPHA = 0.05                      # main padj cutoff
LFC_LEVELS = [0.0, 0.5, 1.0]      # |log2FC| > x on top of padj < 0.05; 0 = main
PLANNED = [15, 30, 60, 120]       # planned Avena points that fall inside 0-120 min
RESULTS = {}                      # numbers quoted in summary.md

# ---- plotting style (reference palette from the dataviz skill) --------------
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
UP, DOWN = "#e34948", "#2a78d6"           # diverging poles: N>KCl red, N<KCl blue
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]  # categorical slots 1-3
SEQ = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "savefig.dpi": 160,
})


def tlabel(t):
    return f"{t} min" if t < 60 else f"{t // 60} h"


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight")
    plt.close(fig)
    print("  wrote", name)


# =============================================================================
# 1. Data, design, provenance check, QC
# =============================================================================
print("== 1. Load data and read the design")
counts = pd.read_csv(os.path.join(DATA, "RootRawCounts_mirror.txt"), sep="\t", index_col=0)
counts = counts.astype(int)

pat = re.compile(r"^(T|C|N)(\d+)-(\d+)$")
rows = []
for c in counts.columns:
    m = pat.match(c)
    if not m:
        raise SystemExit(f"Unexpected column name {c!r}: design does not match T/C/N<time>-<rep>")
    trt = {"T": "T0", "C": "KCl", "N": "N"}[m.group(1)]
    rows.append(dict(sample=c, treatment=trt, time_min=int(m.group(2)), rep=int(m.group(3))))
design = pd.DataFrame(rows).set_index("sample")
design["group"] = [("T0" if r.treatment == "T0" else f"{'C' if r.treatment == 'KCl' else 'N'}{r.time_min}")
                   for r in design.itertuples()]
design["lib_size_M"] = (counts.sum() / 1e6).round(2)

# Cross-check each column against the per-sample GEO files (exact equality).
geo_files = sorted(glob.glob(os.path.join(DATA, "geo_sample_counts", "GSM*.counts.txt.gz")))
gsm = {}
if geo_files:
    geo = {}
    for f in geo_files:
        s = pd.read_csv(f, sep="\t", header=None, index_col=0).iloc[:, 0]
        geo[os.path.basename(f).split("_")[0]] = s[~s.index.str.startswith("__")]
    geo = pd.DataFrame(geo).reindex(counts.index)
    for c in counts.columns:
        hit = [k for k in geo.columns if np.array_equal(geo[k].values, counts[c].values)]
        gsm[c] = hit[0] if len(hit) == 1 else ("AMBIGUOUS" if hit else "NO_MATCH")
    design["GEO_sample"] = pd.Series(gsm)
    n_ok = int(sum(v.startswith("GSM") for v in gsm.values()))
    print(f"  {n_ok}/{counts.shape[1]} columns identical to a GEO GSE97500 per-sample count file")
    RESULTS["geo_columns_matched"] = n_ok
design.to_csv(os.path.join(OUT, "design_samples.csv"))

cell = design.groupby(["treatment", "time_min"]).size().unstack(0).fillna(0).astype(int)
print("  replicates per treatment x time cell:\n" + cell.to_string())
TIMES = sorted(design.loc[design.treatment != "T0", "time_min"].unique())
RESULTS["times"] = [int(t) for t in TIMES]
RESULTS["cells"] = cell.to_dict()
RESULTS["n_genes_in_table"] = int(counts.shape[0])

# Expression filter: >= 10 reads in >= 3 samples (smallest group size).
keep = (counts >= 10).sum(axis=1) >= 3
cts = counts.loc[keep]
RESULTS["n_genes_filtered"] = int(keep.sum())
print(f"  {keep.sum()} of {len(keep)} genes pass >=10 reads in >=3 samples")


def size_factors(x):
    lx = np.log(x.loc[(x > 0).all(axis=1)])
    return np.exp((lx.sub(lx.mean(axis=1), axis=0)).median(axis=0))


sf = size_factors(cts)
logn = np.log2(cts / sf + 1)

# PCA on the 2,000 most variable genes.
top = logn.loc[logn.var(axis=1).sort_values().index[-2000:]]
X = top.sub(top.mean(axis=1), axis=0).values
u, s, vt = np.linalg.svd(X, full_matrices=False)
ve = s ** 2 / (s ** 2).sum()
pcs = pd.DataFrame({"PC1": vt[0] * s[0], "PC2": vt[1] * s[1]}, index=top.columns)
RESULTS["pca_var_explained"] = [round(float(v), 3) for v in ve[:3]]
odd_times = [10, 15]
odd = design.index[design.time_min.isin(odd_times) & (design.treatment != "T0")]
rest = design.index.difference(odd)
RESULTS["pc1_mean_10_15"] = round(float(pcs.loc[odd, "PC1"].mean()), 1)
RESULTS["pc1_mean_other"] = round(float(pcs.loc[rest, "PC1"].mean()), 1)
cpm = cts / cts.sum() * 1e6
expr = cpm[(cpm > 5).sum(axis=1) > 10]
shift = np.log2(expr[odd].mean(axis=1) + 1) - np.log2(expr[rest].mean(axis=1) + 1)
RESULTS["n_genes_2fold_shift_10_15"] = int((shift.abs() > 1).sum())
RESULTS["n_genes_expressed_for_shift"] = int(len(shift))
chl = counts.index.str.startswith("ATCG")
RESULTS["chloroplast_read_frac_10_15"] = float(counts.loc[chl, odd].sum().sum() / counts[odd].sum().sum())
RESULTS["chloroplast_read_frac_other"] = float(counts.loc[chl, rest].sum().sum() / counts[rest].sum().sum())

fig, ax = plt.subplots(figsize=(6.4, 4.6))
for trt, mk in [("T0", "s"), ("KCl", "o"), ("N", "^")]:
    idx = design.index[design.treatment == trt]
    col = np.where(design.loc[idx, "time_min"].isin(odd_times), SERIES[1], SERIES[0])
    ax.scatter(pcs.loc[idx, "PC1"], pcs.loc[idx, "PC2"], c=col, marker=mk, s=46,
               edgecolor=SURFACE, linewidth=1.2, label=trt)
for g in ["C10", "N10", "C15", "N15"]:
    idx = design.index[design.group == g]
    ax.annotate(g.replace("C", "KCl ").replace("N", "N "), pcs.loc[idx].mean(), fontsize=8,
                color=INK2, xytext=(6, 4), textcoords="offset points")
ax.set_xlabel(f"PC1 ({ve[0]:.0%} of variance)")
ax.set_ylabel(f"PC2 ({ve[1]:.0%})")
ax.set_title("Root samples: 10 and 15 min (orange) separate on PC1")
h = [plt.Line2D([], [], marker=m, ls="", color=INK2, label=l) for m, l in
     [("s", "T0"), ("o", "KCl"), ("^", "N")]]
h += [plt.Line2D([], [], marker="o", ls="", color=SERIES[1], label="10 / 15 min"),
      plt.Line2D([], [], marker="o", ls="", color=SERIES[0], label="other times")]
ax.legend(handles=h, loc="upper center", fontsize=8, ncol=5, bbox_to_anchor=(0.5, -0.14))
save(fig, "qc_pca.png")


# =============================================================================
# 2. Spline sanity check: limma-trend style moderated F-test
# =============================================================================
print("== 2. Spline time x treatment model (sanity check vs ~2,174 genes)")


def bspline_basis(x, df, knots_from):
    """Cubic B-spline basis without intercept (all columns 0 at min(x))."""
    from patsy import dmatrix
    inner = np.quantile(np.unique(knots_from), np.linspace(0, 1, df - 3 + 2)[1:-1]) if df > 3 else []
    return np.asarray(dmatrix("bs(x, knots=k, degree=3, include_intercept=False, "
                              "lower_bound=lo, upper_bound=hi) - 1",
                              {"x": x, "k": list(inner), "lo": 0, "hi": max(knots_from)}))


def trigamma_inverse(x):
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.max(-dif / y) < 1e-8:
            break
    return y


def moderated_F(E, Xd, test_cols, amean):
    """limma lmFit + eBayes(trend=TRUE) F-test for a block of coefficients."""
    XtX_inv = np.linalg.pinv(Xd.T @ Xd)
    B = E @ Xd @ XtX_inv
    resid = E - B @ Xd.T
    d = Xd.shape[0] - np.linalg.matrix_rank(Xd)
    s2 = (resid ** 2).sum(axis=1) / d
    # fitFDist with a lowess trend on mean expression (Smyth 2004)
    z = np.log(np.maximum(s2, 1e-12))
    e = z - special.digamma(d / 2) + np.log(d / 2)
    emean = lowess(e, amean, frac=0.3, return_sorted=False)
    evar = np.mean((e - emean) ** 2) * len(e) / (len(e) - 1) - special.polygamma(1, d / 2)
    d0 = 2 * trigamma_inverse(evar) if evar > 0 else np.inf
    s02 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2)) if np.isfinite(d0) else np.exp(emean)
    s2p = (d0 * s02 + d * s2) / (d0 + d) if np.isfinite(d0) else s02
    C = np.zeros((Xd.shape[1], len(test_cols)))
    C[test_cols, np.arange(len(test_cols))] = 1
    V = np.linalg.inv(C.T @ XtX_inv @ C)
    b = B @ C
    F = np.einsum("ij,jk,ik->i", b, V, b) / len(test_cols) / s2p
    df2 = d + d0 if np.isfinite(d0) else 1e6
    p = stats.f.sf(F, len(test_cols), df2)
    return p, d0


def spline_test(samples, df, time_scale):
    lib = cts[samples].sum() * sf[samples] / sf[samples].mean()
    E = np.log2((cts[samples] + 0.5) / (lib + 1) * 1e6).values           # log2-CPM
    amean = E.mean(axis=1)
    t = design.loc[samples, "time_min"].values.astype(float)
    tt = np.log2(t + 1) if time_scale == "log" else t
    basis = bspline_basis(tt, df, tt)
    isN = (design.loc[samples, "treatment"] == "N").values.astype(float)[:, None]
    # T0 is shared: the N-KCl difference is constrained to 0 at t = 0 because the
    # basis has no intercept. Tested block = N x spline(time), df coefficients.
    Xd = np.column_stack([np.ones(len(t)), basis, isN * basis])
    p, d0 = moderated_F(E, Xd, list(range(1 + df, 1 + 2 * df)), amean)
    return pd.Series(p, index=cts.index), d0


spline_rows = []
spline_sets = {}
all_s = list(design.index)
no_odd = list(rest)
for subset_name, samples in [("all 57 samples", all_s), ("without 10/15 min", no_odd)]:
    for scale in ["linear", "log"]:
        for df in [3, 4, 5]:
            p, d0 = spline_test(samples, df, scale)
            q = multipletests(p, method="fdr_bh")[1]
            for fdr in [0.01, 0.05]:
                n = int((q < fdr).sum())
                spline_rows.append(dict(samples=subset_name, time_scale=scale, spline_df=df,
                                        FDR=fdr, n_genes=n, prior_df=round(float(d0), 1)))
            spline_sets[(subset_name, scale, df)] = pd.Series(q, index=cts.index)
spline_tab = pd.DataFrame(spline_rows)
spline_tab.to_csv(os.path.join(OUT, "spline_sanity_check_counts.csv"), index=False)
print(spline_tab.pivot_table(index=["samples", "time_scale", "spline_df"], columns="FDR",
                             values="n_genes").to_string())
RESULTS["spline_grid"] = spline_rows


# =============================================================================
# 3. Per-time-point DE: N vs KCl at each time, pydeseq2, one model ~group
# =============================================================================
print("== 3. Per-time-point DE with pydeseq2")


def run_deseq(samples, tag):
    path = os.path.join(CACHE, f"deseq_{tag}.csv")
    if os.path.exists(path) and not args.recompute:
        print(f"  using cached {os.path.basename(path)}")
        return pd.read_csv(path)
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    inf = DefaultInference(n_cpus=args.cpus)
    meta = design.loc[samples, ["group"]].copy()
    dds = DeseqDataSet(counts=cts[samples].T, metadata=meta, design="~group",
                       refit_cooks=True, inference=inf, quiet=True)
    dds.deseq2()
    out = []
    times = sorted(design.loc[samples].query("treatment != 'T0'").time_min.unique())
    for t in times:
        st = DeseqStats(dds, contrast=["group", f"N{t}", f"C{t}"], alpha=ALPHA,
                        inference=inf, quiet=True)
        st.summary()
        r = st.results_df[["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]].copy()
        r.insert(0, "time_min", t)
        out.append(r.rename_axis("gene").reset_index())
    res = pd.concat(out, ignore_index=True)
    res.to_csv(path, index=False)
    return res


de_main = run_deseq(all_s, "all_samples")
de_sens = run_deseq(no_odd, "without_10_15")


def de_sets(res, lfc, times=None):
    sig = res[(res.padj < ALPHA) & (res.log2FoldChange.abs() > lfc)]
    times = times if times is not None else sorted(res.time_min.unique())
    return {t: set(sig.loc[sig.time_min == t, "gene"]) for t in times}


# Per-gene table of calls, written as the main DE CSV.
de_out = de_main.copy()
de_out["direction"] = np.where(de_out.log2FoldChange > 0, "up_in_N", "down_in_N")
for lfc in LFC_LEVELS:
    de_out[f"sig_padj05_lfc{lfc:g}"] = (de_out.padj < ALPHA) & (de_out.log2FoldChange.abs() > lfc)
de_out.loc[de_out[f"sig_padj05_lfc0"]].sort_values(["time_min", "padj"]).to_csv(
    os.path.join(OUT, "de_genes_per_timepoint.csv"), index=False, float_format="%.4g")
de_out.to_csv(os.path.join(CACHE, "de_all_tests_all_samples.csv"), index=False)

count_rows = []
for name, res in [("all samples", de_main), ("without 10/15 min samples", de_sens)]:
    for lfc in LFC_LEVELS:
        sig = res[(res.padj < ALPHA) & (res.log2FoldChange.abs() > lfc)]
        for t in sorted(res.time_min.unique()):
            s_ = sig[sig.time_min == t]
            count_rows.append(dict(model=name, lfc_threshold=lfc, time_min=t,
                                   up=int((s_.log2FoldChange > 0).sum()),
                                   down=int((s_.log2FoldChange < 0).sum()), total=len(s_)))
count_tab = pd.DataFrame(count_rows)
count_tab.to_csv(os.path.join(OUT, "de_counts_per_timepoint.csv"), index=False)
print(count_tab.pivot_table(index="time_min", columns=["model", "lfc_threshold"],
                            values="total").to_string())
RESULTS["de_counts"] = count_rows


# =============================================================================
# 4. Overlap analysis
# =============================================================================
print("== 4. Overlap analysis")


def cumulative(sets, order):
    seen, rows = set(), []
    for i, t in enumerate(order, 1):
        new = sets[t] - seen
        seen |= sets[t]
        rows.append(dict(n_points=i, time_min=t, set_size=len(sets[t]), new_genes=len(new),
                         cumulative_union=len(seen)))
    df = pd.DataFrame(rows)
    df["frac_of_final_union"] = df.cumulative_union / max(len(seen), 1)
    return df


def jaccard(sets, order):
    J = pd.DataFrame(np.nan, index=order, columns=order)
    for a in order:
        for b in order:
            u_ = sets[a] | sets[b]
            J.loc[a, b] = len(sets[a] & sets[b]) / len(u_) if u_ else np.nan
    return J


def greedy(sets, order):
    """Repeatedly drop the point whose removal loses the fewest genes not
    detected at any other remaining point. Ties -> smaller set dropped first."""
    remaining = list(order)
    total = len(set().union(*sets.values()))
    rows = []
    while len(remaining) > 1:
        loss = {}
        for t in remaining:
            others = set().union(*(sets[o] for o in remaining if o != t))
            loss[t] = len(sets[t] - others)
        t_drop = min(remaining, key=lambda t: (loss[t], len(sets[t]), t))
        remaining.remove(t_drop)
        kept = len(set().union(*(sets[o] for o in remaining)))
        rows.append(dict(step=len(order) - len(remaining), dropped_time_min=t_drop,
                         unique_genes_lost=loss[t_drop],
                         union_remaining=kept, frac_of_full_union=kept / max(total, 1),
                         points_remaining=",".join(map(str, remaining))))
    return pd.DataFrame(rows)


over_rows, greedy_tabs, cum_tabs, jac = [], {}, {}, {}
for scope, times in [("all", TIMES), ("planned", PLANNED)]:
    for lfc in LFC_LEVELS:
        for model, res in [("all samples", de_main), ("without 10/15", de_sens)]:
            tt = [t for t in times if t in set(res.time_min)]
            S = de_sets(res, lfc, tt)
            cum_tabs[(scope, lfc, model)] = cumulative(S, tt)
            greedy_tabs[(scope, lfc, model)] = greedy(S, tt)
            jac[(scope, lfc, model)] = jaccard(S, tt)

jac[("all", 0.0, "all samples")].round(3).rename_axis("time_min").to_csv(
    os.path.join(OUT, "jaccard_matrix_all_timepoints.csv"))
jac[("planned", 0.0, "all samples")].round(3).rename_axis("time_min").to_csv(
    os.path.join(OUT, "jaccard_matrix_planned_timepoints.csv"))

frames = []
for (scope, lfc, model), d in cum_tabs.items():
    frames.append(d.assign(scope=scope, lfc_threshold=lfc, model=model))
pd.concat(frames).to_csv(os.path.join(OUT, "cumulative_union.csv"), index=False, float_format="%.4f")
frames = []
for (scope, lfc, model), d in greedy_tabs.items():
    frames.append(d.assign(scope=scope, lfc_threshold=lfc, model=model))
pd.concat(frames).to_csv(os.path.join(OUT, "greedy_redundancy.csv"), index=False, float_format="%.4f")
RESULTS["cumulative"] = {f"{k[0]}|{k[1]}|{k[2]}": v.to_dict("records") for k, v in cum_tabs.items()}
RESULTS["greedy"] = {f"{k[0]}|{k[1]}|{k[2]}": v.to_dict("records") for k, v in greedy_tabs.items()}

# ---- 4a. bar chart of DE genes per time point, up / down --------------------
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=False)
for ax, lfc in zip(axes, LFC_LEVELS):
    d = count_tab[(count_tab.model == "all samples") & (count_tab.lfc_threshold == lfc)]
    x = np.arange(len(d))
    ax.bar(x - 0.21, d.up, 0.4, color=UP, label="up in N vs KCl")
    ax.bar(x + 0.21, d.down, 0.4, color=DOWN, label="down in N vs KCl")
    for xi, (u_, dn) in enumerate(zip(d.up, d.down)):
        ax.text(xi, max(u_, dn) * 1.02 + 1, f"{u_ + dn}", ha="center", va="bottom", fontsize=7.5, color=INK2)
    ax.set_xticks(x, [tlabel(t) for t in d.time_min], rotation=45, ha="right")
    for lab in ax.get_xticklabels():
        if lab.get_text() in ("10 min", "15 min"):
            lab.set_color(SERIES[1])
    ax.set_title("padj < 0.05" + (f", |log2FC| > {lfc:g}" if lfc else ""))
    ax.grid(axis="x", visible=False)
axes[0].set_ylabel("DE genes (N vs KCl)")
axes[0].legend(loc="upper left", fontsize=8)
fig.suptitle("4a  DE genes per time point (number above bars = total; orange labels = 10/15 min batch)",
             x=0.01, ha="left", fontsize=11, fontweight="bold")
fig.tight_layout()
save(fig, "fig4a_de_counts_per_timepoint.png")


# ---- 4b / 4e. cumulative union ----------------------------------------------
def plot_cumulative(scope, fname, title):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.3))
    d = cum_tabs[(scope, 0.0, "all samples")]
    x = np.arange(len(d))
    a1.bar(x, d.new_genes, 0.6, color=SEQ[2], label="new genes added at this point")
    a1.plot(x, d.cumulative_union, color=SERIES[0], lw=2, marker="o", ms=6,
            markeredgecolor=SURFACE, mew=1.5, label="cumulative distinct DE genes")
    for xi, (c, n) in enumerate(zip(d.cumulative_union, d.new_genes)):
        a1.annotate(f"+{n}", (xi, n), xytext=(0, 3), textcoords="offset points", ha="center",
                    fontsize=7.5, color=INK2)
    a1.set_xticks(x, [f"{i}\n{tlabel(t)}" for i, t in zip(d.n_points, d.time_min)], fontsize=8)
    a1.set_xlabel("time points included (chronological order)")
    a1.set_ylabel("genes (padj < 0.05)")
    a1.set_title("Main call: padj < 0.05")
    a1.grid(axis="x", visible=False)
    a1.legend(loc="upper left", fontsize=8)
    for i, lfc in enumerate(LFC_LEVELS):
        d = cum_tabs[(scope, lfc, "all samples")]
        a2.plot(np.arange(len(d)), d.frac_of_final_union, color=SERIES[i], lw=2, marker="o",
                ms=6, markeredgecolor=SURFACE, mew=1.5,
                label=f"|log2FC| > {lfc:g} (n = {int(d.cumulative_union.iloc[-1])})" if lfc else
                f"padj only (n = {int(d.cumulative_union.iloc[-1])})")
    if scope == "all":
        d = cum_tabs[(scope, 0.0, "without 10/15")]
        a2.plot([TIMES.index(t) for t in d.time_min], d.frac_of_final_union, color=INK2, lw=1.4,
                ls="--", marker="o", ms=4,
                label=f"padj only, 10/15 min removed (n = {int(d.cumulative_union.iloc[-1])})")
    a2.set_xticks(x, [tlabel(t) for t in cum_tabs[(scope, 0.0, 'all samples')].time_min], fontsize=8)
    a2.set_ylim(0, 1.04)
    a2.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    a2.set_ylabel("share of final union")
    a2.set_title("Sensitivity: share of final union reached")
    a2.legend(loc="lower right", fontsize=8)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, fname)


plot_cumulative("all", "fig4b_cumulative_union_all.png",
                "4b  Cumulative union of DE genes, adding time points chronologically (all 9 points)")
plot_cumulative("planned", "fig4e_cumulative_union_planned.png",
                "4e  Cumulative union, planned-design points only (15 min, 30 min, 1 h, 2 h)")


# ---- 4c / 4e. Jaccard heatmap -------------------------------------------------
def plot_jaccard(scope, fname, title):
    J = jac[(scope, 0.0, "all samples")]
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    n = len(J)
    fig, ax = plt.subplots(figsize=(1.0 + 0.62 * n, 0.6 + 0.55 * n))
    M = J.values.copy()
    np.fill_diagonal(M, np.nan)
    vmax = max(0.05, np.nanmax(M))
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=vmax)
    for i in range(n):
        for j in range(n):
            if i == j:
                ax.text(j, i, f"{len(de_sets(de_main, 0.0, [J.index[i]])[J.index[i]])}",
                        ha="center", va="center", fontsize=7, color=INK2, style="italic")
                continue
            v = M[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7.5,
                    color=SURFACE if v > 0.55 * vmax else INK)
    ax.set_xticks(range(n), [tlabel(t) for t in J.columns], rotation=45, ha="right")
    ax.set_yticks(range(n), [tlabel(t) for t in J.index])
    ax.grid(False)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Jaccard index")
    cb.outline.set_visible(False)
    ax.set_title(title, loc="left")
    fig.tight_layout()
    save(fig, fname)


plot_jaccard("all", "fig4c_jaccard_all.png", "4c  Jaccard similarity of DE sets\n(padj < 0.05; diagonal = set size)")
plot_jaccard("planned", "fig4e_jaccard_planned.png", "4e  Jaccard, planned points\n(padj < 0.05; diagonal = set size)")


# ---- 4d / 4e. greedy redundancy ---------------------------------------------
def plot_greedy(fname, title):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [9, 4]})
    for ax, scope in zip(axes, ["all", "planned"]):
        for i, lfc in enumerate(LFC_LEVELS):
            g = greedy_tabs[(scope, lfc, "all samples")]
            full = len(set().union(*de_sets(de_main, lfc, TIMES if scope == "all" else PLANNED).values()))
            y = [1.0] + list(g.frac_of_full_union)
            ax.plot(range(len(y)), y, color=SERIES[i], lw=2, marker="o", ms=6,
                    markeredgecolor=SURFACE, mew=1.5,
                    label=("padj only" if not lfc else f"|log2FC| > {lfc:g}") + f" (n = {full})")
            if lfc == 0.0:
                for k, r in enumerate(g.itertuples(), 1):
                    ax.annotate(f"−{tlabel(r.dropped_time_min)}\n({r.unique_genes_lost})",
                                (k, r.frac_of_full_union), xytext=(0, -26), textcoords="offset points",
                                ha="center", fontsize=7, color=INK2)
        ax.set_ylim(0, 1.05)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
        ax.set_xlabel("time points dropped")
        ax.set_title("All 9 points" if scope == "all" else "Planned points (15, 30, 60, 120 min)")
    axes[0].set_ylabel("share of DE union still detected")
    axes[0].legend(loc="lower left", fontsize=8)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, fname)


plot_greedy("fig4d_greedy_redundancy.png",
            "4d/4e  Greedy removal of the least-informative time point "
            "(labels: point dropped, genes lost, padj-only call)")

for key in [("all", 0.0, "all samples"), ("planned", 0.0, "all samples")]:
    print(f"  greedy {key}:")
    print(greedy_tabs[key][["step", "dropped_time_min", "unique_genes_lost", "union_remaining"]]
          .to_string(index=False))


# =============================================================================
# 5. GO biological process enrichment
# =============================================================================
print("== 5. GO BP enrichment")


def read_obo(path):
    parents, name, ns, obsolete = {}, {}, {}, set()
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = {"is_a": []}
                continue
            if line.startswith("[") and line != "[Term]":
                cur = None
                continue
            if cur is None or not line:
                if cur is not None and "id" in cur:
                    parents[cur["id"]] = cur["is_a"]
                    cur = None
                continue
            k, _, v = line.partition(": ")
            if k == "id":
                cur["id"] = v
            elif k == "name":
                name[cur["id"]] = v
            elif k == "namespace":
                ns[cur["id"]] = v
            elif k == "is_a":
                cur["is_a"].append(v.split(" ! ")[0])
            elif k == "relationship" and v.startswith("part_of "):
                cur["is_a"].append(v.split()[1])
            elif k == "is_obsolete" and v == "true":
                obsolete.add(cur["id"])
    return parents, name, ns, obsolete


parents, goname, gons, obsolete = read_obo(os.path.join(DATA, "go-basic.obo"))
anc_cache = {}


def ancestors(t):
    if t not in anc_cache:
        out = {t}
        for p in parents.get(t, []):
            out |= ancestors(p)
        anc_cache[t] = out
    return anc_cache[t]


slim = pd.read_csv(os.path.join(DATA, "ATH_GO_GOSLIM.txt.gz"), sep="\t", header=None, comment="!",
                   usecols=[0, 3, 5, 7, 9], names=["gene", "relation", "go", "aspect", "evidence"],
                   dtype=str)
slim = slim[(slim.aspect == "P") & ~slim.relation.str.startswith("NOT", na=False)]
slim = slim[slim.gene.str.match(r"^AT[1-5CM]G\d{5}$", na=False)]
gene2go = {}
for g, t in zip(slim.gene, slim.go):
    if t in obsolete or t not in parents:
        continue
    gene2go.setdefault(g, set()).update(a for a in ancestors(t) if gons.get(a) == "biological_process")

universe = set(cts.index) & set(gene2go)
RESULTS["go_universe"] = len(universe)
term2genes = {}
for g in universe:
    for t in gene2go[g]:
        term2genes.setdefault(t, set()).add(g)
term2genes = {t: s for t, s in term2genes.items() if 10 <= len(s) <= 500}


def enrich(genes, label):
    genes = set(genes) & universe
    N, n = len(universe), len(genes)
    rows = []
    for t, members in term2genes.items():
        k = len(members & genes)
        if k < 3:
            continue
        p = stats.hypergeom.sf(k - 1, N, len(members), n)
        rows.append(dict(set=label, GO=t, term=goname.get(t, t), set_size=n, term_size=len(members),
                         overlap=k, fold=(k / n) / (len(members) / N), pvalue=p,
                         genes=";".join(sorted(members & genes)[:40])))
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # BH over all terms tested (terms with < 3 hits have p ~ 1 and are counted)
    p_all = np.concatenate([df.pvalue.values, np.ones(len(term2genes) - len(df))])
    df["padj"] = multipletests(p_all, method="fdr_bh")[1][:len(df)]
    return df.sort_values("pvalue")


S_all = de_sets(de_main, 0.0, TIMES)
union_all = set().union(*S_all.values())
first_time = {}
for t in TIMES:
    for g in S_all[t]:
        first_time.setdefault(g, t)
clusters = {
    "early (first DE at 5-20 min)": {g for g, t in first_time.items() if t <= 20},
    "middle (first DE at 30-45 min)": {g for g, t in first_time.items() if 30 <= t <= 45},
    "late (first DE at 60-120 min)": {g for g, t in first_time.items() if t >= 60},
}
unique_sets = {f"unique to {tlabel(t)}": S_all[t] - set().union(*(S_all[o] for o in TIMES if o != t))
               for t in TIMES}
go_frames = []
for label, genes in list(clusters.items()) + list(unique_sets.items()) + [("all DE genes (union)", union_all)]:
    e = enrich(genes, label)
    if not e.empty:
        go_frames.append(e)
go = pd.concat(go_frames, ignore_index=True)
go.to_csv(os.path.join(OUT, "go_bp_enrichment.csv"), index=False, float_format="%.4g")
RESULTS["cluster_sizes"] = {k: len(v) for k, v in clusters.items()}
RESULTS["unique_sizes"] = {k: len(v) for k, v in unique_sets.items()}
top_go = {}
for label in go.set.unique():
    sub = go[(go.set == label) & (go.padj < 0.05)].head(12)
    top_go[label] = sub[["term", "overlap", "term_size", "fold", "padj"]].round(4).to_dict("records")
    print(f"  {label}: {int((go[go.set == label].padj < 0.05).sum())} BP terms at FDR < 0.05")
    for r in sub.head(6).itertuples():
        print(f"      {r.term} (k={r.overlap}, fold={r.fold:.1f}, padj={r.padj:.2g})")
RESULTS["top_go"] = top_go

# Theme check: expected early uptake/reduction, later amino acid synthesis.
themes = {
    "nitrate transport": ["nitrate transport", "nitrate import"],
    "nitrate assimilation / N reduction": ["nitrate assimilation", "nitrate metabolic process",
                                           "nitrogen compound metabolic process"],
    "amino acid biosynthesis": ["amino acid biosynthetic process", "alpha-amino acid biosynthetic process",
                                "glutamine family amino acid biosynthetic process",
                                "cellular amino acid biosynthetic process"],
    "response to nitrate": ["response to nitrate", "cellular response to nitrate"],
    "pentose-phosphate / glycolysis": ["pentose-phosphate shunt", "glycolytic process",
                                       "NADP metabolic process"],
    "translation / ribosome biogenesis": ["ribosome biogenesis", "translation", "rRNA processing"],
}
theme_rows = []
for label in list(clusters) + list(unique_sets):
    sub = go[go.set == label].set_index("term")
    for theme, terms in themes.items():
        hits = sub.loc[sub.index.intersection(terms)]
        if hits.empty:
            theme_rows.append(dict(set=label, theme=theme, best_term="", overlap=0, fold=np.nan, padj=np.nan))
            continue
        best = hits.sort_values("pvalue").iloc[0]
        theme_rows.append(dict(set=label, theme=theme, best_term=best.name, overlap=int(best.overlap),
                               fold=round(float(best.fold), 2), padj=float(best.padj)))
theme_tab = pd.DataFrame(theme_rows)
theme_tab.to_csv(os.path.join(OUT, "go_theme_check.csv"), index=False, float_format="%.3g")
print(theme_tab[theme_tab.set.str.startswith(("early", "middle", "late"))].to_string(index=False))

# Marker genes: when are they first called?
markers = {"AT1G08090": "NRT2.1", "AT1G12110": "NRT1.1/NPF6.3", "AT1G77760": "NIA1",
           "AT1G37130": "NIA2", "AT2G15620": "NIR1", "AT1G68670": "HHO2",
           "AT5G60770": "NRT2.4", "AT1G13300": "HRS1", "AT3G16857": "ARR1",
           "AT1G66340": "ETR1", "AT5G35630": "GLN2", "AT5G04140": "GLU1/Fd-GOGAT",
           "AT5G53460": "GLT1/NADH-GOGAT", "AT3G47340": "ASN1", "AT5G65010": "ASN2",
           "AT1G24280": "G6PD3", "AT4G37540": "LBD39", "AT1G67110": "CYP735A2",
           "AT5G19040": "IPT3", "AT4G05320": "UBQ10 (control)"}
mk = de_main[de_main.gene.isin(markers)].copy()
mk["symbol"] = mk.gene.map(markers)
mk_tab = mk.pivot_table(index=["gene", "symbol"], columns="time_min", values="log2FoldChange").round(2)
mk_sig = mk.assign(sig=mk.padj < ALPHA).pivot_table(index=["gene", "symbol"], columns="time_min",
                                                   values="sig", aggfunc="first")
mk_out = mk_tab.astype(str) + np.where(mk_sig.fillna(False).astype(bool), "*", "")
mk_out.to_csv(os.path.join(OUT, "marker_genes_log2fc.csv"))
print("  marker log2FC (N vs KCl; * = padj < 0.05):\n" + mk_out.to_string())

with open(os.path.join(OUT, "results.json"), "w") as fh:
    json.dump(RESULTS, fh, indent=1, default=lambda o: o if not isinstance(o, (np.integer, np.floating))
              else o.item())
print("done")
