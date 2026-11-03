#!/usr/bin/env python3
"""
Step 6: Krapp et al. 2011 Plant Physiol 157:1255 (long-term N starvation) vs
Varala et al. 2018 early N resupply, Arabidopsis roots.

Input: the paper's supplementary zip (plphys_v157_3_1255_s1.zip) extracted to
../data/krapp_supp/x/ ; Supplemental Table VII = genome-wide root CATMA
two-colour results, comparisons day2 vs day0, day10 vs day2, day10 vs day0 of
N starvation, 2 biological replicates with dye swap, Bonferroni p-values.
NB: comparisons are against day 0 of the same culture, NOT against an
N-replete control at the same age, so starvation is confounded with 10 days of
growth.  Varala per-time-point results come from the cache written by
varala_timepoint_redundancy.py.

    python3 -I krapp_comparison.py
"""
import glob
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "krapp_supp", "x")
CACHE = os.path.join(HERE, "..", "cache")
ALPHA = 0.05

f = glob.glob(os.path.join(DATA, "*Suppl_Table_VIIfinal.xls"))[0]
raw = pd.read_excel(f, sheet_name=0, header=None, skiprows=24)
cols = ["probe", "catma", "agi", "desc", "_a",
        "i1_02", "i2_02", "lfc_02", "bonf_02", "by_02", "_b",
        "i1_210", "i2_210", "lfc_210", "bonf_210", "by_210", "_c",
        "i1_010", "i2_010", "lfc_010", "bonf_010", "by_010"]
k = raw.iloc[:, :len(cols)].copy()
k.columns = cols
k = k[k.agi.astype(str).str.match(r"^AT[1-5CM]G\d{5}$")].copy()
for c in ["lfc_02", "lfc_210", "lfc_010", "bonf_02", "bonf_210", "bonf_010"]:
    k[c] = pd.to_numeric(k[c], errors="coerce")
# One probe per gene: keep the probe with the smallest 0-10 d p-value.
k = k.sort_values("bonf_010").drop_duplicates("agi").set_index("agi")
print(f"Krapp roots: {len(k)} genes with a CATMA probe")

sig = {c: k[f"bonf_{c}"] < ALPHA for c in ["02", "210", "010"]}
sgn = {c: np.sign(k[f"lfc_{c}"]) for c in ["02", "210", "010"]}


def archetype(g):
    s02, s210, s010 = sig["02"][g], sig["210"][g], sig["010"][g]
    d02, d010 = sgn["02"][g], sgn["010"][g]
    if not (s02 or s210 or s010):
        return "not responsive"
    if s02 and s010 and d02 == d010 and not s210:
        return "early, sustained (by 2 d, still at 10 d)"
    if s02 and s010 and d02 == d010 and s210:
        return "early, still growing 2-10 d"
    if s02 and not s010:
        return "early, transient (by 2 d, gone by 10 d)"
    if not s02 and (s010 or s210):
        return "late (only after 2 d)"
    return "other / mixed"


k["archetype"] = [archetype(g) for g in k.index]
k["direction_0_10d"] = np.where(k.lfc_010 > 0, "up in starvation", "down in starvation")
k.loc[k.archetype == "early, transient (by 2 d, gone by 10 d)", "direction_0_10d"] = np.where(
    k.loc[k.archetype == "early, transient (by 2 d, gone by 10 d)", "lfc_02"] > 0,
    "up in starvation", "down in starvation")
arch_counts = k[k.archetype != "not responsive"].groupby(["archetype", "direction_0_10d"]).size().unstack(fill_value=0)
print(arch_counts.to_string())

# --- Varala resupply genes ---------------------------------------------------
v = pd.read_csv(os.path.join(CACHE, "de_all_tests_all_samples.csv"))
lfc = v.pivot_table(index="gene", columns="time_min", values="log2FoldChange")
padj = v.pivot_table(index="gene", columns="time_min", values="padj").reindex(lfc.index).fillna(1)
times = list(lfc.columns)
common = lfc.index.intersection(k.index)
print(f"genes in both datasets: {len(common)}")

rows = []
for t in times:
    vs = padj.loc[common, t] < ALPHA
    for kc, label in [("02", "starvation 0-2 d"), ("010", "starvation 0-10 d")]:
        ks = sig[kc].reindex(common).fillna(False)
        both = common[vs & ks]
        opp = (np.sign(lfc.loc[both, t]) != sgn[kc][both]).sum()
        # Fisher test of overlap vs genes tested in both
        a, b = len(both), int((vs & ~ks).sum())
        c, d = int((~vs & ks).sum()), int((~vs & ~ks).sum())
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        rho = stats.spearmanr(lfc.loc[both, t], k.loc[both, f"lfc_{kc}"]).statistic if a > 5 else np.nan
        rows.append(dict(varala_time_min=t, krapp_comparison=label, varala_DE=int(vs.sum()),
                         krapp_DE=int(ks.sum()), overlap=a, overlap_odds_ratio=round(odds, 2),
                         overlap_p=p, n_opposite_direction=int(opp),
                         frac_opposite=round(opp / a, 3) if a else np.nan,
                         spearman_lfc_on_overlap=round(rho, 3) if a > 5 else np.nan))
ov = pd.DataFrame(rows)
ov.to_csv(os.path.join(HERE, "krapp_vs_varala_overlap.csv"), index=False, float_format="%.3g")
print(ov.to_string(index=False))

# Archetype membership of Varala early (<=20 min, excluding the 10/15 batch -> 5 & 20)
# and late (60-120 min) resupply genes: which starvation archetypes do they sit in?
early = common[(padj.loc[common, [5, 20]] < ALPHA).any(axis=1)]
late = common[(padj.loc[common, [60, 90, 120]] < ALPHA).any(axis=1)]
tab = pd.DataFrame({
    "all genes in both": k.loc[common, "archetype"].value_counts(),
    "Varala DE at 5 or 20 min": k.loc[early, "archetype"].value_counts(),
    "Varala DE at 60-120 min": k.loc[late, "archetype"].value_counts(),
}).fillna(0).astype(int)
tab.to_csv(os.path.join(HERE, "krapp_archetypes_vs_varala.csv"))
print(tab.to_string())

out = k[k.archetype != "not responsive"][["catma", "desc", "lfc_02", "bonf_02", "lfc_210", "bonf_210",
                                           "lfc_010", "bonf_010", "archetype", "direction_0_10d"]]
out.to_csv(os.path.join(HERE, "krapp_root_starvation_archetypes.csv"), float_format="%.4g")

markers = {"AT1G08090": "NRT2.1", "AT1G12110": "NRT1.1", "AT1G77760": "NIA1", "AT1G37130": "NIA2",
           "AT2G15620": "NIR1", "AT4G37540": "LBD39", "AT5G60770": "NRT2.4", "AT1G13300": "HRS1",
           "AT1G68670": "HHO2", "AT3G47340": "ASN1", "AT5G53460": "GLT1", "AT1G24280": "G6PD3"}
m = k.reindex(list(markers))[["lfc_02", "bonf_02", "lfc_010", "bonf_010", "archetype"]]
m.insert(0, "symbol", [markers[g] for g in m.index])
m["varala_lfc_120min"] = lfc.reindex(m.index)[120].round(2)
print(m.round(3).to_string())
m.to_csv(os.path.join(HERE, "krapp_markers.csv"), float_format="%.3g")
