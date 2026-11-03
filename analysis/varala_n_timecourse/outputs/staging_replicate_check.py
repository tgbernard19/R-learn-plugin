#!/usr/bin/env python3
"""
How many replicates does a reference curve need for staging to work?
Empirical check on Varala et al. 2018 roots (3 reps per time x treatment).

For each time t and replicate r, the held-out "query" is the N-vs-KCl log2
ratio of replicate r at t.  A reference curve is built from k of the other
replicates (k = 1 or 2) at every time point; genes are chosen from the
reference only (top responsive genes), and the query is assigned to the
reference time with the highest Pearson correlation.  The 10/15 min batch is
excluded.

    python3 -I staging_replicate_check.py
"""
import os
from itertools import combinations

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
counts = pd.read_csv(os.path.join(HERE, "..", "data", "RootRawCounts_mirror.txt"), sep="\t", index_col=0)
counts = counts[(counts >= 10).sum(axis=1) >= 3]
lx = np.log(counts[(counts > 0).all(axis=1)])
sf = np.exp(lx.sub(lx.mean(axis=1), axis=0).median(axis=0))
E = np.log2(counts / sf + 1)

TIMES = [5, 20, 30, 45, 60, 90, 120]
REPS = [1, 2, 3]


def col(trt, t, r):
    return f"{trt}{t}-{r}"


rows = []
for n_genes in [200, 500, 2000]:
    for k in [1, 2]:
        for t_q in TIMES:
            for r_q in REPS:
                query = E[col("N", t_q, r_q)] - E[col("C", t_q, r_q)]
                others = [r for r in REPS if r != r_q]
                for ref_reps in combinations(others, k):
                    ref = pd.DataFrame({t: E[[col("N", t, r) for r in ref_reps]].mean(axis=1)
                                        - E[[col("C", t, r) for r in ref_reps]].mean(axis=1)
                                        for t in TIMES})
                    genes = ref.abs().max(axis=1).sort_values().index[-n_genes:]
                    cors = ref.loc[genes].corrwith(query.loc[genes])
                    t_hat = int(cors.idxmax())
                    rows.append(dict(n_genes=n_genes, ref_reps=k, true_time=t_q, query_rep=r_q,
                                     predicted=t_hat, exact=t_hat == t_q,
                                     within_one_step=abs(TIMES.index(t_hat) - TIMES.index(t_q)) <= 1,
                                     best_r=round(float(cors.max()), 3)))
res = pd.DataFrame(rows)
res.to_csv(os.path.join(HERE, "staging_replicate_check.csv"), index=False)
summ = res.groupby(["n_genes", "ref_reps"])[["exact", "within_one_step"]].mean().round(2)
print(summ.to_string())
print(res[res.n_genes == 500].groupby(["ref_reps", "true_time"])[["exact", "within_one_step"]]
      .mean().round(2).unstack(0).to_string())

# Query precision: 1 reference rep vs a query averaged over the other 2 reps
# (or a single query rep), to show replication in the experiment being staged
# matters as much as replication in the reference.
rows = []
for t_q in TIMES:
    for r_ref in REPS:
        ref = pd.DataFrame({t: E[col("N", t, r_ref)] - E[col("C", t, r_ref)] for t in TIMES})
        genes = ref.abs().max(axis=1).sort_values().index[-500:]
        others = [r for r in REPS if r != r_ref]
        for q_reps in [(others[0],), (others[1],), tuple(others)]:
            query = (E[[col("N", t_q, r) for r in q_reps]].mean(axis=1)
                     - E[[col("C", t_q, r) for r in q_reps]].mean(axis=1))
            t_hat = int(ref.loc[genes].corrwith(query.loc[genes]).idxmax())
            rows.append(dict(query_reps=len(q_reps), exact=t_hat == t_q,
                             within_one_step=abs(TIMES.index(t_hat) - TIMES.index(t_q)) <= 1))
q = pd.DataFrame(rows).groupby("query_reps")[["exact", "within_one_step"]].mean().round(2)
print("\n1-rep reference, query = 1 vs 2 reps averaged (500 genes):\n" + q.to_string())
