# Which early N-resupply time points are redundant? Evidence from Varala et al. 2018 (Arabidopsis roots)

Script: `varala_timepoint_redundancy.py` (python ≥3.11, pydeseq2 0.5.4). Inputs come from `../scripts/fetch_data.sh`. All numbers below are written to the CSVs and `results.json` in this folder.

## 0. Data source and design (read from the data, not assumed)

**Source deviation.** Dryad (doi:10.5061/dryad.248g184), GEO and TAIR were all blocked by this environment's network policy. Because of that:

- **Counts.** The root count table is a public GitHub mirror (OceaneCsn/integrative_GRN_N_induction, pinned commit) with the same layout as `RootRawCounts.txt`. **All 57 columns are identical, gene for gene, to the per-sample GEO count files** GSM2570359–GSM2570415, which a second group mirrors independently (IcyFermion/feature_synergy). So the counts are the GEO GSE97500 root data. The sample→GSM map is in `design_samples.csv`.
- **Missing metadata and scripts.** I could not read the GEO sample metadata text or the authors' `DE_Spline.R` / `DEgene_Plots.R`.

**Design**, parsed from the column names (`T0-r`, `C<min>-r` = KCl, `N<min>-r` = N):

| | 0 | 5 | 10 | 15 | 20 | 30 | 45 | 60 | 90 | 120 min |
|---|---|---|---|---|---|---|---|---|---|---|
| T0 (shared, untreated) | 3 | | | | | | | | | |
| KCl mock | | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 |
| N | | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 |

This matches the expected N vs KCl time course. With 3 replicates per cell, per-point tests are feasible (thin but standard), so I did not fall back to the spline peak-time approach. 22,528 of 28,717 genes pass the expression filter (≥10 reads in ≥3 samples).

**QC problem: the 10 and 15 min samples are a separate batch** (`qc_pca.png`).
- **PCA.** PC1 explains 74% of the variance among the top 2,000 genes. It separates all 12 samples at 10 and 15 min (KCl *and* N) from everything else, while 5 and 20 min sit with T0.
- **Size of the shift.** 5,848 of 16,035 expressed genes differ more than 2-fold between these samples and the rest, and their chloroplast-encoded read fraction is 6× lower (0.02% vs 0.14%).
- **Interpretation.** A transient shift like this in both treatments is not plausible biology. It looks like a library/processing batch, and another group using these data dropped the same samples.
- **What it does and doesn't affect.** N vs KCl *within* 10 or 15 min is still a within-batch comparison, but the 10/15 DE sets are not comparable to the other points.
- **How I handled it.** I kept the samples in the main analysis, because the design brief said not to change the design without asking. Every overlap result is repeated without them ("without 10/15").

## 1. Sanity check against the paper's ~2,174 genes (step 2): not reproduced

I could not obtain `DE_Spline.R`, so I cannot summarise how the authors actually called genes. As a stand-in I fit a time × treatment spline model:
- **Model.** log2-CPM, limma-style moderated F-test on the N × B-spline(time) terms, with the N−KCl difference forced to 0 at T0.
- **Grid of choices.** Spline df 3/4/5, linear or log time, FDR 0.01/0.05, with and without 10/15 min (`spline_sanity_check_counts.csv`).

The count ranges from **1,757 to 13,988 genes** depending on those choices:
- The closest settings are linear time, df 3, FDR 0.01 (1,757), and log time, df 3–4, FDR 0.01 with max |fitted log2 diff| > 1 (2,339–2,437).
- That is **not a reproduction**, only the same order of magnitude under some settings. The 8-fold spread across reasonable choices is itself a warning: "number of N-responsive genes" in this dataset depends heavily on the model.

## 2. Per-time-point DE (step 3)

pydeseq2, one model `~group` over all 57 samples, Wald contrast N_t vs KCl_t, padj < 0.05 (`de_counts_per_timepoint.csv`, `de_genes_per_timepoint.csv`).

| min | 5 | 10* | 15* | 20 | 30 | 45 | 60 | 90 | 120 |
|---|---|---|---|---|---|---|---|---|---|
| padj<0.05 | 1,068 | 4,533 | 5,914 | 1,028 | 1,416 | 4,141 | 6,460 | 4,080 | 8,746 |
| + \|log2FC\|>0.5 | 655 | 3,110 | 4,324 | 930 | 1,260 | 3,101 | 4,669 | 3,091 | 6,427 |
| + \|log2FC\|>1 | 80 | 747 | 1,475 | 522 | 629 | 1,072 | 1,737 | 1,326 | 2,485 |

\* batch-affected points.

The sets are large: 15,463 genes, two-thirds of the expressed transcriptome, are DE at some point at padj < 0.05. Small effects are easy to detect at roughly 45 M reads per library.

Two quality flags:
- **10/15 min.** Their counts spike above both neighbours (5 and 20 min).
- **5 min.** Its log2FCs correlate *negatively* with every later point (median ρ −0.32, even for its own DE genes) and track expression level (ρ 0.21 vs ≈0 elsewhere; `lfc_consistency_diagnostics.csv`). The fast primary genes are right at 5 min (NIR1 +1.6, LBD39 +1.8, NIA1 +0.9, NRT1.1 +0.8 log2FC). But most of the 1,068 look like a technical KCl-vs-N library difference.

## 3. Answers

### (1) How much does each added time point contribute?

Adding points in chronological order (`fig4b_cumulative_union_all.png`, `cumulative_union.csv`), new genes per point at padj < 0.05:

| added | 5 | 10 | 15 | 20 | 30 | 45 | 60 | 90 | 120 |
|---|---|---|---|---|---|---|---|---|---|
| new genes | 1,068 | +4,250 | +3,590 | +286 | +502 | +1,422 | +1,839 | +828 | +1,678 |
| cumulative | 1,068 | 5,318 | 8,908 | 9,194 | 9,696 | 11,118 | 12,957 | 13,785 | 15,463 |

- **No plateau by 2 h.** The last point still adds 11% of the union (18% at |log2FC| > 1; 19% without 10/15).
- **20 and 30 min add the least**: 2–3% each.
- **Planned points only** (15, 30, 60, 120 min; `fig4e_cumulative_union_planned.png`): 5,914 → +743 → +3,473 → +2,961. 30 min adds 6% and 2 h still adds 23%. Using 20 min instead of the batch-affected 15 min: 1,181 → +1,000 → +5,794 → +4,072.

### (2) Which time points look redundant?

**20 and 30 min are the most redundant, consistently.**
- **Greedy removal, all 9 points** (`fig4d_greedy_redundancy.png`, `greedy_redundancy.csv`):
  20 (78 genes lost) → 30 (94) → 5 (263) → 45 (287) → 90 (475) → 60 (1,104) → 10 (1,544) → 15 (2,872).
- **Robust across settings.** 20 and 30 min are among the first three dropped in 11 of 12 settings: 3 thresholds × strict or "soft" loss × with or without the 10/15 batch. Under the soft criterion a gene only counts as lost if no remaining point even nominally supports it in the same direction. In the remaining setting, 20 min is dropped fourth.
- **Size of the loss.** Dropping both loses about 1% of the union.

**45 min is the next most redundant.** It is highly correlated with 1 h and 2 h (log2FC ρ 0.82 and 0.78; `fig4c_lfc_correlation_all.png`). 45 min, 1 h, 1.5 h and 2 h report broadly the same response, and it grows over time.

**Within the planned design, 30 min is the most dispensable.**
- **30 min first.** Dropped first in all 6 settings (3 thresholds × strict/soft): 176 genes (1.3%) lost at padj < 0.05, 125 at |log2FC| > 1.
- **1 h and 2 h are not redundant.** Removing 1 h loses 1,297 genes (2,074 in the 20/30/60/120 variant).
- **15 min.** It comes out as the most informative point, but that is inflated by the 10/15 batch. Its stand-in, 20 min, is itself as dispensable as 30 min.

**Read "redundant" carefully.**
- **Low Jaccard everywhere.** Jaccard between points is low (≤ 0.43; `fig4c_jaccard_all.png`), so threshold noise inflates "unique" genes.
- **Small sets look redundant.** 20 and 30 min look redundant mostly because their sets are small (about 1,000–1,400 genes, Jaccard 0.29 with each other), not because a neighbour duplicates them.
- **Unique-gene counts miss timing.** Early points earn their place through *timing*, not unique genes:
  - Uptake/reduction and signalling genes (NIR1, NIA1, LBD39, NRT2.1, HRS1, HHO2, G6PD3) are called at 5–15 min.
  - Assimilation and cytokinin-synthesis genes come later:
    - NADH-GOGAT/GLT1 and CYP735A2 first reach significance at 30 min.
    - GLN2, ASN1 and ASN2 are weakly significant at the batch-affected 15 min point, then consistently significant from 30 min.
  - See `marker_genes_log2fc.csv`.
  - An Avena design without a ≤15 min point could not order these events.

### Biology check (step 5)

GO BP, hypergeometric test, terms propagated through go-basic.obo, background = 17,231 annotated expressed genes (`go_bp_enrichment.csv`, `go_theme_check.csv`).

- **padj-only clusters.** Dominated by ribosome/rRNA biogenesis. The sets are too large for specific signals.
- **|log2FC| > 1 clusters:**
  - The **early cluster** (first DE at 5–20 min, 2,172 genes) is enriched for *response to nitrate* (4.2×, FDR 0.007), *nitrate transmembrane transport* (4.5×, FDR 0.007) and *nitrate import* (3.5×, FDR 0.025).
  - The **late cluster** (60–120 min, 2,050 genes) is enriched for *rRNA modification* (4.3×, FDR 4e-11, almost all up in N), consistent with growth resuming. *Nitrate metabolic process* is borderline (5×, FDR 0.06).
  - **Hypoxia-response terms** are enriched in both clusters, with mixed direction. Nitrate→NO signalling (e.g. HB1/GLB1 up 1.6–2.3 log2FC from 1 h) is a plausible link, but this is speculative.
- **Amino-acid biosynthesis is not enriched in any cluster.** The "later amino-acid synthesis" shift shows only at marker-gene level (ASN1/2, GLT1, GLN2).
- **Verdict.** "Uptake/reduction early" is supported, and "amino-acid synthesis later" is supported only at marker level.
- **Genes unique to a single time point** give few or no significant terms, except at 5 min (respiratory chain, 6 genes) and 2 h (endomembrane transport, cell wall). This fits the idea that most single-point genes are noise rather than distinct biology.

### Optional step 6: Krapp et al. 2011 starvation vs Varala resupply (`krapp_comparison.py`)

**Flag first.** This compares long-term N *starvation* with minutes-to-hours N *resupply*, across labs and platforms (CATMA arrays vs RNA-seq). Any overlap is suggestive only.

**Krapp's design**, from the paper's Supplemental Table VII (supplementary zip):
- Roots after 0, 2 and 10 days of N starvation.
- 2 biological replicates with dye swap; genes called at Bonferroni p < 0.05.
- All comparisons are against day 0 of the same culture. There is no N-replete control at 2 or 10 days, so starvation is confounded with 10 days of growth.

**Starvation archetypes in roots** (`krapp_root_starvation_archetypes.csv`):

| archetype | down in starvation | up in starvation |
|---|---|---|
| early, sustained (by 2 d, still at 10 d) | 133 | 32 |
| early, still growing 2–10 d | 27 | 28 |
| early, transient (gone by 10 d) | 44 | 23 |
| late (only after 2 d) | 556 | 697 |

Most change happens *between 2 and 10 days*.

**Direction: starvation is the mirror image of resupply** (`krapp_vs_varala_overlap.csv`):
- **Overlap is enriched.** At every Varala time point the overlap with Krapp's DE genes is larger than chance (odds ratio 1.4–8).
- **Shared genes move oppositely.** From 20 min on, 73–92% of shared genes go in opposite directions. Log2FC correlation on the overlap runs from ρ −0.38 to −0.71 and strengthens with resupply time: genes repressed by starvation are induced by resupply.
- **5 min doesn't follow** (66% *same* direction for the 0–10 d comparison), consistent with the technical concern about that point.
- **Early resupply genes concentrate in the "early, sustained" class.** Genes DE at 5 or 20 min are enriched about 4× in it (3.3% vs 0.8% of all genes); NIA1, NIR1, LBD39 and G6PD3 are all down by day 2 and still down at day 10.
- **NRT2.4 is a clean low-N marker.** It goes up in starvation (+2.2 log2FC by 10 d) and down on resupply.
- **Time span.** Krapp adds coarse information on the days scale: 2 d and 10 d only, with nothing past 10 d.

### (3) Caveats

- **Species.** Arabidopsis, not Avena. Timing of primary N signalling is broadly conserved, but kinetics and gene sets will differ, and one-to-one orthology for oat is incomplete.
- **Growth system.** Hydroponic/plate-grown seedlings after N starvation, not soil-grown plants. Soil N release is slow and patchy, so responses there will be slower and smaller.
- **Window.** Data stop at 2 h. **Nothing here informs the 4 h, 8 h, 24 h, 2 d … 28 d points.** The response was still growing at 2 h (no plateau), so later points are likely to add many genes.
- **Batch at 10/15 min.** Results at 15 min, one of your planned points, are confounded with a library batch. The 5 min set also looks partly technical.
- **Threshold dependence.** Unique-gene counts and greedy orderings depend on thresholds and on per-point noise. Each time point is a separate harvest, so per-point technical variation is mixed into the "time" signal.
- **Spline cross-check.** The authors' spline result (~2,174 genes) could not be reproduced without their script.
- **Data provenance.** Counts and annotations come from GitHub mirrors. The counts were verified exactly against the GEO per-sample files. The GO table is TAIR's ATH_GO_GOSLIM (2024-10-01) and the ontology is go-basic (2026-03-25); the ontology file's checksum was verified.
