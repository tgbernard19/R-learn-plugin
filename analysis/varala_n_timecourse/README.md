# Varala et al. 2018 root N-resupply: time-point redundancy

Question: how different are the N-responsive root gene sets between time points
after N resupply? Used to judge which early sampling points in a planned Avena
time course are redundant. Results and caveats: `outputs/summary.md`.

```
data/      downloads only (git-ignored; recreate with scripts/fetch_data.sh)
scripts/   fetch_data.sh + list of GEO per-sample files
outputs/   varala_timepoint_redundancy.py, figures (PNG), CSVs, summary.md
cache/     cached pydeseq2 fits and run log (git-ignored)
```

Reproduce (Python 3.11+, `pip install pydeseq2 statsmodels matplotlib patsy`):

```bash
bash scripts/fetch_data.sh
python3 -I outputs/varala_timepoint_redundancy.py --recompute
```
