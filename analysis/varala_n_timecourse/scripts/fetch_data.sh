#!/usr/bin/env bash
# Download inputs for the Varala et al. 2018 time-point redundancy analysis.
#
# The canonical sources (Dryad doi:10.5061/dryad.248g184, GEO GSE97500,
# TAIR, geneontology.org) were not reachable from the environment this was
# built in, so every file below comes from a public GitHub mirror pinned to
# a commit. Integrity checks:
#   * root counts (Dryad-style table, 57 columns) are cross-checked
#     column by column against the per-sample GEO count files (GSM*.counts.txt)
#     mirrored independently by a second group -- done in the analysis script;
#   * go-basic.obo matches its Git-LFS sha256 oid.
# If Dryad/GEO are reachable for you, prefer them (see bottom of file).
set -euo pipefail
DATA="${1:-$(cd "$(dirname "$0")/.." && pwd)/data}"
mkdir -p "$DATA/geo_sample_counts"
LIST="$(cd "$(dirname "$0")" && pwd)/geo_sample_files.txt"
cd "$DATA"

OCEANE=https://raw.githubusercontent.com/OceaneCsn/integrative_GRN_N_induction/f8ba23cd9e405c46473b6488bd97dc5c34ad2558
ICY=https://raw.githubusercontent.com/IcyFermion/feature_synergy/c72781034400dc3ff3dde41868fff0d4ab45f69c
SLIM=https://media.githubusercontent.com/media/mduchoslav/Brassicaceae_orthology/e3666bc961bfd04128f400a6da5e87176ca29cd7
OBO=https://media.githubusercontent.com/media/heart-gen/isograph-brain-aging-benchmarking/aef9f360e4966e67092d3c33fd105d68657ae550

# 1. Root raw counts, Varala 2018 (same layout as Dryad RootRawCounts.txt:
#    T0 x3, KCl 'C' and N 'N' at each time x3)
curl -fsSL -o RootRawCounts_mirror.txt \
  "$OCEANE/data/Arabidopsis_Varala_2018_nitrate_response_roots_with_control.txt"

# 2. Per-sample GEO GSE97500 HTSeq count files (roots + shoots, GSM2570302-415)
while read -r f; do
  [ -s "geo_sample_counts/$f" ] || curl -fsSL -o "geo_sample_counts/$f" "$ICY/data/arabidopsis/GSE97500/raw_counts/$f"
done < "$LIST"

# 3. GO: TAIR gene->GO table (Last Updated 2024-10-01) and go-basic.obo (2026-03-25)
curl -fsSL -o ATH_GO_GOSLIM.txt.gz "$SLIM/functional_annotation/functional_annotation_data/Araport11/ATH_GO_GOSLIM.txt.gz"
curl -fsSL -o go-basic.obo "$OBO/inputs/go_annotations/go-basic.obo"
echo "a77e356737dab39a4f620dce35fc4d6eb531c4b6153af6cacaaa322b49b804bd  go-basic.obo" | sha256sum -c -

# Canonical sources, if reachable:
#   https://datadryad.org/stash/dataset/doi:10.5061/dryad.248g184  (RootRawCounts.txt, DE_Spline.R, DEgene_Plots.R)
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97500
#   https://www.arabidopsis.org/download/ (ATH_GO_GOSLIM.txt), https://current.geneontology.org/ontology/go-basic.obo
