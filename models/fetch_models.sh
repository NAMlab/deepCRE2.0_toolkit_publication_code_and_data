#!/usr/bin/env bash
# Purpose: Download the trained deepCRE HDF5 model archives from GitHub release models-20261126,
#   verify SHA-256 checksums against models/SHA256SUMS, unpack the .h5 weights, and create
#   compatibility symlinks (old_file -> new_id) listed in models/model_id_mapping.tsv.
# Manuscript (Peleke et al., 2026): Provides the 2026 single-species (SSR, 232 files),
#   shuffled-sequence control (SSC, 217 files), and multi-species (MSR, 52 files) models, plus
#   optionally the Peleke et al. (2024) models (P2024, 45 files unpacked to models/peleke2024/)
#   used in Supplementary Figures 5 and 6.
# Usage (from any directory): bash models/fetch_models.sh [SSR] [SSC] [MSR] [P2024]
#   Default groups when omitted: SSR SSC MSR (~2.4 GB); add P2024 (+0.2 GB) for 2024 models.
# How to change or extend: Point URL to a new release tag and update SHA256SUMS /
#   model_id_mapping.tsv when publishing new model weights or additional species/tissues.
set -euo pipefail
cd "$(dirname "$0")"
URL="https://github.com/NAMlab/deepCRE2.0_toolkit_publication_code_and_data/releases/download/models-20261126"
GROUPS_="${*:-SSR SSC MSR}"
for g in $GROUPS_; do
  if [ "$g" = P2024 ]; then t="deepCRE_2024_models_Peleke2024.tar"; else t="deepCRE_TK_2026_models_${g}.tar"; fi
  [ -f "$t" ] || curl -fL -o "$t" "$URL/$t"
  grep " $t\$" SHA256SUMS | sha256sum -c -
  tar -xf "$t" && rm "$t"
done
tail -n +2 model_id_mapping.tsv | while IFS=$'\t' read -r kind code tissue fold old new rest; do
  [ -f "$new" ] && [ "$old" != "$new" ] && ln -sf "$new" "$old"
done
echo "models in $(pwd): $(ls *.h5 2>/dev/null | wc -l) files/links"
