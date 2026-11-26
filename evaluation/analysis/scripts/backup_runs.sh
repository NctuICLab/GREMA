#!/bin/bash
# Pack the run directories of the analyses (outside the repository) into one archive per folder, with a manifest and checksums,
# ready to be kept as a backup or uploaded (e.g. to Zenodo) at the end of a stage.
# usage: backup_runs.sh STAGE_NAME [RUNS_DIR [OUT_DIR]]
#   RUNS_DIR  default /mnt/10T-backup-mingju/GREMA_analysis_runs   (holds external/, runs/, logs/)
#   OUT_DIR   default RUNS_DIR/../GREMA_analysis_backups/STAGE_NAME
# The environment folder .pixi is left out (rebuild it with `pixi install` in evaluation/analysis/env/).
set -eu
STAGE=${1:?usage: backup_runs.sh STAGE_NAME [RUNS_DIR [OUT_DIR]]}
RUNS=${2:-/mnt/10T-backup-mingju/GREMA_analysis_runs}
OUT=${3:-$(dirname "$RUNS")/GREMA_analysis_backups/$STAGE}
mkdir -p "$OUT"
[ -z "$(ls -A "$OUT")" ] || { echo "$OUT is not empty"; exit 1; }
printf 'archive\tfiles\tbytes_uncompressed\n' > "$OUT/MANIFEST.tsv"
for group in external runs logs; do
  [ -d "$RUNS/$group" ] || continue
  for d in "$RUNS/$group"/*; do
    name=$group-$(basename "$d")
    if [ -d "$d" ]; then
      files=$(find "$d" -path '*/.pixi' -prune -o -type f -print | wc -l)
      bytes=$(du -sb --exclude=.pixi "$d" | cut -f1)
      tar --exclude=.pixi -C "$RUNS/$group" -cf - "$(basename "$d")" | zstd -q -T4 -o "$OUT/$name.tar.zst"
    else
      files=1; bytes=$(stat -c %s "$d")
      tar -C "$RUNS/$group" -cf - "$(basename "$d")" | zstd -q -o "$OUT/$name.tar.zst"
    fi
    printf '%s.tar.zst\t%s\t%s\n' "$name" "$files" "$bytes" >> "$OUT/MANIFEST.tsv"
  done
done
( cd "$OUT" && sha256sum *.tar.zst > SHA256SUMS )
echo "stage $STAGE: $(ls "$OUT"/*.tar.zst | wc -l) archives, $(du -sh "$OUT" | cut -f1) in $OUT"
echo "restore: tar --zstd -xf <archive> -C <RUNS_DIR>/<external|runs|logs>"
