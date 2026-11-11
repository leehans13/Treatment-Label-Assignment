#!/usr/bin/env bash
# Run from any directory. Existing complete files are checked before reuse.
set -u
cd "$(dirname "$0")/.." || exit 1
mkdir -p data/cellxgene data/kang data/schmidt
failed=0
fetch() {
  local path="$1" url="$2" size="$3" checksum="${4:-}" actual
  if [[ -f "$path" ]] && [[ $(wc -c < "$path" | tr -d ' ') -gt "$size" ]]; then
    echo "File larger than expected; move it aside and retry: $path" >&2
    failed=1; return
  fi
  if [[ ! -f "$path" ]] || [[ $(wc -c < "$path" | tr -d ' ') != "$size" ]]; then
    curl --fail --location --retry 3 --continue-at - --output "$path" "$url" || { echo "Download failed: $url" >&2; failed=1; return; }
  fi
  actual=$(wc -c < "$path" | tr -d ' ')
  if [[ "$actual" != "$size" ]]; then echo "Incorrect size: $path ($actual, expected $size)" >&2; failed=1; return; fi
  if [[ -n "$checksum" ]]; then
    if command -v md5sum >/dev/null; then actual=$(md5sum "$path" | awk '{print $1}'); else actual=$(md5 -q "$path"); fi
    if [[ "$actual" != "$checksum" ]]; then echo "Checksum mismatch: $path; move it aside and retry" >&2; failed=1; return; fi
  fi
  echo "Verified $path"
}
fetch data/cellxgene/yazar.h5ad https://datasets.cellxgene.cziscience.com/1e44db10-b572-46cc-adae-dcc7acd44ca6.h5ad 4434273970
fetch data/cellxgene/reichart-all.h5ad https://datasets.cellxgene.cziscience.com/6d0b8b8a-2b22-4578-b8ee-8ca4ab893594.h5ad 6890428078
fetch data/cellxgene/reed-global.h5ad https://datasets.cellxgene.cziscience.com/ed8cdda6-9985-4a62-933a-eee60b11cb18.h5ad 8336712858
fetch data/schmidt/CRISPRa-Perturb-seq.zip https://zenodo.org/api/records/5784651/files/CRISPRa-Perturb-seq.zip/content 3055676205 7e8b8f2c12360272a52dd176801093e4
fetch data/kang/pbmcs_ctrl.RDS https://ndownloader.figshare.com/files/40054042 9997173 b1225a20f904ffee0015eb073574b8be
fetch data/kang/pbmcs_stim.RDS https://ndownloader.figshare.com/files/40054048 11091577 37297c31fd7fcac3301dd1ffc752c9fa
exit "$failed"
