#!/bin/sh
set -eu
image=${1:-receptor-signaling:0.1.0}
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
fixture=$(mktemp -d)
trap 'rm -rf "$fixture"' EXIT HUP INT TERM
mkdir "$fixture/input" "$fixture/output"
docker run --rm --network none --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$fixture/input,target=/fixture" \
  --mount "type=bind,source=$root/tests/container_fixture.R,target=/fixture.R,readonly" \
  --entrypoint Rscript "$image" --vanilla /fixture.R /fixture

run_analysis() {
  docker run --rm --network none --read-only --tmpfs /tmp:rw,mode=1777 \
    --user "$(id -u):$(id -g)" \
    --mount "type=bind,source=$fixture/input,target=/input,readonly" \
    --mount "type=bind,source=$fixture/output,target=/output" \
    "$image" /input/config.txt
}
run_analysis
for file in gene_results.csv receptor_effects.csv analysis.rds sessionInfo.txt; do
  test -s "$fixture/output/run/analysis_0001/$file"
done
if run_analysis > "$fixture/retry.log" 2>&1; then
  echo "ERROR: an existing output directory was overwritten" >&2
  exit 1
fi
grep -q 'Output already exists' "$fixture/retry.log"
docker run --rm --network none \
  --mount "type=bind,source=$fixture/output,target=/output,readonly" \
  --entrypoint Rscript "$image" --vanilla -e \
  'r <- readRDS("/output/run/analysis_0001/analysis.rds"); stopifnot(r$results$status == "tested", r$results$q_value < .05)'
echo "Container mount test passed: offline run, read-only inputs/root, host-owned output, significance and overwrite protection."
