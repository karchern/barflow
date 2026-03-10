#!/usr/bin/env bash
set -euo pipefail

base_outdir="/scratch/karcher/barflow_output"
#project_dir="/g/typas/Personal_Folders/Nic/barflow"
project_dir=$(pwd)

if [[ "$#" -lt 1 ]]; then
    echo "Error: RUN_ID is required as the first argument. Just needs to be a unique identifier for this test run (e.g. 'test_run_1')"
    echo "Usage: $0 <RUN_ID>"
    exit 1
fi

RUN_ID="$1"

echo "Starting test run with RUN_ID: ${RUN_ID} with base_outdir ${base_outdir}"

mkdir -p "${base_outdir}/${RUN_ID}"

export NXF_CACHE_DIR="${base_outdir}/${RUN_ID}/.nextflow"   # per-run cache

# nextflow run "${project_dir}/main.nf" \
#   -profile test_run \
#   --test_run true \
#   --samplesheet "${project_dir}/test_data/test_samplesheet.csv" \
#   --good_barcodes_csv "${project_dir}/test_data/good_barcodes_with_locustag.csv" \
#   --outdir "${base_outdir}/${RUN_ID}" \
#   --comparisons "${project_dir}/test_data/comparisons.json" \
#   --singularity true \
#   --conda false \
#   -resume \
#   -work-dir "${base_outdir}/${RUN_ID}/work"

nextflow run "${project_dir}/main.nf" \
  -profile test_run \
  --test_run true \
  --samplesheet "${project_dir}/test_data/test_samplesheet.csv" \
  --sample_librarymap "${project_dir}/test_data/sample_library_map.csv" \
  --good_barcodes_csv "${project_dir}/test_data/good_barcodes_with_locustag.csv" \
  --outdir "${base_outdir}/${RUN_ID}" \
  --comparisons "${project_dir}/test_data/comparisons.json" \
  --singularity true \
  --conda false \
  -work-dir "${base_outdir}/${RUN_ID}/work"