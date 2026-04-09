#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR=$(cd "$(dirname "$0")" && pwd)
ENV_PATH="${PROJECT_DIR}/.conda_cache/nf-test-env"
CFG1="${PROJECT_DIR}/.nf-test.1_2.config"
CFG2="${PROJECT_DIR}/.nf-test.2_2.config"

cleanup() {
    rm -f "${CFG1}" "${CFG2}"
}
trap cleanup EXIT

if [[ ! -d "${ENV_PATH}" ]]; then
    echo "Missing environment at ${ENV_PATH}."
    echo "Create it first with: mamba env create -p .conda_cache/nf-test-env -f envs/nf-test.yaml"
    exit 1
fi

# Optional first argument: test selector path/tag args passed through to both shards.
EXTRA_ARGS=("$@")

run_shard() {
    local shard="$1"
    local shard_id
    shard_id="${shard//\//_}"
    local cfg="${PROJECT_DIR}/.nf-test.${shard_id}.config"

    cat > "${cfg}" <<EOF
config {
    testsDir "tests/nf-test"
    workDir ".nf-test-${shard_id}"
    configFile "tests/nextflow.config"
    profile "test_run"
    withTrace false
}
EOF

    conda run -p "${ENV_PATH}" nf-test test -c "${cfg}" --without-trace --shard "${shard}" "${EXTRA_ARGS[@]}"
}

run_shard "1/2" &
PID1=$!
run_shard "2/2" &
PID2=$!

wait "${PID1}"
wait "${PID2}"

echo "Both shards completed successfully."
