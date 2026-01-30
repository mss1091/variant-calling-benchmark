#!/bin/bash

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFERENCE_DIR="${PROJECT_ROOT}/data/reference"
MINIFORGE_PATH="${HOME}/miniforge3"

source "${MINIFORGE_PATH}/etc/profile.d/conda.sh"
conda activate cosap

export COSAP_LIBRARY_PATH="${REFERENCE_DIR}"
export COSAP_RAMDISK_PATH="/dev/shm"
export COSAP_IN_MEMORY_MODE="false"
export COSAP_THREADS_PER_JOB="2"
