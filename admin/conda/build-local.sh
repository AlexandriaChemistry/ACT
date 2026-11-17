#!/bin/bash
# Builds the conda packages from the last commit of this git repository, with
# rattler-build, into the directory given as the first argument (default
# output-conda in the current directory). Further arguments go to
# rattler-build, for instance --variant mpi=mpich to build only that variant.
set -euo pipefail

recipe_dir=$(cd "$(dirname "$0")" && pwd)
top=$(git -C "${recipe_dir}" rev-parse --show-toplevel)
out=$(mkdir -p "${1:-output-conda}" && cd "${1:-output-conda}" && pwd)
shift || true

# An archive of the last commit, such that untracked files are not copied
work=$(mktemp -d)
trap 'rm -rf "${work}"' EXIT
git -C "${top}" archive --format=tar.gz --prefix=act/ -o "${work}/act.tar.gz" HEAD

ACT_SOURCE="${work}/act.tar.gz" rattler-build build \
    --recipe "${recipe_dir}/recipe.yaml" \
    --variant-config "${recipe_dir}/variants.yaml" \
    --channel conda-forge \
    --output-dir "${out}" \
    "$@"
