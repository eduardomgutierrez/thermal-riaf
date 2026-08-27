#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${project_dir}/radproc/build"

if ! command -v cmake >/dev/null 2>&1; then
    echo "error: CMake is not installed; see README.md" >&2
    exit 2
fi

cmake -S "${project_dir}/radproc" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${build_dir}" --parallel

echo
echo "radproc is ready: ${build_dir}/src/adaf/adaf"
echo "Run a model with:"
echo "  python riaf_pipeline.py examples/external-profile.toml"
