#!/usr/bin/env bash
set -euo pipefail

prefix=${RADSEX_PREFIX:-benchmarks/tools/radsex-cpp}
repo=${RADSEX_REPO:-https://github.com/RomainFeron/RADSex.git}
ref=${RADSEX_REF:-v1.2.0}
jobs=${RADSEX_BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-4}}

mkdir -p "$(dirname "$prefix")"
if [[ ! -d "$prefix/.git" ]]; then
    git clone "$repo" "$prefix"
fi

git -C "$prefix" fetch --tags origin
git -C "$prefix" checkout "$ref"

marker_header="$prefix/src/marker.h"
if [[ -f "$marker_header" ]] && ! grep -q '<cstdint>' "$marker_header"; then
    tmp=$(mktemp)
    {
        echo '#include <cstdint>'
        cat "$marker_header"
    } > "$tmp"
    mv "$tmp" "$marker_header"
fi

if [[ -f "$prefix/Makefile" ]]; then
    make -C "$prefix" -j "$jobs"
elif [[ -f "$prefix/CMakeLists.txt" ]]; then
    cmake -S "$prefix" -B "$prefix/build"
    cmake --build "$prefix/build" --parallel "$jobs"
else
    echo "No supported build system found in $prefix" >&2
    exit 2
fi

for candidate in "$prefix/bin/radsex" "$prefix/radsex" "$prefix/build/radsex"; do
    if [[ -x "$candidate" ]]; then
        echo "$candidate"
        exit 0
    fi
done

echo "RADSex build finished but no radsex executable was found under $prefix" >&2
exit 3
