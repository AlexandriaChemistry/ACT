#!/bin/bash
#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# Copyright (C) 2026
#
# Developers:
#             Mohammad Mehdi Ghahremanpour,
#             Paul J. van Maaren,
#             David van der Spoel (Project leader)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
#
#
# Line and function coverage of ACT by its tests, with the source-based
# coverage of clang. Builds DIR (default build_COVERAGE) with coverage
# instrumentation, configured like build_Release_DOUBLE, runs all tests and
# writes in DIR/coverage:
#   report.txt   Coverage of each file
#   summary.txt  Coverage of each directory of src/act and the files with
#                the most lines that are not covered
#   html/        The source annotated with coverage
# The GROMACS and external code and the tests themselves are left out.
#
# Usage: admin/coverage.sh [-b DIR] [-j N]

set -u

root=$(cd "$(dirname "$0")/.." && pwd)
build=build_COVERAGE
jobs=8
while [ $# -gt 0 ]; do
    case "$1" in
        -b) build=$2; shift ;;
        -j) jobs=$2; shift ;;
        *) echo "Unknown option $1"; exit 2 ;;
    esac
    shift
done
case "$build" in
    /*) ;;
    *) build=$root/$build ;;
esac
ref=$root/build_Release_DOUBLE

cacheValue()
{
    grep "^$2:" "$1/CMakeCache.txt" | head -1 | cut -d= -f2-
}

cmake=$(cacheValue "$ref" CMAKE_COMMAND)
ctest=$(dirname "$cmake")/ctest
if [ ! -f "$build/CMakeCache.txt" ]; then
    flags="-fprofile-instr-generate -fcoverage-mapping"
    "$cmake" -S "$root" -B "$build" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER="$(cacheValue "$ref" CMAKE_C_COMPILER)" \
        -DCMAKE_CXX_COMPILER="$(cacheValue "$ref" CMAKE_CXX_COMPILER)" \
        -DCMAKE_PREFIX_PATH="$(cacheValue "$ref" CMAKE_PREFIX_PATH)" \
        -DACT_DOUBLE="$(cacheValue "$ref" ACT_DOUBLE)" \
        -DACT_MPI="$(cacheValue "$ref" ACT_MPI)" \
        -DCMAKE_C_FLAGS="$flags" -DCMAKE_CXX_FLAGS="$flags" > /dev/null || exit 1
fi
echo "=== Building $build"
"$cmake" --build "$build" -j "$jobs" --target all tests > "$build/coverage-build.log" 2>&1 ||
    { tail -30 "$build/coverage-build.log"; exit 1; }

out=$build/coverage
rm -rf "$out"
mkdir -p "$out/profiles"
echo "=== Running the tests"
(cd "$build" && LLVM_PROFILE_FILE="$out/profiles/%p-%m.profraw" "$ctest" -j "$jobs" > "$out/tests.log" 2>&1)
grep -E "tests passed|\(Failed\)" "$out/tests.log"

xcrun llvm-profdata merge -sparse "$out"/profiles/*.profraw -o "$out/merged.profdata" || exit 1
rm -rf "$out/profiles"
objects=()
for b in "$build"/bin/*-test; do
    objects+=(-object "$b")
done
ignore='(src/gromacs|src/external|testutils|/tests/|googletest|_deps|/usr/|python|config\.h)'
xcrun llvm-cov report "$build/bin/alexandria" "${objects[@]}" -instr-profile="$out/merged.profdata" \
    -ignore-filename-regex="$ignore" > "$out/report.txt"
xcrun llvm-cov show "$build/bin/alexandria" "${objects[@]}" -instr-profile="$out/merged.profdata" \
    -ignore-filename-regex="$ignore" -format=html -output-dir="$out/html" > /dev/null

python3 - "$out/report.txt" > "$out/summary.txt" <<'EOF'
import collections, sys
files = []
for line in open(sys.argv[1]):
    p = line.split()
    if len(p) == 13 and "src/act/" in p[0] and int(p[7]) > 0:
        name = p[0][p[0].index("src/act/") + 8:]
        files.append((name, int(p[4]), int(p[5]), int(p[7]), int(p[8])))
dirs = collections.defaultdict(lambda: [0, 0, 0, 0])
for name, nf, uf, nl, ul in files:
    d = dirs[name.split("/")[0]]
    for k, v in enumerate((nf, uf, nl, ul)):
        d[k] += v
print(f"{'Directory':18s} {'Lines':>7s} {'Not covered':>12s} {'Lines %':>8s} {'Functions %':>12s}")
for name, (nf, uf, nl, ul) in sorted(dirs.items(), key=lambda x: -x[1][3]):
    print(f"{name:18s} {nl:7d} {ul:12d} {100*(nl-ul)/nl:8.1f} {100*(nf-uf)/max(nf, 1):12.1f}")
nf, uf, nl, ul = (sum(d[k] for d in dirs.values()) for k in range(4))
print(f"{'Total':18s} {nl:7d} {ul:12d} {100*(nl-ul)/nl:8.1f} {100*(nf-uf)/max(nf, 1):12.1f}")
print("\nFiles with the most lines not covered")
for name, nf, uf, nl, ul in sorted(files, key=lambda f: -f[4])[:40]:
    print(f"{ul:6d} of {nl:6d} lines not covered  {name}")
EOF
cat "$out/summary.txt" | head -17
echo "Details in $out/summary.txt, $out/report.txt and $out/html/index.html"
