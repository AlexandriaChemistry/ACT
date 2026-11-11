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
# Verification of ACT after a change: builds, runs all tests with and
# without explicit SIMD, optionally without any vectorization or with the
# address sanitizer, and compares a set of training and second virial runs
# with a reference.
#
# Usage: admin/verify.sh [options]
#   -b DIR            Build directory (default build_Release_DOUBLE)
#   -j N              Parallel build and test jobs (default 8)
#   --novec           Also build and test without explicit SIMD and without
#                     vectorization by the compiler, in DIR with _NOVEC
#   --asan            Also build and test the address sanitizer build
#                     build_ASAN_DOUBLE
#   --single          Also build and test the single precision build
#                     build_Release_SINGLE, with and without explicit SIMD
#   --save-reference  Run the regression set and keep it as the reference
#   --regression      Run the regression set and compare it with the reference
#   --tolerance T     Largest relative difference accepted by --regression,
#                     default 0, i.e. identical results
#   --no-tests        Only build and run the regression set
#
# The regression set needs the input files in $ACT_VERIFY_DATA, by default
# $HOME/Alexandria/ACTdata. Its outputs are in DIR/verify/<mode>/<case>, the
# reference in DIR/verify-reference/<mode>/<case>, where mode is simd or
# nosimd. The genetic algorithm amplifies the smallest differences, so the
# trained force fields are only comparable within one build and mode.
#
# Exits with a nonzero status if anything fails.

set -u

root=$(cd "$(dirname "$0")/.." && pwd)
build=build_Release_DOUBLE
jobs=8
novec=0
asan=0
single=0
regression=""
tolerance=0
tests=1
while [ $# -gt 0 ]; do
    case "$1" in
        -b) build=$2; shift ;;
        -j) jobs=$2; shift ;;
        --novec) novec=1 ;;
        --asan) asan=1 ;;
        --single) single=1 ;;
        --save-reference) regression=save ;;
        --regression) regression=compare ;;
        --tolerance) tolerance=$2; shift ;;
        --no-tests) tests=0 ;;
        *) echo "Unknown option $1"; exit 2 ;;
    esac
    shift
done
case "$build" in
    /*) ;;
    *) build=$root/$build ;;
esac
actdata=${ACT_VERIFY_DATA:-$HOME/Alexandria/ACTdata}
failures=()

# The value of a CMake cache variable of a build directory
cacheValue()
{
    grep "^$2:" "$1/CMakeCache.txt" | head -1 | cut -d= -f2-
}

cmake=$(cacheValue "$build" CMAKE_COMMAND)
ctest=$(dirname "$cmake")/ctest
mpiexec=$(cacheValue "$build" MPIEXEC_EXECUTABLE)

# Build a directory, return nonzero on failure
buildDir()
{
    echo "=== Building $1"
    # The test programs are not part of the default target
    "$cmake" --build "$1" -j "$jobs" --target all tests > "$1/verify-build.log" 2>&1 ||
        { tail -30 "$1/verify-build.log"; return 1; }
}

# Run all tests of a build directory, the label says how
runTests()
{
    local dir=$1 label=$2
    echo "=== Tests $label"
    if (cd "$dir" && "$ctest" -j "$jobs" --output-on-failure > "verify-tests-$3.log" 2>&1); then
        grep "tests passed" "$dir/verify-tests-$3.log"
    else
        grep -E "tests passed|\*\*\*Failed|Failed  |Timeout" "$dir/verify-tests-$3.log"
        failures+=("tests $label")
    fi
}

# Run the tests of a build with and without explicit SIMD
testBuild()
{
    local dir=$1 name=$2
    if ! buildDir "$dir"; then
        failures+=("build $name")
        return
    fi
    if [ $tests -eq 1 ]; then
        runTests "$dir" "$name, explicit SIMD if available" simd
        ACT_NO_EXPLICIT_SIMD=1 runTests "$dir" "$name, without explicit SIMD" nosimd
    fi
}

# The regression set, in the current directory, with the binary $1
runCases()
{
    local act=$1 tff=$actdata/Performance/train_ff b2d=$actdata/Performance/b2
    local ga="-seed 219463979 -v 3 -random_init -optimizer HYBRID -n_crossovers 5 -temp 0.025 -anneal_begin 0.1 -anneal_end 0.9 -anneal_global_begin 0.1 -anneal_global_end 0.9 -cp_gen_interval 5 -cp_pop_frac 0.2 -ga_boltz_temp 0.5 -ga_boltz_anneal 0.5 -fc_bound 10000 -ener_boltz_temp 3000 -o ff_inter.xml -g train.log"
    local b2="-T1 273.15 -T2 773.15 -dT 5 -shellMaxIter 100 -shelltoler 1e-6 -maxdimer 2048 -g b2.log -b2 b2"
    run()
    {
        local name=$1; shift
        mkdir -p "$name"
        if ! (cd "$name" && env ${noSimd:+ACT_NO_EXPLICIT_SIMD=1} \
                  "$mpiexec" -n 8 "$act" "$@" > out.txt 2>&1); then
            echo "    $name failed, see $(pwd)/$name/out.txt"
            failures+=("regression run $name")
        fi
    }
    run hf-train train_ff $ga -max_generations 4 -max_iter 10 -pop_size 8 \
        -ff "$tff/4S-lj14-7-pg-hvdw.xml" -mp "$tff/sapt2+aug-cc-pvtz-0.03H-rmin0.9-HF_mp2-hh.xml" \
        -sel "$tff/hhmono_dimer_train-HF.dat" \
        -fit 'charge eta chi zeta vs2a sigma epsilon gamma delta' \
        -fc_allelec 1 -fc_exch 1 -fc_disp 1 -fc_mu 1 -fc_quad 1
    run halgren-train train_ff $ga -max_generations 4 -max_iter 10 -pop_size 8 \
        -ff "$tff/halgren-PG.xml" -mp "$tff/meth-eth-water.xml" -sel "$tff/meth.dat" \
        -fit 'eta chi sigma epsilon' -fc_allelec 1 -fc_exch 1 -fc_disp 1
    run ljpc-intra train_ff $ga -max_generations 2 -max_iter 4 -pop_size 8 \
        -ff "$tff/lj-pc.xml" -mp "$tff/alc-mono.xml" -sel "$tff/mono.dat" \
        -fit 'De b c bondenergy kb bondlength angle kt c1 c2 c3 c4 c5' -fc_epot 1 -fc_force 0.1
    run pol-train train_ff $ga -max_generations 2 -max_iter 2 -pop_size 8 \
        -ff "$tff/alc-pol-wbh.xml" -mp "$tff/meth-eth-water.xml" -sel "$tff/meth.dat" \
        -fit 'zeta charge alpha eta chi sigma epsilon' -fc_allelec 1 -fc_exch 1 -fc_disp 1
    run bk3-b2 b2 -ff "$actdata/LegacyForceFields/BK3.xml" -f "$b2d/water#water.pdb" $b2
    run tip5p-b2 b2 -ff "$actdata/LegacyForceFields/TIP5P.xml" -f "$b2d/water#water.pdb" $b2
    run hf-b2 b2 -ff "$tff/4S-lj14-7-pg-hvdw.xml" -f "$b2d/hydrogen-fluoride#hydrogen-fluoride.pdb" $b2
}

# Run the regression set with and without explicit SIMD, and keep or compare it
regressionSet()
{
    local act=$build/bin/alexandria
    for mode in simd nosimd; do
        echo "=== Regression set, $mode"
        local out=$build/verify/$mode ref=$build/verify-reference/$mode
        rm -rf "$out"
        mkdir -p "$out"
        local nfail=${#failures[@]}
        noSimd=""
        [ $mode = nosimd ] && noSimd=1
        pushd "$out" > /dev/null
        runCases "$act"
        popd > /dev/null
        if [ "$regression" = save ] && [ ${#failures[@]} -gt $nfail ]; then
            echo "    not kept as the reference, since runs failed"
        elif [ "$regression" = save ]; then
            rm -rf "$ref"
            mkdir -p "$(dirname "$ref")"
            cp -R "$out" "$ref"
            echo "    kept as the reference"
        elif ! python3 "$root/admin/verify_compare.py" --tolerance "$tolerance" "$ref" "$out"; then
            failures+=("regression $mode")
        fi
    done
}

testBuild "$build" "$(basename "$build")"

if [ $novec -eq 1 ]; then
    nv=${build}_NOVEC
    if [ ! -f "$nv/CMakeCache.txt" ]; then
        cxx=$(cacheValue "$build" CMAKE_CXX_COMPILER)
        if "$cxx" --version 2>&1 | grep -qi clang; then
            noVecFlags="-fno-vectorize -fno-slp-vectorize"
        else
            noVecFlags="-fno-tree-vectorize"
        fi
        echo "=== Configuring $nv"
        "$cmake" -S "$root" -B "$nv" \
            -DCMAKE_BUILD_TYPE="$(cacheValue "$build" CMAKE_BUILD_TYPE)" \
            -DCMAKE_C_COMPILER="$(cacheValue "$build" CMAKE_C_COMPILER)" \
            -DCMAKE_CXX_COMPILER="$cxx" \
            -DCMAKE_PREFIX_PATH="$(cacheValue "$build" CMAKE_PREFIX_PATH)" \
            -DACT_DOUBLE="$(cacheValue "$build" ACT_DOUBLE)" \
            -DACT_MPI="$(cacheValue "$build" ACT_MPI)" \
            -DACT_OPENMP="$(cacheValue "$build" ACT_OPENMP)" \
            -DACT_SIMD=NONE \
            -DCMAKE_C_FLAGS="$noVecFlags" -DCMAKE_CXX_FLAGS="$noVecFlags" \
            > "$root/verify-configure.log" 2>&1 ||
            { tail -30 "$root/verify-configure.log"; failures+=("configure $nv"); }
        mv -f "$root/verify-configure.log" "$nv/" 2>/dev/null
    fi
    [ -f "$nv/CMakeCache.txt" ] && testBuild "$nv" "$(basename "$nv"), no vectorization"
fi

if [ $asan -eq 1 ]; then
    testBuild "$root/build_ASAN_DOUBLE" "build_ASAN_DOUBLE"
fi

if [ $single -eq 1 ]; then
    testBuild "$root/build_Release_SINGLE" "build_Release_SINGLE"
fi

if [ -n "$regression" ]; then
    regressionSet
fi

if [ ${#failures[@]} -gt 0 ]; then
    echo "=== FAILED:"
    printf '    %s\n' "${failures[@]}"
    exit 1
fi
echo "=== All passed"
