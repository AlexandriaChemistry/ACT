#!/bin/bash
# Build script of the conda package alexandria-python, see recipe.yaml: the
# Python package act with the module act._core, which has the classes of ACT.
set -euxo pipefail

case "${target_platform}" in
    osx-arm64|linux-aarch64) SIMD=ARM_NEON ;;
    linux-64)                SIMD=AVX2 ;;
    *)                       SIMD=NONE ;;
esac

# The libraries of ACT are static here and go into the module, such that the
# module does not depend on how the libraries of the package alexandria were
# built, and loads nothing of a GPU when it is imported. It uses MPI as the
# libraries do, so the package is for one MPI library, like alexandria.
cmake ${CMAKE_ARGS} --compile-no-warning-as-error -S . -B build-python \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DPython3_EXECUTABLE="${PREFIX}/bin/python" \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_TESTING=OFF \
    -DACT_MPI=ON \
    -DACT_SIMD=${SIMD} \
    -DACT_GPU=NONE \
    -DACT_OPENMP=OFF \
    -DACT_DOUBLE=$([[ "${ACT_PRECISION:-single}" == double ]] && echo ON || echo OFF) \
    -DACT_PYTHON=ON \
    -DACT_INSTALL_PYTHON_PACKAGE=ON \
    -DACT_BUILD_TARBALL=ON

# Only the module and what it needs, not the programs
cmake --build build-python --parallel ${CPU_COUNT} --target act_python_core

if [[ "${target_platform}" == linux-* ]]; then
    bash "${RECIPE_DIR}/glibc-check.sh" build-python/python/act/_core*.so
fi

# Only the Python package: the rest is in the package alexandria
cmake --install build-python --component python

if [[ "${CONDA_BUILD_CROSS_COMPILATION:-0}" != "1" ]]; then
    # The tests of the classes with the installed package and the program of
    # the package alexandria. For UCX_TLS see build.sh.
    SITE=$("${PREFIX}/bin/python" -c "import sysconfig; print(sysconfig.get_path('purelib'))")
    TESTS=$(pwd)/python/tests
    TESTDATA=$(pwd)/src/act/testdata
    (cd "$(mktemp -d)" && \
     UCX_TLS=self,sm,tcp \
     ACT_PYTHON_PATH="${SITE}" \
     ACT_ALEXANDRIA="${PREFIX}/bin/alexandria" \
     ACT_TESTDATA="${TESTDATA}" \
     ACT_REQUIRE_CLASSES=1 \
     "${PREFIX}/bin/python" "${TESTS}/test_core.py" -v)
fi
