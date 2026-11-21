#!/bin/bash
# Build script of the conda package alexandria, see recipe.yaml
set -euxo pipefail

case "${target_platform}" in
    osx-arm64|linux-aarch64) SIMD=ARM_NEON ;;
    linux-64)                SIMD=AVX2 ;;
    *)                       SIMD=NONE ;;
esac

GPU=(-DACT_GPU=NONE)
if [[ "${ACT_GPU_VARIANT:-none}" == sycl_* ]]; then
    GPU=(-DACT_GPU=SYCL -DAdaptiveCpp_DIR="${PREFIX}/lib/cmake/AdaptiveCpp")
    if [[ "${target_platform}" == linux-* ]]; then
        # acpp compiles and links libactGpu with the clang of AdaptiveCpp,
        # which has no sysroot of its own and would take the headers and the
        # C library of the build machine (glibc 2.39 on the GitHub runners),
        # so that the package does not run on systems with an older glibc.
        # Give it the sysroot and the GCC (libstdc++) of the conda compiler.
        GCC_DIR=$(dirname "$(${CXX} -print-libgcc-file-name)")
        GPU+=("-DACPP_EXTRA_ARGS=--sysroot=$(${CXX} -print-sysroot) --gcc-install-dir=${GCC_DIR}")
    fi
fi

# The Python package act is not installed: it would go to the site-packages
# of the Python of the build, while this package does not depend on a Python
# version. It is in the package alexandria-python, see build-python.sh.
# Warnings are errors in the ACT build, which is good for development, but a
# newer compiler or library should not stop the package from building.
cmake ${CMAKE_ARGS} --compile-no-warning-as-error -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DPython3_EXECUTABLE="${PREFIX}/bin/python" \
    -DBUILD_SHARED_LIBS=ON \
    -DACT_MPI=ON \
    -DACT_SIMD=${SIMD} \
    "${GPU[@]}" \
    -DACT_OPENMP=OFF \
    -DACT_DOUBLE=OFF \
    -DACT_INSTALL_PYTHON_PACKAGE=OFF \
    -DACT_BUILD_TARBALL=ON

cmake --build build --parallel ${CPU_COUNT}

if [[ "${target_platform}" == linux-* ]]; then
    # No library may need a newer glibc than that of the sysroot, which the
    # package requires (__glibc)
    bash "${RECIPE_DIR}/glibc-check.sh" build/lib/*.so* build/bin/alexandria
fi

if [[ "${CONDA_BUILD_CROSS_COMPILATION:-0}" != "1" ]]; then
    # The test programs are not part of the default target
    cmake --build build --parallel ${CPU_COUNT} --target tests
    # MPICH of conda-forge uses UCX, which may fail to open the network
    # devices of a build machine, for instance on the GitHub runners. The tests
    # run on one machine, which needs only shared memory and TCP.
    (cd build && UCX_TLS=self,sm,tcp ctest --output-on-failure --parallel ${CPU_COUNT})
fi

cmake --install build
