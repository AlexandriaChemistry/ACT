#!/bin/bash

# Shell script for macOS and Linux
# Can use variables such as ARCH (architecture) and PREFIX (future environment prefix)

# Installing openbabel
echo "Installing openbabel..."

cd ${SRC_DIR}/openbabel
git checkout alexandria
mkdir ${SRC_DIR}/openbabel/build

[ ${ARCH} == "64" ] && DLIB_SUFFIX="-DLIB_SUFFIX=64" || DLIB_SUFFIX=""
# echo ${DLIB_SUFFIX}
cmake \
	-DBUILD_SHARED=ON \
	-DBUILD_GUI=OFF \
	-DCMAKE_PREFIX_PATH=${BUILD_PREFIX} \
	-DRUN_SWIG=ON \
	-DPYTHON_BINDINGS=ON \
	-DOPTIMIZE_NATIVE=OFF \
	-DWITH_COORDGEN=OF \
	-DWITH_MAEPARSER=OFF \
	-DCMAKE_CXX_COMPILER=${CXX_FOR_BUILD} \
	-DCMAKE_C_COMPILER=${CC_FOR_BUILD} \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_INSTALL_PREFIX=${PREFIX} \
	${DLIB_SUFFIX} \
	-S ${SRC_DIR}/openbabel -B ${SRC_DIR}/openbabel/build
# -DCMAKE_C_FLAGS=${CFLAGS} \
# -DCMAKE_CXX_FLAGS=${CXXFLAGS} \

(cd ${SRC_DIR}/openbabel/build && make -j ${CPU_COUNT} install)


# Installing ACT
echo "Installing ACT..."
mkdir ${SRC_DIR}/ACT/build

cmake \
	-DGMX_DOUBLE=ON \
	-DMPIEXEC=${BUILD_PREFIX}/bin/mpirun \
	-DMPIEXEC_NUMPROC_FLAG='-n' \
	-DGMX_EXTERNAL_BLAS=OFF \
	-DGMX_EXTERNAL_LAPACK=OFF \
	-DGMX_X11=OFF \
	-DBUILD_SHARED_LIBS=OFF \
	-DGMX_OPENMP=OFF \
	-DGMX_MPI=ON \
	-DCMAKE_INSTALL_PREFIX=${PREFIX} \
	-DCMAKE_CXX_COMPILER=${CXX_FOR_BUILD} \
	-DCMAKE_C_COMPILER=${CC_FOR_BUILD} \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_PREFIX_PATH=${BUILD_PREFIX} \
	-DGMX_BUILD_MANUAL=OFF \
	-DGMX_COMPACT_DOXYGEN=ON \
	-DREGRESSIONTEST_DOWNLOAD=OFF \
	-DGMX_DEFAULT_SUFFIX=OFF \
	-S ${SRC_DIR}/ACT -B ${SRC_DIR}/ACT/build

(cd ${SRC_DIR}/ACT/build && make -j ${CPU_COUNT} install tests)
