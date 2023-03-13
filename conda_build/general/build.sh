#!/bin/bash

# Shell script for macOS and Linux
# Can use variables such as ARCH (architecture) and PREFIX (future environment prefix)

# Installing openbabel
echo "Installing openbabel..."

cd ${SRC_DIR}/openbabel
git checkout alexandria
mkdir ${SRC_DIR}/openbabel/build
cd ${SRC_DIR}/openbabel/build

[ ${ARCH} == "64" ] && DLIB_SUFFIX="-DLIB_SUFFIX=64" || DLIB_SUFFIX=""
# echo ${DLIB_SUFFIX}
cmake -S${SRC_DIR}/openbabel -B${SRC_DIR}/openbabel/build \
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
	..
# -DCMAKE_C_FLAGS=${CFLAGS} \
# -DCMAKE_CXX_FLAGS=${CXXFLAGS} \

(cd ${SRC_DIR}/openbabel/build && make -j ${CPU_COUNT} install)
