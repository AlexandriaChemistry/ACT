#!/bin/bash

# Shell script for macOS and Linux
# Can use variables such as ARCH (architecture) and PREFIX (build environment prefix)

# uncompress build
tar -xf ${SRC_DIR}/act_osx-64.tar.gz -C ${SRC_DIR}

# bin
mkdir -p ${PREFIX}/bin
cp ${SRC_DIR}/bin/* ${PREFIX}/bin

# include
mkdir -p ${PREFIX}/include
cp -R ${SRC_DIR}/include/* ${PREFIX}/include

# share
mkdir -p ${PREFIX}/share
cp -R ${SRC_DIR}/share/* ${PREFIX}/share

# lib64
mkdir -p ${PREFIX}/lib64
cp -R ${SRC_DIR}/lib64/* ${PREFIX}/lib64

# lib
mkdir -p ${PREFIX}/lib
rsync -a --exclude python* ${SRC_DIR}/lib/ ${PREFIX}/lib/

# Python wrappers/interfaces
mkdir -p ${SP_DIR}
cp -R ${SRC_DIR}/lib/python*/site-packages/* ${SP_DIR}
