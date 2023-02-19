#!/bin/bash

# Shell script for macOS and Linux
# Can use variables such as ARCH (architecture) and PREFIX (build environment prefix)

tar -xf act_darwin64.tar.gz

ls

# Python wrappers/interfaces
mkdir -p ${SP_DIR}
cp -R ${SRC_DIR}/lib/python3.11/site-packages/act ${SP_DIR}/act
cp -R ${SRC_DIR}/lib/python3.11/site-packages/openbabel ${SP_DIR}/openbabel
