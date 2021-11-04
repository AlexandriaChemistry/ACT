#!/usr/bin/env python3

import os, shutil, glob, sys

DEST      = os.environ["HOME"] + "/tools/openbabel-install"
CXX       = shutil.which("clang++")
CC        = shutil.which("clang")
# SW        = os.environ["HOME"] + "/opt/miniconda3"
SW        = "/opt/anaconda3/envs/ACT"
FLAGS = ( "-DBUILD_SHARED=ON -DCMAKE_PREFIX_PATH=%s -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON -DOPTIMIZE_NATIVE=OFF -DWITH_COORDGEN=OFF  -DWITH_MAEPARSER=OFF -DCMAKE_CXX_COMPILER=%s -DCMAKE_C_COMPILER=%s -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%s" % ( SW, CXX, CC, DEST ))

ob = "openbabel"
if not os.path.exists(ob):
    os.system("git clone https://github.com/dspoel/openbabel")
    os.chdir(ob)
    os.system("git checkout alexandria")
    os.chdir("..")

do_build = True
hack = True
if os.path.exists(ob):
    if do_build:
        bdir = "openbabel-build"
        os.makedirs(bdir, exist_ok=True)
        os.chdir(bdir)
        if os.path.exists("CMakeCache.txt"):
            os.unlink("CMakeCache.txt")
        os.system("cmake %s ../%s >& cmake.log" % ( FLAGS, ob))
        os.system("make -j 8 install >& make.log")
        os.chdir("..")
    oblib = "tools/openbabel-install/lib/openbabel"
    if hack:
        if os.path.exists(oblib):
            print('Hacking...')
            os.chdir(oblib)
            for obversion in glob.glob("*"):
                os.chdir(obversion)
                for lib in glob.glob("*.so"):
                    for rpath in [ "libz.1.dylib", "libcairo.2.dylib" ]:
                        os.system("install_name_tool -change @rpath/%s %s/lib/%s %s" % ( rpath, SW, rpath, lib) )
                os.chdir("..")
        else:
            print(f"Not hacking... Path not found: {oblib}")
else:
    print("You started building from the wrong directory")
