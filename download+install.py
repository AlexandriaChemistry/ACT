#!/usr/bin/env python3

import os, shutil, argparse, subprocess, sys, glob

debug = True

def get_host():
    '''Utility to get a string describing the build system'''
    if "HOST" in os.environ:
        HOST = os.environ["HOST"]
    elif "SNIC_SITE" in os.environ:
        HOST = os.environ["SNIC_SITE"]
    else:
        HOST = sys.platform
    return HOST

def get_compilers(ostype: str):
    '''Utility to return the C and C++ compiler respectively'''
    if ostype=="darwin":
        CXX = shutil.which("clang++")
        CC  = shutil.which("clang")
    else:
        CXX = None
        for cxxtry in [ "mpicxx", "cxx", "ccc", "g++", "gcc" ]:
            CXX = shutil.which(cxxtry)
            if CXX:
                break
        for cctry in [ "mpicc", "cc", "gcc" ]:
            CC  = shutil.which(cctry)
            if CC:
                break
    if not CC or not CXX:
        sys.exit("Cannot find compilers")
    return CC, CXX

def get_mpirun():
    for mrun in [ "srun", "mpirun" ]:
        mpirun = shutil.which("srun")
        if mpirun:
            return mpirun
    return None

def install_openbabel(anonymous, clone, destination, prefix_path, build_type, CXX, CC, ncores, HOST):
    # Get working directory
    pwd = os.getcwd()
    # Do we need to clone the repo?
    if clone:
        if anonymous:
            os.system("git clone https://github.com/dspoel/openbabel.git")
        else:
            os.system("git clone git@github.com:dspoel/openbabel.git")
    # See if we succeeded, if so let's go
    swdir = "openbabel"
    if not os.path.isdir(swdir):
        sys.exit("No directory %s. You may want to use the 'clone' flag" % swdir)
    os.chdir(swdir)
    
    # To to start the build
    bdir = "build"
    os.makedirs(bdir, exist_ok=True)
    if os.path.isdir(bdir):
        os.chdir(bdir)
        # Get the correct branch
        os.system("git checkout alexandria")
        FLAGS = ( "-DBUILD_SHARED=ON -DCMAKE_PREFIX_PATH=%s -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON -DOPTIMIZE_NATIVE=OFF -DWITH_COORDGEN=OFF  -DWITH_MAEPARSER=OFF -DCMAKE_CXX_COMPILER=%s -DCMAKE_C_COMPILER=%s -DCMAKE_BUILD_TYPE=%s -DCMAKE_INSTALL_PREFIX=%s" % ( prefix_path, CXX, CC, build_type, destination ))
        # Get rid of history, if any and run cmake and make
        if os.path.exists("CMakeCache.txt"):
            os.unlink("CMakeCache.txt")
        os.system("cmake %s .. &> cmake.log" % FLAGS)
        os.system("make -j %d install &> make.log" % ncores)
        
        # Check the result
        obabel = ("%s/bin/obabel" % destination)
        if not os.path.exists(obabel):
            sys.exit("Could not build %s, check cmake.log and make.log in %s" % ( obabel, os.getcwd() ) )
        else:
            print("Succesfully installed OpenBabel in %s/bin. Please add to your PATH variable" % destination)
        # Some hacking for MacOS only
        if HOST == "darwin":
            oblib = ( "%s/lib/openbabel" % destination )
            if os.path.exists(oblib):
                os.chdir(oblib)
                for obversion in glob.glob("*"):
                    os.chdir(obversion)
                    for lib in glob.glob("*.so"):
                        for rpath in [ "libz.1.dylib", "libcairo.2.dylib" ]:
                            os.system("install_name_tool -change @rpath/%s %s/lib/%s %s" % ( rpath, SW, rpath, lib) )
                    os.chdir("..")
    # Go back to where we came from
    os.chdir(pwd)

def install_gmx(args, CXX, CC, HOST):
    # Get working directory
    pwd = os.getcwd()
    act = args.branch == "main" or args.branch == "david"
    if act:
        if args.clone:
            if args.anonymous:
                os.system("git clone https://github.com/dspoel/ACT.git")
            else:
                os.system("git clone git@github.com:dspoel/ACT.git")
        swdir = "ACT"
    else:
        os.makedirs(args.branch, exist_ok=True)
        os.chdir(args.branch)
        if args.clone:
            if args.anonymous:
                os.system("git clone https://gitlab.com:/gromacs/gromacs.git")
            else:
                os.system("git clone git@gitlab.com:gromacs/gromacs.git")
        swdir = "gromacs"

    if not os.path.isdir(swdir):
        print("No directory %s. Maybe you want to use the -clone flag" % swdir)
        exit(1)
    os.chdir(swdir)
    if ((act and args.branch != "main") or 
        (not act and args.branch != "master")):
        os.system("git checkout --track origin/%s -b %s" % ( args.branch, args.branch ) )
        os.system("git pull")

    extra_dirs = []
    LBFLAGS    = ""
    mpirun     = get_mpirun()
    if HOST == "darwin":
        # MacOS machines
        LAPACK   = ( "%s/lib/liblapack.dylib" % args.prefix)
        BLAS     = ( "%s/lib/libblas.dylib" % args.prefix)
        LBFLAGS  = ( "-DGMX_BLAS_USER=%s -DGMX_LAPACK_USER=%s" % ( BLAS, LAPACK ) )
    elif HOST.find("nsc") >= 0:
        LAPACK = "/software/sse/easybuild/prefix/software/ScaLAPACK/2.0.2-gompi-2018a-OpenBLAS-0.2.20/lib/libscalapack.a" 
        BLAS   = "/software/sse/easybuild/prefix/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/libopenblas.so.0"
        LBFLAGS = ( "-DGMX_BLAS_USER=%s -DGMX_LAPACK_USER=%s" % ( BLAS, LAPACK ) )
    elif HOST.find("hpc2n") >= 0:
#        LAPACK = "/usr/lib/lapack/liblapack.so.3"
        LAPACK = "/hpc2n/eb/software/ScaLAPACK/2.1.0-gompi-2020b-bf/lib/libscalapack.so"
        BLAS   = "/hpc2n/eb/software/OpenBLAS/0.3.12-GCC-10.2.0/lib/libopenblas.so"
        LBFLAGS = ( "-DGMX_BLAS_USER=%s -DGMX_LAPACK_USER=%s" % ( BLAS, LAPACK ) )
    elif HOST.find("csb") >= 0:
        extra_dirs = []
        cinc = "CPLUS_INCLUDE_PATH"
        if cinc in os.environ:
            for inc in os.environ[cinc].split(":"):
                # Add the directory excluding the /include part
                extra_dirs.append(inc[:-8])
    else:
        sys.exit("Don't know how to commpile on host %s" % HOST)

    # The prefix path where to look for libraries
    PPATH = args.prefix
    for ed in [ args.destination ] + extra_dirs:
        PPATH = PPATH + ";" + ed

    # Construct cmake flags
    FLAGS = ("-DMPIEXEC=%s -DMPIEXEC_NUMPROC_FLAG='-n' -DGMX_X11=OFF -DGMX_LOAD_PLUGINS=OFF -DBUILD_SHARED_LIBS=OFF -DGMX_OPENMP=OFF -DGMX_MPI=ON -DGMX_GPU=OFF -DCMAKE_INSTALL_PREFIX=%s -DCMAKE_CXX_COMPILER=%s -DCMAKE_C_COMPILER=%s -DCMAKE_BUILD_TYPE=%s -DCMAKE_PREFIX_PATH='%s' -DGMX_BUILD_MANUAL=OFF -DGMX_COMPACT_DOXYGEN=ON -DREGRESSIONTEST_DOWNLOAD=OFF -DGMX_DEFAULT_SUFFIX=OFF -DGMX_LIBXML2=ON -DGMX_EXTERNAL_BLAS=ON -DGMX_EXTERNAL_LAPACK=ON %s" % ( mpirun, args.destination, CXX, CC, args.build, PPATH, LBFLAGS ) )
    if args.cln:
        FLAGS += " -DGMX_CLN=ON"
    # Make correct build directory
    bdir = "build_" + args.build
    if not args.single:
        FLAGS += " -DGMX_DOUBLE=ON"
        bdir = bdir + "_DOUBLE"
    os.makedirs(bdir, exist_ok=True)
    if debug:
        print("install_gmx: FLAGS = %s" % FLAGS)
    os.chdir(bdir)
    # Clean up old CMake stuff
    cmc = "CMakeCache.txt"
    if os.path.exists(cmc):
        os.unlink(cmc)
    # Run cmake and make!
    os.system("cmake %s .. &> cmake.log" % FLAGS)
    os.system("make -j %d install tests &> make.log" % args.ncores)
    # Check the result
    alex = ("%s/bin/alexandria" % args.destination)
    if not os.path.exists(alex):
        sys.exit("Could not build '%s', check cmake.log and make.log in %s" % ( alex, os.getcwd() ) )
    else:
        print("Succesfully installed ACT in %s/bin. Please add to your PATH variable" % args.destination)
    
    # Go back to where we came from
    os.chdir(pwd)
    
def parseArguments():
    parser = argparse.ArgumentParser()
    branch = "main"  
    parser.add_argument("-b", "--branch", help="Branch to clone from either ACT or GROMACS depending on branch name. Default is ACT with branch "+branch,             type=str, default=branch)
    parser.add_argument("-f", "--flags",  help="Additional compilation flags",type=str, default="")
    parser.add_argument("-anon", "--anonymous",  help="If you do not have an account at gitlab.com (gromacs) or github.com (ACT, OpenBabel) select this option.", action="store_true")
    parser.add_argument("-clone", "--clone", help="Clone git repository",     action="store_true")
    parser.add_argument("-ncores","--ncores", help="Number of cores for compiling", type=int, default=8)
    parser.add_argument("-bt", "--build", help="Build type for cmake. Typically Release or Debug, but Profile and ASAN (Adress Sanitizer) are available as well.", type=str, default="Release")
    parser.add_argument("-single", "--single", help="Single precision build (will not work well)", action="store_true")
    parser.add_argument("-cln", "--cln", help="Use the Class Library for Numbers", action="store_true")
    parser.add_argument("-ob", "--openbabel", help="Install the correct openbabel version", action="store_true")
    home = "HOME"
    if home in os.environ:
        dest = os.environ["HOME"]
    else:
        dest = "."
    dest += "/tools"
    parser.add_argument("-dest", "--destination", help="Where to install the software. Default: "+dest , type=str, default=dest)
    conda = shutil.which("conda")
    if conda and len(conda) > 10:
        prefix = conda[:-10]
    else:
        prefix = ""
    parser.add_argument("-prefix", "--prefix", help="Directory where libraries are installed, sometimes you may need to explicitly use -prefix ''. Default: "+prefix, type=str, default=prefix)
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args    = parseArguments()
    
    HOST    = get_host()
    CC, CXX = get_compilers(HOST)
    SW      = args.prefix
    if args.openbabel:
        install_openbabel(args.anonymous, args.clone, args.destination, args.prefix,
                          args.build, CXX, CC, args.ncores, HOST)
    
    install_gmx(args, CXX, CC, HOST)
