#!/usr/bin/env python3

import os, shutil, argparse

def parseArguments():
    parser = argparse.ArgumentParser()
    branch = "main"  
    parser.add_argument("-b", "--branch", help="Branch to clone from either ACT or GROMACS depending on branch name. Default is ACT with branch "+branch,             type=str, default=branch)
    parser.add_argument("-f", "--flags",  help="Additional compilation flags",type=str, default="")
    parser.add_argument("-u", "--user",   help="Account name at gerrit.gromacs.org", type=str, default="dspoel")
    parser.add_argument("-clone", "--clone", help="Clone git repository",     action="store_true")
    parser.add_argument("-ncores","--ncores", help="Number of cores for compiling", type=int, default=8)
    parser.add_argument("-bt", "--build", help="Build type for cmake. Typicall Release or Debug", type=str, default="Release")
    parser.add_argument("-sngl", "--single", help="Single precision build (will not work well)", action="store_true")
    parser.add_argument("-cln", "--cln", help="Use the Class Library for Numbers", action="store_true")
  
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parseArguments()  
    act  = args.branch == "main" or args.branch == "david"
    if act:
        if args.clone:
            os.system("git clone git@github.com:dspoel/ACT.git")
        swdir = "ACT"
    else:
        os.makedirs(args.branch, exist_ok=True)
        os.chdir(args.branch)
        if args.clone:
            os.system("git clone git@gitlab.com:gromacs/gromacs.git")
        swdir = "gromacs"
    if not os.path.isdir(swdir):
        print("No directory %s. Maybe you want to use the -clone flag" % swdir)
        exit(1)
    os.chdir(swdir)
    if ((act and args.branch != "main") or 
        (not act and args.branch != "master")):
        os.system("git checkout --track origin/%s -b %s" % (args.branch, args.branch ) )
        os.system("git pull")

    if "HOST" in os.environ:
        HOST = os.environ["HOST"]
    elif "SNIC_SITE" in os.environ:
        HOST = os.environ["SNIC_SITE"]
    else:
        print("Don't know how to work on %s" % os.environ["HOST"] )
        exit(1)
    extra_dirs = []
    ROOT = os.environ["HOME"]
    DEST   = ( "%s/%s-%s/" % ( ROOT, swdir, args.branch ) )

    if HOST == "BirdMac":
        LAPACK = "/usr/lib/liblapack.dylib"
        BLAS   = "/usr/lib/libblas.dylib"
        LBFLAGS = ( "-DGMX_BLAS_USER=%s -DGMX_LAPACK_USER=%s" % ( BLAS, LAPACK ) )
    elif HOST.find("tetralith") >= 0:
        ROOT   = ROOT + "/wd"
        LAPACK = "/software/sse/easybuild/prefix/software/ScaLAPACK/2.0.2-gompi-2018a-OpenBLAS-0.2.20/lib/libscalapack.a" 
        BLAS   = "/software/sse/easybuild/prefix/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/libopenblas.so.0"
        LBFLAGS = ( "-DGMX_BLAS_USER=%s -DGMX_LAPACK_USER=%s" % ( BLAS, LAPACK ) )
    elif HOST.find("csb") >= 0:
        extra_dirs = []
        for libs in [ "LIBXML2", "OPENBLAS" ]:
            if libs in os.environ:
                extra_dirs.append(os.environ[libs])
        LBFLAGS    = ""

    PPATH  = ( "%s/GG/openbabel-alexandria/install" % ROOT )
    for ed in extra_dirs:
        PPATH = PPATH + ";" + ed

    CXX    = shutil.which("mpicxx")
    CC     = shutil.which("mpicc")
    gmxdouble = ""
    if not args.single:
        gmxdouble = "-DGMX_DOUBLE=ON"
    FLAGS = ("-DMPIEXEC=/usr/bin/srun -DMPIEXEC_NUMPROC_FLAG='-n' -DGMX_X11=OFF -DGMX_LOAD_PLUGINS=OFF -DBUILD_SHARED_LIBS=OFF -DGMX_OPENMP=OFF -DGMX_MPI=ON -DGMX_GPU=OFF -DCMAKE_INSTALL_PREFIX=%s -DCMAKE_CXX_COMPILER=%s -DCMAKE_C_COMPILER=%s -DCMAKE_BUILD_TYPE=%s -DCMAKE_PREFIX_PATH='%s' -DGMX_BUILD_MANUAL=OFF -DGMX_COMPACT_DOXYGEN=ON -DREGRESSIONTEST_DOWNLOAD=OFF %s -DGMX_DEFAULT_SUFFIX=OFF -DGMX_LIBXML2=ON -DGMX_EXTERNAL_BLAS=ON -DGMX_EXTERNAL_LAPACK=ON %s" % ( DEST, CXX, CC, args.build, PPATH, gmxdouble, LBFLAGS ) )
    if args.cln:
        FLAGS += " -DGMX_CLN=ON")

    bdir = "build_" + args.build
    if not args.single:
        bdir = bdir + "_DOUBLE"
    os.makedirs(bdir, exist_ok=True)
    print("FLAGS: %s" % FLAGS)
    os.chdir(bdir)
    os.system("cmake %s .. >& cmake.log" % FLAGS)
    os.system("make -j %d install tests >& make.log" % args.ncores)

