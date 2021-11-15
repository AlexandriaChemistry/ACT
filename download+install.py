#!/usr/bin/env python3

import os, shutil, argparse, subprocess, sys, glob
from subprocess import Popen, PIPE

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
    if "darwin" in ostype:
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
        mpirun = shutil.which(mrun)
        if mpirun:
            return mpirun
    return None

def get_prefix(myprefix):
    cpp  = "CMAKE_PREFIX_PATH"
    cinc = "CPLUS_INCLUDE_PATH"
    prefix = []
    if cpp in os.environ:
        prefix = os.environ[cpp].split(":")
    elif cinc in os.environ:
        for inc in os.environ[cinc].split(":"):
            # Add the directory excluding the /include part
            prefix.append(inc[:-8])
    else:
        conda = shutil.which("conda")
        if conda and len(conda) > 10:
            prefix = [ conda[:-10] ]
    if len(prefix) == 0:
        return ""
    else:
        pp = prefix[0]
        if myprefix and len(myprefix) > 0:
            pp += ";" + myprefix
        for p in range(1, len(prefix)):
            pp += ";" + prefix[p]
        return pp
    
def find_lib(prefix, libname, verbose):
    if verbose:
        print('Entering find_lib...')
        print(f'prefix: {prefix}')
        print(f'libname: {libname}')
    for pp in prefix.split(";"):
        libtry = pp+"/lib/"+libname
        if os.path.exists(libtry):
            if verbose:
                print('Leaving find_lib...')
            return libtry
    libpath = "LIBRARY_PATH"
    if libpath in os.environ:
        for pp in os.environ[libpath].split(":"):
            libtry = pp+"/"+libname
            if os.path.exists(libtry):
                return libtry
    sys.exit("Cannot find %s" % libname)

def run_command(command, flags):
    cmdlist = [ shutil.which(command) ] + flags.split()
    with Popen(cmdlist, stdout=PIPE, stderr=PIPE, stdin=PIPE) as p:
        output, error = p.communicate()
        with open(command+".out", "w") as outf:
            for line in output.decode('utf-8').splitlines():
                outf.write("%s\n" % line)
        with open(command+".err", "w") as outf:
            for line in error.decode('utf-8').splitlines():
                outf.write("%s\n" % line)

def check_for_openbabel(destination):
    obabel = ("%s/bin/obabel" % destination)
    if not os.path.exists(obabel):
        sys.exit("Could not build %s, check cmake.log and make.log in %s" % ( obabel, os.getcwd() ) )
    return True
    
def install_openbabel(anonymous, clone, destination, prefix, build_type, CXX, CC, ncores, HOST, verbose):
    # Get working directory
    pwd = os.getcwd()
    swdir = "openbabel"
    # Do we need to clone the repo?
    if not os.path.isdir(swdir):
        if clone:
            if anonymous:
                os.system("git clone https://github.com/dspoel/openbabel.git")
            else:
                os.system("git clone git@github.com:dspoel/openbabel.git")
        # See if we succeeded, if so let's go
        if not os.path.isdir(swdir):
            sys.exit("No directory %s. You may want to use the '-clone_OB' flag" % swdir)
        os.chdir(swdir)
    else:
        print("Will not clone again since directory %s exists. Will git pull instead." % swdir)
        os.chdir(swdir)
        run_command("git", "pull")
    
    # To to start the build
    bdir = "build"
    os.makedirs(bdir, exist_ok=True)
    if os.path.isdir(bdir):
        os.chdir(bdir)
        # Get the correct branch
        os.system("git checkout alexandria")
        FLAGS = ( "-DBUILD_SHARED=ON -DBUILD_GUI=OFF -DCMAKE_PREFIX_PATH=%s -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON -DOPTIMIZE_NATIVE=OFF -DWITH_COORDGEN=OFF  -DWITH_MAEPARSER=OFF -DCMAKE_CXX_COMPILER=%s -DCMAKE_C_COMPILER=%s -DCMAKE_BUILD_TYPE=%s -DCMAKE_INSTALL_PREFIX=%s" % ( prefix, CXX, CC, build_type, destination ))
        if verbose:
            print("install_openbabel cmake FLAGS: %s" % FLAGS)
        # Get rid of history, if any and run cmake and make
        if os.path.exists("CMakeCache.txt"):
            os.unlink("CMakeCache.txt")
        run_command("cmake", FLAGS + " ..")
        run_command("make", (" -j %d install" % ncores))
        
        # Check the result
        if check_for_openbabel(destination):
            print("Succesfully installed OpenBabel in %s/bin. Please add to your PATH variable" % destination)
        # Some hacking for MacOS only
        if "darwin" in HOST:
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

def conda_prefix_provided(prefix: str):
    if prefix is not None:
        for path in prefix.split(';'):
            for ele in path.split('/'):
                if ele in ['anaconda3', 'miniconda3']:
                    return True
        return False
    else:
        return False

def install_act(args, CXX, CC, HOST, prefix):
    # Get working directory
    pwd = os.getcwd()
    act = args.branch == "main" or args.branch == "david"
    if act:
        swdir = "ACT"
        
        if args.clone_ACT:
            if os.path.isdir(swdir):
                print("Will not clone since directory %s exists. Will gitt pull instead")
                os.chdir(swdir)
                run_command("git", "pull")
            else:
                if args.anonymous:
                    os.system("git clone https://github.com/dspoel/ACT.git")
                else:
                    os.system("git clone git@github.com:dspoel/ACT.git")
                if not os.path.isdir(swdir):
                    sys.exit("No directory %s. Maybe you want to use the -clone_ACT flag" % swdir)
                os.chdir(swdir)
        else:
            os.chdir(swdir)
    else:
        os.makedirs(args.branch, exist_ok=True)
        os.chdir(args.branch)
        if args.clone_ACT:
            if args.anonymous:
                os.system("git clone https://gitlab.com:/gromacs/gromacs.git")
            else:
                os.system("git clone git@gitlab.com:gromacs/gromacs.git")
        swdir = "gromacs"
        os.chdir(swdir)

    if ((act and args.branch != "main") or 
        (not act and args.branch != "master")):
        os.system("git checkout --track origin/%s -b %s" % ( args.branch, args.branch ) )
        os.system("git pull")

    # The prefix path where to look for libraries
    # Delete anaconda from prefix if provided by args
    if conda_prefix_provided(args.prefix):
        split_prefix = [ele for ele in prefix.split(';')
                        if ('anaconda3' not in ele.split('/') and 'miniconda3' not in ele.split('/'))]
        prefix = ';'.join(split_prefix)
        prefix += f';{args.prefix}'
    prefix = prefix.lstrip(';')
    PPATH  = prefix + ";" + args.destination
    LAPACK = None
    BLAS   = None
    if "darwin" in HOST:
        # MacOS machines
        LAPACK   = find_lib(PPATH, "liblapack.dylib", args.verbose)
        BLAS     = find_lib(PPATH, "libblas.dylib", args.verbose)
    elif HOST.find("nsc") >= 0 or HOST.find("hpc2n") >= 0:
        LAPACK = "" #find_lib(PPATH, "libscalapack.a", args.verbose)
        BLAS   = "" #find_lib(PPATH, "libopenblas.so", args.verbose)
    else:
        print("No specific knowledge on, how to commpile on host %s, trying anyway." % HOST)
        
    # Construct cmake flags
    FLAGS = ("-DMPIEXEC=%s -DMPIEXEC_NUMPROC_FLAG='-n' -DGMX_X11=OFF -DBUILD_SHARED_LIBS=OFF -DGMX_OPENMP=OFF -DGMX_MPI=ON -DCMAKE_INSTALL_PREFIX=%s -DCMAKE_CXX_COMPILER=%s -DCMAKE_C_COMPILER=%s -DCMAKE_BUILD_TYPE=%s -DCMAKE_PREFIX_PATH='%s' -DGMX_BUILD_MANUAL=OFF -DGMX_COMPACT_DOXYGEN=ON -DREGRESSIONTEST_DOWNLOAD=OFF -DGMX_DEFAULT_SUFFIX=OFF -DGMX_LIBXML2=ON" % ( get_mpirun(), args.destination, CXX, CC, args.build, PPATH ) )
    if args.cln:
        FLAGS += " -DGMX_CLN=ON"
    if LAPACK:
        FLAGS += "  -DGMX_EXTERNAL_LAPACK=ON -DGMX_LAPACK_USER=" + LAPACK
    if BLAS:
        FLAGS += "  -DGMX_EXTERNAL_BLAS=ON -DGMX_BLAS_USER=" + BLAS

    # Make correct build directory
    bdir = "build_" + args.build
    if not args.single:
        FLAGS += " -DGMX_DOUBLE=ON"
        bdir = bdir + "_DOUBLE"
    if args.verbose:
        print("Will create %s" % bdir)
    os.makedirs(bdir, exist_ok=True)
    os.chdir(bdir)
    if args.verbose:
        print("install_act cmake FLAGS: %s in directory %s" % ( FLAGS, os.getcwd() ) )
    # Clean up old CMake stuff
    cmc = "CMakeCache.txt"
    if os.path.exists(cmc):
        os.unlink(cmc)
    # Run cmake and make!
    run_command("cmake", FLAGS + " ..")
    run_command("make", (" -j %d install tests" % args.ncores))
    # Check the result
    alex = ("%s/bin/alexandria" % args.destination)
    if not os.path.exists(alex):
        sys.exit("Could not build '%s', check cmake.* and make.* in %s" % ( alex, os.getcwd() ) )
    else:
        print("Succesfully installed ACT in %s/bin. Please add this dir to your PATH variable" % args.destination)
    
    # Go back to where we came from
    os.chdir(pwd)
    
def parseArguments():
    parser = argparse.ArgumentParser()
    branch = "main"  
    parser.add_argument("-b", "--branch", help="Branch to clone from either ACT or GROMACS depending on branch name. Default is ACT with branch "+branch,             type=str, default=branch)
    parser.add_argument("-f", "--flags",  help="Additional compilation flags",type=str, default="")
    parser.add_argument("-anon", "--anonymous",  help="If you do not have an account at gitlab.com (gromacs) or github.com (ACT, OpenBabel) select this option.", action="store_true")
    parser.add_argument("-clone_OB", "--clone_OB", help="Clone openbabel git repository", action="store_true")
    parser.add_argument("-clone_ACT", "--clone_ACT", help="Clone ACT git repository", action="store_true")
    parser.add_argument("-ncores","--ncores", help="Number of cores for compiling", type=int, default=8)
    parser.add_argument("-bt", "--build", help="Build type for cmake. Typically Release or Debug, but Profile and ASAN (Adress Sanitizer) are available as well.", type=str, default="Release")
    parser.add_argument("-single", "--single", help="Single precision build (will not work well)", action="store_true")
    parser.add_argument("-cln", "--cln", help="Use the Class Library for Numbers (optional)", action="store_true")
    parser.add_argument("-ob", "--openbabel", help="Install the correct openbabel version", action="store_true")
    home = "HOME"
    if home in os.environ:
        dest = os.environ["HOME"]
    else:
        dest = "."
    dest += "/tools"
    parser.add_argument("-dest", "--destination", help="Where to install the software. Default: "+dest , type=str, default=dest)
    parser.add_argument("-prefix", "--prefix", help="Directory where libraries are installed, default autodetect", type=str, default=None)
    parser.add_argument("-v", "--verbose", help="Print more stuff on the terminal", action="store_true")
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args    = parseArguments()
    
    HOST    = get_host()
    CC, CXX = get_compilers(HOST)
    SW      = get_prefix(args.prefix)
    if args.openbabel or args.clone_OB:
        install_openbabel(args.anonymous, args.clone_OB, args.destination, SW,
                          args.build, CXX, CC, args.ncores, HOST, args.verbose)
    # First check whether our OB is installed where it should be
    if check_for_openbabel(args.destination):
        install_act(args, CXX, CC, HOST, SW)
