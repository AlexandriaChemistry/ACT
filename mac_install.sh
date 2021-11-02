#!/bin/bash

# Installs all dependencies for ACT and then ACT itself
# Should be executed from the home directory.
# Pre-requisites: conda, wget, git, and clang.
# Anaconda should be located at /opt/anaconda3

# Get user name
USERNAME=$(whoami)

# Create the <tools> directory
echo "Creating the <tools> directory"
mkdir tools

# Update conda
echo "Updating conda"
conda update conda -y

# Create conda environment
echo "Creating the conda environment"
conda create -n ACT python -y

# Activate the environment
echo "Activating the conda environment"
conda activate ACT

# Install dependencies with conda
echo "Installing dependencies with conda"
conda install -c anaconda cmake sqlite libxml2 graphviz gmp -y
conda install -c conda-forge openmpi openmpi-mpicc openmpi-mpixx eigen fftw liblapack libblas doxygen  -y

# Install the patched OpenBabel locally under /Users/tools/openbabel-install
echo "Installing patched OpenBabel"
git clone https://github.com/dspoel/openbabel
mkdir openbabel-build
cd openbabel-build
cmake ../openbabel -DCMAKE_INSTALL_PREFIX=~/tools/openbabel-install
make && make install
cd ..

# Install the Class Library for Numbers locally under
echo "Installing Class Library for Numbers"
wget https://www.ginac.de/CLN/cln-1.3.6.tar.bz2
tar -xf cln-1.3.6.tar.bz2
cd cln-1.3.6
./configure --prefix='/Users/${USERNAME}/tools/cln-install'
# make
# make check
# make install
make && make install
cd ..

# Run the python script for cloning and installation of ACT
echo "Fetching and running download_install_mac.py"
wget https://jcodingstuff.github.io/docs/download_install_mac.py
chmod 755 download_install_mac.py
python download_install_mac.py -clone -cln
rm -f download_install_mac.py
