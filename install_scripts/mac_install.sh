#!/bin/bash

# Installs all dependencies for ACT and then ACT itself
# Should be executed from the home directory.
# Pre-requisites: conda
# Anaconda should be located at /opt/anaconda3

# Check if conda exists
if ! command -v conda
then
  echo "Please install Anaconda3 under /opt"
  echo "See https://www.anaconda.com/products/individual"
  exit 1
else
  echo "conda is installed, continuing run..."
fi

# cd to home directory
cd ~

# Create the <tools> directory
echo "If <tools> directory does not exist, creating it..."
mkdir -p tools

# Update conda
echo "Updating conda..."
conda update conda -y

# Create conda environment
echo "Creating the conda environment..."
conda create -n ACT python -y

# Activate the environment
echo "Activating the conda environment..."
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate ACT

# Install missing commands with conda
if ! command -v wget
then
  echo "Missing wget, installing..."
  conda install -c anaconda wget -y
else
  echo "wget requirement is satisfied..."
fi

if ! command -v git
then
  echo "Missing git, installing..."
  conda install -c anaconda git -y
else
  echo "git requirement is satisfied..."
fi

# Clang gets installed automatically when mpi is installed
# if ! command -v clang || ! command -v clang++
# then
#   echo "Missing clang/clang++, installing..."
#   conda install -c conda-forge clang clangxx clang-tools -y
# else
#   echo "clang/clang++ requirement is satisfied..."
# fi

# Install ACT dependencies with conda
echo "Installing ACT dependencies with conda..."
conda install -c anaconda cmake perl sqlite libxml2 cairo graphviz gmp -y
conda install -c conda-forge openmpi openmpi-mpicc openmpi-mpicxx boost eigen fftw liblapack libblas doxygen  -y

# Install the patched OpenBabel locally under /Users/tools/openbabel-install
echo "Installing patched OpenBabel..."
wget https://jcodingstuff.github.io/docs/build_openbabel_mac.py
chmod 755 build_openbabel_mac.py
./build_openbabel_mac.py
rm -f build_openbabel_mac.py
rm -rf openbabel-build
# git clone https://github.com/dspoel/openbabel
# cd openbabel
# git checkout alexandria
# cd ..
# mkdir openbabel-build
# cd openbabel-build
# cmake ../openbabel -DCMAKE_INSTALL_PREFIX=~/tools/openbabel-install
# make && make install
# cd ..
# rm -rf openbabel-build

# Install the Class Library for Numbers locally under
echo "Installing Class Library for Numbers..."
wget https://www.ginac.de/CLN/cln-1.3.6.tar.bz2
tar -xf cln-1.3.6.tar.bz2
rm -f cln-1.3.6.tar.bz2
cd cln-1.3.6
./configure --prefix=/Users/copernico/tools/cln-install
make
make check
make install
# make && make install
cd ..

# Run the python script for cloning and installation of ACT
echo "Fetching and running download_install_mac.py..."
wget https://jcodingstuff.github.io/docs/download_install_mac.py
chmod 755 download_install_mac.py
./download_install_mac.py -clone -cln
rm -f download_install_mac.py
