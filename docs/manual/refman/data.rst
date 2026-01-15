.. _sec-data:
        
*************
Training Data
*************

===================
Using existing data
===================
A recent review from the ACT developers~ :cite:p:`Kriz2023a` discusses the different available quantum chemistry data sets. Some of these can be used in the ACT, for instance the coupled-cluster dimer data set due to Donchev {\em et al}~ :cite:p:`Donchev2021a`, the Non-covalent interaction atlas from the Czech group led by \u{R}ez{\'a}\u{c}~ :cite:p:`Rezac2020a_NCIA` and some more.
The SAPT dataset on protein side-chain analogs and backbone analogs by Burns~ :cite:p:`Burns2017_biofragment` can in principle be used in the ACT as well.

In principle, the ANI-1 dataset~ :cite:p:`Smith2017a` can be used to provide off-equilibrium energies of small compounds, but in the first version it was limited to compound containing the elements C, H, N, O only.

==================
Alexandria Library
==================
The Alexandria Library contains energies at optimized conformations of about 5000 compounds, thermochemistry, electric multipoles and electrostatic potential at grid points around molecules~ :cite:p:`Ghahremanpour2017a,Ghahremanpour2018a`.
It is possible to train electrostatic models using data from the library.

================
Donchev data set
================
To use the dataset due to Donchev {\em et al} :cite:p:`Donchev2021a`, we provide a script that reads the comma-separated value file provided by those authors. Please refer to their article for download information. 
Then refer to the built-in help in the script for more guidance by executing::

  donchev2molprop -h

and investigating the output.

==============================
Non-covalent interaction atlas
==============================
The Non-covalent interaction atlas~ :cite:p:`Rezac2020a_NCIA` can be used in a similar manner. At the time of writing it can be \href{http://www.nciatlas.org}{found here}.
To convert the data to ACT files, start by::

  ncia2molprop -h

and investigate your options.

========
ACT data
========
Quantum chemical data used in ACT-related publications will be uploaded to a sharing site.

.. _sapt:

====================
Generating SAPT data
====================
A recent review described the available data sets available for machine learning or force field training~ :cite:p:`Kriz2023a`.
There is however a lack of certain data, or data sets are incomplete. For this reason there is a git repository `SaptACT`_ that provides script to run SAPT calculations using the Psi4 software~ :cite:p:`psi4` including tools to convert the output to ACT compatible inputs (see below). 
The steps needed to prepare data for training in ACT are as follows:

.. _SaptACT: https://github.com/AlexandriaChemistry/SaptACT

* Clone the SaptACT_ repository using *git*
* Prepare selection of compounds, e.g. water, methanol, ethanol, and make sure that monomer structures for these compounds are present in the *xyz/monomers* catalog.
* Generate a dimer selection file using the *gen_dimers.py* script in the SaptACT repository::

    ./gen_dimers.py -sel water methanol ethanol -o alcohol.dat

  which will generate a file containing::

    water#water
    water#methanol
    water#ethanol
    methanol#methanol
    methanol#ethanol
    ethanol#ethanol

  where the # symbol separates the monomers. Note that the order on the command line determines the order in the selection file.
* Perform SAPT calculations using Psi4 by generating randomly oriented dimers at a series of distances defined by scaling the shortest distance between any pair of atoms in the generated orientations. In this case, six distances varying from 0.9 to 1.4 times the sum of the van der Waals radii of the nearest atoms will be generated.::

     ./run_calcs.py -dimers alchol.dat -ndist 6 -mindist 0.9 -maxdist 1.4 -norient 10

  Note, that in the above command 6$\times$10 calculations are started for each dimer in {\tt alcohol.dat}, that is 360 calculations in total. You will obviously need a compute cluster to do these calculations, in particular if you select a more accurate SAPT level of theory~ :cite:p:`Parker2014a` than the default (sapt2+/aug-cc-pvdz). 

.. _qmdata:

===============================
Generating single molecule data
===============================
\label{qmdata}
Single molecules off-equilibrium energies and forces are needed to parameterize the intramolecular potential functions. Scripts are available that will take a structure of a monomer, perform a 50 ps MD simulation in the gas phase at elevated temperature using the GAFF force field~ :cite:p:`Wang2004a` and the GROMACS software~ :cite:p:`Pronk2013a`.  
Then,  conformations are extracted from the simulation trajectory and these are subjected to quantum chemistry calculations using Psi4~ :cite:p:`psi4`.

===============================
Conversion to ACT molprop files
===============================
The ACT uses the \href{https://en.wikipedia.org/wiki/XML{eXtensible Markup Language} to store both data for training and force field files. Once your SAPT calculations are finished you need to perform the following steps to generate the ACT input:

* Convert the Psi4 outputs to compact and human-readable \href{https://en.wikipedia.org/wiki/JSON}{json} files using another script in SaptACT::

     ./generate_json
    
  which will traverse the output directories, find Psi4 outputs, and if there is not already a {\tt results.json} file, will generate it. Please familiarize yourself with the content of these files.
* When the json files have been generated, you are ready to convert them to a molprop file::

    ./write_molprop.py -o molprop.xml
    
  please inspect the output file from this script to verify that your calculations are present.

Both these scripts have more flags that can be useful, please investigate those using the {\tt -h} option.

