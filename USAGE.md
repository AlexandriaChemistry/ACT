Using the Alexandria Chemistry Toolkit
======================================
To use the ACT, you first need to [install it](INSTALL.md). Once you have finished that, please try the command

```alexandria -h```

and

```alexandria help commands```

Note that one can get detailed help for the alexandria modules using the ```-h``` flag, e.g.:

```alexandria gentop -h```

Try it!

Optimizing your first force field parameters
--------------------------------------------
You will need to download the Alexandria Library 
from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1170597.svg)](https://doi.org/10.5281/zenodo.1170597), unpack it, and set an environment variable 
```AlexandriaLib``` pointing to the root directory of the library. For instance,if you unpacked the Alexandria Library in your home catalogue:

```export AlexandriaLib=$HOME/AlexandriaLib```

Then, copy the examples directory from the ACT over to a working directory, for instance if the ACT source is in your home catalogue:

```cp -r $HOME/ACT/examples .; cd examples```

In the ```examples``` catalogue you will find a force field file, ```ACS-pg.dat```, and three subdirectories. Feel free to explore the files and scripts, but to run the example, please take the following steps in the ```XML``` directory:

```./extract_xml.py; make -j 8 alcohol.dat```

Then, head over to the ```TUNE_EEM``` catalogue and run the example optimization using:

```./run_alcohol.py```

the adventure has begun. Inspect the different ```.xvg``` output files using a graphing program, and the ```tune_eem.log``` file using a text viewer.

Running optimizations using parallel processing
-----------------------------------------------
The optimization process is heavy in computer time. Therefore, it has been parallellized using the MPI library. Since this is q prerequisite for compiling the ACT, you likely have it installed if you got this far. To enable it, just add ```mpirun -np 2``` before the ```alexandria``` command in the example script. As far as we have tested, the code is quite efficient down to 4-5 molecules per core. In the above example there are 10 compounds in the training set, such that it is not worthwhile to use more than 2-3 cores, but do experiment with the number of cores.

