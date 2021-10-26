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

In the ```examples``` catalogue you will find a force field file, ```ACS-pg.dat```, and three subdirectories. Feel free to explore the files and scripts, but to run the example, please take the following steps:

```cd XML; ./extract_xml.py; make -j 8 alcohol.dat; cd ..```

Then, head over to the ```TUNE_EEM``` catalogue and run the example optimization using:

```cd TUNE_EEM; ./run_alcohol.py```

the adventure has begun. Inspect the different ```.xvg``` output files using a graphing program, and the ```tune_eem.log``` file using a text viewer.
