Installing ACT
==============
The Alexandria Chemistry Toolkit (ACT) relies on a number of libraries. Even though we tried to keep it to a minimum,
some more or less standard libraries are needed. ACT should compile fine on any UNIX (including MacOs) or Linux machine (but no warranty!).

Prerequisites
-------------
+ A patched version of OpenBabel is needed that can be found at
https://github.com/dspoel/openbabel. The oficial OpenBabel can be found [here](https://github.com/openbabel).

+ Some version of a library that supports the message passing interface (MPI) for parallel programming. A popular version is [Open MPI](https://open-mpi.org).

+ The [SQLite](https://www.sqlite.org/index.html) database engine is needed to process experimental data as well as quantum chemistry data from the Alexandria Library, available from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1170597.svg)](https://doi.org/10.5281/zenodo.1170597).

+ The [Class Library for Numbers](https://www.ginac.de/CLN/) is used in an optional part of the code and can be omitted.

+ The [cmake](https://cmake.org) tools are needed for compiling the code, along with an up-to-date C++ compiler supporting [C++14](https://en.wikipedia.org/wiki/C++14) at least.

+ For linear algebra operations we use [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/).

+ The [LibXml2](http://xmlsoft.org) is needed for processing the [XML](https://en.wikipedia.org/wiki/XML) data files used by the ACT.

+ You will also need [Python])https://www.python.org), version 3, and a number of Python libraries.

Install ACT on your computer
----------------------------
The easiest way to get going is to fetch the [download+install.py](download+install.py) script to your working directory of choice. Then, start by executing
```./download+install.py -h```
Let's say you have a four core machine and you want to install first time around, then this command should do the trick:
```./download+install.py -clone -ncore 4

