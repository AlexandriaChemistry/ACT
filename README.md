Alexandria Chemistry Toolkit
============================

The Alexandria Chemistry Toolkit is a software package that will derive a complete force field for molecular simulation from datasets of quantum chemistry calculations. 

The paper describing this software is available at [Digital Discovery 4 (2025) 1925](https://doi.org/10.1039/D5DD00178A).

The ACT manual, containing installation instructions, user and reference documentation as well as developer information is [on-line here](https://alexandriachemistry.github.io/ACT/).

## Installation

ACT is available as the conda package `alexandria`, for Macs with Apple silicon and for Linux on x86-64 processors with AVX2 and on 64-bit ARM processors, with support for NVIDIA GPUs on Linux. With [Miniforge](https://conda-forge.org/download/) installed:

```
conda create -n act -c alexandriachemistry -c conda-forge alexandria
conda activate act
```

For other platforms, double precision or AMD GPUs, build ACT from source. The [installation chapter](https://alexandriachemistry.github.io/ACT/userman/installation.html) of the manual describes both ways.

ACT is being developed by

+ [David van der Spoel](https://github.com/dspoel) (team leader),
+ [Mohammad Ghahremanpour](https://github.com/mmghahremanpour) and
+ [Paul J. van Maaren](https://github.com/maaren)

with contributions from

+ [Julian Marrades](https://github.com/jCodingStuff),
+ [Marie-Madeleine Walz](https://github.com/MMW1),
+ [Julian K&ouml;ngeter](https://github.com/JulianKoengeter),
+ [A. Najla Hosseini](https://github.com/najla23),
+ [Kristian Kriz](https://github.com/kkriz26),
+ [Alfred T. Nordman](https://github.com/pastaalfredo)
+ [M&aring;ns Rosenbaum](https://github.com/MansRosenbaum) and
+ [Oskar Tegby](https://github.com/OskarTegby)

In addition there are many past contributions to the [GROMACS](https://www.gromacs.org) part of the code by developers from that software package.


