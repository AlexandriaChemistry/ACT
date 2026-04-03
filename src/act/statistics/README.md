# statistics

This directory provides basic statistical analysis utilities used across the ACT toolkit for validating force-field parameters and summarising property datasets.

## Contents

- **statistics.h / statistics.cpp** — `Statistics` class: accumulates a series of numerical observations and computes summary statistics including mean, standard deviation, root-mean-square deviation (RMSD), minimum, maximum, and linear-regression metrics (R², slope, intercept). Used throughout ACT to report agreement between computed and reference molecular properties.
