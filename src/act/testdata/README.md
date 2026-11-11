# Test data

The input files of the tests of ACT. The tests find them with
`gmx::test::TestFileManager::getInputFilePath()`, which looks in the
directory of the test first and then here (`GMX_TESTSIMULATIONDATABASE_DIR`
in `src/buildinfo.h.cmakein`), so a test refers to a file as, for instance,
`forcefields/ACS-g1.xml`. The reference data of the tests stay in the
`refdata` directories next to the tests.

| Directory     | Contents |
|---------------|----------|
| `forcefields` | Force fields (`alexandria_chemistry` XML files), and the CSV input of `gen_ff`. `getForceField("name")` reads `forcefields/name.xml`. |
| `molprop`     | Molprop files (`molecules` XML files) and their selection files (`.dat`), and `make_subsets.py`, which made the small molprop files of the `train_ff` tests. |
| `molecules`   | Structures (SDF, PDB, XYZ) and Gaussian output files. |
| `other`       | Gene pools of the GA tests, a table, and GROMACS files that no test uses. |

## Renamed files

The files used to be in the test directories, where some different files
had the same name. They got a number:

| Now | Was |
|-----|-----|
| `forcefields/ACS-g1.xml` | `src/act/alexandria/tests/ACS-g.xml` |
| `forcefields/ACS-g2.xml` | `src/act/alexandria/tests/ACS-g2.xml` (unchanged) |
| `forcefields/ACS-g3.xml` | `src/act/ga/tests/ACS-g.xml` |
| `forcefields/ACS-g4.xml` | `src/act/qgen/tests/ACS-g.xml` |
| `forcefields/ACS-pg1.xml` | `src/act/alexandria/tests/ACS-pg.xml` |
| `forcefields/ACS-pg2.xml` | `src/act/qgen/tests/ACS-pg.xml` |
| `forcefields/ACM-g1.xml` | `src/act/forcefield/tests/ACM-g.xml` |
| `forcefields/ACM-g2.xml` | `src/act/qgen/tests/ACM-g.xml` |
| `forcefields/ACM-pg1.xml` | `src/act/forcefield/tests/ACM-pg.xml` |
| `forcefields/ACM-pg2.xml` | `src/act/qgen/tests/ACM-pg.xml` |
| `forcefields/water_hyperpol1.xml` | `src/act/alexandria/tests/water_hyperpol.xml` |
| `forcefields/water_hyperpol2.xml` | `src/act/forces/tests/water_hyperpol.xml` |

`ACS-pol.xml` and `water_dimer.pdb` were in two directories with the same
contents and are here once. The molecules of `src/act/alexandria/tests/mols`
and the files of `src/act/alexandria/tests/trainff` are in the directories
above, the GROMACS files of `src/testutils/simulationdatabase` in `other`.
