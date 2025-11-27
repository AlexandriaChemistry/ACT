#include <cstdio>
#include <string>
#include <vector>
#include "act/alexandria/topology.h"
#include "gromacs/fileio/pdbio.h"

namespace alexandria
{

//! Function to write a pdb file.
void pdbWriter(FILE                           *out, 
               const char                     *title,
               const std::vector<ActAtom>     &atoms,
               const std::vector<gmx::RVec>   &x,
               const std::vector<std::string> &residueNames,
               int                             ePBC,
               const matrix                    box,
               char                            chain,
               int                             model_nr,
               const std::vector<int>         &index,
               gmx_conect                      conect,
               bool                            renumberAtoms=false);

} // namespace
