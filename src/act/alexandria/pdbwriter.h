#include <cstdio>
#include <string>
#include <vector>

#include "gromacs/fileio/pdbio.h"

namespace alexandria
{

class ActAtom;

/*! \brief Make unique names for the particles in each residue
 * The name is the atom name followed by the number of the atom or
 * virtual site in its residue, such as C1 or H2, the numbering starts
 * at one in every residue. A shell gets the number of the atom
 * preceding it. The OpenMM force field and the PDB file use the same
 * names, such that OpenMM can find the bonds from them.
 * \param[in] atoms The particles
 * \return The names
 */
std::vector<std::string> particleNames(const std::vector<ActAtom> &atoms);

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
