




#ifndef FUNCTIONTYPE_H
#define FUNCTIONTYPE_H

#include <string>

namespace alexandria
{

enum class FunctionType

{

    // Potential functions for describing InteractionType::BOND

    //! Harmonic function
    HARMONIC;
    
    //! Morse function
    MORSE;

    // Potential functions for describing InteractionType::ANGLE

    //!
    UREY_BRADLEY;

    //!
    COSINE;


    // Potential functions for describing InteractionType::LINEAR_ANGLES
    LINEAR_ANGLE;


    // Potential functions for describing InteractionType::IMPROPER_DIHEDRAL

    IDIHS;

    // Potential functions for describing InteractionType::PROPER_DIHEDRAL

    //!
    PDIHS;

    //!
    FOURDIHS;

    // Potential functions for describing InteractionType::CHARGEDISTRIBUTION

    //! Dirac delat function
    DIRAC_DELTA;

    //! Gaussain distribution function
    GAUSSIAN;

    //! Slater distribution function
    SLATER;

    // Potential functions for describing InteractionType::VDW

    //! Lennard-Jones function
    LENNARD_JONES

    //! Buckingham function
    BUCKINGHAM;

    //! Wang-Buckingham function
    WANG_BUCKINGHAM;
};

const std::string &functionTypeToString(FunctionType fType);

FunctionType stringToFfunctionType(const std::string &name);




} // namespace alexandria

#endif

