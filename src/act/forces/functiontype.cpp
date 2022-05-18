
#include "functiontype.h"

#include <map>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{


typedef struct
{
    std::string name, description;
} FuncDescr;

std::map<FunctionType, FuncDescr> eftNames = {
    { FunctionType::MORSE,              { "Morse", "Morse potential function" } },
    { FunctionType::HARMONIC,           { "Harmonic", "Harmonic potential function" } },
    { FunctionType::UREY_BRADLEY,       { "Urey_Bradley", "UREY_BRADLEY for angles between three atoms" } },
    { FunctionType::COSINE,   			 { "Cosine", "Cosine function for angles between three atoms" } },
    { FunctionType::LINEAR_ANGLE, 		 { "LINEAR_ANGLE", "Linear angles" } },
    { FunctionType::IDIHS,              { "IDIHS", "Out of plane dihedral angles" } },
    { FunctionType::PDIHS, 		       { "PDIHS", "Proper dihedrals" } },
    { FunctionType::FOURDIHS,           { "FOURDIHS", "Proper dihedral" } },
    { FunctionType::DIRAC_DELTA,        { "Dirac", "Coulomb between two point charges" } },
    { FunctionType::GAUSSIAN,           { "Gaussian", "Coulomb between two Gaussian-type charges" } },
    { FunctionType::SLATER,	          { "Slater", "Coulomb between two Slater-type charges" } },
    { FunctionType::LENNARD_JONES, 		 { "Lennard Jones", "VDW functions" } },
    { FunctionType::BUCKINGHAM,    		 { "Buckingham", "VDW functions" } },
    { FunctionType::WANG_BUCKINGHAM, 	 { "Wang_Buckingham", "Buffered Buckingham potential" } }
};

const std::string &functionTypeToString(FunctionType fType)
{
    auto en = eftNames.find(fType);
    if (en == eftNames.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No string corresponding to function type %d",
                                                       static_cast<int>(fType)).c_str()));
    }
    return en->second.name;
}

const std::string &functionTypeToDescription(FunctionType fType)
{
    auto en = eftNames.find(fType);
    if (en == eftNames.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No string corresponding to function type %d",
                                                       static_cast<int>(fType)).c_str()));
    }
    return en->second.description;
}

FunctionType stringToFunctionType(const std::string &name)
{
    for (auto &eft : eftNames)
    {
        if (name == eft.second.name)
        {
            return eft.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such function type %s",
                                                       name.c_str()).c_str()));
    return FunctionType::MORSE;
}

int functionTypeToNparams(FunctionType fType)
{
    switch (fType)
    {
    case FunctionType::FOURDIHS:
        return 4;
    case FunctionType::MORSE:
    case FunctionType::BUCKINGHAM:
    case FunctionType::WANG_BUCKINGHAM:
        return 3;
    case FunctionType::HARMONIC:
    case FunctionType::LINEAR_ANGLE:
    case FunctionType::IDIHS:
    case FunctionType::LENNARD_JONES:
        return 2;
    case FunctionType::GAUSSIAN:
    case FunctionType::SLATER:
        return 1;
    case FunctionType::DIRAC_DELTA:
        return 0;
    }
    return 0;
}
} // namespace
