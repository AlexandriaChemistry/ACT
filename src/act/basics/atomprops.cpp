#include "atomprops.h"

#include "act/utility/stringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

namespace alexandria
{

std::map<std::string, AtomProp> readAtomProps()
{
    const char *actdata = getenv("ACTDATA");
    if (nullptr == actdata)
    {
        GMX_THROW(gmx::InvalidInputError("Environment variable ACTDATA is not set"));
    }
    std::string     props = gmx::formatString("%s/atomprops.csv", actdata);
    gmx::TextReader tr(props);
    std::map<std::string, AtomProp> table;
    std::string              tmp;
    while (tr.readLine(&tmp))
    {
        auto ptr = split(tmp, '|');
        if (ptr.size() == 4)
        {
            AtomProp ap(ptr[1], 
                        my_atoi(ptr[2].c_str(), "atomnumber"),
                        my_atof(ptr[3].c_str(), "mass"));
            table.insert({ ptr[0], ap });
        }
    }
    return table;
}

} // namespace
