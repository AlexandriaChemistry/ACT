
#include "dataset.h"

    
std::map<iMolSelect, const char *> MolSelect_Names = {
    { iMolSelect::Train,   "Train"   },
    { iMolSelect::Test,    "Test"    },
    { iMolSelect::Ignore,  "Ignore"  }
};

const std::map<iMolSelect, const char *> &iMolSelectNames()
{
    return MolSelect_Names;
}

const char *iMolSelectName(iMolSelect ims)
{
    return MolSelect_Names[ims];
}

bool name2molselect(const std::string &name, iMolSelect *ims)
{
    for (auto iter = MolSelect_Names.begin(); iter != MolSelect_Names.end(); ++iter)
    {
        if (name.compare(iter->second) == 0)
        {
            *ims = iter->first;
            return true;
        }
    }
    return false;
}
