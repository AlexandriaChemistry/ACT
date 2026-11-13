/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023-2026
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "allmols.h"

#include <algorithm>
#include <string>

#include "act/basics/libraryfile.h"
#include "act/utility/jsonreader.h"
#include "act/utility/stringutil.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

/*! \brief A member of a molecule of alexandria.json as text, empty if it has none
 * \param[in] mol  The molecule
 * \param[in] name The member
 * \return The string, or the number as text
 */
static std::string jsonText(const JsonValue &mol, const std::string &name)
{
    const JsonValue *v = mol.find(name);
    if (!v || v->isNull())
    {
        return std::string();
    }
    if (v->type() == JsonValue::Type::Number)
    {
        return gmx::formatString("%.17g", v->asNumber());
    }
    return v->asString();
}

/*! \brief A member of a molecule of alexandria.json as a number, 0 if it has none
 * As my_atoi and my_atof of an empty column of the former alexandria.csv.
 * \param[in] mol  The molecule
 * \param[in] name The member
 * \return The number
 */
static double jsonNumber(const JsonValue &mol, const std::string &name)
{
    const JsonValue *v = mol.find(name);
    if (!v || v->isNull())
    {
        return 0;
    }
    if (v->type() == JsonValue::Type::String)
    {
        return my_atof(v->asString(), name);
    }
    return v->asNumber();
}

AlexandriaMol::AlexandriaMol(const JsonValue &mol)
{
    iupac          = jsonText(mol, "iupac");
    formula        = jsonText(mol, "formula");
    charge         = static_cast<int>(jsonNumber(mol, "charge"));
    mult           = static_cast<int>(jsonNumber(mol, "multiplicity"));
    mass           = jsonNumber(mol, "weight");
    cas            = jsonText(mol, "cas");
    csid           = static_cast<int>(jsonNumber(mol, "csid"));
    pubid          = static_cast<int>(jsonNumber(mol, "pubchem"));
    inchi          = jsonText(mol, "inchi");
    inchikey       = jsonText(mol, "inchikey");
    pointGroup     = jsonText(mol, "point_group");
    symmetryNumber = static_cast<int>(jsonNumber(mol, "symmetry_number"));
    auto list = [&mol](const std::string &name)
    {
        std::vector<std::string> items;
        const JsonValue *v = mol.find(name);
        if (v && v->type() == JsonValue::Type::Array)
        {
            for (const auto &e : v->asArray())
            {
                items.push_back(e.asString());
            }
        }
        return items;
    };
    for (auto ss : list("classification"))
    {
        classid.insert(ss);
        std::replace(ss.begin(), ss.end(), ' ', '-');
        classid.insert(ss);
    }
    // Synonyms must include the iupac as well in case we are
    // looking that one up later.
    synonyms.insert(iupac);
    for (auto ss : list("synonyms"))
    {
        if (ss.empty())
        {
            continue;
        }
        synonyms.insert(ss);
        std::replace(ss.begin(), ss.end(), ' ', '-');
        synonyms.insert(ss);
    }
    // Add the filename without its extension as well
    const auto filename = jsonText(mol, "filename");
    if (filename.size() > 4)
    {
        synonyms.insert(filename.substr(0, filename.size()-4));
    }
}

AlexandriaMols::AlexandriaMols()
    : AlexandriaMols(findLibrary("alexandria.json", false))
{
}

AlexandriaMols::AlexandriaMols(const std::string &fileName)
{
    if (fileName.empty())
    {
        return;
    }
    // The molecules of the Alexandria library in UTF-8, see src/act/python/mol_json.py
    const auto doc = readJsonFile(fileName);
    for (const auto &mol : doc.asArray())
    {
        const AlexandriaMol am(mol);
        if (!am.inchi.empty())
        {
            mols_.insert( { am.inchi, am } );
            for(const auto &syn : am.synonyms)
            {
                nameToInChi_.insert( { syn, am.inchi } );
            }
        }
    }
}

void AlexandriaMols::dump(FILE *fp)
{
    if (fp)
    {
        fprintf(fp, "There are %zu molecules and %zu synonyms\n", mols_.size(), nameToInChi_.size());
        for(const auto &n2i : nameToInChi_)
        {
            fprintf(fp, "%s %s\n", n2i.first.c_str(), n2i.second.c_str());
        }
    }
}
      
const AlexandriaMol *AlexandriaMols::findInChi(const std::string &inchi) const
{
    if (inchi.empty())
    {
        return nullptr;
    }
    const auto mptr = mols_.find(inchi);
    if (mols_.end() == mptr)
    {
        return nullptr;
    }
    return &mptr->second;
}

const AlexandriaMol *AlexandriaMols::findMol(const std::string &name) const
{
    const auto n2i = nameToInChi_.find(name);
    if (n2i == nameToInChi_.end())
    {
        return nullptr;
    }
    const auto mptr = mols_.find(n2i->second);
    if (mols_.end() == mptr)
    {
        return nullptr;
    }
    return &mptr->second;
}

}
