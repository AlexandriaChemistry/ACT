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
#include "gromacs/utility/textreader.h"

namespace alexandria
{

AlexandriaMol::AlexandriaMol(const std::vector<std::string> &line)
{
    if (line.size() >= 20)
    {
        iupac    = line[2];
        formula  = line[3];
        charge   = my_atoi(line[4], "charge");
        mult     = my_atoi(line[5], "multiplicity");
        mass     = my_atof(line[6], "mass");
        cas      = line[7];
        csid     = my_atoi(line[8], "chemspider");
        pubid    = my_atoi(line[9], "pubchem");
        inchi    = line[10];
        inchikey = line[11];
        pointGroup     = line[15];
        symmetryNumber = my_atoi(line[16], "symmetry number");
        for(auto ss : split(line[12], ';'))
        {
            classid.insert(ss);
            std::replace(ss.begin(), ss.end(), ' ', '-');
            classid.insert(ss);
        }
        // Synonyms must include the iupac as well in case we are
        // looking that one up later.
        synonyms.insert(iupac);
        for(auto ss : split(line[13], ';'))
        {
            if (ss.empty())
            {
                continue;
            }
            synonyms.insert(ss);
            std::replace(ss.begin(), ss.end(), ' ', '-');
            synonyms.insert(ss);
        }
        // Add the filename as well
        const auto fn = line[14].substr(0, line[14].size()-4);
        if (!fn.empty())
        {
            synonyms.insert(fn);
        }
    }
}

/*! \brief The columns of alexandria.csv of a molecule of alexandria.json
 * The names of the members are those of COLUMNS in src/act/python/mol_csv.py,
 * which writes the JSON file. Numbers are written such that they are read
 * back as they are, lists are separated by ';', and missing members are empty.
 * \param[in] mol The molecule
 * \return The columns
 */
static std::vector<std::string> jsonColumns(const JsonValue &mol)
{
    static const std::vector<std::string> columns = {
        "alexandria", "parentid", "iupac", "formula", "charge",
        "multiplicity", "weight", "cas", "csid", "pubchem", "inchi",
        "inchikey", "classification", "synonyms", "filename",
        "point_group", "symmetry_number", "rotatable_bonds",
        "html_name", "latex_name", "author", "timestamp"
    };
    std::vector<std::string> line;
    for (const auto &name : columns)
    {
        const JsonValue *v = mol.find(name);
        std::string      word;
        if (v)
        {
            switch (v->type())
            {
            case JsonValue::Type::String:
                word = v->asString();
                break;
            case JsonValue::Type::Number:
                word = gmx::formatString("%.17g", v->asNumber());
                break;
            case JsonValue::Type::Array:
                for (const auto &e : v->asArray())
                {
                    word += (word.empty() ? "" : ";") + e.asString();
                }
                break;
            default:
                break;
            }
        }
        line.push_back(std::move(word));
    }
    return line;
}

void AlexandriaMols::add(const std::vector<std::string> &words)
{
    if (words.size() >= 20)
    {
        const AlexandriaMol am(words);
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

AlexandriaMols::AlexandriaMols()
    // The JSON file in UTF-8 of src/act/python/mol_csv2json, otherwise the
    // csv file, whichever the first library directory with either has
    : AlexandriaMols(findLibraryFirst({ "alexandria.json", "alexandria.csv" }, false))
{
}

AlexandriaMols::AlexandriaMols(const std::string &alexandria)
{
    if (alexandria.empty())
    {
        return;
    }
    if (gmx::endsWith(alexandria, ".json"))
    {
        const auto doc = readJsonFile(alexandria);
        for (const auto &mol : doc.asArray())
        {
            add(jsonColumns(mol));
        }
        return;
    }
    gmx::TextReader tr(alexandria);
    std::string line;
    while (tr.readLine(&line))
    {
        add(gmx::splitDelimitedString(line, '|'));
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
