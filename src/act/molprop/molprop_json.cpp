/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
/*! \internal \file
 * \brief
 * Implements functions that describe the content of molprop files in a
 * JsonTree.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <map>
#include <set>
#include <string>
#include <vector>

#include "act/basics/msg_handler.h"
#include "act/molprop/experiment.h"
#include "act/molprop/fragment.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molpropobservable.h"
#include "act/utility/jsontree.h"
#include "molprop_json.h"

namespace alexandria
{

/*! \brief Count the values of each property of a molecule
 * \param[in]    mp     The molecule
 * \param[inout] counts The number of values by the name of the property
 */
static void countProperties(const MolProp &mp, std::map<std::string, int> *counts)
{
    for (const auto &exp : mp.experimentConst())
    {
        for (const auto &prop : exp.propertiesConst())
        {
            (*counts)[mpo_name(prop.first)] += static_cast<int>(prop.second.size());
        }
    }
}

/*! \brief Make an object with a number for each name
 * \param[in] key    The key of the object
 * \param[in] counts The numbers
 * \return the tree, an object also if there are no numbers
 */
static JsonTree countsToJson(const std::string &key, const std::map<std::string, int> &counts)
{
    JsonTree jc = JsonTree::object(key);
    for (const auto &count : counts)
    {
        jc.addObject(count.first, count.second);
    }
    return jc;
}

JsonTree molPropToJson(const MolProp &mp, const std::string &key)
{
    JsonTree jm(key);
    jm.addObject("name", mp.getMolname());
    jm.addObject("iupac", mp.getIupac());
    jm.addObject("formula", mp.formula());

    int      charge = 0;
    int      natom  = 0;
    JsonTree jfrags = JsonTree::array("fragments");
    for (const auto &frag : mp.fragments())
    {
        JsonTree jf("fragment");
        jf.addObject("inchi", frag.inchi());
        jf.addObject("iupac", frag.iupac());
        jf.addObject("formula", frag.formula());
        jf.addObject("charge", frag.charge());
        jf.addObject("multiplicity", frag.multiplicity());
        jf.addObject("natom", static_cast<int>(frag.atoms().size()));
        jfrags.addObject(jf);
        charge += frag.charge();
        natom  += static_cast<int>(frag.atoms().size());
    }
    if (natom == 0 && !mp.experimentConst().empty())
    {
        // Without fragments, the atoms of the first experiment
        natom = mp.experimentConst()[0].NAtom();
    }
    jm.addObject("natom", natom);
    jm.addObject("charge", charge);
    jm.addObject(jfrags);

    jm.addObject("experiments", static_cast<int>(mp.experimentConst().size()));
    std::map<std::string, int> jobtypes;
    std::set<std::string>      methods;
    for (const auto &exp : mp.experimentConst())
    {
        jobtypes[jobType2string(exp.getJobtype())] += 1;
        if (!exp.getMethod().empty() || !exp.getBasisset().empty())
        {
            methods.insert(exp.getMethod() + "/" + exp.getBasisset());
        }
    }
    jm.addObject(countsToJson("jobtypes", jobtypes));
    jm.addArray("methods", std::vector<std::string>(methods.begin(), methods.end()));
    std::map<std::string, int> counts;
    countProperties(mp, &counts);
    jm.addObject(countsToJson("properties", counts));
    return jm;
}

JsonTree molPropSummaryToJson(const std::vector<MolProp> &mps, const std::string &key)
{
    int                        monomers = 0, dimers = 0, larger = 0, experiments = 0;
    std::map<std::string, int> counts;
    for (const auto &mp : mps)
    {
        const auto nfrag = mp.fragments().size();
        if (nfrag <= 1)
        {
            monomers += 1;
        }
        else if (nfrag == 2)
        {
            dimers += 1;
        }
        else
        {
            larger += 1;
        }
        experiments += static_cast<int>(mp.experimentConst().size());
        countProperties(mp, &counts);
    }
    JsonTree js(key);
    js.addObject("molecules", static_cast<int>(mps.size()));
    js.addObject("monomers", monomers);
    js.addObject("dimers", dimers);
    js.addObject("larger", larger);
    js.addObject("experiments", experiments);
    js.addObject(countsToJson("properties", counts));
    return js;
}

JsonTree molPropValidationToJson(const MolPropValidation               &validation,
                                 const std::vector<ACTRecordedMessage> &messages,
                                 const std::string                     &key)
{
    JsonTree jv(key);
    jv.addBool("checked", true);
    jv.addBool("ok", validation.nproblem == 0);
    jv.addObject("problems", static_cast<int>(validation.nproblem));
    jv.addArray("invalid", std::vector<std::string>(validation.invalid.begin(), validation.invalid.end()));
    JsonTree jmsgs = JsonTree::array("messages");
    for (const auto &message : messages)
    {
        JsonTree jmsg("message");
        jmsg.addObject("level", actStatusName(message.level));
        jmsg.addObject("kind", actMessage(message.actm));
        jmsg.addObject("text", message.text);
        jmsgs.addObject(jmsg);
    }
    jv.addObject(jmsgs);
    return jv;
}

} // namespace alexandria
