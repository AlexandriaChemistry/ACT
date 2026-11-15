/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020,2025,2026
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "act/molprop/molprop_util.h"

#include <algorithm>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <utility>

#include "act/molprop/molprop_xml.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace alexandria
{

void generate_index(std::vector<MolProp> *mp)
{
    int index = 1;
    for (auto mpi = mp->begin(); mpi < mp->end(); ++mpi)
    {
        // We have to set some kind of index here to distinguish
        // the compounds.
        //! \todo: implement alexandria ID here.
        mpi->setIndex(index++);
    }
}

/*! \brief
 * Merge multiple molprops for molecules into one, e.g. one molprop for water
 * one for methane etc.
 *
 * \param[in]     msg_handler For storing warnings etc.
 * \param[inout]  mp          The vector of MolProps
 * \ingroup group_molprop
 */
static void MergeDoubleMolprops(MsgHandler                       *msg_handler,
                                std::vector<alexandria::MolProp> *mp)
{
    alexandria::MolPropIterator mpi, mmm[2];
    std::string                 molname[2];
    std::string                 form[2];
    int   i, ndouble = 0;
    bool  bDouble;
    int   cur   = 0;
#define prev (1-cur)
    gmx::TextWriter *tw = nullptr;
    if (msg_handler && msg_handler->debug())
    {
        tw = msg_handler->twDebug();
    }
    i = 0;
    for (mpi = mp->begin(); (mpi < mp->end()); )
    {
        bDouble = false;
        mpi->Dump(tw);
        mmm[cur]     = mpi;
        molname[cur] = mpi->getMolname();
        form[cur]    = mpi->formula();
        if (i > 0)
        {
            if (molname[prev] == molname[cur])
            {
                const auto warn = mmm[prev]->sameCompound(&(*mmm[cur]));
                if (msg_handler && !warn.empty())
                {
                    if (msg_handler)
                    {
                        for (const auto &w : warn)
                        {
                            msg_handler->msg(ACTStatus::Warning, w);
                        }
                    }
                }
                else
                {
                    const auto mergeWarn = mmm[prev]->Merge(&(*mmm[cur]));
                    if (msg_handler && !mergeWarn.empty())
                    {
                        for (const auto &w : mergeWarn)
                        {
                            msg_handler->msg(ACTStatus::Warning, w);
                        }
                    }
                    mpi    = mp->erase(mmm[cur]);

                    bDouble = true;
                    ndouble++;
                    if (mpi != mp->end())
                    {
                        mmm[cur]      = mpi;
                        molname[cur]  = mmm[cur]->getMolname();
                        form[cur]     = mmm[cur]->formula();
                    }
                }
            }
        }
        if (!bDouble)
        {
            cur = prev;
            mpi++;
            i++;
        }
    }
    for (mpi = mp->begin(); (mpi < mp->end()); mpi++)
    {
        mpi->Dump(tw);
    }
    if (msg_handler)
    {
        msg_handler->write(gmx::formatString("There were %d double entries, leaving %d after merging.\n",
                                             ndouble, (int)mp->size()));
    }
}

void merge_xml(MsgHandler                             *msg_handler,
               const gmx::ArrayRef<const std::string> &filens,
               std::vector<alexandria::MolProp>       *mpout)
{
    for (auto &fn : filens)
    {
        std::vector<alexandria::MolProp> mp;
        if (!gmx_fexist(fn))
        {
            if (msg_handler)
            {
                msg_handler->msg(ACTStatus::Warning, gmx::formatString("Molprop file named '%s' does not exist, but continuing anyway.", fn.c_str()));
            }
            continue;
        }
        MolPropRead(msg_handler, fn.c_str(), &mp);
        for (auto &mpi : mp)
        {
            mpout->push_back(std::move(mpi));
        }
    }
    auto mpBefore = mpout->size();
    gmx::TextWriter *tw = nullptr;
    if (msg_handler && msg_handler->debug())
    {
        tw = msg_handler->twDebug();

        tw->writeStringFormatted("mpout->size() = %u mpout->max_size() = %u\n",
                                (unsigned int)mpout->size(), (unsigned int)mpout->max_size());
        for (auto mpi = mpout->begin(); mpi < mpout->end(); ++mpi)
        {
            mpi->Dump(tw);
        }
    }
    std::sort(mpout->begin(), mpout->end(),
              [](const MolProp &ma, const MolProp &mb) { return ma.getMolname() < mb.getMolname(); });
    MergeDoubleMolprops(msg_handler, mpout);
    printf("There were %zu molecule entries before merging, %zu after.\n",
           mpBefore, mpout->size());
}

/*! \brief The first experiment of a molecule that has atoms
 * After MolProp::validate all experiments have their atoms in its order.
 * \param[in] mp The molecule
 * \return The experiment, nullptr if none
 */
static const Experiment *referenceExperiment(const MolProp &mp)
{
    for (const auto &ex : mp.experimentConst())
    {
        if (ex.NAtom() > 0)
        {
            return &ex;
        }
    }
    return nullptr;
}

/*! \brief The key of the compound of a fragment: its InChI and formula
 * A fragment with another identifier than an InChI, such as "1", does not
 * identify its compound, and validateMolProps reports it.
 * \param[in] frag The fragment
 * \return The key, empty if the fragment has no InChI
 */
static std::string compoundKey(const Fragment &frag)
{
    if (frag.inchi().rfind("InChI=", 0) != 0)
    {
        return "";
    }
    return frag.inchi() + " " + frag.formula();
}

/*! \brief The order of the atoms of a fragment, as names and types
 * \param[in] atoms The atoms of the molecule
 * \param[in] frag  The fragment
 * \return The atoms of the fragment in their order, empty if an atom is out of range
 */
static std::string fragmentOrder(const std::vector<CalcAtom> &atoms, const Fragment &frag)
{
    std::string order;
    for (int a : frag.atoms())
    {
        if (a < 0 || static_cast<size_t>(a) >= atoms.size())
        {
            return "";
        }
        order += gmx::formatString("%s%s:%s", order.empty() ? "" : " ",
                                   atoms[a].getName().c_str(), atoms[a].getObtype().c_str());
    }
    return order;
}

MolPropValidation validateMolProps(MsgHandler           *msghandler,
                                   std::vector<MolProp> *mps,
                                   std::vector<MolProp> *references,
                                   ACTStatus             level)
{
    MolPropValidation result;
    for (auto &mp : *mps)
    {
        if (!mp.validate(msghandler))
        {
            result.nproblem += 1;
            result.invalid.insert(mp.getMolname());
        }
        // Without an InChI the compound of a fragment is unknown, and its
        // order and charges cannot be compared with other molecules
        for (const auto &frag : mp.fragments())
        {
            if (compoundKey(frag).empty())
            {
                result.nproblem += 1;
                result.invalid.insert(mp.getMolname());
                msghandler->msg(level, ACTMessage::NoInChI,
                                gmx::formatString("Molecule %s: fragment %s has the identifier '%s', which is not an InChI.",
                                                  mp.getMolname().c_str(),
                                                  frag.iupac().empty() ? frag.formula().c_str() : frag.iupac().c_str(),
                                                  frag.inchi().c_str()));
            }
        }
    }
    // Where an order of the atoms of a compound occurs
    struct Occurrence
    {
        std::string molecule, datafile, order;
        //! The number of molecules with this order, and the first of them
        size_t      count = 0, first = 0;
        bool        reference = false;
    };
    // For each compound, by InChI and formula, its orders
    std::map<std::string, std::map<std::string, Occurrence>> compounds;
    size_t index = 0;
    const auto collect = [&compounds, &index](const MolProp &mp, bool reference)
    {
        const auto ref = referenceExperiment(mp);
        if (!ref)
        {
            return;
        }
        for (const auto &frag : mp.fragments())
        {
            const auto key   = compoundKey(frag);
            const auto order = key.empty() ? std::string() : fragmentOrder(ref->calcAtomConst(), frag);
            if (order.empty())
            {
                continue;
            }
            auto &occ = compounds[key][order];
            if (occ.count == 0 && !occ.reference)
            {
                occ = { mp.getMolname(), ref->getDatafile(), order, 0, index, reference };
            }
            occ.count     += reference ? 0 : 1;
            occ.reference  = occ.reference || reference;
            index         += 1;
        }
    };
    if (references)
    {
        for (auto &mp : *references)
        {
            if (mp.validate(msghandler))
            {
                collect(mp, true);
            }
        }
    }
    for (const auto &mp : *mps)
    {
        collect(mp, false);
    }
    // The reference order of each compound: that of the references, or else
    // the most common one, the first of those if several are
    std::map<std::string, const Occurrence *> chosen;
    for (const auto &c : compounds)
    {
        const Occurrence *best = nullptr;
        for (const auto &o : c.second)
        {
            const auto &occ = o.second;
            if (!best ||
                (occ.reference && !best->reference) ||
                (occ.reference == best->reference &&
                 (occ.count > best->count || (occ.count == best->count && occ.first < best->first))))
            {
                best = &occ;
            }
        }
        chosen[c.first] = best;
    }
    for (const auto &mp : *mps)
    {
        const auto ref = referenceExperiment(mp);
        if (!ref)
        {
            continue;
        }
        for (const auto &frag : mp.fragments())
        {
            const auto key   = compoundKey(frag);
            const auto order = key.empty() ? std::string() : fragmentOrder(ref->calcAtomConst(), frag);
            if (order.empty())
            {
                continue;
            }
            const auto *best = chosen[key];
            if (best->order != order)
            {
                result.nproblem += 1;
                result.invalid.insert(mp.getMolname());
                msghandler->msg(level, ACTMessage::InconsistentCompoundOrder,
                                gmx::formatString("Molecule %s (data file %s): compound %s has its atoms in the order %s, "
                                                  "but in the reference, molecule %s (data file %s%s), in the order %s.",
                                                  mp.getMolname().c_str(), ref->getDatafile().c_str(),
                                                  frag.iupac().empty() ? frag.inchi().c_str() : frag.iupac().c_str(),
                                                  order.c_str(), best->molecule.c_str(), best->datafile.c_str(),
                                                  best->reference ? " of the charges" : "", best->order.c_str()));
            }
        }
    }
    // Dimers A#B with their fragments in the order B, A. A file may have
    // both A#B and B#A, for instance with another hydrogen bond donor.
    for (const auto &mp : *mps)
    {
        const auto &name = mp.getMolname();
        const auto  hash = name.find('#');
        if (hash == std::string::npos || name.find('#', hash + 1) != std::string::npos)
        {
            continue;
        }
        const auto a = name.substr(0, hash);
        const auto b = name.substr(hash + 1);
        const auto &frags = mp.fragments();
        if (a != b && frags.size() == 2 && frags[0].iupac() == b && frags[1].iupac() == a)
        {
            result.nproblem += 1;
            result.invalid.insert(name);
            const auto ref = referenceExperiment(mp);
            msghandler->msg(level, ACTMessage::SwappedDimer,
                            gmx::formatString("Dimer %s (data file %s) has its fragments in the order %s, %s.",
                                              name.c_str(), ref ? ref->getDatafile().c_str() : "none",
                                              b.c_str(), a.c_str()));
        }
    }
    return result;
}

MolProp extractFragment(const MolProp &mp, size_t f)
{
    MolProp     out;
    const auto &frag = mp.fragments()[f];
    const std::string name = mp.fragments().size() == 1 ? mp.getMolname() :
        (frag.iupac().empty() ? frag.formula() : frag.iupac());
    out.SetMolname(name);
    out.SetIupac(frag.iupac().empty() ? name : frag.iupac());
    out.SetInchi(frag.inchi());
    const auto *ref = referenceExperiment(mp);
    if (!ref)
    {
        return out;
    }
    // The new index of each atom of the fragment
    std::map<int, int> index;
    std::vector<int>   atoms;
    Experiment         exper(ref->getProgram(), ref->getMethod(), ref->getBasisset(),
                             ref->getReference(), ref->getConformation(),
                             ref->getDatafile(), JobType::OPT);
    for (int a : frag.atoms())
    {
        index[a] = atoms.size();
        auto ca  = ref->calcAtomConst()[a];
        ca.setAtomid(atoms.size());
        exper.AddAtom(ca);
        atoms.push_back(atoms.size());
    }
    auto newFrag = frag;
    newFrag.setAtoms(atoms);
    out.addFragment(newFrag);
    for (const auto &b : mp.bondsConst())
    {
        const auto ai = index.find(b.aI());
        const auto aj = index.find(b.aJ());
        if (ai != index.end() && aj != index.end())
        {
            int    i, j;
            double bondorder;
            b.get(&i, &j, &bondorder);
            out.AddBond(Bond(ai->second, aj->second, bondorder));
        }
    }
    out.AddExperiment(std::move(exper));
    out.generateComposition();
    return out;
}

std::vector<MolProp> compoundReferences(MsgHandler                         *msghandler,
                                        const std::vector<const MolProp *> &molecules,
                                        const std::vector<MolProp>         &all,
                                        const std::vector<MolProp>         *references)
{
    // The compounds needed, in the order they are found, with the first
    // molecule and fragment they are in
    std::vector<std::pair<std::string, std::pair<const MolProp *, size_t>>> needed;
    std::set<std::string> keys;
    for (const auto *mp : molecules)
    {
        for (size_t f = 0; f < mp->fragments().size(); f++)
        {
            const auto key = compoundKey(mp->fragments()[f]);
            if (!key.empty() && keys.insert(key).second)
            {
                needed.push_back({ key, { mp, f } });
            }
        }
    }
    // The first monomer of each compound in a set of molecules
    const auto monomers = [](const std::vector<MolProp> &mps)
    {
        std::map<std::string, const MolProp *> found;
        for (const auto &mp : mps)
        {
            if (mp.fragments().size() == 1 && referenceExperiment(mp))
            {
                found.insert({ compoundKey(mp.fragments()[0]), &mp });
            }
        }
        return found;
    };
    const auto fromReferences = references ? monomers(*references) : std::map<std::string, const MolProp *>();
    const auto fromAll        = monomers(all);
    std::vector<MolProp> result;
    for (const auto &n : needed)
    {
        const MolProp *mp = nullptr;
        size_t         f  = 0;
        std::string    source;
        if (const auto r = fromReferences.find(n.first); r != fromReferences.end())
        {
            mp     = r->second;
            source = "the monomer in the charge file";
        }
        else if (const auto a = fromAll.find(n.first); a != fromAll.end())
        {
            mp     = a->second;
            source = "the monomer in the molecules";
        }
        else
        {
            mp     = n.second.first;
            f      = n.second.second;
            source = gmx::formatString("fragment %zu of %s", f + 1, mp->getMolname().c_str());
        }
        result.push_back(extractFragment(*mp, f));
        if (msghandler)
        {
            const auto *ref = referenceExperiment(*mp);
            msghandler->msg(ACTStatus::Info,
                            gmx::formatString("The charges of compound %s are those of %s (data file %s).",
                                              result.back().getMolname().c_str(), source.c_str(),
                                              ref ? ref->getDatafile().c_str() : ""));
        }
    }
    return result;
}

} // namespace alexandria
