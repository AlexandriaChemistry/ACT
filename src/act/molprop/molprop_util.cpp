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

#include "act/molprop/molprop_xml.h"
#include "act/topology/actmol.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/futil.h"
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

} // namespace alexandria
