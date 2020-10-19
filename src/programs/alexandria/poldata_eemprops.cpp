/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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

#include "poldata_eemprops.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <vector>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/stringutil.h"

#include "gmx_simple_comm.h"
#include "plistwrapper.h"
#include "stringutil.h"

namespace alexandria
{

RowZetaQ::RowZetaQ(int row, double zeta, double q)

    :

      row_(row),
      zeta_(zeta),
      q_(q),
      zetaRef_(zeta)


{
    zindex_ = -1;
    char buf[256];
#if HAVE_LIBCLN
    row_ = std::min(row_, SLATER_MAX_CLN);
    if (row_ < row && debug)
    {
        fprintf(debug, "Reducing row from %d to %d\n", row, row_);
    }
    snprintf(buf, sizeof(buf), "Row (%d) is bigger than Slater Max (%d)", row_, SLATER_MAX_CLN);
    GMX_RELEASE_ASSERT(row_ <= SLATER_MAX_CLN, buf);
#else
    row_ = std::min(row_, SLATER_MAX);
    if (row_ < row && debug)
    {
        fprintf(debug, "Reducing row from %d to %d\n", row, row_);
    }
    snprintf(buf, sizeof(buf), "Row (%d) is bigger than Slater Max (%d)", row_, SLATER_MAX);
    GMX_RELEASE_ASSERT(row_ <= SLATER_MAX, buf);
#endif
    fixedQ_ = (q != 0);
}

CommunicationStatus RowZetaQ::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, row_);
        gmx_send_double(cr, dest, zeta_);
        gmx_send_double(cr, dest, q_);
        gmx_send_double(cr, dest, zetaRef_);
        gmx_send_int(cr, dest, zindex_);
        gmx_send_int(cr, dest, fixedQ_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent RowZetaQ %d %g %g %g %d %d, status %s\n",
                    row_, zeta_, q_, zetaRef_, zindex_, fixedQ_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus RowZetaQ::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        row_     = gmx_recv_int(cr, src);
        zeta_    = gmx_recv_double(cr, src);
        q_       = gmx_recv_double(cr, src);
        zetaRef_ = gmx_recv_double(cr, src);
        zindex_  = gmx_recv_int(cr, src);
        fixedQ_  = gmx_recv_int(cr, src);

        if (nullptr != debug)
        {
            fprintf(debug, "Received RowZetaQ %d %g %g %g %d %d, status %s\n",
                    row_, zeta_, q_, zetaRef_, zindex_, fixedQ_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

Eemprops::Eemprops(const std::string        &name,
                   const std::string        &rowstr,
                   const std::string        &zetastr,
                   const std::string        &zeta_sigma,
                   const std::string        &qstr,
                   double                    qref,
                   double                    J0,
                   double                    J0_sigma,
                   double                    chi0,
                   double                    chi0_sigma)
    :
      name_(name),
      rowstr_(rowstr),
      zetastr_(zetastr),
      zeta_sigma_(zeta_sigma),
      qstr_(qstr),
      qref_(qref),
      J0_(J0),
      J0_sigma_(J0_sigma),
      chi0_(chi0),
      chi0_sigma_(chi0_sigma)
{
    setRowZetaQ(rowstr, zetastr, qstr);
}

CommunicationStatus Eemprops::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &name_);
        gmx_send_str(cr, dest, &rowstr_);
        gmx_send_str(cr, dest, &zetastr_);
        gmx_send_str(cr, dest, &zeta_sigma_);
        gmx_send_str(cr, dest, &qstr_);
        gmx_send_double(cr, dest, qref_);
        gmx_send_double(cr, dest, J0_);
        gmx_send_double(cr, dest, J0_sigma_);
        gmx_send_double(cr, dest, chi0_);
        gmx_send_double(cr, dest, chi0_sigma_);
        gmx_send_int(cr, dest, rzq_.size());

        for (auto &rzq : rzq_)
        {
            cs = rzq.Send(cr, dest);
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Eemprops %s %s %s %s %g %g, status %s\n",
                    name_.c_str(), rowstr_.c_str(),
                    zetastr_.c_str(), qstr_.c_str(), J0_, chi0_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Eemprops::Receive(const t_commrec *cr, int src)
{
    size_t              nrzq;
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &name_);
        gmx_recv_str(cr, src, &rowstr_);
        gmx_recv_str(cr, src, &zetastr_);
        gmx_recv_str(cr, src, &zeta_sigma_);
        gmx_recv_str(cr, src, &qstr_);
        qref_       = gmx_recv_double(cr, src);
        J0_         = gmx_recv_double(cr, src);
        J0_sigma_   = gmx_recv_double(cr, src);
        chi0_       = gmx_recv_double(cr, src);
        chi0_sigma_ = gmx_recv_double(cr, src);
        nrzq        = gmx_recv_int(cr, src);

        rzq_.clear();
        for (size_t n = 0; n < nrzq; n++)
        {
            RowZetaQ rzq;
            cs = rzq.Receive(cr, src);
            if (CS_OK == cs)
            {
                rzq_.push_back(rzq);
            }
        }
    }
    return cs;
}

void Eemprops::setRowZetaQ(const std::string &rowstr,
                           const std::string &zetastr,
                           const std::string &qstr)
{
    rzq_.clear();
    std::vector<std::string>  sz, sq, sr;
    sr = gmx::splitString(rowstr);
    sz = gmx::splitString(zetastr);
    sq = gmx::splitString(qstr);
    size_t nn = std::min(sz.size(), std::min(sq.size(), sr.size()));

    std::string buf = 
        gmx::formatString("%d zeta (%s), %d q (%s) and %d row (%s) values for %s\n",
                          static_cast<int>(sz.size()), zetastr.c_str(),
                          static_cast<int>(sq.size()), qstr.c_str(),
                          static_cast<int>(sr.size()), rowstr.c_str(),
                          getName());
    GMX_RELEASE_ASSERT(sz.size() == sr.size() && sz.size() == sq.size(), buf.c_str());

    for (size_t n = 0; n < nn; n++)
    {
        RowZetaQ rzq(atoi(sr[n].c_str()), my_atof(sz[n].c_str(), "zeta"), my_atof(sq[n].c_str(), "q"));
        rzq_.push_back(rzq);
    }
}

} // namespace alexandria
