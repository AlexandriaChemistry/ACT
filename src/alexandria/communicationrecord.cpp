/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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
#include "communicationrecord.h"

#include "mpi.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{
CommunicationRecord::CommunicationRecord()
{
    cr_   = init_commrec();
    rank_ = cr_->nodeid;
    size_ = cr_->nnodes;
}

void CommunicationRecord::init(int nmiddleman)
{
    nmiddlemen_ = nmiddleman;
    if (nmiddlemen_ > size_-1 || nmiddlemen_ < 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Cannot handle %d middlemen (individuals) with %d cores/threads",
                  nmiddlemen_, size_).c_str()));
    }
    cr_->nmiddlemen = nmiddlemen_;
    if (nmiddlemen_ == 0)
    {
        nhelper_per_middleman_ = size_-1;
    }
    else
    {
        nhelper_per_middleman_ = (size_-1) / nmiddlemen_ - 1;
        
        // We are picky. Each individual needs the same number of helpers
        if (!((nhelper_per_middleman_ == 0 && 
               size_ == 1+nmiddlemen_) ||
              (nhelper_per_middleman_ > 0 && 
               size_ % (1+nhelper_per_middleman_) == 1)))
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("The number of cores/threads (%d) should be the product of the number of helpers (%d) and the number of individuals (%d) plus 1 for the overlord", size_, nhelper_per_middleman_, nmiddlemen_).c_str()));
        }
    }
    cr_->nhelper_per_middleman = nhelper_per_middleman_;
    // Select the node type etc.
    if (rank_ == 0)
    {
        nt_ = NodeType::Master;
        for (int i = 1; i < size_; i++)
        {
            helpers_.push_back(i);
            // Not updating superior_ from the default value. If it is
            // used for communication the program will crash, since
            // it is an error to do so.
            
            // Not updating ordinal_ for the same reason.
            if (nmiddlemen_ > 0)
            {
                if ((i - 1) % (1+nhelper_per_middleman_) == 0)
                {
                    middlemen_.push_back(i);
                }
            }
        }
        if (nmiddlemen_ != static_cast<int>(middlemen_.size()))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Inconsistency in number of middlemen. Expected %d but found %zu.", nmiddlemen_, middlemen_.size()).c_str()));
        }                   
    }
    else
    {
        if (nmiddlemen_ == 0)
        {
            nt_       = NodeType::Helper;
            superior_ = 0;
            // Not updating ordinal_, see above.
        }
        else
        {
            if ((rank_ - 1) % (1+nhelper_per_middleman_) == 0)
            {
                nt_       = NodeType::MiddleMan;
                superior_ = 0;
                ordinal_  = (rank_ - 1) / (1+nhelper_per_middleman_);
                for (int i = rank_ + 1; i <= nhelper_per_middleman_; i++)
                {
                    helpers_.push_back(i);
                }
            }
            else
            {
                nt_       = NodeType::Helper;
                superior_ = (rank_-1)/(1+nhelper_per_middleman_) + 1;
            }
        }
    }
    
    // Create ACT communicators for helpers and middlemen
    int color = MPI_UNDEFINED;
    int key   = 0;
    if (rank_ > 0)
    {
        color = 0;
        key   = rank_-1;
    }
    MPI_Comm_split(cr_->mpi_comm_mysim, color, key, &cr_->mpi_act_not_master);
    // Default value
    cr_->mpi_act_helpers = MPI_COMM_WORLD;
    if (!MASTER(cr_))
    {
        int nonMasterSize;
        MPI_Comm_size(cr_->mpi_act_not_master, &nonMasterSize);
        int myrank;
        MPI_Comm_rank(cr_->mpi_act_not_master, &myrank);
        if (nmiddlemen_ > 0)
        {
            // Split the non-masters into nmiddlemen X nhelpers
            MPI_Comm_split(cr_->mpi_act_not_master, myrank / nmiddlemen_,
                           myrank % nmiddlemen_, &cr_->mpi_act_helpers);
        }
    }
}

CommunicationRecord::~CommunicationRecord()
{
    if (cr_)
    {
        done_commrec(cr_);
    }
}

/*************************************************
 *           LOW LEVEL ROUTINES                  *
 *************************************************/
void CommunicationRecord::send(int dest, const void *buf, int bufsize) const
{
    int         tag = 0;
    MPI_Status  status;
    MPI_Request req;

    if (MPI_Isend(buf, bufsize, MPI_BYTE, RANK(cr, dest), tag,
                  cr_->mpi_comm_mygroup, &req) != MPI_SUCCESS)
    {
        gmx_comm("MPI_Isend Failed");
    }
    if (MPI_Wait(&req, &status) != MPI_SUCCESS)
    {
        gmx_comm("MPI_Wait failed");
    }
}

void CommunicationRecord::recv( int src, void *buf, int bufsize) const
{
    int         tag = 0;
    MPI_Status  status;
    MPI_Request req;

    if (MPI_Irecv(buf, bufsize, MPI_BYTE, src, tag,
                  cr_->mpi_comm_mygroup, &req) != MPI_SUCCESS)
    {
        gmx_comm("MPI_Irecv Failed");
    }
    if (MPI_Wait(&req, &status) != MPI_SUCCESS)
    {
        gmx_comm("MPI_Wait failed");
    }
}

void CommunicationRecord::send_str(int dest, const std::string *str) const
{
    int len = str->size();

    if (nullptr != debug)
    {
        fprintf(debug, "Sending string '%s' to %d\n", str->c_str(), dest);
    }
    send(dest, &len, sizeof(len));
    if (!str->empty())
    {
        send(dest, (void *)str->data(), len);
    }
}

void CommunicationRecord::recv_str(int src, std::string *str) const
{
    int   len;

    recv(src, &len, sizeof(len));
    if (len == 0)
    {
        str->clear();
    }
    else
    {
        str->resize(len);
        recv(src, (void*)str->data(), len);
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Received string '%s' from %d\n", str->c_str(), src);
    }
}

void CommunicationRecord::send_double(int dest, double d) const
{
    if (nullptr != debug)
    {
        fprintf(debug, "Sending double '%g' to %d\n", d, dest);
    }
    send(dest, &d, sizeof(d));
}

double CommunicationRecord::recv_double(int src) const
{
    double d;

    recv(src, &d, sizeof(d));
    if (nullptr != debug)
    {
        fprintf(debug, "Received double '%g' from %d\n", d, src);
    }

    return d;
}

void CommunicationRecord::send_int(int dest, int d) const
{
    if (nullptr != debug)
    {
        fprintf(debug, "Sending int '%d' to %d\n", d, dest);
    }
    send(dest, &d, sizeof(d));
}

int CommunicationRecord::recv_int(int src) const
{
    int d;

    recv(src, &d, sizeof(d));
    if (nullptr != debug)
    {
        fprintf(debug, "Received int '%d' from %d\n", d, src);
    }

    return d;
}

void CommunicationRecord::send_double_vector(int dest,
                                             const std::vector<double> *d) const
{
    send_int(dest, d->size());
    if (d->size() > 0)
    {
        send(dest, static_cast<const void *>(d->data()),
             d->size()*sizeof(double));
    }
}

void CommunicationRecord::recv_double_vector(int src,
                                             std::vector<double> *d) const
{
    int len = recv_int(src);
    d->resize(len);
    if (len > 0)
    {
        recv(src, static_cast<void *>(d->data()), len*sizeof(double));
    }
}

void CommunicationRecord::sumd_helpers(int    nr, 
                                       double r[]) const
{
    int size = 0;
    if (MPI_Comm_size(cr_->mpi_act_helpers, &size) && size > 1)
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM,
                      cr_->mpi_act_helpers);
    }
}

void CommunicationRecord::sumi_helpers(int nr, 
                                       int r[]) const
{
    int size = 0;
    if (MPI_Comm_size(cr_->mpi_act_helpers, &size) && size > 1)
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM,
                      cr_->mpi_act_helpers);
    }
}

} // namespace alexandria
