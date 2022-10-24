/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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
#include "communicationrecord.h"

#include <map>

#include "mpi.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

std::map<CommunicationStatus, const char*> csToString = {
    { CommunicationStatus::OK,        "Communication OK"        },
    { CommunicationStatus::DONE,      "Communication Done"      },
    { CommunicationStatus::ERROR,     "Communication Error"     },
    { CommunicationStatus::SEND_DATA, "Communication sent data" },
    { CommunicationStatus::RECV_DATA, "Communication OK"        }
};

std::map<bool, const char *> boolName = {
    { false, "false" },
    { true, "true" }
};

const char *cs_name(CommunicationStatus cs)
{
    return csToString[cs];
};

std::map<NodeType, const char *> ntToString = {
    { NodeType::Master,    "Master"    },
    { NodeType::MiddleMan, "MiddleMan" },
    { NodeType::Helper,    "Helper"    }  
};

CommunicationRecord::CommunicationRecord()
{
    cr_            = init_commrec();
    mpi_act_world_ = MPI_COMM_WORLD;
    (void) MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    (void) MPI_Comm_size(MPI_COMM_WORLD, &size_);
    if (rank_ == 0)
    {
        nt_ = NodeType::Master;
    }
}

void CommunicationRecord::check_init() const
{
    if (!initCalled_)
    {
        GMX_THROW(gmx::InternalError("ComunicationRecord::init has not been called"));
    }
}

void CommunicationRecord::print(FILE *fp)
{
    if (!fp)
    {
        return;
    }
    std::string strToPrint = gmx::formatString("rank: %d/%d nodetype: %s superior: %d\n",
                                               rank_, size_, ntToString[nt_], superior_);
    if (isMaster())
    {
        strToPrint.append(gmx::formatString("nmiddlemen: %d nhelper_per_middleman: %d\n", nmiddlemen_, 
                                            nhelper_per_middleman_));
    }
    if (isMasterOrMiddleMan())
    {
        strToPrint.append(gmx::formatString("ordinal: %d nhelper: %d\n",
                                            ordinal_, nhelper_per_middleman_));
    }
    if (!helpers_.empty())
    {
        strToPrint.append("helpers:");
        for (auto &h : helpers_)
        {
            strToPrint.append(gmx::formatString(" %3d", h));
        }
        strToPrint.append("\n");
    }
    if (!middlemen_.empty())
    {
        strToPrint.append("middlemen:");
        for (auto &m : middlemen_)
        {
            strToPrint.append(gmx::formatString(" %3d", m));
        }
        strToPrint.append("\n");
    }
    fprintf(fp, "%s", strToPrint.c_str());
}

CommunicationStatus CommunicationRecord::init(int nmiddleman)
{
    if (initCalled_)
    {
        return CommunicationStatus::TWOINIT;
    }
    nmiddlemen_ = nmiddleman;
    if (nmiddlemen_ > size_ || nmiddlemen_ < 1)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Cannot handle %d middlemen (individuals) with %d cores/threads",
                  nmiddlemen_, size_).c_str()));
    }
    if (nmiddlemen_ == 1)  // Only MASTER + HELPERS
    {
        nhelper_per_middleman_ = size_-1;
    }
    else
    {
        nhelper_per_middleman_ = size_ / nmiddlemen_ - 1;
        
        // We are picky. Each individual needs the same number of helpers
        if (nmiddlemen_ * (1+nhelper_per_middleman_) != size_)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("The number of cores/threads (%d) should be the product of the number of workers (per individual) (%d) and the number of individuals (%d)", size_, nhelper_per_middleman_ + 1, nmiddlemen_).c_str()));
        }
    }
    // Select the node type etc.
    if (rank_ == 0)  // If I am the MASTER
    {
        // Set ordinal to 0, as the first middleman
        ordinal_ = 0;
        // middlemen_.push_back(0);
        if (nmiddlemen_ == 1) // If only master and helpers
        {
            for (int i = 1; i < size_; i++)
            {
                helpers_.push_back(i);
            }
        }
        else  // If there are helpers
        {
            // Fill in my helpers
            for (int i = 0; i < nhelper_per_middleman_; i++)
            {
                helpers_.push_back(i+1);
            }
            // Fill in the middlemen
            for (int i = 1+nhelper_per_middleman_; i < size_; i += 1+nhelper_per_middleman_)
            {
                middlemen_.push_back(i);
            }
        }
        if (nmiddlemen_ - 1 != static_cast<int>(middlemen_.size()))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Inconsistency in number of middlemen. Expected %d but found %zu.", nmiddlemen_, middlemen_.size()).c_str()));
        }                   
    }
    else  // If I am NOT the MASTER
    {
        if (nmiddlemen_ == 1)
        {
            // No middlemen, then I am a helper reporting to the master
            nt_       = NodeType::Helper;
            superior_ = 0;
            // Not updating ordinal_
        }
        else
        {
            if (rank_ % (1+nhelper_per_middleman_) == 0)
            {
                nt_       = NodeType::MiddleMan;
                superior_ = 0;
                ordinal_  = rank_ / (1+nhelper_per_middleman_);
                for (int i = 0; i < nhelper_per_middleman_; i++)
                {
                    helpers_.push_back(rank_ + i + 1);
                }
            }
            else
            {
                nt_       = NodeType::Helper;
                superior_ = (rank_/(1+nhelper_per_middleman_)) * (1+nhelper_per_middleman_);
            }
        }
    }
    
    // Create ACT communicators for helpers and middlemen
    // Default value
    mpi_act_helpers_ = mpi_act_world_;
    // Split the non-masters into nmiddlemen X nhelpers
    if (nhelper_per_middleman_ > 0)
    {
        MPI_Comm_split(mpi_act_world_, rank_ / (nhelper_per_middleman_ + 1),
                       rank_ % (nhelper_per_middleman_ + 1), &mpi_act_helpers_);
    }
    print(stdout);
    initCalled_ = true;
    return CommunicationStatus::OK;
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

static void check_return(const char *msg, int returnvalue)
{
    switch (returnvalue)
    {
    case MPI_SUCCESS:
        return;
    case MPI_ERR_COMM:
        GMX_THROW(gmx::InternalError(gmx::formatString("Invalid communicator. %s.", msg).c_str()));
        break;
    case MPI_ERR_ARG:
        GMX_THROW(gmx::InternalError(gmx::formatString("Invalid argument. %s.", msg).c_str()));
        break;
    }
}

void CommunicationRecord::send(int dest, const void *buf, int bufsize) const
{
    int         tag = 0;
    MPI_Status  status;
    MPI_Request req;

    check_init();
    check_return("MPI_Isend Failed",
                 MPI_Isend(buf, bufsize, MPI_BYTE, RANK(cr, dest), tag,
                           mpi_act_world_, &req));
    check_return("MPI_Wait failed", MPI_Wait(&req, &status));
}

void CommunicationRecord::recv( int src, void *buf, int bufsize) const
{
    int         tag = 0;
    MPI_Status  status;
    MPI_Request req;

    check_init();
    check_return ("MPI_Irecv Failed",
                  MPI_Irecv(buf, bufsize, MPI_BYTE, src, tag,
                            mpi_act_world_, &req));
    check_return("MPI_Wait failed", MPI_Wait(&req, &status));
}

void CommunicationRecord::bcast_str(std::string *str,
                                    bool helpers_only) const
{
    check_init();
    int ssize = str->size();
    MPI_Bcast((void *)&ssize, 1, MPI_INT, superior(),
              helpers_only ? mpi_act_helpers_ : mpi_act_world_);
    if (0 != rank_)
    {
        str->resize(ssize);
    }
    MPI_Bcast((void *)str->data(), ssize, MPI_BYTE, superior(), 
              helpers_only ? mpi_act_helpers_ :  mpi_act_world_);
}
    
void CommunicationRecord::bcast_int(int *i,
                                    bool helpers_only) const
{
    check_init();
    MPI_Bcast((void *)i, 1, MPI_INT, superior(),
              helpers_only ? mpi_act_helpers_ :  mpi_act_world_);
}
    
void CommunicationRecord::bcast_bool(bool *b,
                                     bool helpers_only) const
{
    check_init();
    int d = *b ? 1 : 0;
    MPI_Bcast((void *)&d, 1, MPI_INT, superior(),
              helpers_only ? mpi_act_helpers_ :  mpi_act_world_);
    *b = d;
}
    
void CommunicationRecord::bcast_double(double *d,
                                       bool helpers_only) const
{
    check_init();
    MPI_Bcast((void *)d, 1, MPI_DOUBLE, superior(),
              helpers_only ? mpi_act_helpers_ :  mpi_act_world_);
}
   
void CommunicationRecord::bcast_double_vector(std::vector<double> *d,
                                              bool helpers_only) const
{
    check_init();
    int ssize = d->size(); 
    MPI_Bcast((void *)&ssize, 1, MPI_INT, superior(),
              helpers_only ? mpi_act_helpers_ :  mpi_act_world_);
    if (0 != rank_)
    {
        d->resize(ssize);
    }

    MPI_Bcast((void *)d->data(), ssize, MPI_DOUBLE, 0, mpi_act_world_);
}

void CommunicationRecord::send_str(int dest, const std::string *str) const
{
    check_init();
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
    check_init();
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
    check_init();
    if (nullptr != debug)
    {
        fprintf(debug, "Sending double '%g' to %d\n", d, dest);
    }
    send(dest, &d, sizeof(d));
}

double CommunicationRecord::recv_double(int src) const
{
    check_init();
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
    check_init();
    if (nullptr != debug)
    {
        fprintf(debug, "Sending int '%d' to %d\n", d, dest);
    }
    send(dest, &d, sizeof(d));
}

void CommunicationRecord::send_bool(int dest, bool b) const
{
    check_init();
    if (nullptr != debug)
    {
        fprintf(debug, "Sending bool '%s' to %d\n", boolName[b], dest);
    }
    int d = b ? 1 : 0;
    send(dest, &d, sizeof(d));
}

int CommunicationRecord::recv_int(int src) const
{
    check_init();
    int d;

    recv(src, &d, sizeof(d));
    if (nullptr != debug)
    {
        fprintf(debug, "Received int '%d' from %d\n", d, src);
    }

    return d;
}

bool CommunicationRecord::recv_bool(int src) const
{
    check_init();
    int d;

    recv(src, &d, sizeof(d));
    bool b = d;
    if (nullptr != debug)
    {
        fprintf(debug, "Received bool '%s' from %d\n", boolName[b], src);
    }

    return b;
}

void CommunicationRecord::send_ff_middleman_mode(int dest, TuneFFMiddlemanMode mode) const
{
    check_init();
    send_int(dest, static_cast<int>(mode));
}

TuneFFMiddlemanMode CommunicationRecord::recv_ff_middleman_mode(int src) const
{
    check_init();
    return static_cast<TuneFFMiddlemanMode>( recv_int(src) );
}

std::map<iMolSelect, int> imsInt = 
    {
        { iMolSelect::Train, 3 },
        { iMolSelect::Test,  7 },
        { iMolSelect::Ignore, 13 }
    };
    
void CommunicationRecord::send_iMolSelect(int dest, iMolSelect ims) const
{
    if (nullptr != debug)
    {
        fprintf(debug, "Sending iMolSelect '%s' to %d\n",
                iMolSelectName(ims), dest);
    }
    int d = imsInt[ims];
    send(dest, &d, sizeof(d));
}   

iMolSelect CommunicationRecord::recv_iMolSelect(int src) const
{
    int d;

    recv(src, &d, sizeof(d));
    for(auto &ims : imsInt)
    {
        if (ims.second == d)
        {
            if (nullptr != debug)
            {
                fprintf(debug, "Received iMolSelect '%s' from %d\n",
                        iMolSelectName(ims.first), src);
            }
            return ims.first;
        }
    }
    GMX_THROW(gmx::InternalError(gmx::formatString("Received unknown iMolSelect, int code %d", d).c_str()));
    return iMolSelect::Ignore;
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
    if (nhelper_per_middleman_ == 0)
    {
        return;
    }
    int size = 0;
    if (MPI_SUCCESS == MPI_Comm_size(mpi_act_helpers_, &size) && size > 1)
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM,
                      mpi_act_helpers_);
    }
}

void CommunicationRecord::sumi_helpers(int nr, 
                                       int r[]) const
{
    if (nhelper_per_middleman_ == 0)
    {
        return;
    }
    int size = 0;
    if (MPI_SUCCESS == MPI_Comm_size(mpi_act_helpers_, &size) && size > 1)
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM,
                      mpi_act_helpers_);
    }
}

#define ACT_SEND_DATA 19823
#define ACT_SEND_DONE -666
CommunicationStatus CommunicationRecord::bcast_data() const
{
    int acd = ACT_SEND_DATA;
    bcast_int(&acd);

    if (ACT_SEND_DATA == acd)
    {
        return CommunicationStatus::OK;
    }
    else
    {
        return CommunicationStatus::ERROR;
    }
}

CommunicationStatus CommunicationRecord::send_data(int dest) const
{
    send_int(dest, ACT_SEND_DATA);

    return CommunicationStatus::SEND_DATA;
}

CommunicationStatus CommunicationRecord::send_done(int dest) const
{
    send_int(dest, ACT_SEND_DONE);

    return CommunicationStatus::OK;
}

CommunicationStatus CommunicationRecord::bcast_done() const
{
    int acd = ACT_SEND_DONE;
    bcast_int(&acd);

    if (ACT_SEND_DONE == acd)
    {
        return CommunicationStatus::OK;
    }
    else
    {
        return CommunicationStatus::ERROR;
    }
}

CommunicationStatus CommunicationRecord::recv_data(int src) const
{
    int kk = recv_int(src);

    switch (kk)
    {
    case ACT_SEND_DATA:
        return CommunicationStatus::RECV_DATA;
    case ACT_SEND_DONE:
        return CommunicationStatus::DONE;
    default:
        GMX_THROW(gmx::InternalError(gmx::formatString("Received %d from src %d in gmx_recv_data. Was expecting either %d or %d\n.", kk, src, (int)ACT_SEND_DATA, (int)ACT_SEND_DONE).c_str()));
    }
    return CommunicationStatus::ERROR;
}

#undef ACT_SEND_DATA
#undef ACT_SEND_DONE

} // namespace alexandria
