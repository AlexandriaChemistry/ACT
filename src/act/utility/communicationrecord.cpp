/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2026
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
#include <set>

#include "mpi.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

//! Map CommunicationStatus tp string
std::map<CommunicationStatus, const char*> csToString = {
    { CommunicationStatus::OK,        "Communication OK"        },
    { CommunicationStatus::DONE,      "Communication Done"      },
    { CommunicationStatus::ERROR,     "Communication Error"     },
    { CommunicationStatus::SEND_DATA, "Communication sent data" },
    { CommunicationStatus::RECV_DATA, "Communication OK"        },
    { CommunicationStatus::TWOINIT,   "Double initiation"       }
};

//! Map bool to string
std::map<bool, const char *> boolName = {
    { false, "false" },
    { true, "true" }
};

const char *cs_name(CommunicationStatus cs)
{
    return csToString[cs];
};

//! \brief Map NodeType to string
std::map<NodeType, const char *> ntToString = {
    { NodeType::Master,    "Master"    },
    { NodeType::MiddleMan, "MiddleMan" },
    { NodeType::Helper,    "Helper"    }  
};

//! \brief Evaluate return value and crash if needed.
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

CommunicationRecord::CommunicationRecord(MsgHandler *msghandler)
{
    cr_            = init_commrec();
    check_return("MPI_Comm_dup", MPI_Comm_dup(MPI_COMM_WORLD, &mpi_act_world_));
    check_return("MPI_Comm_rank", MPI_Comm_rank(mpi_act_world_, &rank_));
    check_return("MPI_Comm_size", MPI_Comm_size(mpi_act_world_, &size_));
    if (rank_ == 0)
    {
        nt_ = NodeType::Master;
    }
    msg_handler_ = msghandler;
}

void CommunicationRecord::check_init_done() const
{
    if (!initCalled_)
    {
        GMX_THROW(gmx::InternalError("ComunicationRecord::init has not been called"));
    }
    if (doneCalled_)
    {
        GMX_THROW(gmx::InternalError("Comunication attempted after done routine was called"));
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
    // Split the non-masters into nmiddlemen X nhelpers
    if (nhelper_per_middleman_ > 0)
    {
        int colour = rank_ / (nhelper_per_middleman_ + 1);
        int key    = rank_ % (nhelper_per_middleman_ + 1);
        check_return("MPI_Comm_split", MPI_Comm_split(mpi_act_world_, colour, key, &mpi_act_helpers_));
    }
    if (MPI_COMM_NULL == mpi_act_helpers_)
    {
        check_return("MPI_Comm_dup", MPI_Comm_dup(mpi_act_world_, &mpi_act_helpers_));
    }
    if (isMaster())
    {
        print(stdout);
    }
    initCalled_ = true;
    return CommunicationStatus::OK;
}

CommunicationRecord::~CommunicationRecord()
{
    if (cr_)
    {
        done_commrec(cr_);
    }
    if (mpi_act_helpers_ != mpi_act_world_ && mpi_act_helpers_ != MPI_COMM_NULL)
    {
        check_return("MPI_Comm_free", MPI_Comm_free(&mpi_act_helpers_));
    }
    if (MPI_COMM_NULL != mpi_act_world_)
    {
        check_return("MPI_Comm_free", MPI_Comm_free(&mpi_act_world_));
    }
}

MPI_Comm CommunicationRecord::create_column_comm(int column, int superior) const
{
    MPI_Comm my_comm = MPI_COMM_NULL;
    if (nhelper_per_middleman() > 0)
    {
        std::set<int> ranks;
        ranks.insert(superior);
        for(int i = 0; i < size(); i++)
        {
            if (i % (1 + nhelper_per_middleman()) == column)
            {
                ranks.insert(i);
            }
        }
        MPI_Group world_group, mygroup;
        check_return("MPI_Comm_group", MPI_Comm_group(MPI_COMM_WORLD, &world_group));
        std::vector<int> iranks;
        std::string      cgroup;
        for(int i : ranks)
        {
            iranks.push_back(i);
            cgroup += gmx::formatString(" %d", i);
        }
        if (msg_handler_)
        {
            msg_handler_->writeDebug(gmx::formatString("Rank: %d cgroup %s\n", rank(), cgroup.c_str()));
        }
        check_return("MPI_Group_incl", MPI_Group_incl(world_group, iranks.size(), iranks.data(), &mygroup));
        check_return("MPI_Comm_create_group", MPI_Comm_create_group(MPI_COMM_WORLD, mygroup, 0, &my_comm));
    }
    if (my_comm == MPI_COMM_NULL)
    {
        check_return("MPI_Comm_dup", MPI_Comm_dup(comm_world(), &my_comm));
    }
    return my_comm;
}


/*************************************************
 *           LOW LEVEL ROUTINES                  *
 *************************************************/

void CommunicationRecord::send_low(int dest, const void *buf, int bufsize) const
{
    int         tag = 0;
    MPI_Status  status;
    MPI_Request req;

    check_init_done();
    check_return("MPI_Isend Failed",
                 MPI_Isend(buf, bufsize, MPI_BYTE, RANK(cr, dest), tag,
                           mpi_act_world_, &req));
    check_return("MPI_Wait failed", MPI_Wait(&req, &status));
}

void CommunicationRecord::recv_low( int src, void *buf, int bufsize) const
{
    int         tag = 0;
    MPI_Status  status;
    MPI_Request req;

    check_init_done();
    check_return ("MPI_Irecv Failed",
                  MPI_Irecv(buf, bufsize, MPI_BYTE, src, tag,
                            mpi_act_world_, &req));
    check_return("MPI_Wait failed", MPI_Wait(&req, &status));
}

//! \brief Broadcase a string
template <> void CommunicationRecord::bcast<std::string>(std::string *str, MPI_Comm comm, int root) const
{
    check_init_done();
    int ssize = str->size();
    check_return("MPI_Bcast", MPI_Bcast((void *)&ssize, 1, MPI_INT, root, comm));
    if (root != rank_)
    {
        str->resize(ssize);
    }
    check_return("MPI_Bcast", MPI_Bcast((void *)str->data(), ssize, MPI_BYTE, root, comm));
}
    
//! \brief Broadcase an int
template <> void CommunicationRecord::bcast<int>(int *i, MPI_Comm comm, int root) const
{
    check_init_done();
    check_return("MPI_Bcast", MPI_Bcast((void *)i, 1, MPI_INT, root, comm));
}

//! \brief Broadcase an unsigned int
template <> void CommunicationRecord::bcast<unsigned int>(unsigned int *i, MPI_Comm comm, int root) const
{
    check_init_done();
    check_return("MPI_Bcast", MPI_Bcast((void *)i, 1, MPI_UNSIGNED, root, comm));
}

//! \brief Broadcase a bool
template <> void CommunicationRecord::bcast<bool>(bool *b, MPI_Comm comm, int root) const
{
    check_init_done();
    int d = *b ? 1 : 0;
    check_return("MPI_Bcast", MPI_Bcast((void *)&d, 1, MPI_INT, root, comm));
    *b = d;
}

//! \brief Broadcase a double
template <> void CommunicationRecord::bcast<double>(double *d, MPI_Comm comm, int root) const
{
    check_init_done();
    check_return("MPI_Bcast", MPI_Bcast((void *)d, 1, MPI_DOUBLE, root, comm));
}

//! \brief Broadcase a size_t
template <> void CommunicationRecord::bcast<size_t>(size_t *d, MPI_Comm comm, int root) const
{
    check_init_done();
    check_return("MPI_Bcast", MPI_Bcast((void *)d, 1, MPI_UNSIGNED_LONG, root, comm));
}

//! \brief Broadcase a vector of doubles
template <> void CommunicationRecord::bcast<std::vector<double>>(std::vector<double> *d,
                                                                 MPI_Comm             comm,
                                                                 int                  root) const
{
    check_init_done();
    int ssize = d->size(); 
    check_return("MPI_Bcast", MPI_Bcast((void *)&ssize, 1, MPI_INT, root, comm));
    if (root != rank_)
    {
        d->resize(ssize);
    }

    check_return("MPI_Bcast", MPI_Bcast((void *)d->data(), ssize, MPI_DOUBLE, root, comm));
}

//! \brief Send a string
template <> void CommunicationRecord::send<std::string>(int dest, const std::string &str) const
{
    check_init_done();
    int len = str.size();

    if (nullptr != msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Sending string '%s' to %d\n", str.c_str(), dest));
    }
    send_low(dest, &len, sizeof(len));
    if (!str.empty())
    {
        send_low(dest, (void *)str.data(), len);
    }
}

//! \brief Receive a string
template <> void CommunicationRecord::recv<std::string>(int src, std::string *str) const
{
    check_init_done();
    int   len;

    recv_low(src, &len, sizeof(len));
    if (len == 0)
    {
        str->clear();
    }
    else
    {
        str->resize(len);
        recv_low(src, (void*)str->data(), len);
    }
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Received string '%s' from %d\n", str->c_str(), src));
    }
}

//! \brief Send a double
template <> void CommunicationRecord::send<double>(int dest, const double &d) const
{
    check_init_done();
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Sending double '%g' to %d\n", d, dest));
    }
    send_low(dest, &d, sizeof(d));
}

//! \brief Receive a double
template <> void CommunicationRecord::recv<double>(int src, double *dd) const
{
    check_init_done();

    recv_low(src, dd, sizeof(*dd));
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Received double '%g' from %d\n", *dd, src));
    }
}

//! \brief Send an int
template <> void CommunicationRecord::send<int>(int dest, const int &d) const
{
    check_init_done();
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Sending int '%d' to %d\n", d, dest));
    }
    send_low(dest, &d, sizeof(d));
}

//! \brief Send an unsigned int
template <> void CommunicationRecord::send<unsigned int>(int dest, const unsigned int &d) const
{
    check_init_done();
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Sending unsigned int '%u' to %d\n", d, dest));
    }
    send_low(dest, &d, sizeof(d));
}

//! \brief Send a size_t
template <> void CommunicationRecord::send<size_t>(int dest, const size_t &d) const
{
    check_init_done();
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Sending size_t '%zu' to %d\n", d, dest));
    }
    send_low(dest, &d, sizeof(d));
}

//! \brief Send a bool
template <> void CommunicationRecord::send<bool>(int dest, const bool &b) const
{
    check_init_done();
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Sending bool '%s' to %d\n", boolName[b], dest));
    }
    int d = b ? 1 : 0;
    send_low(dest, &d, sizeof(d));
}

//! \brief Receive an int
template <> void CommunicationRecord::recv<int>(int src, int *t) const
{
    check_init_done();
    recv_low(src, t, sizeof(*t));
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Received int '%d' from %d\n", *t, src));
    }
}

//! \brief Receive an unsigned int
template <> void CommunicationRecord::recv<unsigned int>(int src, unsigned int *t) const
{
    check_init_done();
    recv_low(src, t, sizeof(*t));
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Received unsigned int '%u' from %d\n", *t, src));
    }
}

//! \brief Receive a size_t
template <> void CommunicationRecord::recv<size_t>(int src, size_t *t) const
{
    check_init_done();
    recv_low(src, t, sizeof(*t));
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Received size_t '%zu' from %d\n", *t, src));
    }
}

//! \brief Receive a bool
template <> void CommunicationRecord::recv<bool>(int src, bool *b) const
{
    check_init_done();
    int d;

    recv_low(src, &d, sizeof(d));
    *b = d;
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Received bool '%s' from %d\n", boolName[*b], src));
    }
}

//! \brief Send a TrainFFMiddlemanMode
template <> void CommunicationRecord::send<TrainFFMiddlemanMode>(int dest, const TrainFFMiddlemanMode &mode) const
{
    check_init_done();
    send(dest, static_cast<int>(mode));
}

//! \brief Receive a TrainFFMiddlemanMode
template <> void CommunicationRecord::recv<TrainFFMiddlemanMode>(int src, TrainFFMiddlemanMode *t) const
{
    check_init_done();
    int i;
    recv(src, &i);
    *t = static_cast<TrainFFMiddlemanMode>(i);
}

//! \brief iMolSelect to int
std::map<iMolSelect, int> imsInt = 
    {
        { iMolSelect::Train, 3 },
        { iMolSelect::Test,  7 },
        { iMolSelect::Ignore, 13 }
    };
    
//! \brief Send a iMolSelect
template <> void CommunicationRecord::send<iMolSelect>(int dest, const iMolSelect &ims) const
{
    if (msg_handler_)
    {
        msg_handler_->writeDebug(gmx::formatString("Sending iMolSelect '%s' to %d\n",
                                                   iMolSelectName(ims).c_str(), dest));
    }
    int d = imsInt[ims];
    send_low(dest, &d, sizeof(d));
}   

//! \brief Receive a iMolSelect
template <> void CommunicationRecord::recv<iMolSelect>(int src, iMolSelect *t) const
{
    int d;

    *t = iMolSelect::Ignore;
    recv_low(src, &d, sizeof(d));
    for(auto &ims : imsInt)
    {
        if (ims.second == d)
        {
            if (msg_handler_)
            {
                msg_handler_->writeDebug(gmx::formatString("Received iMolSelect '%s' from %d\n",
                                                           iMolSelectName(ims.first).c_str(), src));
            }
            *t = ims.first;
            return;
        }
    }
    GMX_THROW(gmx::InternalError(gmx::formatString("Received unknown iMolSelect, int code %d", d).c_str()));
}

//! \brief Send a vector of doubles
template <> void CommunicationRecord::send<std::vector<double>>(int dest,
                                                                const std::vector<double> &d) const
{
    send(dest, static_cast<int>(d.size()));
    if (!d.empty())
    {
        send_low(dest, static_cast<const void *>(d.data()),
                 d.size()*sizeof(double));
    }
}

//! \brief Receive a vector of doubles
template <> void CommunicationRecord::recv<std::vector<double>>(int src,
                                                                std::vector<double> *d) const
{
    int len;
    recv(src, &len);
    d->resize(len);
    if (len > 0)
    {
        recv_low(src, static_cast<void *>(d->data()), len*sizeof(double));
    }
}

//! \brief Send a vector of RVec (e.g. coordinates)
template <> void CommunicationRecord::send<std::vector<gmx::RVec>>(int dest,
                                                                   const std::vector<gmx::RVec> &d) const
{
    send(dest, static_cast<int>(d.size()));
    if (!d.empty())
    {
        send_low(dest, static_cast<const void *>(d.data()),
                 d.size()*sizeof(gmx::RVec));
    }
}

//! \brief Receive a vector of RVecs.
template <> void CommunicationRecord::recv<std::vector<gmx::RVec>>(int src,
                                                                   std::vector<gmx::RVec> *d) const
{
    int len;
    recv(src, &len);
    d->resize(len);
    if (len > 0)
    {
        recv_low(src, static_cast<void *>(d->data()), len*sizeof(gmx::RVec));
    }
}

//! \brief Send a vector of ints
template <> void CommunicationRecord::send<std::vector<int>>(int dest,
                                                             const std::vector<int> &d) const
{
    send(dest, static_cast<int>(d.size()));
    if (!d.empty())
    {
        send_low(dest, static_cast<const void *>(d.data()),
                 d.size()*sizeof(int));
    }
}

//! \brief Receive a vector of ints
template <> void CommunicationRecord::recv<std::vector<int>>(int src,
                                                             std::vector<int> *d) const
{
    int len;
    recv(src, &len);
    d->resize(len);
    if (len > 0)
    {
        recv_low(src, static_cast<void *>(d->data()), len*sizeof(int));
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
    check_return("MPI_Comm_size", MPI_Comm_size(mpi_act_helpers_, &size));
    if (size > 1)
    {
        check_return("MPI_Allreduce", MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM,
                                                    mpi_act_helpers_));
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
    check_return("MPI_Comm_size", MPI_Comm_size(mpi_act_helpers_, &size));
    if (size > 1)
    {
        check_return("MPI_Allreduce", MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM,
                                                    mpi_act_helpers_));
    }
}

#define ACT_SEND_DATA 19823
#define ACT_SEND_DONE -666
CommunicationStatus CommunicationRecord::bcast_data(MPI_Comm comm, int root) const
{
    int acd = ACT_SEND_DATA;
    bcast(&acd, comm, root);

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
    send(dest, ACT_SEND_DATA);

    return CommunicationStatus::SEND_DATA;
}

CommunicationStatus CommunicationRecord::send_done(int dest) const
{
    send(dest, ACT_SEND_DONE);

    return CommunicationStatus::OK;
}

CommunicationStatus CommunicationRecord::bcast_done(MPI_Comm comm, int root) const
{
    int acd = ACT_SEND_DONE;
    bcast(&acd, comm, root);

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
    int kk;
    recv(src, &kk);

    switch (kk)
    {
    case ACT_SEND_DATA:
        return CommunicationStatus::RECV_DATA;
    case ACT_SEND_DONE:
        return CommunicationStatus::DONE;
    default: // throws
        GMX_THROW(gmx::InternalError(gmx::formatString("Received %d from src %d in recv_data. Was expecting either %d or %d\n.", kk, src, (int)ACT_SEND_DATA, (int)ACT_SEND_DONE).c_str()));
    }
    return CommunicationStatus::ERROR;
}

#undef ACT_SEND_DATA
#undef ACT_SEND_DONE

} // namespace alexandria
