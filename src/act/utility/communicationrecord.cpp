/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Stub implementation for non-MPI builds.
 */
#include "communicationrecord.h"

#include <map>
#include <set>

#include "gromacs/utility/exceptions.h"

namespace alexandria
{

//! Map CommunicationStatus to string
std::map<CommunicationStatus, const std::string> csToString = {
    { CommunicationStatus::OK,        "Communication OK"        },
    { CommunicationStatus::DONE,      "Communication Done"      },
    { CommunicationStatus::ERROR,     "Communication Error"     },
    { CommunicationStatus::SEND_DATA, "Communication sent data" },
    { CommunicationStatus::RECV_DATA, "Communication OK"        },
    { CommunicationStatus::TWOINIT,   "Double initiation"       }
};

const std::string &cs_name(CommunicationStatus cs)
{
    return csToString[cs];
}

// Constructors and destructors
CommunicationRecord::CommunicationRecord(MsgHandler *msghandler) : msg_handler_(msghandler)
{
}

CommunicationRecord::~CommunicationRecord()
{
}

CommunicationStatus CommunicationRecord::init(int nmiddleman)
{
    (void)nmiddleman;
    if (initCalled_)
    {
        return CommunicationStatus::TWOINIT;
    }
    initCalled_ = true;
    rank_       = 0;
    size_       = 1;
    nt_         = NodeType::Master;
    return CommunicationStatus::OK;
}

void CommunicationRecord::send_low(int dest, const void *buf, int bufsize) const
{
    (void)dest;
    (void)buf;
    (void)bufsize;
}

void CommunicationRecord::recv_low(int src, void *buf, int bufsize) const
{
    (void)src;
    (void)buf;
    (void)bufsize;
}

void CommunicationRecord::print(FILE *fp)
{
    (void)fp;
}

void CommunicationRecord::check_init_done() const
{
    if (!initCalled_)
    {
        GMX_THROW(gmx::FileIOError("CommunicationRecord::check_init_done: init has not been called"));
    }
    if (doneCalled_)
    {
        GMX_THROW(gmx::FileIOError("CommunicationRecord::check_init_done: done has been called"));
    }
}

MPI_Comm CommunicationRecord::create_column_comm(int column, int superior) const
{
    (void)column;
    (void)superior;
    return MPI_COMM_NULL;
}

CommunicationStatus CommunicationRecord::bcast_data(MPI_Comm comm, int root) const
{
    (void)comm;
    (void)root;
    return CommunicationStatus::OK;
}

CommunicationStatus CommunicationRecord::send_data(int dest) const
{
    (void)dest;
    return CommunicationStatus::OK;
}

CommunicationStatus CommunicationRecord::bcast_done(MPI_Comm comm, int root) const
{
    (void)comm;
    (void)root;
    return CommunicationStatus::OK;
}

CommunicationStatus CommunicationRecord::send_done(int dest) const
{
    (void)dest;
    return CommunicationStatus::OK;
}

CommunicationStatus CommunicationRecord::recv_data(int src) const
{
    (void)src;
    return CommunicationStatus::OK;
}

void CommunicationRecord::sumd_helpers(int nr, double r[]) const
{
    (void)nr;
    (void)r;
}

void CommunicationRecord::sumi_helpers(int nr, int r[]) const
{
    (void)nr;
    (void)r;
}

// Template implementations for send, recv, and bcast
template<typename T>
void CommunicationRecord::send(int dest, const T &t) const
{
    (void)dest;
    (void)t;
}

template<typename T>
void CommunicationRecord::recv(int src, T *t) const
{
    (void)src;
    (void)t;
}

template<typename T>
void CommunicationRecord::bcast(T *t, MPI_Comm comm, int root) const
{
    (void)t;
    (void)comm;
    (void)root;
}

// Explicit template instantiations
template void CommunicationRecord::send<int>(int, const int &) const;
template void CommunicationRecord::send<double>(int, const double &) const;
template void CommunicationRecord::send<bool>(int, const bool &) const;
template void CommunicationRecord::send<unsigned long>(int, const unsigned long &) const;
template void CommunicationRecord::send<unsigned int>(int, const unsigned int &) const;
template void CommunicationRecord::send<std::string>(int, const std::string &) const;
template void CommunicationRecord::send<std::vector<double>>(int, const std::vector<double> &) const;
template void CommunicationRecord::send<std::vector<int>>(int, const std::vector<int> &) const;
template void CommunicationRecord::send<std::vector<gmx::BasicVector<double>>>(int, const std::vector<gmx::BasicVector<double>> &) const;
template void CommunicationRecord::send<TrainFFMiddlemanMode>(int, const TrainFFMiddlemanMode &) const;
template void CommunicationRecord::send<iMolSelect>(int, const iMolSelect &) const;

template void CommunicationRecord::recv<int>(int, int *) const;
template void CommunicationRecord::recv<double>(int, double *) const;
template void CommunicationRecord::recv<bool>(int, bool *) const;
template void CommunicationRecord::recv<unsigned long>(int, unsigned long *) const;
template void CommunicationRecord::recv<unsigned int>(int, unsigned int *) const;
template void CommunicationRecord::recv<std::string>(int, std::string *) const;
template void CommunicationRecord::recv<std::vector<double>>(int, std::vector<double> *) const;
template void CommunicationRecord::recv<std::vector<int>>(int, std::vector<int> *) const;
template void CommunicationRecord::recv<std::vector<gmx::BasicVector<double>>>(int, std::vector<gmx::BasicVector<double>> *) const;
template void CommunicationRecord::recv<TrainFFMiddlemanMode>(int, TrainFFMiddlemanMode *) const;
template void CommunicationRecord::recv<iMolSelect>(int, iMolSelect *) const;

template void CommunicationRecord::bcast<int>(int *, MPI_Comm, int) const;
template void CommunicationRecord::bcast<double>(double *, MPI_Comm, int) const;
template void CommunicationRecord::bcast<bool>(bool *, MPI_Comm, int) const;
template void CommunicationRecord::bcast<unsigned long>(unsigned long *, MPI_Comm, int) const;
template void CommunicationRecord::bcast<unsigned int>(unsigned int *, MPI_Comm, int) const;
template void CommunicationRecord::bcast<std::string>(std::string *, MPI_Comm, int) const;
template void CommunicationRecord::bcast<std::vector<double>>(std::vector<double> *, MPI_Comm, int) const;
template void CommunicationRecord::bcast<std::vector<int>>(std::vector<int> *, MPI_Comm, int) const;
template void CommunicationRecord::bcast<std::vector<gmx::BasicVector<double>>>(std::vector<gmx::BasicVector<double>> *, MPI_Comm, int) const;
template void CommunicationRecord::bcast<TrainFFMiddlemanMode>(TrainFFMiddlemanMode *, MPI_Comm, int) const;
template void CommunicationRecord::bcast<iMolSelect>(iMolSelect *, MPI_Comm, int) const;

} // namespace alexandria
