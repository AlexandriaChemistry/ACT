/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2024
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
#ifndef ACT_COMMUNICATIONRECORD_H
#define ACT_COMMUNICATIONRECORD_H

#include <string>
#include <vector>

#include "act/basics/dataset.h"
#include "gromacs/mdtypes/commrec.h"

namespace alexandria
{

/*! \brief Distinguish the node types in use in ACT
 */
enum class NodeType {
    //! Master, also known as the Overlord
    Master,
    //! The process that mediates contacts and does some calculations too
    MiddleMan,
    //! The helper process, only working on low-level calculations
    Helper
};

/*! \brief
 * Enumerated type holding the result status of communication operations
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum CommunicationStatus {
    OK        = 6666,
    DONE      = 6789,
    ERROR     = 7777,
    SEND_DATA = 8888,
    RECV_DATA = 9999,
    TWOINIT   = 4567
};

/*!
 * \brief Operation modes of the middleman in train_ff
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum class TrainFFMiddlemanMode {
    FITNESS  = 1999,  // Compute fitness and send fitness back
    MUTATION = 2022   // Mutate, compute fitness, send parameters and fitness back
};

//! \return A string corresponding to a CommunicationStatus
const char *cs_name(CommunicationStatus cs);

class CommunicationRecord
{
private:
    //! The GROMACS data structure
    t_commrec        *cr_                    = nullptr;
    //! My NodeType
    NodeType          nt_                    = NodeType::Helper;
    //! MPI Communicator for the whole system
    MPI_Comm          mpi_act_world_         = MPI_COMM_NULL;
    /*! \brief MPI Communicator for my local helpers.
     * This corresponds to the rows in the matrix of processors
     */
    MPI_Comm          mpi_act_helpers_       = MPI_COMM_NULL;
    //! My MPI rank
    int               rank_                  = 0;
    //! Total number of MPI ranks
    int               size_                  = 0;
    //! Total number of middlemen (includes the MASTER node)
    /*! \brief The number of middle men in a 2D parallellization
     * Of all the size_ nodes, the middlemen each have
     * nhelper_per_middleman_ = size_/nmiddlemen_ - 1
     * helpers. 
     */
    int               nmiddlemen_            = 0;
    //! Number of helpers per middleman
    int               nhelper_per_middleman_ = 0;
    //! My source for receiving instructions and reporting back to
    int               superior_              = 0;
    //! If I am a middleman, this is my ordinal in the ranks of middlemen
    int               ordinal_               = -1;
    //! My helpers (if any)
    std::vector<int>  helpers_;
    //! All the middlemen (not including the MASTER), as if anyone should care
    std::vector<int>  middlemen_;
    //! Check whether init has been called
    bool              initCalled_            = false;
    //! Check whether done has been called
    bool              doneCalled_            = false;
    /*************************************************
     *           LOW LEVEL ROUTINES                  *
     *************************************************/

    /*! Send a buffer to another processor.
     * \param[in] dest    The destination processor
     * \param[in] buf     Pointer to the byte-buffer
     * \param[in] bufsize The size of the buffer in bytes
     */
    void send_low(int dest, const void *buf, int bufsize) const;

    /*! Receive a buffer from another processor.
     * \param[in] src     The source processor
     * \param[in] buf     Pointer to the byte-buffer
     * \param[in] bufsize The size of the buffer in bytes
     */
    void recv_low(int src, void *buf, int bufsize) const;
    /*! \brief Print the values to a file
     * \param[in] fp The file pointer to print to
     */
    void print(FILE *fp);
    //! Check whether init has been called and done has not been called and throw if not
    void check_init_done() const;
    
public:
    //! \brief Constructor
    CommunicationRecord();

    //! \brief Destructor is needed to get rid of commrec
    ~CommunicationRecord();

    /*! \brief Initiate the internal data once and for all.
     * This routine should be called exactly once. 
     * Only constant data can be extracted from this.
     * \param[in] nmiddlemen The number of middlemen. Knowing the total 
     *                       number of cores is then sufficient to derive
     *                       the rest.
     * \return Outcome of the initiation.
     */
    CommunicationStatus init(int nmiddleman);
    
    /*! \brief Tell the communication record that we are done.
     * Any additional communication will cause a throw.
     */
    void done()
    {
        doneCalled_ = true;
    }
    
    //! \return my MPI rank
    int rank() const { return rank_; }

    //! \return number of MPI ranks
    int size() const { return size_; }
    
    //! \return whether this is a parallel calculation
    bool isParallelCalc() const { return nhelper_per_middleman_ > 0; }

    //! \return whether there is more than one processor
    bool isParallel() const { return size_ > 1; }
    
    //! \return my helper nodes, or empty vector if none
    const std::vector<int> &helpers() const { return helpers_; }

    //! \return all the middlemen except the master, or empty vector if none
    const std::vector<int> &middlemen() const { return middlemen_; }
    
    //! \return communicator for sending to helpers
    MPI_Comm send_helpers() const { return mpi_act_helpers_; }

    //! \return communicator for all cores
    MPI_Comm comm_world() const { return mpi_act_world_; }
    
    /*! \brief Create special communicator for one column in the
     * the node matrix. It includes the master node, allowing to
     * broadcast. The communicator should be freed by the caller.
     * \param[in] column The column number
     * \param[in] superior The master or middleman (if wanted). Default master.
     * \return communicator
     */
    MPI_Comm create_column_comm(int column, int superior = 0) const;

    //! \return communicator for subsystem
    // MPI_Comm comm_subsystem() const { return mpi_act_subsystem_; }

    //! \return my superior node, or -1 if none
    int superior() const { return superior_; }

    //! \return my superior node, or -1 if none
    int superiorOrMeIfMiddleman() const { return 0; } //if (isMiddleMan()) { return rank_; } else { return superior_; } }

    //! \return number of middlemen
    int nmiddlemen() const { return nmiddlemen_; }
    
    //! \return number of helpers per middlemen
    int nhelper_per_middleman() const { return nhelper_per_middleman_; }
    
    //! \return my ordinal in the rank of middlemen or -1 if I am not one
    int middleManOrdinal() const { return ordinal_; }

    //! \return my NodeType
    NodeType nodeType() const { return nt_; }
    
    //! \return whether I am the master or a middleman
    bool isMasterOrMiddleMan() const { return NodeType::Helper != nt_; }

    //! \return whether I am the master
    bool isMaster() const { return NodeType::Master == nt_; }

    //! \return whether I am a middleman
    bool isMiddleMan() const { return NodeType::MiddleMan == nt_; }

    //! \return whether I am the master
    bool isHelper() const { return NodeType::Helper == nt_; }

    //! \return the GROMACS communication record
    const t_commrec *commrec() const { return cr_; }
    
    /*************************************************
     *           LOW LEVEL ROUTINES                  *
     *************************************************/
    /*! Broadcast data to helpers or all processors from the master.
     * \param[inout] str  Pointer to the string
     * \param[in]    comm MPI communicator
     * \param[in]    root Who is the root of this communicatione
     */
    template<typename T>
    void bcast(T *t, MPI_Comm comm, int root=0) const;

    /*! Send data to another processor.
     * \param[in] dest The destination processor
     * \param[in] t    Pointer to the data
     */
    template<typename T>
    void send(int dest, const T &t) const;
    
    /*! Receive data from another processor.
     * \param[in] src The destination processor
     * \param[in] t   Pointer to the data
     */
    template<typename T>
    void recv(int src, T *t) const;
    
    /*! \brief Calculate the global sum of an array of doubles
     * but only on my act helpers.
     * \param[in]    nr The number of doubles
     * \param[inout] r  The list of doubles input and output
     */
    void sumd_helpers(int    nr, 
                      double r[]) const;
    
    /*! \brief Calculate the global sum of an array of ints 
     * but only on my act helper group
     * \param[in] nr The number of ints
     * \param[in] r  The ints.
     */
    void sumi_helpers(int nr, 
                      int r[]) const;

    /*********************************************************
     * Routines to initiate and finalize data transmissions
     *********************************************************/
    /*! \brief Initiate broadcasting data to a processor
     * \param[in] comm The MPI communicator
     * \param[in] root Who is the root for this communication
     * \return CommunicationStatus::OK, if OK
     */
    CommunicationStatus  bcast_data(MPI_Comm comm, int root=0) const;
    
    /*! \brief Initiate sending data to a processor
     * \param[in] dest The destination processor
     * \return CommunicationStatus::SEND_DATA, if OK
     */
    CommunicationStatus send_data(int dest) const;
    
    /*! \brief Finish broadcasting data to a processor
     * \param[in] comm The MPI communicator
     * \param[in] root Who is the root for this communication
     * \return CommunicationStatus::OK, if OK
     */
    CommunicationStatus  bcast_done(MPI_Comm comm, int root=0) const;
    
    /*! \brief Finalize sending data to a processor
     * \param[in] dest The destination processor
     * \return CommunicationStatus::OK, if OK
     */
    CommunicationStatus send_done(int dest) const;

    /*! \brief Initiate or finalize receiving data from a processor
     * \param[in] src The source processor
     * \return CommunicationStatus::RECV_DATA, if OK
     */
    CommunicationStatus recv_data(int src) const;

};

} // namespace alexandria

#endif
