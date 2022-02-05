#ifndef ACT_COMMUNICATIONRECORD_H
#define ACT_COMMUNICATIONRECORD_H

#include <string>
#include <vector>

#include "ga/Dataset.h"
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
    RECV_DATA = 9999
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
    //! MPI Communicator for every processor except the master
    MPI_Comm          mpi_act_not_master_    = MPI_COMM_NULL;
    //! MPI Communicator for my local helpers
    MPI_Comm          mpi_act_helpers_       = MPI_COMM_NULL;
    //! My MPI rank
    int               rank_                  = 0;
    //! Total number of MPI ranks
    int               size_                  = 0;
    //! Total number of middlemen
    /*! \brief The number of middle men in a 2D parallellization
     * Of all the size_ nodes, the middlemen each have
     * nhelper_per_middleman_ = size_/nmiddlemen_
     * helpers. 
     */
    int               nmiddlemen_            = 0;
    //! Number of helpers per middleman
    int               nhelper_per_middleman_ = 0;
    //! My source for receiving instructions and reporting back to
    int               superior_              = -1;
    //! If I am a middleman, this is my ordinal in the ranks of middlemen
    int               ordinal_               = -1;
    //! My helpers (if any)
    std::vector<int>  helpers_;
    //! All the middlemen, as if anyone should care
    std::vector<int>  middlemen_;
    
    /*************************************************
     *           LOW LEVEL ROUTINES                  *
     *************************************************/

    /*! Send a buffer to another processor.
     * \param[in] dest    The destination processor
     * \param[in] buf     Pointer to the byte-buffer
     * \param[in] bufsize The size of the buffer in bytes
     */
    void send(int dest, const void *buf, int bufsize) const;

    /*! Receive a buffer from another processor.
     * \param[in] src     The source processor
     * \param[in] buf     Pointer to the byte-buffer
     * \param[in] bufsize The size of the buffer in bytes
     */
    void recv(int src, void *buf, int bufsize) const;
    /*! \brief Print the values to a file
     * \param[in] fp The file pointer to print to
     */
    void print(FILE *fp);
    
public:
    //! \brief Constructor
    CommunicationRecord();

    //! \brief Destructor is needed to get rid of commrec
    ~CommunicationRecord();

    /*! \brief Initiate the internal data  once and for all.
     * Only constant data can be extracted from this.
     * \param[in] nmiddlemen The number of middlemen. Knowing the total 
     *                       number of cores is then sufficient to derive
     *                       the rest.
     */
    void init(int nmiddleman);
    
    //! \return my MPI rank
    int rank() const { return rank_; }

    //! \return number of MPI ranks
    int size() const { return size_; }
    
    //! \return whether this is a parallel calculation
    bool isParallel() const { return size_ > 0; }
    
    //! \return my helper nodes, or empty vector if none
    const std::vector<int> &helpers() const { return helpers_; }

    //! \return all the middlemen, or empty vector if none
    const std::vector<int> &middlemen() const { return middlemen_; }

    //! \return my superior node, or -1 if none
    int superior() const { return superior_; }

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

    /*! Send a string to another processor.
     * \param[in] dest The destination processor
     * \param[in] str  Pointer to the string
     */
    void send_str(int dest, const std::string *str) const;
    
    /*! Receive a string from another processor.
     * \param[in] src The destination processor
     * \param[in] str Pointer to the string
     */
    void recv_str(int src, std::string *str) const;
    
    /*! Send a double to another processor.
     * \param[in] dest The destination processor
     * \param[in] d    The value of the double
     */
    void send_double(int dest, double d) const;
    
    /*! \brief Receive a double
     * \param[in] src The source processor
     * \return The double received
     */
    double recv_double(int src) const;
    
    /*! Send a double vector to another processor.
     * \param[in] dest The destination processor
     * \param[in] d    Pointer to vector of the doubles
     */
    void send_double_vector(int dest, const std::vector<double> *d) const;
    
    /*! Receive a double vector from another processor.
     * \param[in] dest The source processor
     * \param[in] d    Pointer to vector of the doubles
     */
    void recv_double_vector(int src, std::vector<double> *d) const;
    
    /*! Send an int to another processor.
     * \param[in] dest The destination processor
     * \param[in] d    The value of the integer
     */
    void send_int(int dest, int d) const;
    
    /*! \brief Receive an int
     * \param[in] src The source processor
     * \return The int received
     */
    int recv_int(int src) const;
 
    /*! Send an iMolSelect to another processor.
     * \param[in] dest The destination processor
     * \param[in] ims  The value of the iMolSelect
     */
    void send_iMolSelect(int dest, iMolSelect ims) const;
    
    /*! \brief Receive an iMolSelect
     * \param[in] src The source processor
     * \return The iMolSelect received
     * \throws if something is wrong
     */
    iMolSelect recv_iMolSelect(int src) const;
 
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
    /*! \brief Initiate sending data to a processor
     * \param[in] dest The destination processor
     * \return CommunicationStatus::SEND_DATA, if OK
     */
    CommunicationStatus send_data(int dest) const;
    
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
