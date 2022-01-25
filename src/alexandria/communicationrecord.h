#ifndef ACT_COMMUNICATIONRECORD_H
#define ACT_COMMUNICATIONRECORD_H

#include <vector>

#include "gromacs/mdtypes/commrec.h"

namespace alexandria
{

enum class NodeType {
    Master, MiddleMan, Helper
};

class CommunicationRecord
{
private:
    //! The GROMACS data structure
    t_commrec        *cr_                    = nullptr;
    //! My NodeType
    NodeType          nt_                    = NodeType::Helper;
    //! My MPI rank
    int               rank_                  = 0;
    //! Total number of MPI ranks
    int               size_                  = 0;
    //! Total number of middlemen
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
    
    //! \return whether I am the master
    bool isMaster() const { return NodeType::Master == nt_; }

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
};

} // namespace alexandria

#endif
