#include <vector>

#include "act/qgen/qgen_acm.h"
#include "act/molprop/fragment.h"

namespace alexandria
{

    /*! Class to control charge generation for fragments
     */
    class FragmentHandler
    {
    private:
        //! We need one ACM structure for each fragment
        std::vector<QgenAcm>               QgenAcm_;
        //! A complete topology for each fragment is needed to compute energies
        std::vector<Topology>              topologies_;
        //! And a supporting atoms structure too
        std::vector<std::vector<ActAtom> > FragAtoms_;
        //! And a vector of bonds
        std::vector<std::vector<Bond> >    bonds_;
        //! Array denoting where the atoms start in the global
        std::vector<size_t>                atomStart_;
        //! Total number of atoms
        size_t                             natoms_ = 0;
    public:
        /*! Constructor
         * \param[in] pd            Force field data
         * \param[in] coordinates   The atomic coordinates
         * \param[in] atoms         The atoms
         * \param[in] bonds         The bonds
         * \param[in] fragments     The fragmentation information
         * \param[in] shellRenumber Info on renumbering atoms because of shells
         */
        FragmentHandler(const Poldata                    *pd,
                        const gmx::HostVector<gmx::RVec> &coordinates,
                        const std::vector<ActAtom>       &atoms,
                        const std::vector<Bond>          &bonds,
                        const std::vector<Fragment>      *fragments,
                        const std::vector<int>           &shellRenumber);

        /*! \brief Fetch charges for all atoms
         * \param[out] qq Vector that will be reinitialized at correct lenght
         */
        void fetchCharges(std::vector<double> *qq);
        
        //! \return the atomStart_ vector, containing one extra index more than the number of fragments.
        const std::vector<size_t> atomStart() const { return atomStart_; }

        //! \return the vector of Topology structures        
        const std::vector<Topology> topologies() const { return topologies_; }

        /*! \brief Generate charges for all fragments
         * \param[in]  fp      Debug file pointer, may be nullptr
         * \param[in]  molname Molecule name for printing
         * \param[in]  x       Atomic coordinates
         * \param[in]  pd      Force field file
         * \param[out] atoms   Atoms to store the charges in
         * \return Error code or OK if all is fine.
         */
        eQgen generateCharges(FILE                             *fp,
                              const std::string                &molname,
                              const gmx::HostVector<gmx::RVec> &x,
                              const Poldata                    *pd,
                              std::vector<ActAtom>             *atoms);
    };
} // namespace alexandria
