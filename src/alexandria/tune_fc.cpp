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

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <map>
#include <random>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"

#include "alex_modules.h"
#include "communication.h"
#include "gentop_core.h"
#include "gmx_simple_comm.h"
#include "molgen.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "molselect.h"
#include "mymol.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "regression.h"
#include "tuning_utility.h"
#include "tune_fc_utils.h"

/*! \brief Write a csv file containing molecule names and bond energy
 *
 * Writes the whole bond energy matrix.
 */
static void dump_csv(const std::vector<std::string>        &ctest,
                     const std::vector<alexandria::MyMol>  &mm,
                     const std::vector<int>                &ntest,
                     const std::vector<double>             &Edissoc,
                     const MatrixWrapper                   &a,
                     const double                           x[])
{
    FILE *csv = gmx_ffopen("tune_fc.csv", "w");
    fprintf(csv, ",");
    for (auto j : ctest)
    {
        fprintf(csv, "%s,", j.c_str());
    }
    fprintf(csv, "\n");
    int i = 0;
    for (auto &mymol : mm)
    {
        fprintf(csv, "%s,", mymol.getMolname().c_str());
        for (size_t j = 0; (j < ctest.size()); j++)
        {
            fprintf(csv, "%g,", a.get(j, i));
        }
        fprintf(csv, "%.3f\n", x[i]);
        i++;
    }
    fprintf(csv, "Total,");
    for (auto j : ntest)
    {
        fprintf(csv, "%d,", j);
    }
    fprintf(csv, "\n");
    fprintf(csv, "Edissoc,");
    for (auto j : Edissoc)
    {
        fprintf(csv, "%.3f,", j);
    }
    fprintf(csv, "\n");
    fclose(csv);
}

namespace alexandria
{

typedef struct MolEnergyEntry
{
    // Job type
    int    jobType;
    // Reference energy
    double EReference;
    // Alexandria energy
    double EAlexandria;
} MolEnergyEntry;

class MolEnergy
{
    public:
        MolEnergy() : Force2_(0.0) {}

        void setForce2(double Force2) { Force2_ = Force2; }

        double force2() const { return Force2_; }

        void setTerms(real term[F_NRE]);

        real term(int ftype) const { return term_[ftype]; }

        void addE(int jobType, double Ereference, double Ealexandria);

        void clear();

        const std::vector<MolEnergyEntry> &entries() { return MolEnergyEntry_; }
    private:
        // Sum of forces squared at the OPT geometry
        double                      Force2_;
        // Array of energies at all geometries
        std::vector<MolEnergyEntry> MolEnergyEntry_;
        // Array of energy components at the OPT geometry
        real                        term_[F_NRE] = { 0 };
};

void MolEnergy::addE(int jobType, double Ereference, double Ealexandria)
{
    MolEnergyEntry_.push_back({ jobType, Ereference, Ealexandria });
}

void MolEnergy::clear()
{
    Force2_ = 0;
    for (int i = 0; i < F_NRE; i++)
    {
        term_[i] = 0;
    }
    MolEnergyEntry_.clear();
}

void MolEnergy::setTerms(real term[F_NRE])
{
    for (int i = 0; i < F_NRE; i++)
    {
        term_[i] = term[i];
    }
}

class Optimization : public MolGen, Bayes
{
    using param_type = std::vector<double>;

    private:
        std::map<InteractionType, ForceConstants> ForceConstants_;
        std::vector<PoldataUpdate>                poldataUpdates_;
        bool                        bDissoc_          = true;
        real                        w_dhf_            = 1;
        const char                 *lot_              = nullptr;
        bool                        calcAll_          = false;
        // Map from molid to MolEnergy
        std::map<int, MolEnergy>    MolEnergyMap_;
        std::string                 outputFile_;
    public:


        /*! \brief
         *
         * Constructor
         */
        Optimization() {};

        /*! \brief
         * Add command line options.
         * \param[inout] pargs Command line options
         */
        void add_pargs(std::vector<t_pargs> *pargs);

        //! \brief To be called after processing options
        void optionsFinished(const std::string &outputFile)
        {
            MolGen::optionsFinished();
            setBoxConstraint(bConstrain());
            outputFile_ = outputFile;
        }

        void saveState();
        
        void setCalcAll(bool calcAll) { calcAll_ = calcAll; }

        //! \brief Return the level of theory used as reference
        const char *lot() const { return lot_; }

        /*! \brief
         *
         * Check whether molecules are supported by the force field.
         * Check whether all the atomtypes in the molecules are present
         * in the force field file. If not the molecules are removed.
         *
         * Routine also divides molecules over processors.
         */
        void checkSupport(FILE *fp);

        /*! \brief
         *
         * Fill parameter vector using ForceConstants which
         * is built based on Poldata. Will generate parameters
         * within the bounds set or free, depening in the input in gentop.dat file.
         * \param[in] bRandom Whether or not to generate random parameters 
         */
        void polData2TuneFc(bool bRandom);

        /*! \brief
         *
         * Copy the optimized parameters back to Poldata
         *
         * \param[in] changed Boolean list stating whether a parameter has changed
         */
        void toPoldata(const std::vector<bool> &changed);

        /*! \brief
         * Broadcast changes in Poldata to the
         * slaves when in parallel.
         */
        void broadcastPoldataUpdate();

        /*! \brief
         *
         * Compute the dissociation energies for all the bonds.
         * Given all the bonds and the enthalpies of formation of all
         * molecules, we can approximate the dissociation enthalpy (D0 in the
         * Morse potential by least squares fitting the D0 to reproduce the
         * molecular energy (Delta H formation of molecule - Delta H formation of
         * the atoms). This is a crude approximation since all other energy
         * terms in the force field are ignored, however the dissociation
         * energy is the largest contribution to the molecular energy.
         */
        void getDissociationEnergy(FILE *fplog);

        /*! \brief
         * Initialize the optimization algorithm.
         * \param[in] fplog Log file for dumping information.
         */
        void InitOpt(FILE *fplog, bool bRandom);

        /*! \brief
         * Do the actual optimization.
         * \param[in] fplog  FILE pointer for logging
         * \param[in] oenv   Output environment for managing xvg files etc.
         * \param[in] xvgconv Output file monitoring parameters
         * \param[in] xvgepot Output file monitoring penalty function
         * \return true if better parameters were found.
         */
        bool optRun(FILE                   *fplog,
                    const gmx_output_env_t *oenv,
                    const char             *xvgconv,
                    const char             *xvgepot);

        /*! \brief
         * Print the results of the optimization per molecule
         * as well as statistics over the minimum energy
         * structures.
         * \param[in] fp     File pointer to write to
         * \param[in] title  Text string to write
         * \param[in] xvg    Filename for correlation plot for DHformation at
         *                   minimum energy structures.
         * \param[in] HF_xvg Filename for energy correlation per compound
         * \param[in] oenv   GROMACS file management structure.
         */
        void printResults(FILE                   *fp,
                          const char             *title,
                          const char             *xvg,
                          const char             *HF_xvg,
                          const gmx_output_env_t *oenv);

        /*! \brief
         * Print data for all the molecules.
         * \param[in] fp     File pointer to write to
         * \param[in] bForce Print forces and energy terms at optimal geometry
         * \param[in] bMtop  Print the molecular topolog7
         */
        void printMolecules(FILE *fp,
                            bool  bForce,
                            bool  bMtop);
        /*! \brief
         * Compute the deviation from the reference energies etc.
         */
        virtual double calcDeviation();

        /*! \brief
         * Send the information from this processor to destination
         * \param[in] cr   Communication data structure
         * \param[in] dest Destination processor
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receive information from another processor to this
         * \param[in] cr   Communication data structure
         * \param[in] src Source processor
         */
        CommunicationStatus Receive(t_commrec *cr, int src);

        //! Distribute all information over all processors
        void broadcast();
};

CommunicationStatus Optimization::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, ForceConstants_.size());
        for (auto &fc : ForceConstants_)
        {
            std::string itype(interactionTypeToString(fc.first));
            gmx_send_str(cr, dest, &itype);
            fc.second.Send(cr, dest);
        }
    }
    return cs;
}

CommunicationStatus Optimization::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        int nfc = gmx_recv_int(cr, src);

        for (int n = 0; (CS_OK == cs) && (n < nfc); n++)
        {
            std::string itype;
            gmx_recv_str(cr, src, &itype);
            alexandria::ForceConstants fc;
            cs = fc.Receive(cr, src);
            if (CS_OK == cs)
            {
                ForceConstants_.insert({stringToInteractionType(itype), fc});
            }
        }
    }
    return cs;
}

void Optimization::add_pargs(std::vector<t_pargs> *pargs)
{
    t_pargs pa[] =
    {
        { "-dissoc",  FALSE, etBOOL, {&bDissoc_},
          "Derive dissociation energy from the enthalpy of formation. If not chosen, the dissociation energy will be read from the gentop.dat file." },
        { "-weight_dhf", FALSE, etREAL, {&w_dhf_},
          "Fitting weight of the minimum energy structure, representing the enthalpy of formation relative to high energy structures." }
    };
    for (size_t i = 0; i < sizeof(pa)/sizeof(pa[0]); i++)
    {
        pargs->push_back(pa[i]);
    }
    addOptions(pargs, etuneFC);
    Bayes::add_pargs(pargs);
}

void Optimization::broadcastPoldataUpdate()
{
    t_commrec *cr = commrec();
    if (!PAR(cr))
    {
        return;
    }

    if (MASTER(cr))
    {
        for (int i = 1; i < cr->nnodes; i++)
        {
            gmx_send_int(cr, i, static_cast<int>(poldataUpdates_.size()));
            for (auto &p : poldataUpdates_)
            {
                p.Send(cr, i);
            }
        }
    }
    else
    {
        int n = gmx_recv_int(cr, 0);
        poldataUpdates_.resize(n);
        for (int i = 0; i < n; i++)
        {
            poldataUpdates_[i].Receive(cr, 0);
        }
        if (debug)
        {
            fprintf(debug, "broadcastPoldataUpdate: received %d updates\n",
                    n);
        }
    }
}

void Optimization::broadcast()
{
    t_commrec *cr = commrec();
    if (!PAR(cr))
    {
        return;
    }
    const int src = 0;
    if (MASTER(cr))
    {
        for (int dest = 1; dest < cr->nnodes; dest++)
        {
            Send(cr, dest);
        }
    }
    else
    {
        if (nullptr != debug)
        {
            fprintf(debug, "Updating ForceConstants on node %d\n",
                    cr->nodeid);
        }
        Receive(cr, src);
    }
    poldata()->broadcast(cr);
}

void Optimization::checkSupport(FILE *fp)
{
    auto ntotal  = mymols().size();
    auto nlocal  = 0;

    for (auto mymol =  mymols().begin(); mymol <  mymols().end(); )
    {
        if (mymol->eSupp_ != eSupport::Local)
        {
            mymol++;
            continue;
        }

        bool bSupport = true;
        for (auto &bt : iopt()) 
        {
            if (bSupport && bt.second)
            {
                auto iType = bt.first;
                auto fs    = poldata()->findForces(iType);
                int ft     = fs->fType();
                bSupport   = (mymol->ltop_ != nullptr);
                for (int i = 0; (i < mymol->ltop_->idef.il[ft].nr) && bSupport;
                     i += interaction_function[ft].nratoms+1)
                {
                    // Loop starts from 1 because the first value is the function type
                    std::vector<std::string> atoms;
                    auto myatoms = mymol->atoms();
                    for(int j = 1; j <= interaction_function[ft].nratoms && bSupport; j++)
                    {
                        std::string aa;
                        int         ai = mymol->ltop_->idef.il[ft].iatoms[i+j];
                        if (!poldata()->atypeToBtype(*myatoms->atomtype[ai], &aa))
                        {
                            bSupport = false;
                        }
                        else
                        {
                            atoms.push_back(aa);
                        }
                    }
                    if (bSupport)
                    {
                        // TODO: Check whether true is correct
                        Identifier bondId(atoms, CanSwap::Yes);
                        bSupport = fs->parameterExists(bondId);
                    }
                }
            }
        }
        if (!bSupport)
        {
            fprintf(stderr, "No force field support for %s\n",
                    mymol->getMolname().c_str());
            mymol =  mymols().erase(mymol);
        }
        else
        {
            mymol++;
            nlocal++;
        }
    }
    if (PAR(commrec()))
    {
        gmx_sumi(1, &nlocal, commrec());
    }
    if (nullptr != fp)
    {
        fprintf(fp, "%d out of %zu molecules have support in the force field.\n",
                nlocal, ntotal);
    }
}

void Optimization::polData2TuneFc(bool bRandom)
{
    for (auto &fc : ForceConstants_)
    {
        auto fs = poldata()->findForces(fc.first);
        for (auto &b : fc.second.bondNamesConst())
        {
            ForceFieldParameterMap *ptr = fs->findParameters(b.first);
            for (auto &p : *ptr)
            { 
                Bayes::addParam("test",
                                p.second.value(), p.second.minimum(),
                                p.second.maximum(), bRandom);
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "poldata2TuneFc: Copied %zu parameters.\n",
                Bayes::nParam());
    }
}

void Optimization::saveState()
{
    writePoldata(outputFile_, poldata(), false);
}

void Optimization::toPoldata(const std::vector<bool> &changed)
{
    size_t n   = 0;
    poldataUpdates_.clear();
    auto param = Bayes::getParam();
    for (auto &fc : ForceConstants_)
    {
        const auto iType = fc.first;
        for (auto &b : fc.second.bondNames())
        {
            std::vector<double> parameterValues;
            bool                bondChanged = false;
            for (size_t p = 0; p < b.second.nParams(); p++)
            {
                bondChanged = bondChanged || changed[n];
                parameterValues.push_back(param[n++]);
            }
            if (bondChanged)
            {
                b.second.setParameterValues(parameterValues);
                poldataUpdates_.push_back(PoldataUpdate(iType, b.first,
                                                        parameterValues));
            }
        }
    }
    GMX_RELEASE_ASSERT(n == param.size(), "Number of parameters set should be equal to the length of the parameter array");
}

void Optimization::getDissociationEnergy(FILE *fplog)
{
    std::vector<double>         rhs;
    std::vector<int>            ntest;
    std::vector<std::string>    ctest;

    int nD   = ForceConstants_[InteractionType::BONDS].nbad();
    int nMol = mymols().size();

    if ((0 == nD) || (0 == nMol))
    {
        gmx_fatal(FARGS, "Number of variables is %d and number of molecules is %d",
                  nD, nMol);
    }

    MatrixWrapper a(nD, nMol);
    MatrixWrapper a_copy(nD, nMol);
    ntest.resize(nD, 0);
    ctest.resize(nD);

    fprintf(fplog, "There are %d different bondtypes to optimize the heat of formation\n", nD);
    fprintf(fplog, "There are %d (experimental) reference heat of formation.\n", nMol);

    auto fs  = poldata()->findForces(InteractionType::BONDS);
    auto j   = 0;

    for (auto &mymol :  mymols())
    {
        auto myatoms = mymol.atoms();
        for (auto &b : mymol.bondsConst())
        {
            const char *atypeI = *myatoms->atomtype[b.getAi()];
            const char *atypeJ = *myatoms->atomtype[b.getAj()];
            std::string btypeI, btypeJ;
            if (poldata()->atypeToBtype(atypeI, &btypeI) &&
                poldata()->atypeToBtype(atypeJ, &btypeJ))
            {
                Identifier bondId({btypeI, btypeJ}, b.getBondOrder(), CanSwap::Yes);
                auto f   = fs->findParameterTypeConst(bondId, "Dm");
                auto gt  = f.index();
                auto gti = ForceConstants_[InteractionType::BONDS].reverseIndex(gt);
                a.set(gti, j, a.get(gti, j) + 1);
                a_copy.set(gti, j, a.get(gti, j));
                ntest[gti]++;
                if (ctest[gti].empty())
                {
                    ctest[gti].assign(bondId.id());
                }
            }
            else
            {
                gmx_fatal(FARGS, "No parameters for bond in the force field, atoms %s-%s mol %s",
                          atypeI, atypeJ,
                          mymol.getIupac().c_str());
            }
        }
        rhs.push_back(-mymol.Emol_);
    }

    char buf[STRLEN];
    snprintf(buf, sizeof(buf), "Inconsistency in number of energies nMol %d != #rhs %zu", nMol, rhs.size());
    GMX_RELEASE_ASSERT(static_cast<int>(rhs.size()) == nMol, buf);

    auto nzero = std::count_if(ntest.begin(), ntest.end(), [](const int n)
                               {
                                   return n == 0;
                               });

    GMX_RELEASE_ASSERT(nzero == 0, "Inconsistency in the number of bonds in poldata and ForceConstants_");

    std::vector<double> Edissoc(nD);
    a.solve(rhs, &Edissoc);
    if (debug)
    {
        dump_csv(ctest,  mymols(), ntest, Edissoc, a_copy, rhs.data());
    }
    for (size_t i = 0; i < ctest.size(); i++)
    {
        if (fplog)
        {
            fprintf(fplog, "Optimized dissociation energy for %8s with %4d copies to %g\n",
                    ctest[i].c_str(), ntest[i], Edissoc[i]);
        }
    }

    int i = 0;
    for (auto &b : ForceConstants_[InteractionType::BONDS].bondNames())
    {
        auto fs = poldata()->findForces(InteractionType::BONDS);
        for(auto &fp : *(fs->findParameters(b.first)))
        {
            if (fp.second.mutability() == Mutability::Free ||
                fp.second.mutability() == Mutability::Bounded)
            {
                if (fp.first == "De")
                {
                    fp.second.setValue(std::max(100.0, Edissoc[i++]));
                }
            }
        }
    }
}

void Optimization::InitOpt(FILE *fplog, bool bRandom)
{
    for (auto fs : poldata()->forcesConst())
    {
        if (optimize(fs.first))
        {
            ForceConstants fc(fs.second.fType(), fs.first, true);
            fc.analyzeIdef(mymols(), poldata());
            fc.makeReverseIndex();
            fc.dump(fplog);
            ForceConstants_.insert({fs.first, std::move(fc)});
        }
    }
    if (bDissoc_)
    {
        if (ForceConstants_[InteractionType::BONDS].nbad() <= mymols().size())
        {
            getDissociationEnergy(fplog);
        }
        else
        {
            printf("\n"
                   "WARNING: %zu molecule(s) is (are) not enough to calculate dissociation\n"
                   "         energy for %zu bond type(s) using linear regression. Default\n"
                   "         values from gentop.dat will be used as the initial guess.\n"
                   "         Recomendation is to add more molecules having the same bond types.\n\n",
                   mymols().size(), ForceConstants_[InteractionType::BONDS].nbad());
        }
    }
    
    polData2TuneFc(bRandom);
}

double Optimization::calcDeviation()
{
    if (!calcAll_)
    {
        if (PAR(commrec()))
        {
            broadcastPoldataUpdate();
        }
    }
    for (auto &pUpd : poldataUpdates_)
    {
        pUpd.execute(poldata());
    }
    poldataUpdates_.clear();
    if (calcAll_)
    {
        Bayes::printParameters(debug);
    }
    for (auto &mem : MolEnergyMap_)
    {
        mem.second.clear();
    }
    // First compute all the energies and store them
    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupport::Local) ||
            (calcAll_ && (mymol.eSupp_ == eSupport::Remote)))
        {
            int      nSP    = 0, nOpt = 0;
            int      natoms = mymol.mtop_->natoms;
            gmx_bool bpolar = (mymol.shellfc_ != nullptr);
            double   optHF;
            if (mymol.getOptHF(&optHF))
            {
                for (const auto &fc : ForceConstants_)
                {
                    if (fc.second.nbad() > 0)
                    {
                        mymol.UpdateIdef(poldata(), fc.first);
                    }
                }

                mymol.f_.resizeWithPadding(natoms);
                mymol.optf_.resizeWithPadding(natoms);
                int  molid          = mymol.getIndex();
                auto molEnergyEntry = MolEnergyMap_.find(molid);
                if (molEnergyEntry == MolEnergyMap_.end())
                {
                    // Create new (empty) MolEnergy and insert it into the map
                    MolEnergy me;
                    molEnergyEntry = 
                        MolEnergyMap_.insert(MolEnergyMap_.begin(),
                                             std::pair<int, MolEnergy>(molid, std::move(me)));
                }
                molEnergyEntry->second.clear();
                // Now loop over experiments!
                for (auto &ei : mymol.experimentConst())
                {
                    auto jtype = ei.getJobtype();
                    // Exclude other experimental data points!
                    if (jtype == JOB_OPT || jtype == JOB_SP)
                    {
                        double spHF;
                        if (!ei.getHF(&spHF))
                        {
                            continue;
                        }

                        FILE *dbcopy = debug;
                        debug  = nullptr;
                        mymol.changeCoordinate(ei, bpolar);
                        double rmsf;
                        mymol.computeForces(debug, commrec(), &rmsf);
                        debug  = dbcopy;

                        double deltaEref  = spHF - optHF;
                        if (deltaEref < 0)
                        {
                            gmx_fatal(FARGS, "Energy mismatch for %s. spHF = %g (%s) optHF = %g",
                                      mymol.getMolname().c_str(),
                                      spHF,
                                      ei.getDatafile().c_str(),
                                      optHF);
                        }
                        double deltaEalex = mymol.enerd_->term[F_EPOT] - mymol.Emol_;
                        molEnergyEntry->second.addE(jtype, deltaEref, deltaEalex);

                        // Force is added only for the opt geometry
                        double OptForce2 = 0.0;
                        if (jtype == JOB_OPT)
                        {
                            for (int j = 0; j < natoms; j++)
                            {
                                OptForce2 += iprod(mymol.f_[j], mymol.f_[j]);
                                copy_rvec(mymol.f_[j], mymol.optf_[j]);
                            }
                            OptForce2 /= natoms;
                            molEnergyEntry->second.setForce2(OptForce2);
                            molEnergyEntry->second.setTerms(mymol.enerd_->term);
                            nOpt += 1;
                        }
                        else
                        {
                            nSP += 1;
                        }
                        if (nullptr != debug)
                        {
                            int angleType = poldata()->findForces(InteractionType::ANGLES)->fType();
                            int pdihType  = poldata()->findForces(InteractionType::PROPER_DIHEDRALS)->fType();
                            int idihType  = poldata()->findForces(InteractionType::IMPROPER_DIHEDRALS)->fType();
                            int vdwType   = poldata()->findForces(InteractionType::VDW)->fType();
                            fprintf(debug, "spHF: %g  optHF: %g  deltaRef: %g  deltaEalex: %g\n",
                                    spHF, optHF, deltaEref, deltaEalex);
                            fprintf(debug, "%s Chi2 %g Morse %g  "
                                    "%s %g Langle %g %s %g %s %g Coul %g VdW %g POL %g  Force %g\n",
                                    mymol.getMolname().c_str(),
                                    gmx::square(deltaEalex-deltaEref),
                                    mymol.enerd_->term[F_MORSE],
                                    interaction_function[angleType].name,
                                    mymol.enerd_->term[angleType],
                                    mymol.enerd_->term[F_LINEAR_ANGLES],
                                    interaction_function[pdihType].name,
                                    mymol.enerd_->term[pdihType],
                                    interaction_function[idihType].name,
                                    mymol.enerd_->term[idihType],
                                    mymol.enerd_->term[F_COUL_SR],
                                    mymol.enerd_->term[vdwType],
                                    mymol.enerd_->term[F_POLARIZATION],
                                    std::sqrt(OptForce2));
                        }
                    }
                }
            }
            else
            {
                gmx_fatal(FARGS, "There is no optimized structure for %s\n",
                          mymol.getMolname().c_str());
            }
            if (debug)
            {
                fprintf(debug, "%d: %s nOpt %d nSP %d\n",
                        commrec()->nodeid,
                        mymol.getMolname().c_str(),
                        nOpt, nSP);
            }
        }
    }
    // Now compute the deviation for the fitting or otherwise
    resetChiSquared();
    double nCalc = 0;
    double ePot2 = 0;
    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupport::Local) ||
            (calcAll_ && (mymol.eSupp_ == eSupport::Remote)))
        {
            int  molid          = mymol.getIndex();
            auto molEnergyEntry = MolEnergyMap_.find(molid);
            if (molEnergyEntry != MolEnergyMap_.end())
            {
                increaseChiSquared(ermsForce2, 1, molEnergyEntry->second.force2());
                if (debug)
                {
                    fprintf(debug, "%d: %s molid: %d nEntries: %zu\n",
                            commrec()->nodeid,
                            mymol.getMolname().c_str(),
                            molid,
                            molEnergyEntry->second.entries().size());
                }
                for (const auto &d : molEnergyEntry->second.entries())
                {
                    ePot2 += gmx::square(d.EReference-d.EAlexandria);
                    nCalc += 1;
                }
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "%d: ePot2 = %g nCalc = %g\n",
                commrec()->nodeid, ePot2, nCalc);
    }
    if (!calcAll_ && PAR(commrec()))
    {
        gmx_sumd(1, &ePot2, commrec());
        gmx_sumd(1, &nCalc, commrec());
    }
    double chi2 = ePot2/nCalc;
    if (debug)
    {
        fprintf(debug, "%d: ePot2 = %g nCalc = %g chi2 = %g\n",
                commrec()->nodeid, ePot2, nCalc, chi2);
    }
    setChiSquared(ermsEPOT, 1, chi2);
    setChiSquared(ermsTOT, 1, chi2);
    printChiSquared(debug);
    
    return chiSquared(ermsTOT);
}

bool Optimization::optRun(FILE                   *fplog,
                          const gmx_output_env_t *oenv,
                          const char             *xvgconv,
                          const char             *xvgepot)
{
    bool bMinimum = false;
    if (MASTER(commrec()))
    {
        if (PAR(commrec()))
        {
            // Tell the slave nodes how many times they have
            // to run calcDeviation.
            int niter = 1+Bayes::numberObjectiveFunctionCalls();
            for (int dest = 1; dest < commrec()->nnodes; dest++)
            {
                gmx_send_int(commrec(), dest, niter);
            }
        }
        double chi2_min = calcDeviation();
        if (fplog)
        {
            fprintf(fplog, "Initial chi2 %g\n", chi2_min);
        }
        {
            std::vector<std::string> paramClass;
            Bayes::setOutputFiles(xvgconv, paramClass, xvgepot, oenv);
        }
        double chi2 = Bayes::MCMC(fplog);
        if (fplog)
        {
            fprintf(fplog, "Final chi2 %g\n", chi2_min);
        }
        if (chi2 < chi2_min)
        {
            chi2_min = chi2;
            bMinimum = true;

            setCalcAll(true);
            if (fplog)
            {
                auto pmean  = Bayes::getPmean();
                auto psigma = Bayes::getPsigma();
                auto best   = Bayes::getBestParam();
                // This call copies data to poldata as well.
                if (!best.empty())
                {
                    std::vector<bool> changed;
                    changed.resize(best.size(), true);
                    toPoldata(changed);
                    double chi2 = calcDeviation();
                    fprintf(fplog, "\nLowest RMSD value during optimization: %g.\n",
                            std::sqrt(chi2));
                    fprintf(fplog, "Parameters after the optimization:\n");
                    fprintf(fplog, "%-5s  %10s  %10s  %10s\n", "Index",
                            "Average", "Std. Dev.", "Optimum");
                    for (size_t k = 0; k < Bayes::nParam(); k++)
                    {
                        fprintf(fplog, "%5zu  %10g  %10g  %10g\n",
                                k, pmean[k], psigma[k], best[k]);
                    }
                }
            }
            printChiSquared(fplog);
        }
    }
    else
    {
        /* S L A V E   N O D E S */
        int niter = gmx_recv_int(commrec(), 0);
        for (int n = 0; n < niter; n++)
        {
            calcAll_ = false;
            (void) calcDeviation();
        }
    }
    return bMinimum;
}

void Optimization::printMolecules(FILE *fp,
                                  bool  bForce,
                                  bool  bMtop)
{
    int j, k;
    for (auto &mi : mymols())
    {
        int nSP = 0, nOpt = 0;
        for (auto &ei : mi.experimentConst())
        {
            auto jtype = ei.getJobtype();
            if (jtype == JOB_SP)
            {
                nSP += 1;
            }
            else if (jtype == JOB_OPT)
            {
                nOpt += 1;
            }
        }
        if (nOpt != 1)
        {
            fprintf(stderr, "Number of optimized conformations is %d. Check your input.\n", nOpt);
        }
        fprintf(fp, "%s natoms: %d Opt conformations: %d SP conformations: %d\n",
                mi.getMolname().c_str(),
                mi.mtop_->natoms,
                nOpt,
                nSP);
        auto myatoms = mi.atoms();
        for (j = 0; j < myatoms->nr; j++)
        {
            fprintf(fp, "  %-5s  %-5s  q = %10g",
                    *(myatoms->atomname[j]),
                    *(myatoms->atomtype[j]),
                    myatoms->atom[j].q);
            if (bForce)
            {
                fprintf(fp, "   f = %8.3f  %8.3f  %8.3f",
                        mi.optf_[j][XX], mi.optf_[j][YY], mi.optf_[j][ZZ]);
            }
            fprintf(fp, "\n");
        }

        if (bMtop)
        {
            pr_mtop(fp, 0, mi.getMolname().c_str(), mi.mtop_, true, false);
        }
    }
    if (bForce)
    {
        for (const auto &mi : mymols())
        {
            int  molid     = mi.getIndex();
            auto molEnergy = MolEnergyMap_.find(molid);
            if (molEnergy != MolEnergyMap_.end())
            {
                fprintf(fp, "%-20s", mi.getMolname().c_str());
                for (k = 0; k < F_NRE; k++)
                {
                    real term = molEnergy->second.term(k);
                    if (term != 0 ||
                        (mi.mtop_->moltype[0].ilist[k].size() > 0))
                    {
                        fprintf(fp, " %s: %.2f", interaction_function[k].name,
                                mi.enerd_->term[k]);
                    }
                }
                fprintf(fp, "\n");
            }
        }
    }
}

void Optimization::printResults(FILE                   *fp,
                                const char             *title,
                                const char             *hform_xvg,
                                const char             *HF_xvg,
                                const gmx_output_env_t *oenv)
{
    FILE       *xfp = nullptr, *hfp = nullptr;
    gmx_stats_t gst;

    gst = gmx_stats_init();
    {
        const char *title = "Enthalpy of Formation";
        const char *yaxis = "Alexandria (kJ/mol)";
        if (nullptr != hform_xvg)
        {
            xfp = xvgropen(hform_xvg, title, "Experiment (kJ/mol)", yaxis, oenv);
        }
        if (nullptr != HF_xvg)
        {
            hfp = xvgropen(HF_xvg, title,
                           "Experiment + \\f{12}DD\\f{4}E(B3LYP/aug-cc-pVTZ) (kJ/mol)",
                           yaxis, oenv);
            xvgr_view(hfp, 0.15, 0.15, 0.75, 0.85, oenv);
        }
    }
    fprintf(fp, "%s\n", title);
    fprintf(fp, "Fit of energy at different conformations to y = ax\n");
    fprintf(fp, "Nr.   %-30s %10s %10s %7s %7s %7s %7s %7s %4s\n",
            "Molecule", "DHf@298K", "Emol@0K", "rmsF@0K",
            "rms E", "MSE E", "a", "r2", "N");

    int    imol          = 0;
    int    nconformation = 0;
    for (auto mi = mymols().begin(); mi < mymols().end(); mi++, imol++)
    {
        int  molid      = mi->getIndex();
        auto molEnergy  = MolEnergyMap_.find(molid);
        if (molEnergy != MolEnergyMap_.end())
        {
            if (nullptr != hfp)
            {
                fprintf(hfp, "@ s%d legend \"%s\"\n", imol,
                        mi->getMolname().c_str());
                fprintf(hfp, "@type xy\n");
            }
            gmx_stats_t gmol = gmx_stats_init();
            for (const auto &entry : molEnergy->second.entries())
            {
                real deltaE = entry.EAlexandria - entry.EReference;
                if (nullptr != xfp && entry.jobType == JOB_OPT)
                {
                    fprintf(xfp, "%10g  %10g\n", mi->Hform_, mi->Hform_ + deltaE);
                }
                real Hexper = mi->Hform_ + entry.EReference;
                real Halex  = mi->Hform_ + entry.EAlexandria;
                if (nullptr != hfp)
                {
                    fprintf(hfp, "%10g  %10g\n", Hexper, Halex);
                }
                gmx_stats_add_point(gmol, Hexper, Halex, 0, 0);
                if (entry.jobType == JOB_OPT)
                {
                    gmx_stats_add_point(gst, Hexper, Halex, 0, 0);
                }
            }
            if (nullptr != hfp)
            {
                fprintf(hfp, "&\n");
            }

            real a, chi2, da, rmsd, Rfit, Rdata, mse, mae;
            int  N;
            // Note we are ignoring the return value for these functions
            (void) gmx_stats_get_a(gmol, 0, &a, &da, &chi2, &Rfit);
            (void) gmx_stats_get_rmsd(gmol, &rmsd);
            (void) gmx_stats_get_mse_mae(gmol, &mse, &mae);
            (void) gmx_stats_get_npoints(gmol, &N);
            nconformation += N;
            gmx_stats_get_corr_coeff(gmol, &Rdata);
            fprintf(fp, "%-5d %-30s %10g %10g %7.1f %7.3f %7.3f %7.3f %6.1f%% %4d\n",
                    imol,
                    mi->getMolname().c_str(),
                    mi->Hform_,
                    mi->Emol_,
                    std::sqrt(molEnergy->second.force2()), rmsd,
                    mse, a, Rdata*100, N);
            gmx_stats_free(gmol);
        }
    }
    fprintf(fp, "RMSD from target energies for %zu compounds and %d conformation is %g.\n",
            mymols().size(), nconformation, std::sqrt(chiSquared(ermsTOT)));
    fprintf(fp, "\n");
    if (nullptr != hform_xvg)
    {
        xvgrclose(xfp);
        do_view(oenv, hform_xvg, nullptr);
    }
    if (nullptr != HF_xvg)
    {
        xvgrclose(hfp);
        do_view(oenv, HF_xvg, nullptr);
    }

    real a, da, chi2, rmsd, Rfit, Rdata;
    int  N;
    gmx_stats_get_a(gst, 1, &a, &da, &chi2, &Rfit);
    gmx_stats_get_npoints(gst, &N);
    gmx_stats_get_corr_coeff(gst, &Rdata);
    gmx_stats_get_rmsd(gst, &rmsd);
    fprintf(fp, "Regression analysis fit to y = ax of %d minimum energy structures:\n", N);
    fprintf(fp, "a = %.3f  r2 = %.1f%%  rmsd = %.2f kJ/mol\n",
            a, Rdata*100, rmsd);
    gmx_stats_free(gst);
    fflush(fp);
}

} // namespace alexandria

int alex_tune_fc(int argc, char *argv[])
{
    const char           *desc[] = {
        "tune_fc read a series of molecules and corresponding experimental",
        "heats of formation from a file, and tunes parameters in an algorithm",
        "until the experimental energies are reproduced by the force field.[PAR]",
        "Minima and maxima for the parameters can be set, these are however",
        "not strictly enforced, but rather they are penalized with a harmonic",
        "function, for which the force constant can be set explicitly.[PAR]",
        "At every reinit step parameters are changed by a random amount within",
        "the fraction set by step size, and within the boundaries given",
        "by the minima and maxima. If the [TT]-random[tt] flag is",
        "given a completely random set of parameters is generated at the start",
        "of each run. At reinit steps however, the parameters are only changed",
        "slightly, in order to speed-up local search but not global search."
        "In other words, complete random starts are done only at the beginning of each",
        "run, and only when explicitly requested.[PAR]",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored."
    };

    t_filenm              fnm[] = {
        { efDAT, "-f",     "allmols",    ffREAD  },
        { efDAT, "-d",     "gentop",     ffOPTRD },
        { efDAT, "-o",     "tune_fc",    ffWRITE },
        { efDAT, "-sel",   "molselect",  ffREAD  },
        { efXVG, "-table", "table",      ffOPTRD },
        { efLOG, "-g",     "tune_fc",    ffWRITE },
        { efXVG, "-x",     "hform-corr", ffWRITE },
        { efXVG, "-hf",    "hf-corr",    ffWRITE },
        { efXVG, "-conv",  "param-conv", ffWRITE },
        { efXVG, "-epot",  "param-epot", ffWRITE }
    };

    const int             NFILE         = asize(fnm);

    int                   reinit        = 0;
    int                   compress      = 0;
    int                   nmultisim     = 0;
    bool                  bRandom       = false;
    gmx_bool              bZPE          = false;
    gmx_bool              bZero         = true;
    gmx_bool              bTestPar      = false;
    gmx_bool              bForceOutput  = true;
    
    static const char          *select_types[]   = {nullptr, "Train", "Test", "Ignore", "Unknown", nullptr};
    
    t_pargs               pa[]          = {
        { "-multi",   FALSE, etINT, {&nmultisim},
          "Do optimization in multiple simulation" },
        { "-reinit",  FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A value of 0 means this is never done at all." },
        { "-zpe",     FALSE, etBOOL, {&bZPE},
          "Consider zero-point energy from thermochemistry calculations in order to calculate the reference enthalpy of the molecule. If set, the zero point energy will be subtracted from the target energy when optimizing the force field model." },
        { "-testpar", FALSE, etBOOL, {&bTestPar},
          "Test the parallel execution gives the same result." },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML file" },
        { "-random", FALSE, etBOOL, {&bRandom},
          "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step and before each subsequent run." },
        { "-force_output", FALSE, etBOOL, {&bForceOutput},
          "Write output even if no new minimum is found" },
        { "-select", FALSE, etENUM, {select_types},
          "Select type for making the dataset for training or testing." }
    };

    FILE                 *fplog;
    gmx_output_env_t     *oenv;
    MolSelect             gms;

    std::vector<t_pargs>  pargs;
    for (size_t i = 0; i < asize(pa); i++)
    {
        pargs.push_back(pa[i]);
    }
    alexandria::Optimization opt;
    opt.add_pargs(&pargs);

    if (!parse_common_args(&argc, 
                           argv, 
                           PCA_CAN_VIEW, 
                           NFILE, 
                           fnm,
                           pargs.size(), 
                           pargs.data(),
                           asize(desc), 
                           desc, 0, 
                           nullptr, 
                           &oenv))
    {
        return 0;
    }
    opt.optionsFinished(opt2fn("-o", NFILE, fnm));

    if (MASTER(opt.commrec()))
    {
        fplog = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");
        print_header(fplog, pargs);
    }
    else
    {
        fplog = nullptr;
    }

    if (MASTER(opt.commrec()))
    {
        gms.read(opt2fn_null("-sel", NFILE, fnm));
    }

    const char *tabfn = opt2fn_null("-table", NFILE, fnm);

    iMolSelect select_type = name2molselect(select_types[0]);
    opt.Read(fplog ? fplog : (debug ? debug : nullptr),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             bZero,
             gms,
             bZPE,
             true,
             tabfn,
             select_type);

    opt.checkSupport(fplog);

    if (MASTER(opt.commrec()))
    {
        opt.InitOpt(fplog, bRandom);
        opt.printMolecules(fplog, false, false);
    }

    opt.broadcast();

    if (bTestPar)
    {
        opt.setCalcAll(false);
        auto chi2 = opt.calcDeviation();
        if (MASTER(opt.commrec()))
        {
            fprintf(fplog, "chi2 = %g\n", chi2);
            opt.printResults(fplog, (char *)"Before optimization - parallel test",
                             nullptr, nullptr, oenv);
        }
    }
    if (MASTER(opt.commrec()))
    {
        opt.setCalcAll(true);
        auto chi2 = opt.calcDeviation();
        fprintf(fplog, "chi2 = %g\n", chi2);
        opt.printResults(fplog, (char *)"Before optimization",
                         nullptr, nullptr, oenv);
        opt.setCalcAll(false);
    }
    if (bTestPar)
    {
        gmx_ffclose(fplog);
        return 0;
    }
    // Optimize the parameters to minimize the penalty function.
    bool bMinimum = opt.optRun(fplog,
                               oenv,
                               opt2fn("-conv", NFILE, fnm),
                               opt2fn("-epot", NFILE, fnm));

    if (MASTER(opt.commrec()))
    {
        if (bMinimum || bForceOutput)
        {
            // Now print the output.
            opt.printMolecules(fplog, true, false);
            opt.printResults(fplog, (char *)"After optimization",
                             opt2fn("-x", NFILE, fnm),
                             opt2fn("-hf", NFILE, fnm),
                             oenv);
            //writePoldata(opt2fn("-o", NFILE, fnm), opt.poldata(), compress);
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
        gmx_ffclose(fplog);
    }

    return 0;
}
