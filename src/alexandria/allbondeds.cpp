/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "allbondeds.h"

#include "gromacs/fileio/xvgr.h"

#include "tuning_utility.h"
#include "act/utility/units.h"
    
static void round_numbers(real *av, real *sig, int power10)
{
    *av  = ((int)(*av*power10))/(1.0*power10);
    *sig = ((int)(*sig*1.5*power10))/(1.0*power10);
}

namespace alexandria
{

void OneBonded::addPoint(double x)
{
    auto N  = lsq_.get_npoints();
    auto ok = lsq_.add_point(N, x, 0, 0);
    
    if (eStats::OK != ok)
    {
        fprintf(stderr, "Problem adding a point %s\n", gmx_stats_message(ok));
    }
}

void OneBonded::writeHistogram(const char             *fn_prefix,
                               const char             *xaxis,
                               const gmx_output_env_t *oenv,
                               double                  binwidth)
{
    int     nbins      = 0;
    int     normalized = 0;
    std::vector<double> x, y;
    eStats  estats     = lsq_.make_histogram(binwidth, &nbins,
                                             eHisto::Y, normalized, &x, &y);
    if (eStats::OK != estats)
    {
        fprintf(stderr, "Could not make a histogram for %s because of %s\n",
                id_.id().c_str(), gmx_stats_message(estats));
        return;
    }
    auto Nsample = lsq_.get_npoints();
    
    std::string  title = gmx::formatString("%s N = %d", id_.id().c_str(),
                                           static_cast<int>(Nsample));
    std::string  fn    = gmx::formatString("%s_%s.xvg", fn_prefix, id_.id().c_str());
    FILE        *fp    = xvgropen(fn.c_str(), title.c_str(),
                                  xaxis, "N", oenv);
    for (int j = 0; (j < nbins); j++)
    {
        fprintf(fp, "%g  %g\n", x[j], y[j]);
    }
    xvgrclose(fp);
}
    
eStats OneBonded::getAverageSigmaN(real   *average,
                                   real   *sigma,
                                   size_t *N)
{
    eStats ok = lsq_.get_average(average);
    if (eStats::OK == ok)
    {
        ok = lsq_.get_sigma(sigma);
    }
    if (eStats::OK == ok)
    {
        *N = lsq_.get_npoints();
    }
    return ok;
}

void AllBondeds::addOptions(std::vector<t_pargs> *pargs)
{
    t_pargs mypargs[] = 
        {
            { "-De",    FALSE, etREAL, {&De_},
              "Dissociation energy (kJ/mol). Ignore if the [TT]-dissoc[tt] option is used." },
            { "-beta",    FALSE, etREAL, {&beta_},
              "Steepness of the Morse potential (1/nm)" },
            { "-kt",    FALSE, etREAL, {&kt_},
              "Angle force constant (kJ/mol/rad^2)" },
            { "-klin",  FALSE, etREAL, {&klin_},
              "Linear angle force constant (kJ/mol/nm^2)" },
            { "-kp",    FALSE, etREAL, {&kp_},
              "Dihedral angle force constant (kJ/mol/rad^2)" },
            { "-kimp",    FALSE, etREAL, {&kimp_},
              "Improper dihedral angle force constant (kJ/mol/rad^2)" },
            { "-kub",   FALSE, etREAL, {&kub_},
              "Urey_Bradley force constant" },
            { "-bond_tol",   FALSE, etREAL, {&bond_tol_},
              "Tolerance for warning about large sigma in bond-lengths (pm)" },
            { "-angle_tol",   FALSE, etREAL, {&angle_tol_},
              "Tolerance for harmonic and linear angles" },
            { "-factor", FALSE, etREAL, {&factor_},
              "Scale factor to set minimum and maximum values of parameters. Min will be set to value*factor and max to value/factor, assuming factor < 1." },
            { "-bspacing", FALSE, etREAL, {&bspacing_},
              "Spacing for bond histograms in pm" },
            { "-aspacing", FALSE, etREAL, {&aspacing_},
              "Spacing for angle histograms in degrees" },
            { "-laspacing", FALSE, etREAL, {&laspacing_},
              "Spacing for linear angle histograms (no unit)" },
            { "-dspacing", FALSE, etREAL, {&dspacing_},
              "Spacing for dihedral angle histograms in degrees" },
        };
    doAddOptions(pargs, sizeof(mypargs)/sizeof(mypargs[0]), mypargs);
}

void AllBondeds::addBonded(FILE                           *fplog, 
                           InteractionType                 iType,
                           const MyMol                    &mmi,
                           const Identifier               &bondId,
                           const std::vector<int>         &atomid)
{
    auto x = mmi.x();
    // We need to check for linear angles here before we
    // add this to the tables.
    double refValue = 0;
    switch(iType)
    {
    case InteractionType::ANGLES:
    case InteractionType::LINEAR_ANGLES:
        {
            rvec rij, rkj;
            rvec_sub(x[atomid[0]], x[atomid[1]], rij);
            rvec_sub(x[atomid[2]], x[atomid[1]], rkj);
            refValue = RAD2DEG*gmx_angle(rij, rkj);
            if ( (refValue > 175) || (refValue < 5) )
            {
                iType = InteractionType::LINEAR_ANGLES;
            }
        }
        break;
    case InteractionType::BONDS:
        {
            rvec rij;
            rvec_sub(x[atomid[0]], x[atomid[1]], rij);
            refValue = norm(rij)*1000;
        }
        break;
    case InteractionType::PROPER_DIHEDRALS:
        {
            rvec  r_ij, r_kj, r_kl, mm, nn;
            int   t1, t2, t3;
            t_pbc pbc;
            matrix box = {{ 0 }};
            set_pbc(&pbc, epbcNONE, box);
            refValue = RAD2DEG*dih_angle(x[atomid[0]], 
                                         x[atomid[1]],
                                         x[atomid[2]],
                                         x[atomid[3]],
                                         &pbc, r_ij, r_kj, r_kl, mm, nn,
                                         &t1, &t2, &t3);
            if (refValue < 0)
            {
                refValue += 360;
            }
        }
        break;
    case InteractionType::IMPROPER_DIHEDRALS:
        {
            rvec  r_ij, r_kj, r_kl, mm, nn;
            int   t1, t2, t3;
            t_pbc pbc;
            matrix box = {{ 0 }};
            set_pbc(&pbc, epbcNONE, box);
            refValue = RAD2DEG*dih_angle(x[atomid[0]], 
                                         x[atomid[1]],
                                         x[atomid[2]],
                                         x[atomid[3]],
                                         &pbc, r_ij, r_kj, r_kl, mm, nn,
                                         &t1, &t2, &t3);
            if (refValue < 0)
            {
                refValue += 360;
            }
            if (InteractionType::IMPROPER_DIHEDRALS == iType)
            {
                while (refValue > 170)
                {
                    refValue -= 180;
                }
            }
        }
        break;
    case InteractionType::COULOMB:
    case InteractionType::VDW:
        break;
    default:
        {
            gmx_fatal(FARGS, "Help!");
        }
    }
    // Look up the interaction type
    if (bondeds_.find(iType) == bondeds_.end())
    {
        std::vector<OneBonded> ob;
        bondeds_.insert(std::pair<InteractionType, std::vector<OneBonded> >(iType, std::move(ob)));
    }
    // Look up the bondId in the data structure
    auto ob = std::find_if(bondeds_[iType].begin(), bondeds_[iType].end(),
                           [bondId](const OneBonded &b){ return bondId == b.id(); });
    if (bondeds_[iType].end() == ob)
    {
        bondeds_[iType].push_back(OneBonded(bondId));
        ob = bondeds_[iType].end()-1;
    }
    ob->addPoint(refValue);
    
    if (nullptr != fplog)
    {
        fprintf(fplog, "%s %s-%s %g\n",
                mmi.getMolname().c_str(), interactionTypeToString(iType).c_str(),
                bondId.id().c_str(), refValue);
    }
}

void AllBondeds::writeHistogram(const gmx_output_env_t *oenv)
{
    std::map<InteractionType, const char *> xaxis = {
        { InteractionType::BONDS,              "Distance (pm)" },
        { InteractionType::ANGLES,             "Angle (deg.)" },
        { InteractionType::LINEAR_ANGLES,      "Linear Angle (deg.)" },
        { InteractionType::PROPER_DIHEDRALS,   "Dihedral angle (deg.)" },
        { InteractionType::IMPROPER_DIHEDRALS, "Improper angle (deg.)" }
    };
    std::map<InteractionType, real> spacing = {
        { InteractionType::BONDS,              bspacing_ },
        { InteractionType::ANGLES,             aspacing_ },
        { InteractionType::LINEAR_ANGLES,      laspacing_ },
        { InteractionType::PROPER_DIHEDRALS,   dspacing_ },
        { InteractionType::IMPROPER_DIHEDRALS, dspacing_ }
    };
    for (auto &mm : bondeds_)
    {
        for (auto &nn : mm.second)
        {
            nn.writeHistogram(interactionTypeToString(mm.first).c_str(),
                              xaxis[mm.first],
                              oenv, spacing[mm.first]);
        }
    }
}
    
static real calc_r13(const Poldata    *pd,
                     const Identifier &bondId,
                     const real        angle)
{
    double r12 = 0, r23 = 0, r13 = 0;
    auto atoms      = bondId.atoms();
    auto bondOrders = bondId.bondOrders();
    Identifier aij ({atoms[0], atoms[1] }, { bondOrders[0] }, CanSwap::Yes);
    Identifier akj ({atoms[2], atoms[1] }, { bondOrders[1] }, CanSwap::Yes);

    std::string type("bondlength");
    auto fs  = pd->findForcesConst(InteractionType::BONDS);
    if (fs.parameterExists(aij) && fs.parameterExists(akj))
    {
        auto bij = fs.findParameterTypeConst(aij, type);
        auto bkj = fs.findParameterTypeConst(akj, type);

        r12 = convertToGromacs(bij.value(), bij.unit());
        r23 = convertToGromacs(bkj.value(), bkj.unit());

        r13 = std::sqrt((r12*r12) + (r23*r23) - (2*r12*r23*std::cos(DEG2RAD*angle)));

        return r13;
    }
    else
    {
        fprintf(stderr, "Cannot find one of %s or %s\n", aij.id().c_str(), akj.id().c_str());
    }
    return 0.0;
}


/*! \brief Compute geometry parameters for linear angles
 * \param[in]  pd       Force field
 * \param[in]  bondId   Identifier corresponding to angle, i.e. three atoms
 * \param[out] alin     The constant determining the center of the bond
 * \param[out] sigmalin The uncertainty in alin.
 */
static void calc_linear_angle_a(const Poldata    *pd,
                                const Identifier &bondId,
                                double           *alin,
                                double           *sigmalin)
{
    auto atoms      = bondId.atoms();
    auto bondOrders = bondId.bondOrders();
    auto bij = Identifier({atoms[0], atoms[1]}, {bondOrders[0]}, CanSwap::Yes);
    auto bjk = Identifier({atoms[1], atoms[2]}, {bondOrders[1]}, CanSwap::Yes);
    auto fs  = pd->findForcesConst(InteractionType::BONDS);
    if (!fs.parameterExists(bij) || !fs.parameterExists(bjk))
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find bond %s or %s in force field", bij.id().c_str(), bjk.id().c_str()).c_str()));
    }
    std::string blen("bondlength");
    auto pij = fs.findParameterTypeConst(bij, blen);
    auto pjk = fs.findParameterTypeConst(bjk, blen);
    *alin = pjk.value()/(pij.value()+pjk.value());
    *sigmalin = std::sqrt(gmx::square(pij.uncertainty())+
                          gmx::square(pjk.uncertainty()));
}

void AllBondeds::updatePoldata(FILE    *fp,
                               Poldata *pd)
{
    auto bType = InteractionType::BONDS;
    auto fs    = pd->findForces(bType);
    fs->eraseParameter();
    auto fType = fs->fType();
    for(auto &i : bondeds_[bType])
    {
        size_t N = 0;
        real   av = 0, sig = 0;
        i.getAverageSigmaN(&av, &sig, &N);
        auto bondId = i.id();
        // Note that the order of parameters is important!
        // Rounding the numbers to 1/10 pm and 1/10 degree
        round_numbers(&av, &sig, 10);
        switch (fType)
        {
        case F_MORSE:
            {
                fs->addParameter(bondId, "bondlength",
                                 ForceFieldParameter("pm", av, sig, N, av*factor_, av/factor_, Mutability::Bounded, false, true));
                fs->addParameter(bondId, "De",
                                 ForceFieldParameter("kJ/mol", De_, 0, 1, De_*factor_, De_/factor_, Mutability::Bounded, false, true));
                fs->addParameter(bondId, "beta",
                                 ForceFieldParameter("1/nm", beta_, 0, 1, beta_*factor_, beta_/factor_, Mutability::Bounded, false, true));
                real D0_ = 0;
                fs->addParameter(bondId, "D0",
                                 ForceFieldParameter("kJ/mol", D0_, 0, 1, -400, 400, Mutability::Bounded, false, true));
            }
            break;
        case F_BONDS:
            {
                fs->addParameter(bondId, "bondlength",
                                 ForceFieldParameter("pm", av, sig, N, av*factor_, av/factor_, Mutability::Bounded, false, true));
                fs->addParameter(bondId, "kb",
                                 ForceFieldParameter("kJ/mol/nm2", kb_, 0, 1, kb_*factor_, kb_/factor_, Mutability::Bounded, false, true));
            }
            break;
        default:
            gmx_fatal(FARGS, "Don't know what to do for ftype %s", interaction_function[fType].name);
        }
                    
        fprintf(fp, "bond-%s len %g sigma %g (pm) N = %d%s\n",
                i.id().id().c_str(), av, sig, static_cast<int>(N), (sig > bond_tol_) ? " WARNING" : "");
    }

    for(auto &bb : bondeds_)
    {
        auto iType = bb.first;
        auto fs    = pd->findForces(iType);
        fType      = fs->fType();
        if (iType != bType &&
            iType != InteractionType::VDW &&
            iType != InteractionType::COULOMB)
        {
            fs->eraseParameter();
        }
        for (auto &i : bb.second)
        {
            size_t N = 0;
            real   av = 0, sig = 0;
            i.getAverageSigmaN(&av, &sig, &N);
            auto bondId = i.id();
            switch (iType)
            {
            case InteractionType::BONDS:
                // Done earlier!
                break;
            case InteractionType::ANGLES:
                {
                    round_numbers(&av, &sig, 10);
                    fs->addParameter(bondId, "angle",
                                     ForceFieldParameter("degree", av, sig, N, av*factor_, 
                                                         std::min(180.0, av/factor_), Mutability::Bounded, false, true));
                    fs->addParameter(bondId, "kt",
                                     ForceFieldParameter("kJ/mol/rad2", kt_, 0, 1, kt_*factor_, kt_/factor_, Mutability::Bounded, false, false));
                    fprintf(fp, "angle-%s angle %g sigma %g (deg) N = %d%s\n",
                            bondId.id().c_str(), av, sig, static_cast<int>(N), (sig > angle_tol_) ? " WARNING" : "");
                    if (fType == F_UREY_BRADLEY)
                    {
                        std::string unit("pm");
                        double r13 = convertFromGromacs(calc_r13(pd, bondId, av), unit);
                        fs->addParameter(bondId, "r13", 
                                         ForceFieldParameter(unit, r13, 0, N, r13*factor_, r13/factor_,
                                                             Mutability::Bounded, false, true));
                        fs->addParameter(bondId, "kub", 
                                         ForceFieldParameter("kJ/mol/nm2", kub_, 0, 1, kub_*factor_, kub_/factor_, Mutability::Bounded, false, false));
                    }
                }
                break;
            case InteractionType::LINEAR_ANGLES:
                {
                    double alin, sigmalin, myfactor = 0.99;
                    
                    calc_linear_angle_a(pd, bondId, &alin, &sigmalin);
                    fs->addParameter(bondId, "a",
                                     ForceFieldParameter("", alin, sigmalin, N, alin*myfactor, alin/myfactor, Mutability::Bounded, false, true));
                    fs->addParameter(bondId, "klin", 
                                     ForceFieldParameter("kJ/mol/nm2", klin_, 0, 1, klin_*factor_, klin_/factor_, Mutability::Bounded, false, true));
                    
                    fprintf(fp, "linear_angle-%s angle %g sigma %g N = %d%s\n",
                            bondId.id().c_str(), av, sig, static_cast<int>(N), (sig > angle_tol_) ? " WARNING" : "");
                
                    std::string unit("pm");
                    double r13 = convertFromGromacs(calc_r13(pd, bondId, av), unit);
                    fs->addParameter(bondId, "r13lin", 
                                     ForceFieldParameter(unit, r13, 0, N, r13*factor_, r13/factor_,
                                                         Mutability::Bounded, false, true));
                    fs->addParameter(bondId, "kublin", 
                                     ForceFieldParameter("kJ/mol/nm2", kub_, 0, 1, kub_*factor_, kub_/factor_, Mutability::Bounded, false, false));
                }
                break;
            case InteractionType::PROPER_DIHEDRALS:
                {
                    switch (fType)
                    {
                    case F_FOURDIHS:
                        {
                            double val = 1;
                            std::vector<std::string> cname = { "c0", "c1", "c2", "c3" };
                            for(auto &c : cname)
                            {
                                fs->addParameter(bondId, c,
                                                 ForceFieldParameter("kJ/mol", val, 0, 1, val*factor_, val/factor_, Mutability::Bounded, false, false));
                            }
                        }
                        break;
                    case F_PDIHS:
                        {
                            round_numbers(&av, &sig, 10);
                            fs->addParameter(bondId, "angle", 
                                             ForceFieldParameter("degree", av, sig, N, av*factor_, av/factor_, Mutability::Bounded, false, true));
                            fs->addParameter(bondId, "kp", 
                                             ForceFieldParameter("kJ/mol", kp_, 0, 1, kp_*factor_, kp_/factor_, Mutability::Bounded, false, true));
                            fs->addParameter(bondId, "mult", 
                                             ForceFieldParameter("", 3, 0, 1, 3, 3, Mutability::Fixed, true, true));
                        }
                        break;
                    default:
                        GMX_THROW(gmx::InternalError(gmx::formatString("Unsupported dihedral type %s",
                                                                       interaction_function[fType].name).c_str()));
                    }
                    fprintf(fp, "dihedral-%s angle %g sigma %g (deg)\n",
                            bondId.id().c_str(), av, sig);
                }
                break;
            case InteractionType::IMPROPER_DIHEDRALS:
                {
                    round_numbers(&av, &sig, 10);
                    // TODO make this a parameter
                    if (fabs(av) > 4)
                    {
                        fprintf(stderr, "Warning: large improper dihedral %g for %s\n",
                                av, bondId.id().c_str());
                    }
                    fs->addParameter(bondId, "phi", 
                                     ForceFieldParameter("degree", 0, 0, N, 0, 0, Mutability::Fixed, false, false));
                    fs->addParameter(bondId, "kimp", 
                                     ForceFieldParameter("kJ/mol", kimp_, 0, 1, kimp_*factor_, kimp_/factor_, Mutability::Bounded, false, true));
                    
                    fprintf(fp, "improper-%s angle %g sigma %g (deg)\n",
                            bondId.id().c_str(), av, sig);
                }
                break;
            default:
                break;
            }
        }
    }
}

void AllBondeds::extractGeometries(FILE                       *fp,
                                   const std::vector<MolProp> &mp,
                                   std::vector<MyMol>         *mymols,
                                   const Poldata              &pd,
                                   const MolSelect            &gms)
                                   
{
    for (auto mpi = mp.begin(); mpi < mp.end(); mpi++)
    {
        iMolSelect imol;
        if (gms.status(mpi->getIupac(), &imol) && 
            (imol == iMolSelect::Train || imol == iMolSelect::Test))
        {
            alexandria::MyMol mmi;
            int               i;
            mmi.Merge(&(*mpi));
            if (mmi.getMolname().size() == 0)
            {
                fprintf(fp, "Empty molname for molecule with formula %s\n",
                        mmi.formula().c_str());
                continue;
            }
            auto imm = mmi.GenerateTopology(fp, &pd,
                                            missingParameters::Generate);
            if (immStatus::OK != imm)
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "Could not make topology for %s, reason %s\n",
                            mmi.getMolname().c_str(),
                            immsg(imm) );
                }
                continue;
            }
            
            auto myatoms = mmi.atomsConst();
#define ATP(ii) (*myatoms.atomtype[ii])
            for (i = 0; i < myatoms.nr; i++)
            {
                std::string btpi;
                if (!pd.atypeToBtype(*myatoms.atomtype[i], &btpi))
                {
                    if (nullptr != debug)
                    {
                        fprintf(debug, "No bond-type support for atom %s in %s\n",
                                *myatoms.atomtype[i], mmi.getMolname().c_str());
                    }
                    break;
                }
            }
            if ((myatoms.nr <= 0) || (i < myatoms.nr))
            {
                if (nullptr != debug)
                {
                    fprintf(debug, "You may need to check the number of atoms for %s\n",
                            mmi.getMolname().c_str());
                }
                continue;
            }
            auto top = mmi.topology();
            for(const auto &entry : top->entries())
            {
                for (auto topentry : entry.second)
                {
                    addBonded(fp, entry.first, mmi, topentry->ids()[0], topentry->atomIndices());
                }
            }
            mymols->push_back(mmi);
        }
    }
}
    
void AllBondeds::writeSummary(FILE *fp)
{
    for(auto &bb : bondeds_)
    {
        fprintf(fp, "Extracted %zu %s\n",
                bb.second.size(), interactionTypeToString(bb.first).c_str());
    }
}

} // namespace alexandria

