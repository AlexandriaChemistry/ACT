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

#include "actpre.h"

#include "babel_io.h"

#include "config.h"

#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <fstream>
#include <iostream>

#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "molprop.h"
#include "molprop_util.h"
#include "mymol.h"
#include "poldata.h"
#include "utility/stringutil.h"
#include "utility/units.h"

// Include Open Babel classes for OBMol and OBConversion
#if HAVE_LIBOPENBABEL3
// Hack to make this compile!
#undef ANGSTROM
#ifdef HAVE_SYS_TIME_H
#define KOKO HAVE_SYS_TIME_H
#undef HAVE_SYS_TIME_H
#endif
#include <openbabel/atom.h>
#include <openbabel/babelconfig.h>
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/data_utilities.h>
#include <openbabel/elements.h>
#include <openbabel/forcefield.h>
#include <openbabel/generic.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/residue.h>
#include <openbabel/typer.h>
#include <openbabel/math/vector3.h>

#ifdef KOKO
#ifndef HAVE_SYS_TIME_H
#define HAVE_SYS_TIME_H KOKO
#endif
#undef KOKO
#endif

using namespace alexandria;

static inline double A2PM(double a) {return a*1.0e+2; }                /* Angstrom to pm */

static inline double NM_cubed_to_A_cubed(double a) {return a*1.0e+3; } /* nm^3 to A^3 */

enum einformat{
    einfGaussian    = 0,
    einfNotGaussian = 1,
    einfNR
};

BabelFile::BabelFile(BabelFileType      ftype,
                     const std::string &ext,
                     const std::string &InFormat)
    :
      ftype_(ftype),
      ext_(ext),
      InFormat_(InFormat)
{}

BabelFiles::BabelFiles ()
{
    bfiles_.push_back(BabelFile(ebftPDB,  ".pdb",  "pdb"));
    bfiles_.push_back(BabelFile(ebftXYZ,  ".xyz",  "xyz"));
    bfiles_.push_back(BabelFile(ebftSDF,  ".sdf",  "sdf"));
    bfiles_.push_back(BabelFile(ebftMOL,  ".mol",  "mol"));
    bfiles_.push_back(BabelFile(ebftMOL2, ".mol2", "mol2"));
    bfiles_.push_back(BabelFile(ebftG09,  ".log",  "g03"));
}

BabelFileIterator BabelFiles::findBabelFile(const std::string &fn)
{
    auto              extension = fn.substr(fn.find_last_of("."));
    BabelFileIterator fb        = bfiles_.begin(), fe = bfiles_.end();
    return std::find_if(fb, fe, [extension](const BabelFile &bf)
                        {
                            return (extension == bf.ext());
                        });
}

static void merge_electrostatic_potential(alexandria::MolProp                             *mpt,
                                          std::vector<alexandria::ElectrostaticPotential> &espv,
                                          int                                              natom,
                                          int                                              maxPotential)
{
    maxPotential = std::max(0, std::min(maxPotential, 100));
    int npot   = espv.size() - natom;
    int maxpot = (npot * maxPotential)/100;
    int mod    = npot / maxpot;
    int i      = 0;
    for (auto esi = espv.begin(); esi < espv.end(); esi++, i++)
    {
        if ((i < natom) || (((i-natom) % mod) == 0))
        {
            mpt->LastExperiment()->AddPotential(*esi);
        }
    }
}

static bool isGzipFile(const std::string &fileName,
                       std::string       *strippedFileName)
{
    auto gzPosition = fileName.find(".gz");
    if (gzPosition == fileName.size()-3)
    {
        strippedFileName->assign(fileName.substr(0, gzPosition));
        return true;
    }
    return false;
}

static OpenBabel::OBConversion *readBabel(const std::string &g09,
                                          OpenBabel::OBMol  *mol,
                                          einformat         *inputformat)
{
    std::string fileName(g09);
    if (!gmx_fexist(fileName.c_str()))
    {
        fileName += ".gz";
        printf("fileName2: %s\n", fileName.c_str());
        if (!gmx_fexist(fileName.c_str()))
        {
            return nullptr;
        }
    }
    
    std::ifstream g09f;
    std::string   strippedFileName;
    bool          isGzip = isGzipFile(fileName, &strippedFileName);
    
    g09f.open(fileName, std::ios::in);
    
    if (!g09f.is_open())
    {
        gmx_fatal(FARGS, "Cannot open file %s for reading", g09.c_str());
    }
    
    OpenBabel::OBConversion *conv       = new OpenBabel::OBConversion(&g09f, &std::cout); // Read from g09f
    auto                     babelfiles = BabelFiles();
    const char              *informat   = nullptr;
    if (isGzip)
    {
        informat = babelfiles.findBabelFile(strippedFileName)->informat().c_str();
    }
    else
    {
        informat = babelfiles.findBabelFile(g09)->informat().c_str();
    }

    if (strcmp (informat, "g03") == 0 || strcmp (informat, "g09") == 0)
    {
        *inputformat = einfGaussian;
    }

    if (conv->SetInFormat(informat, isGzip))
    {
        bool read_ok = false;
        try
        {
            read_ok = conv->Read(mol, &g09f);
        }
        catch (const std::exception &ex)
        {
            gmx::printFatalErrorMessage(stderr, ex);
        }

        if (read_ok)
        {
            g09f.close();
            return conv; // exit with success
        }
        else
        {
            fprintf(stderr, "Could not read input file %s with OpenBabel3.\n", g09.c_str());
        }
    }
    else
    {
        fprintf(stderr, "Input file %s has incomprehensible format.\n", g09.c_str());
    }
    g09f.close();
    return nullptr;
}

static void checkBondOrders(alexandria::MolProp *mpt)
{
    // We cannot be sure that this is the right one
    auto exper = mpt->LastExperiment();
    std::vector<std::string> atomName;
    for(auto &ca : exper->calcAtomConst())
    {
        atomName.push_back(ca.getObtype());
    }
    std::map<std::string, std::string> my_pairs =
        { { "no", "on" }, { "cm", "om" }, { "p5", "om" }, { "s6", "om" }, { "s4", "om" }, { "s3", "om" }, { "py", "om" }, { "cz", "n" } };
    for(const auto &my_pair : my_pairs)
    {
        for(size_t i = 0; i < atomName.size(); i++)
        {
            if (atomName[i] == my_pair.first)
            {
                std::vector<size_t> bondIndex;
                auto allBonds = mpt->bonds();
                for(size_t b = 0; b < (*allBonds).size(); b++)
                {
                    size_t ai = (*allBonds)[b].aI();
                    size_t aj = (*allBonds)[b].aJ();
                    if ((ai == i && atomName[aj] == my_pair.second) ||
                        (aj == i && atomName[ai] == my_pair.second))
                    {
                        bondIndex.push_back(b);
                    }
                }
                // The first atom has to have at least two bonds to the
                // second atom in the map.
                if (bondIndex.size() >= 2)
                {
                    for (size_t j = 0; j < bondIndex.size(); j++)
                    {
                        (*allBonds)[bondIndex[j]].setBondOrder(0, 1.5);
                    }
                }
            }
        }
    }
}

static bool getBondsFromOpenBabel(OpenBabel::OBMol    *mol,
                                  alexandria::MolProp *mpt,
                                  const char          *inputFile,
                                  bool                 changeAromaticBondOrders)
{
    // Bonds
    auto OBbi = mol->BeginBonds();
    if (OBbi != mol->EndBonds())
    {
        for (auto OBb = mol->BeginBond(OBbi); (nullptr != OBb); OBb = mol->NextBond(OBbi))
        {
            double bo = OBb->GetBondOrder();
            if (changeAromaticBondOrders && OBb->IsAromatic())
            {
                bo = 1.5;
            }
            alexandria::Bond ab(OBb->GetBeginAtom()->GetIndex(),
                                OBb->GetEndAtom()->GetIndex(),
                                bo);
            mpt->AddBond(ab);
        }
        if (changeAromaticBondOrders)
        {
            checkBondOrders(mpt);
        }
        return true;
    }
    else
    {
        fprintf(stderr, "No bond is found for %s\n", inputFile);
        return false;
    }
}

bool readBabel(const char          *g09,
               alexandria::MolProp *mpt,
               const char          *molnm,
               const char          *iupac,
               const char          *conformation,
               const char          *basisset,
               int                  maxPotential,
               int                  nsymm,
               const char          *jobType,
               double              *qtot,
               bool                 addHydrogen)
{
    std::string                formula;
    std::string                attr;
    std::string                value;
    einformat                  inputformat = einfNotGaussian;
    const char                *reference   = "Ghahremanpour2022a";
    const char                *mymol       = "AMM";
    const char                *myprogram   = "ACT2022";
    const char                *mybasis     = "";


    /* Variables to read a Gaussian log file */
    char                      *g09ptr;
    OpenBabel::OBMol           mol;
    OpenBabel::OBAtomIterator  OBai;
    OpenBabel::OBBondIterator  OBbi;
    OpenBabel::OBPairData     *OBpd;
    OpenBabel::OBPcharge      *OBpc       = new OpenBabel::OBPcharge();
    OpenBabel::OBVectorData   *dipole     = new OpenBabel::OBVectorData;
    OpenBabel::OBMatrixData   *quadrupole = new OpenBabel::OBMatrixData;
    OpenBabel::OBMatrixData   *pol_tensor = new OpenBabel::OBMatrixData;
    OpenBabel::OBFreeGrid     *esp;
    std::vector<alexandria::ElectrostaticPotential> espv;

    alexandria::JobType         jobtype = alexandria::string2jobType(jobType);
    OpenBabel::OBConversion    *conv    = readBabel(g09, &mol, &inputformat);
    if (nullptr == conv)
    {
        fprintf(stderr, "Failed reading %s\n", g09);
        return false;
    }
    delete conv;
    conv = new OpenBabel::OBConversion(&std::cin, &std::cout);

    // Chemical Categories
    if (conv->SetOutFormat("fpt"))
    {
        std::vector<std::string> excluded_categories = {
            ">", "C_ONS_bond", "Rotatable_bond", "Conjugated_double_bond", "Conjugated_triple_bond",
            "Chiral_center_specified", "Cis_double_bond", "Bridged_rings", "Conjugated_tripple_bond",
            "Trans_double_bond"
        };

        conv->AddOption("f", OpenBabel::OBConversion::OUTOPTIONS, "FP4");
        conv->AddOption("s");
        conv->Convert();
        OpenBabel::OBMol         mol_copy   = mol;  // We need a copy here because WriteString removes the H.
        std::vector<std::string> categories = gmx::splitString(conv->WriteString(&mol_copy, false));
        for (const auto &category : categories)
        {
            size_t j;
            for (j = 0; j < excluded_categories.size(); j++)
            {
                if (excluded_categories[j] == category)
                {
                    break;
                }
            }
            if (j == excluded_categories.size())
            {
                std::string dup = category;
                std::replace_if(dup.begin(), dup.end(), [](const char c) {return c == '_'; }, ' ');
                mpt->AddCategory(dup);
            }
        }
    }

    // Basis Set
    std::string basis;
    OBpd = (OpenBabel::OBPairData *)mol.GetData("basis");
    if ((nullptr != basisset) && (strlen(basisset) > 0))
    {
        basis.assign(basisset);
    }
    else if (nullptr != OBpd)
    {
        basis = OBpd->GetValue();
        size_t p = basis.find(" (5D, 7F)");
        if (p != basis.npos)
        {
            basis.erase(p, basis.npos);
        }
    }
    else
    {
        basis.assign(mybasis);
    }

    // QM Program
    std::string program(myprogram);
    OBpd = (OpenBabel::OBPairData *)mol.GetData("program");
    if (nullptr != OBpd)
    {
        program.assign(OBpd->GetValue());
    }

    // Method
    std::string method;
    OBpd = (OpenBabel::OBPairData *)mol.GetData("method");
    if (nullptr != OBpd)
    {
        method.assign(OBpd->GetValue());
    }
    g09ptr = (char *) strrchr(g09, '/');
    if (nullptr == g09ptr)
    {
        g09ptr = (char *)g09;
    }
    else
    {
        g09ptr++;
        if (strlen(g09ptr) == 0)
        {
            g09ptr = (char *)g09;
        }
    }

    alexandria::Experiment exp(program, method, basis, reference, conformation, g09ptr, jobtype);
    mpt->AddExperiment(exp);
    mpt->SetMass(mol.GetMolWt());
    mpt->SetMultiplicity(mol.GetTotalSpinMultiplicity());
    mpt->SetFormula(mol.GetFormula());
    // We don't just set this here, since the user may override the value
    // However, it seems that OB does not extract this correctly from
    // the input, unless it is a Gaussian log file.
    *qtot = mol.GetTotalCharge();

    if (nullptr != molnm)
    {
        mpt->SetMolname(molnm);
    }
    else
    {
        mpt->SetMolname(mymol);
    }

    if (nullptr != iupac)
    {
        mpt->SetIupac(iupac);
    }
    else
    {
        mpt->SetIupac(mymol);
    }
    
    // Units for conversion
    std::string energyUnit("kcal/mol");
    std::string entropyUnit("cal/mol K");
    std::string qm_type("electronic");

    // Thermochemistry
    if (inputformat == einfGaussian)
    {
        double              temperature = 0, DeltaHf0 = 0, DeltaHfT = 0;
        double              DeltaGfT    = 0, DeltaSfT = 0, S0T      = 0;
        double              CVT         = 0, CPT      = 0, ZPE      = 0;
        std::vector<double> Scomponents;
        if (extract_thermochemistry(mol, false, &nsymm,
                                    0, 0.0,
                                    &temperature,
                                    &DeltaHf0,
                                    &DeltaHfT,
                                    &DeltaGfT,
                                    &DeltaSfT,
                                    &S0T,
                                    &CVT,
                                    &CPT,
                                    Scomponents,
                                    &ZPE))
        {
            {
                auto mpo = MolPropObservable::DHFORM;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, 0, ePhase::GAS,
                                                           alexandria::convertToGromacs(DeltaHf0, energyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            {
                auto mpo = MolPropObservable::DHFORM;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, temperature, ePhase::GAS,
                                                           alexandria::convertToGromacs(DeltaHfT, energyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            {
                auto mpo = MolPropObservable::DGFORM;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, temperature, ePhase::GAS,
                                                           alexandria::convertToGromacs(DeltaGfT, energyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            {
                auto mpo = MolPropObservable::DSFORM;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, temperature, ePhase::GAS,
                                                           alexandria::convertToGromacs(DeltaSfT, entropyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            {
                auto mpo = MolPropObservable::ENTROPY;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, temperature, ePhase::GAS,
                                                           alexandria::convertToGromacs(S0T, entropyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            {
                auto mpo = MolPropObservable::CP;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, temperature, ePhase::GAS,
                                                           alexandria::convertToGromacs(CPT, entropyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
            std::vector<MolPropObservable> mpos = { 
                MolPropObservable::STRANS, 
                MolPropObservable::SROT,
                MolPropObservable::SVIB
            };
            for (size_t i = 0; (i < mpos.size()); i++)
            {
                auto me = new alexandria::MolecularEnergy(mpos[i], qm_type, temperature, ePhase::GAS,
                                                          alexandria::convertToGromacs(Scomponents[i], entropyUnit), 0);
                mpt->LastExperiment()->addProperty(mpos[i], me);
            }
            {
                auto mpo = MolPropObservable::ZPE;
                auto me  = new alexandria::MolecularEnergy(mpo, qm_type, 0, ePhase::GAS,
                                                           alexandria::convertToGromacs(ZPE, energyUnit), 0);
                mpt->LastExperiment()->addProperty(mpo, me);
            }
        }
    }

    // HF Eenergy
    { 
        auto mpo = MolPropObservable::HF;
        auto me  = new alexandria::MolecularEnergy(mpo, qm_type, 0, ePhase::GAS, 
                                                   alexandria::convertToGromacs(mol.GetEnergy(), energyUnit), 0);
        mpt->LastExperiment()->addProperty(mpo, me);
    }

    if (addHydrogen)
    {
        mol.AddHydrogens();
    }
    // Atoms
    const std::string forcefield("alexandria");
    auto       *ff         = OpenBabel::OBForceField::FindForceField(forcefield);
    if (ff && (ff->Setup(mol)))
    {
        ff->GetAtomTypes(mol);
        FOR_ATOMS_OF_MOL (atom, mol)
        {
            OpenBabel::OBPairData *type = (OpenBabel::OBPairData*) atom->GetData("FFAtomType");
            if (nullptr == type)
            {
                fprintf(stderr, "Error: OpenBabel cannot find the %s atomtype for atom %d\n",
                        forcefield.c_str(), static_cast<int>(atom->GetIdx()));
                return false;
            }
            if (nullptr != debug)
            {
                fprintf(debug, "atom %d gafftype %s OBtype %s\n", atom->GetIdx(), type->GetValue().c_str(), atom->GetType());
            }

            alexandria::CalcAtom ca(OpenBabel::OBElements::GetSymbol(atom->GetAtomicNum()),
                                    type->GetValue(), atom->GetIdx());

            ca.SetUnit("pm");
            ca.SetCoords(A2PM(atom->x()), A2PM(atom->y()), A2PM(atom->z()));
            auto myres = atom->GetResidue();
            ca.SetResidue(myres->GetName(), myres->GetNum());
            ca.SetChain(myres->GetChainNum(), myres->GetChain());
            if (inputformat == einfGaussian)
            {
                for (const auto &cs : qTypes())
                {
                    if (cs.first == qType::Mulliken ||
                        cs.first == qType::ESP ||
                        cs.first == qType::Hirshfeld ||
                        cs.first == qType::CM5)
                    {
                        std::string qstr = cs.second;
                        qstr.append(" charges");
                        OBpd = (OpenBabel::OBPairData *) mol.GetData(qstr.c_str());
                        if (nullptr != OBpd)
                        {
                            OBpc = (OpenBabel::OBPcharge *) mol.GetData(qstr.c_str());
                            if (OBpc && !OBpc->GetPartialCharge().empty())
                            {
                                ca.AddCharge(stringToQtype(cs.second),
                                             OBpc->GetPartialCharge()[atom->GetIdx()-1]);
                            }
                            else
                            {
                                fprintf(stderr, "Inconsistency reading %s from %s", qstr.c_str(), g09);
                                return false;
                                
                            }
                        }
                    }
                }
            }
            mpt->LastExperiment()->AddAtom(ca);
        }
    }
    else
    {
        fprintf(stderr, "Error: OpenBabel cannot read the '%s' force field\n",
                forcefield.c_str());
        return false;
    }

    // Bonds
    getBondsFromOpenBabel(&mol, mpt, g09, forcefield.compare("alexandria") == 0);

    // Dipole
    dipole = (OpenBabel::OBVectorData *) mol.GetData("Dipole Moment");
    if (nullptr != dipole)
    {
        OpenBabel::vector3            v3 = dipole->GetData();
        auto dp = new alexandria::MolecularDipole(qm_type,
                                                  0.0,
                                                  v3.GetX(),
                                                  v3.GetY(),
                                                  v3.GetZ(),
                                                  v3.length(),
                                                  0.0);
        mpt->LastExperiment()->addProperty(MolPropObservable::DIPOLE, dp);
    }
    
    // Quadrupole
    quadrupole = (OpenBabel::OBMatrixData *) mol.GetData("Traceless Quadrupole Moment");
    if (nullptr != quadrupole)
    {
        OpenBabel::matrix3x3            m3 = quadrupole->GetData();
        double                          mm[9];
        m3.GetArray(mm);
        auto mq = new alexandria::MolecularQuadrupole(qm_type, 0.0, mm[0], mm[4], mm[8],
                                                      mm[1], mm[2], mm[5]);
        mpt->LastExperiment()->addProperty(MolPropObservable::QUADRUPOLE, mq);
    }

    // Polarizability
    pol_tensor = (OpenBabel::OBMatrixData *) mol.GetData("Exact polarizability");
    if (nullptr != pol_tensor)
    {
        OpenBabel::matrix3x3 m3 = pol_tensor->GetData();
        double               mm[9], alpha, fac;
        int                  i;
        m3.GetArray(mm);
        fac = NM_cubed_to_A_cubed(pow(alexandria::convertToGromacs(1, "Bohr"), 3));
        for (i = 0; i < 9; i++)
        {
            mm[i] *= fac;
        }
        alpha = (mm[0]+mm[4]+mm[8])/3.0;

        auto mdp = new alexandria::MolecularPolarizability(qm_type, 0.0, mm[0], mm[4], mm[8],
                                                           mm[1], mm[2], mm[5], alpha, 0);
        mpt->LastExperiment()->addProperty(MolPropObservable::POLARIZABILITY, mdp);
    }

    // Electrostatic potential
    esp = (OpenBabel::OBFreeGrid *) mol.GetData("Electrostatic Potential");
    if (nullptr != esp && maxPotential > 0)
    {
        OpenBabel::OBFreeGridPoint        *fgp;
        OpenBabel::OBFreeGridPointIterator fgpi;
        std::string                        xyz_unit("pm");
        std::string                        V_unit("Hartree/e");
        int                                espid = 0;

        fgpi = esp->BeginPoints();
        for (fgp = esp->BeginPoint(fgpi); (nullptr != fgp); fgp = esp->NextPoint(fgpi))
        {
            alexandria::ElectrostaticPotential ep(xyz_unit, V_unit, ++espid,
                                                  A2PM(fgp->GetX()),
                                                  A2PM(fgp->GetY()),
                                                  A2PM(fgp->GetZ()),
                                                  fgp->GetV());
            espv.push_back(ep);
        }
        merge_electrostatic_potential(mpt, espv, mol.NumAtoms(), maxPotential);
    }
    return true;
}

bool SetMolpropAtomTypesAndBonds(alexandria::MolProp *mmm)
{
    OpenBabel::OBMol mol;
    mol.BeginModify();
    mol.SetTotalCharge(mmm->totalCharge());
    mol.SetTotalSpinMultiplicity(mmm->getMultiplicity());
    auto ei  = mmm->experiment()->begin();
    mol.ReserveAtoms(ei->NAtom());
    int  idx = 0;
    for (auto &ca : ei->calcAtomConst())
    {
        int  atomicNum = OpenBabel::OBElements::GetAtomicNum(ca.getName().c_str());
        auto atom      = mol.NewAtom(idx++);
        atom->SetVector(ca.getX(), ca.getY(), ca.getZ());
        atom->SetAtomicNum(atomicNum);
    }
    mol.ConnectTheDots();
    mol.PerceiveBondOrders();
    mol.EndModify();
    mmm->SetFormula(mol.GetFormula());
    mmm->SetMass(mol.GetMolWt());
    const char *forcefield = "alexandria";
    auto        pFF        = OpenBabel::OBForceField::FindForceField(forcefield);
    if (!pFF)
    {
        fprintf(stderr, "could not find forcefield %s\n", forcefield);
        return false;
    }
    if (debug)
    {
        pFF->SetLogFile(&std::cout);
        pFF->SetLogLevel(OBFF_LOGLVL_HIGH);
    }
    if (!pFF->Setup(mol))
    {
        fprintf(stderr, "could not setup force field %s.\n", forcefield);
        return false;
    }
    pFF->GetAtomTypes(mol);
    auto ca = ei->calcAtom()->begin();
    for (OpenBabel::OBMolAtomIter a(mol); a; ++a, ++ca)
    {
        auto type = static_cast<OpenBabel::OBPairData*>(a->GetData("FFAtomType"));
        ca->setObtype(type->GetValue());
    }
    // Bonds
    getBondsFromOpenBabel(&mol, mmm, "psi4 input file", true);

    return true;
}

#endif
