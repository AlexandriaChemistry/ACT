#include "actmol_util.h"
#include "act/alexandria/babel_io.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/atype_mapping.h"
#include "act/alexandria/actmol.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

void initACTMol(const char          *molname, 
                const ForceField    *pd,
                ForceComputer       *fcomp,
                std::vector<ACTMol> *mps)
    {
        int           maxpot   = 100;
        int           nsymm    = 0;
        const char   *conf     = (char *)"minimum";
        std::string   method, basis;
        const char   *jobtype  = (char *)"Opt";
        
        std::string   dataName = gmx::test::TestFileManager::getInputFilePath(molname);
        std::vector<alexandria::MolProp> molprops;
        double        qtot     = 0;
        // Charge gen params
        auto alg = ChargeGenerationAlgorithm::NONE;
        std::vector<double> qcustom;
        matrix box;
        bool readOK = readBabel(pd, dataName.c_str(), &molprops, molname, molname,
                                conf, &method, &basis,
                                maxpot, nsymm, jobtype, &qtot, false, box);
        EXPECT_TRUE(readOK);
        if (readOK)
        {
            for(auto &molprop : molprops)
            {
                ACTMol mm;
                mm.Merge(&molprop);
                auto imm = mm.GenerateTopology(stdout, pd,
                                               missingParameters::Error);
                EXPECT_TRUE(immStatus::OK ==imm);
                std::map<MolPropObservable, iqmType> iqmMap = 
                    {
                        { MolPropObservable::DELTAE0,           iqmType::QM },
                        { MolPropObservable::POTENTIAL,         iqmType::QM },
                        { MolPropObservable::INTERACTIONENERGY, iqmType::QM },
                        { MolPropObservable::DIPOLE,            iqmType::QM },
                        { MolPropObservable::QUADRUPOLE,        iqmType::QM },
                        { MolPropObservable::OCTUPOLE,          iqmType::QM },
                        { MolPropObservable::HEXADECAPOLE,      iqmType::QM },
                        { MolPropObservable::POLARIZABILITY,    iqmType::QM },
                        { MolPropObservable::CHARGE,            iqmType::QM }
                    };
                imm = mm.getExpProps(pd, iqmMap, 0);
                std::vector<gmx::RVec> forces(mm.atomsConst().size());
                std::vector<gmx::RVec> coords = mm.xOriginal();
                mm.GenerateCharges(pd, fcomp, alg, qType::Calc, qcustom, &coords, &forces);
                mps->push_back(mm);
            }
        }
    }    
}
