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
                t_inputrec          *inputrec,
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
        // Needed for GenerateCharges
        CommunicationRecord cr;
        // Generate charges and topology
        fill_inputrec(inputrec);
        // Charge gen params
        auto alg = ChargeGenerationAlgorithm::NONE;
        std::vector<double> qcustom;
        bool qSymm = false;
        matrix box;
        bool readOK = readBabel(pd, dataName.c_str(), &molprops, molname, molname,
                                conf, &method, &basis,
                                maxpot, nsymm, jobtype, &qtot, false, box);
        EXPECT_TRUE(readOK);
        if (readOK)
        {
            gmx::MDLogger mdlog {};
            for(auto &molprop : molprops)
            {
                ACTMol mm;
                mm.Merge(&molprop);
                auto imm = mm.GenerateTopology(stdout, pd,
                                               missingParameters::Error, true);
                EXPECT_TRUE(immStatus::OK ==imm);
                mm.setInputrec(inputrec);
                mm.symmetrizeCharges(pd, qSymm, nullptr);
                std::vector<gmx::RVec> forces(mm.atomsConst().size());
                std::vector<gmx::RVec> coords = mm.xOriginal();
                mm.GenerateCharges(pd, fcomp, mdlog, &cr, alg, qType::Calc, qcustom, &coords, &forces);
                mps->push_back(mm);
            }
        }
    }    
}
