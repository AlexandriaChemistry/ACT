#!/usr/bin/env python3

import os

filelist = { 
    "1-butanol.pdb": "Butanol",
    "1-nitroethyne.pdb": "NitroEthyne",
    "guanidine.pdb": "Guanidine",
    "guanidinium.sdf": "Guanidinium",
    "dimethyl-carbonate.pdb": "DimethylCarbonate",
    "1--ethoxyethylphosphoryl-oxyethane-3-oep.log.pdb": "1EthoxyethylPhosphorylOxyethane",
    "1-amino-1-hydroxyguanidine-3-oep.log.pdb": "1Amino1Hydroxyguanidine",
    "1-1-dimethylguanidinium.sdf": "Dimethylguanidinium",
    "acetate-3-oep.log.pdb": "Acetate",
    "acetic-acid-3-oep.log.pdb": "AceticAcid",
    "acetone-3-oep.log.pdb": "Acetone",
    "diethyl-sulfate-3-oep.log.pdb": "DiethylSulfate",
    "disulfur-monoxide-3-oep.log.pdb": "DisulfurMonoxide",
    "ethyl-sulfate-3-oep.log.pdb": "EthylSulfate",
    "glutamate-3-oep.log.pdb": "Glutamate",
    "glutamic-acid-3-oep.log.pdb": "GlutamicAcid",
    "histidine-hdhe-3-oep.log.pdb": "HistidineHdHe",
    "nitromethane-3-oep.log.pdb": "Nitromethane",
    "sulfur-dioxide-3-oep.log.pdb": "SulfurDioxide",
    "thiazyl-fluoride-3-oep.log.pdb": "ThiazylFluoride",
    "water-3-oep.log.pdb": "Water",
    "aniline.sdf": "Aniline",
    "ethane-12-diamine.sdf": "EthaneDiamine",
    "123-trimethyl-imidazolium.sdf": "TrimethylImidazolium",
    "12-dimethyl-imidazole.sdf": "DimethylImidazole",
    "propa-1-2-dienylidenephosphane.sdf": "PropaDienelydinephosphane",
    "3H-diphosphole.sdf": "3H_Diphosphole",
    "acetonitrile.sdf": "Acetonitrile",
    "N--ethenyl-N-hydroxyethanimidamide.sdf": "N__ethenyl_N_hydroxyethanimidamide",
    "pyridine.sdf": "Pyridine",
    "fluoranthene.sdf": "Fluoranthene",
    "1-buten-3-yne.sdf": "1_Buten_3_yne",
    "cyclopropene.sdf": "Cyclopropene",
    "cyclobutene.sdf": "Cyclobutene",
    "formamide.sdf": "Formamide",
    "N-methylmethanimine.sdf": "N_MethylMethanimine",
    "ammonia.sdf": "Ammonia",
    "125-thiadiazole.sdf": "125_Thiadiazole",
    "acetaldehyde.sdf": "Acetaldehyde",
    "carbon-dioxide.sdf": "CarbonDioxide",
    "methyl-acetate.sdf": "MethylAcetate",
    "furan.sdf": "Furan",
    "dimethylether.sdf": "Dimethylether",
    "thiophene.sdf": "Thiophene",
    "isoxazole.sdf": "Isoxazole",
    "pyrazole.sdf": "Pyrazole",
    "phosphorus-nitride.sdf": "PhosphorusNitride",
    "oxophosphane.sdf": "Oxophosphane",
    "phosphine.sdf": "Phosphine",
    "phosphoethanoamine.sdf": "Phosphoethanoamine",
    "dihydrogen-phosphate.sdf": "DihydrogenPhosphate",
    "pentafluorophosphorane.sdf": "Pentafluorophosphorane", 
    "pentachlorophosphorane.sdf": "Pentachlorophosphorane",
    "formylphosphonic-acid.sdf": "FormylphosphonicAcid",
    "phosphinine.sdf": "Phosphinine",
    "sulfur-monoxide.sdf": "SulfurMonoxide",
    "sulfate.sdf": "Sulfate",
    "methanethiol.sdf": "Methanethiol", 
    "hydrogen-sulfide.sdf": "HydrogenSulfide",
    "hydrogen-chloride.sdf": "HydrogenChloride",
    "dimethyl-sulfide.sdf": "DimethylSulfide",
    "1-methylsulfinylethene.sdf": "1MethylSulfinylethene",
    "ethylsulfonylformaldehyde.sdf": "EthylSulfonylFormaldehyde",
    "thiazirene.sdf": "Thiazirene",
    "uracil.sdf": "Uracil"
}

mytests  = [ "AtomtypeTest", "BondtypeTest" ]
myfunc   = [ "testAtype", "testBtype" ]
typetest = "atom_bond_test_include.cpp"
with open(typetest, "w") as outf:
    for test in range(len(mytests)):
        for fl in filelist:
            outf.write("TEST_F (%s, %s)\n" % (mytests[test], filelist[fl]))
            outf.write("{\n    %s(\"%s\");\n}\n\n" % ( myfunc[0], fl))
