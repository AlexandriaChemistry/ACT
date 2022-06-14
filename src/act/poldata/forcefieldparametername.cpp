#include "forcefieldparametername.h"

namespace alexandria
{

const char *lj_name[ljNR] = { "c6_ij", "c12_ij" };

const char *wbh_name[wbhNR] = { "sigma_ij", "epsilon_ij", "gamma_ij" };
    
const char *coul_name[coulNR] = { "zeta_i", "zeta_j" };

const char *bond_name[bondNR] = { "kb", "bondlength" };

const char *angle_name[angleNR] = { "kt", "angle" };

const char *ub_name[ubNR] = { "kt", "angle", "r13", "kub" };

const char *ps_names[psNR] = { "angle", "rij0", "rjk0" };

const char *pol_name[polNR] = { "alpha", "kshell" };

const char *morse_name[morseNR] = { "beta", "De", "D0", "bondlength" };
    
const char *linang_name[linangNR] = { "a", "klin" };
    
const char *idih_name[idihNR] = { "kimp" }; 

const char *fdih_name[fdihNR] = { "c0", "c1", "c2", "c3" };

} // namespace alexandria
