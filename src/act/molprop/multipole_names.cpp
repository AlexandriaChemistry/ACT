#include "act/molprop/multipole_names.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

std::map<MolPropObservable, std::vector<std::string> >        multipole_names;
std::map<MolPropObservable, std::map<std::string, int> >      multipole_map;
std::map<MolPropObservable, std::map<int, std::string> >      multipole_rmap;
std::map<MolPropObservable, std::map<std::vector<int>, int> > multipole_index;

std::map<size_t, MolPropObservable> sizeMpo = {
    { 1, MolPropObservable::DIPOLE       },
    { 2, MolPropObservable::QUADRUPOLE   },
    { 3, MolPropObservable::OCTUPOLE     },
    { 4, MolPropObservable::HEXADECAPOLE }
};

static void mk_names()
{
    if (multipole_names.size() == 0)
    {
        for (auto &mpo : mpoMultiPoles)
        {
            std::vector<std::string> dummy;
            multipole_names.insert({mpo, dummy});
            std::map<std::string, int> mymap;
            multipole_map.insert({mpo, mymap});
            std::map<int, std::string> myrmap;
            multipole_rmap.insert({mpo, myrmap});
            std::map<std::vector<int>, int> myindex;
            multipole_index.insert({mpo, myindex});
        }
        int qindex = 0;
        int oindex = 0;
        int hindex = 0;
        std::array<std::string, 3> xyz = { "x", "y", "z" };
        for (int m = 0; m < DIM; m++)
        {
            auto dmpo = MolPropObservable::DIPOLE;
            multipole_names[dmpo].push_back(xyz[m]);
            multipole_index[dmpo].insert({ { m }, m});
            multipole_rmap[dmpo].insert({m, xyz[m]});
            multipole_map[dmpo].insert({xyz[m], m});
            for (int n = m; n < DIM; n++)
            {
                std::string qentry = xyz[m]+xyz[n];
                auto qmpo = MolPropObservable::QUADRUPOLE;
                multipole_names[qmpo].push_back(qentry);
                multipole_index[qmpo].insert({ { m, n }, qindex});
                multipole_rmap[qmpo].insert({qindex, qentry});
                multipole_map[qmpo].insert({qentry, qindex++});
                for (int o = n; o < DIM; o++)
                {
                    std::string oentry = xyz[m]+xyz[n]+xyz[o];
                    auto ompo = MolPropObservable::OCTUPOLE;
                    multipole_names[ompo].push_back(oentry);
                    multipole_index[ompo].insert({ { m, n, o }, oindex});
                    multipole_rmap[ompo].insert({oindex, oentry});
                    multipole_map[ompo].insert({oentry, oindex++});
                    for (int p = o; p < DIM; p++)
                    {
                        std::string hentry(xyz[m]+xyz[n]+xyz[o]+xyz[p]);
                        auto hmpo = MolPropObservable::HEXADECAPOLE;
                        multipole_names[hmpo].push_back(hentry);
                        multipole_index[hmpo].insert({ { m, n, o, p }, hindex});
                        multipole_rmap[hmpo].insert({hindex, hentry});
                        multipole_map[hmpo].insert({hentry, hindex++});
                    }
                }
            }
        }
    }
}

bool multipoleIndex(const std::string &myid, int *index)
{
    mk_names();
    if (sizeMpo.find(myid.size()) != sizeMpo.end())
    {
        auto mpo = sizeMpo[myid.size()];
        auto dm  = multipole_map[mpo].find(myid);
        if (multipole_map[mpo].end() != dm)
        {
            *index = dm->second;
            return true;
        }
    }
    return false;
}

const std::vector<std::string> &multipoleNames(MolPropObservable mpo)
{
    mk_names();
    return multipole_names[mpo];
}

const std::string &multipoleName(const std::vector<int> &m)
{
    mk_names();
    if (sizeMpo.find(m.size()) != sizeMpo.end())
    {
        auto mpo   = sizeMpo[m.size()];
        int  index = multipoleIndex(m);
        return multipole_rmap[mpo][index];
    }
    GMX_THROW(gmx::InternalError(gmx::formatString("Don't know how to handle multipoles of order %zu", m.size()).c_str()));
}

int multipoleIndex(const std::vector<int> &m)
{
    mk_names();
    return multipole_index[sizeMpo[m.size()]][m];
}

std::vector<std::string> formatMultipole(MolPropObservable          mpo,
                                         const std::vector<double> &values)
{
    std::vector<std::string>  fmp;
    if (multipole_names[mpo].size() != values.size())
    {
        GMX_THROW(gmx::InvalidInputError("Order does not match values"));
    }
    size_t delta = 4;
    for(size_t i = 0; i < values.size(); i+= delta)
    {
        std::string fmpi;
        for(size_t j = i; j < std::min(values.size(), i+delta) ; j++)
        {
            fmpi += gmx::formatString(" %4s: %8.3f", multipole_names[mpo][j].c_str(), values[j]);
        }
        fmp.push_back(fmpi);
    }
    return fmp;
}

void printMultipole(FILE                      *fp,
                    MolPropObservable          mpo,
                    const std::vector<double> &values)
{
    for(auto &f : formatMultipole(mpo, values))
    {
        fprintf(fp, "%s\n", f.c_str());
    }
}

} // namespace alexandria
