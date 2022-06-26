#include "fragmenthandler.h"

#include <vector>

#include "act/qgen/qgen_acm.h"
#include "act/molprop/fragment.h"
#include "act/utility/stringutil.h"
#include "gromacs/topology/atoms.h"

namespace alexandria
{    

FragmentHandler::FragmentHandler(const Poldata               *pd,
                                 const std::vector<ActAtom>  &atoms,
                                 const std::vector<Bond>     &bonds,
                                 const std::vector<Fragment> *fragments,
                                 const std::vector<int>      &shellRenumber)
{
    GMX_RELEASE_ASSERT(fragments != nullptr,
                       "Empty fragments passed. Wazzuppwitdat?");
    GMX_RELEASE_ASSERT(fragments->size() > 0, "No fragments. Huh?");
    FragAtoms_.resize(fragments->size());
    bonds_.resize(fragments->size());
    natoms_           = 0;
    size_t  ff        = 0;
    atomStart_.push_back(0);
    for(auto f = fragments->begin(); f < fragments->end(); ++f)
    {
        if (f->atoms().size() == 0)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("No atoms in fragment %zu with formula %s", ff, f->formula().c_str()).c_str()));
        }
        // Count the number of atoms
        std::vector<int> toAdd;
        for(auto &a : f->atoms())
        {
            auto anew = a;
            if (!shellRenumber.empty())
            {
                GMX_RELEASE_ASSERT(a < static_cast<int>(shellRenumber.size()), "Atom number out of range");
                anew = shellRenumber[a];
            }
            // We add the new atom index, but relative to the first atom
            // in the compound.
            toAdd.push_back(anew-atomStart_[ff]);
            if (pd->polarizable())
            {
                auto fa = pd->findParticleType(atoms[anew].ffType());
                if (fa->hasInteractionType(InteractionType::POLARIZATION))
                {
                    // TODO remove assumption that shell is next to the atom in order
                    toAdd.push_back(anew+1-atomStart_[ff]);
                }
            }
        }
        // The next fragment, if there is one, starts after this one!
        atomStart_.push_back(atomStart_[ff]+toAdd.size());

        int j = 0;
        for(auto &i : toAdd)
        {
            auto fa      = pd->findParticleType(atoms[i+atomStart_[ff]].ffType());
            int  anumber = my_atoi(fa->optionValue("atomnumber").c_str(), "atomic number");
            ActAtom newat(atoms[i+atomStart_[ff]].name(),
                          fa->optionValue("element"),
                          atoms[i+atomStart_[ff]].ffType(),
                          atoms[i+atomStart_[ff]].pType(),
                          anumber,
                          atoms[i+atomStart_[ff]].mass(),
                          atoms[i+atomStart_[ff]].charge());
            FragAtoms_[ff].push_back(newat);
            j++;
        }
        QgenAcm_.push_back(QgenAcm(pd, FragAtoms_[ff], f->charge()));
        int offset = 0;
        for(int k = 0; k < static_cast<int>(ff); k++)
        {
            offset += (*fragments)[k].atoms().size();
        }
        for(const auto &b : bonds)
        {
            int ai = b.aI();
            int aj = b.aJ();
            if (std::find(f->atoms().begin(), f->atoms().end(),
                          ai) != f->atoms().end() &&
                std::find(f->atoms().begin(), f->atoms().end(),
                          aj) != f->atoms().end())
            {
                if (!shellRenumber.empty())
                {
                    //ai = shellRenumber[ai];
                    //aj = shellRenumber[aj];
                }
                // Bonds should be numbered from the start of the atom.
                // Now the shells should not be taken into account, since
                // the ACM code will do it.
                Bond bb(ai - offset, aj - offset, b.bondOrder());
                bonds_[ff].push_back(bb);
            }
        }
        natoms_ += toAdd.size();
        ff      += 1;
    }
    if (debug)
    {
        fprintf(debug, "FragmentHandler: atoms.size() %lu natoms %zu nbonds %zu nfragments %zu\n",
                atoms.size(), natoms_, bonds.size(), FragAtoms_.size());
    }
}

void FragmentHandler::fetchCharges(std::vector<double> *qq)
{
    qq->resize(natoms_, 0);
    size_t ff = 0;
    for (auto &fa : FragAtoms_)
    {
        for (size_t a = 0; a < fa.size(); a++)
        {
            (*qq)[atomStart_[ff] + a] = QgenAcm_[ff].getQ(a);
        }
        ff += 1;
    }
}

eQgen FragmentHandler::generateCharges(FILE                             *fp,
                                       const std::string                &molname,
                                       const gmx::HostVector<gmx::RVec> &x,
                                       const Poldata                    *pd,
                                       std::vector<ActAtom>             *atoms)
{
    auto   eqgen = eQgen::OK;
    size_t ff    = 0;
    for (auto &fa : FragAtoms_)
    {
        // TODO only copy the coordinates if there is more than one fragment.
        gmx::HostVector<gmx::RVec> xx;
        xx.resizeWithPadding(fa.size());
        for(size_t a = 0; a < fa.size(); a++)
        {
            copy_rvec(x[atomStart_[ff]+a], xx[a]);
        }
        eqgen = QgenAcm_[ff].generateCharges(fp, molname, pd, 
                                             &fa, xx, bonds_[ff]);
        if (eQgen::OK != eqgen)
        {
            break;
        }
        for(size_t a = 0; a < fa.size(); a++)
        {
            (*atoms)[atomStart_[ff]+a].setCharge(fa[a].charge());
        }
        ff += 1;
    }
    return eqgen;
}

} // namespace alexandria
