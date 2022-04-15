#include "fragmenthandler.h"

#include <vector>

#include "act/qgen/qgen_acm.h"
#include "act/molprop/fragment.h"
#include "gromacs/topology/atoms.h"

namespace alexandria
{    

FragmentHandler::FragmentHandler(const Poldata               *pd,
                                 const t_atoms               *atoms,
                                 const std::vector<Bond>     &bonds,
                                 const std::vector<Fragment> *fragments,
                                 const std::vector<int>      &shellRenumber)
{
    GMX_RELEASE_ASSERT(fragments != nullptr,
                       "Empty fragments passed. Wazzuppwitdat?");

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
                GMX_RELEASE_ASSERT(a < shellRenumber.size(), "Atom number out of range");
                anew = shellRenumber[a];
            }
            // We add the new atom index, but relative to the first atom
            // in the compound.
            toAdd.push_back(anew-atomStart_[ff]);
            if (pd->polarizable())
            {
                auto fa = pd->findParticleType(*atoms->atomtype[anew]);
                if (fa->hasInteractionType(InteractionType::POLARIZATION))
                {
                    // TODO remove assumption that shell is next to the atom in order
                    toAdd.push_back(anew+1-atomStart_[ff]);
                }
            }
        }
        // The next fragment, if there is one, starts after this one!
        atomStart_.push_back(atomStart_[ff]+toAdd.size());
        
        add_t_atoms(&FragAtoms_[ff], toAdd.size(), 1);
        // The stupid routine above will not handle the atomtype array
        snew(FragAtoms_[ff].atomtype, toAdd.size());
        int j = 0;
        for(auto &i : toAdd)
        {
            FragAtoms_[ff].atom[j]     = atoms->atom[i+atomStart_[ff]];
            FragAtoms_[ff].atomtype[j] = atoms->atomtype[i+atomStart_[ff]];
            j++;
        }
        QgenAcm_.push_back(QgenAcm(pd, &FragAtoms_[ff], f->charge()));
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
                auto offset = 0;
                if (ff > 0)
                {
                    offset = (*fragments)[ff-1].atoms().size();
                }
                Bond bb(ai - offset, aj - offset, b.bondOrder());
                bonds_[ff].push_back(bb);
            }
        }
        natoms_ += toAdd.size();
        ff      += 1;
    }
    if (debug)
    {
        fprintf(debug, "FragmentHandler: atoms->nr %d natoms %zu nbonds %zu nfragments %zu\n",
                atoms->nr, natoms_, bonds.size(), FragAtoms_.size());
    }
}

void FragmentHandler::fetchCharges(std::vector<double> *qq)
{
    qq->resize(natoms_, 0);
    size_t ff = 0;
    for (auto &fa : FragAtoms_)
    {
        for (int a = 0; a < fa.nr; a++)
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
                                       t_atoms                          *atoms)
{
    auto   eqgen = eQgen::OK;
    size_t ff    = 0;
    for (auto &fa : FragAtoms_)
    {
        gmx::HostVector<gmx::RVec> xx;
        xx.resizeWithPadding(fa.nr);
        for(int a = 0; a < fa.nr; a++)
        {
            copy_rvec(x[atomStart_[ff]+a], xx[a]);
        }
        eqgen = QgenAcm_[ff].generateCharges(fp, molname, pd, 
                                             &fa, xx, bonds_[ff]);
        if (eQgen::OK != eqgen)
        {
            break;
        }
        for(int a = 0; a < fa.nr; a++)
        {
            atoms->atom[atomStart_[ff]+a].q = fa.atom[a].q;
        }
        ff += 1;
    }
    return eqgen;
}

} // namespace alexandria
