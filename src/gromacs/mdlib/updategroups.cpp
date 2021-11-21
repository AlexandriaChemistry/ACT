/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/* \internal \file
 *
 * \brief Implements functions for generating update groups
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 */

#include "actpre.h"

#include "updategroups.h"

#include <cmath>

#include <unordered_map>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

/*! \brief Struct for returning an atom range */
struct AtomIndexExtremes
{
    int minAtom; //!< The minimum atom index
    int maxAtom; //!< The maximum atom index
};

/*! \brief Returns a list of update groups for \p moltype */
static RangePartitioning
makeUpdateGroups(gmx_unused const gmx_moltype_t            &moltype,
                 gmx_unused gmx::ArrayRef<const t_iparams>  iparams)
{
    RangePartitioning groups;

    return groups;
}

std::vector<RangePartitioning> makeUpdateGroups(const gmx_mtop_t &mtop)
{
    std::vector<RangePartitioning> updateGroups;

    bool                           systemSatisfiesCriteria = true;
    for (const gmx_moltype_t &moltype : mtop.moltype)
    {
        updateGroups.push_back(makeUpdateGroups(moltype, mtop.ffparams.iparams));

        if (updateGroups.back().numBlocks() == 0)
        {
            systemSatisfiesCriteria = false;
        }
    }

    if (!systemSatisfiesCriteria)
    {
        updateGroups.clear();
    }

    return updateGroups;
}

/*! \brief When possible, computes the maximum radius of constrained atom in an update group
 *
 * Supports groups with 2 or 3 atoms where all partner atoms are connected to
 * each other by angle potentials. The temperature is used to compute a radius
 * that is not exceeded with a chance of 10^-9. Note that this computation
 * assumes there are no other strong forces working on these angular
 * degrees of freedom.
 * The return value is -1 when all partners are not connected to each other
 * by one angle potential, when a potential is perturbed or when an angle
 * could reach more than 180 degrees.
 */
template <int numPartnerAtoms> static real
constraintGroupRadius(const gmx_moltype_t                     &moltype,
                      gmx::ArrayRef<const t_iparams>           iparams,
                      const int                                centralAtom,
                      const t_blocka                          &at2con,
                      const std::unordered_multimap<int, int> &angleIndices,
                      const real                               constraintLength,
                      const real                               temperature)
{
    const int numConstraints = at2con.index[centralAtom + 1] - at2con.index[centralAtom];
    GMX_RELEASE_ASSERT(numConstraints == numPartnerAtoms, "We expect as many constraints as partner atoms here");

    std::array<int, numPartnerAtoms> partnerAtoms;
    for (int i = 0; i < numPartnerAtoms; i++)
    {
        const int ind = at2con.a[at2con.index[centralAtom] + i]*3;
        if (ind >= moltype.ilist[F_CONSTR].size())
        {
            /* This is a flexible constraint, we don't optimize for that */
            return -1;
        }
        const int a1    = moltype.ilist[F_CONSTR].iatoms[ind + 1];
        const int a2    = moltype.ilist[F_CONSTR].iatoms[ind + 2];
        partnerAtoms[i] = (a1 == centralAtom ? a2 : a1);
    }

    const InteractionList            &angles      = moltype.ilist[F_ANGLES];
    auto                              range       = angleIndices.equal_range(centralAtom);
    int                               angleType   = -1;
    std::array<int, numPartnerAtoms>  numAngles   = { 0 };
    bool                              areSameType = true;
    for (auto it = range.first; it != range.second; ++it)
    {
        /* Check if the outer atoms in the angle are both partner atoms */
        int numAtomsFound = 0;
        for (int ind = it->second + 1; ind < it->second + 4; ind += 2)
        {
            for (const int &partnerAtom : partnerAtoms)
            {
                if (angles.iatoms[ind] == partnerAtom)
                {
                    numAtomsFound++;
                }
            }
        }
        if (numAtomsFound == 2)
        {
            /* Check that the angle potentials have the same type */
            if (angleType == -1)
            {
                angleType = angles.iatoms[it->second];
            }
            else if (angles.iatoms[it->second] != angleType)
            {
                areSameType = false;
            }
            /* Count the number of angle interactions per atoms */
            for (int ind = it->second + 1; ind < it->second + 4; ind += 2)
            {
                for (size_t i = 0; i < partnerAtoms.size(); i++)
                {
                    if (angles.iatoms[ind] == partnerAtoms[i])
                    {
                        numAngles[i]++;
                    }
                }
            }
        }
    }

    bool criteriaSatisfied = areSameType;
    for (int i = 0; i < numPartnerAtoms; i++)
    {
        if (numAngles[i] != numPartnerAtoms - 1)
        {
            criteriaSatisfied = false;
        }
    }

    /* We don't bother optimizing the perturbed angle case */
    const t_iparams &angleParams = iparams[angleType];
    if (criteriaSatisfied &&
        angleParams.harmonic.rB == angleParams.harmonic.rA &&
        angleParams.harmonic.krB == angleParams.harmonic.krA)
    {
        /* Set number of stddevs such that change of exceeding < 10^-9 */
        constexpr real c_numSigma = 6.0;
        /* Compute the maximally stretched angle */
        const real     eqAngle  = angleParams.harmonic.rA*DEG2RAD;
        const real     fc       = angleParams.harmonic.krA;
        const real     maxAngle = eqAngle + c_numSigma*BOLTZ*temperature/((numPartnerAtoms - 1)*fc);
        if (maxAngle >= M_PI)
        {
            return -1;
        }

        if (numPartnerAtoms == 2)
        {
            /* With two atoms constrainted to a cental atom we have a triangle
             * with two equal sides because the constraint type is equal.
             * Return the distance from the COG to the farthest two corners,
             * i.e. the partner atoms.
             */
            real distMidPartner = std::sin(0.5*maxAngle)*constraintLength;
            real distCentralMid = std::cos(0.5*maxAngle)*constraintLength;
            real distCogMid     = distCentralMid*numPartnerAtoms/(numPartnerAtoms + 1);
            real distCogPartner = std::sqrt(gmx::square(distMidPartner) + gmx::square(distCogMid));

            return distCogPartner;
        }
        else if (numPartnerAtoms == 3)
        {
            /* With three atoms constrainted to a cental atom we have the
             * largest COG-atom distance when one partner atom (the outlier)
             * moves out by stretching its two angles by an equal amount.
             * The math here gets a bit more involved, but it is still
             * rather straightforward geometry.
             * We first compute distances in the plane of the three partners.
             * Here we have two "equilibrium" partners and one outlier.
             * We make use of the "Mid" point between the two "Eq" partners.
             * We project the central atom on this plane.
             * Then we compute the distance of the central atom to the plane.
             * The distance of the COG to the ourlier is returned.
             */
            real halfDistEqPartner    = std::sin(0.5*eqAngle)*constraintLength;
            real distPartnerOutlier   = 2*std::sin(0.5*maxAngle)*constraintLength;
            real distMidOutlier       = std::sqrt(gmx::square(distPartnerOutlier) - gmx::square(halfDistEqPartner));
            real distMidCenterInPlane = 0.5*(distMidOutlier - gmx::square(halfDistEqPartner)/distMidOutlier);
            real distCentralToPlane   = std::sqrt(gmx::square(constraintLength) - gmx::square(halfDistEqPartner) - gmx::square(distMidCenterInPlane));
            real distCogOutlierH      = distCentralToPlane/(numPartnerAtoms + 1);
            real distCogOutlierP      = distMidOutlier - (distMidOutlier + distMidCenterInPlane)/(numPartnerAtoms + 1);
            real distCogOutlier       = std::sqrt(gmx::square(distCogOutlierH) + gmx::square(distCogOutlierP));

            return distCogOutlier;
        }
        else
        {
            GMX_RELEASE_ASSERT(false, "Only 2 or 3 constraints are supported here");
        }
    }

    return -1;
}

/*! \brief Returns the maximum update group radius for \p moltype */
real
computeMaxUpdateGroupRadius(const gmx_mtop_t                       &mtop,
                            gmx::ArrayRef<const RangePartitioning>  updateGroups,
                            gmx_unused real                                    temperature)
{
    if (updateGroups.empty())
    {
        return 0;
    }

    GMX_RELEASE_ASSERT(static_cast<size_t>(updateGroups.size()) == mtop.moltype.size(),
                       "We need one update group entry per moleculetype");

    real maxRadius = 1;

    return maxRadius;
}

} // namespace gmx
