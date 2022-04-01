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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Marrades <julian.marrades@hotmail.es>
 */

#include "molhandler.h"

namespace alexandria
{

double MolHandler::computeHessian(      MyMol               *mol,
                                  const t_commrec           *crtmp,
                                  const std::vector<int>    &atomIndex,
                                        MatrixWrapper       *hessian,
                                        std::vector<double> *forceZero)
{
    std::vector<double> fneg;
    fneg.resize(DIM*atomIndex.size(), 0.0);
    
    double shellForceRMS = 0;
    (void) mol->calculateEnergy(crtmp, &shellForceRMS);
    double epot0  = mol->enerd_->term[F_EPOT];

    // Store central force vector
    forceZero->clear();
    for(auto &atom : atomIndex)
    {
        for(int m = 0; m < DIM; m++)
        {
            forceZero->push_back(mol->f_[atom][m]);
        }
    }
    double    stepSize     = 1e-6; // 1 pm
    for(size_t ai = 0; ai < atomIndex.size(); ai++)
    {
        auto atomI = atomIndex[ai];
        for(int atomXYZ = 0; atomXYZ < DIM; atomXYZ++)
        {
            double xyzRef = mol->state_->x[atomI][atomXYZ];
            for(int delta  = 0; delta <= 1; delta++)
            {
                mol->state_->x[atomI][atomXYZ] = xyzRef + (2*delta-1)*stepSize;
                (void) mol->calculateEnergy(crtmp, &shellForceRMS);
                if (delta == 0)
                {
                    for(size_t aj = 0; aj < atomIndex.size(); aj++)
                    {
                        int atomJ = atomIndex[aj];
                        for (int d = 0; d < DIM; d++)
                        {
                            fneg[aj*DIM+d] = mol->f_[atomJ][d];
                        }
                    }
                }
            }
            mol->state_->x[atomI][atomXYZ] = xyzRef;  // Restore
            {
                int col = ai*DIM+atomXYZ;
                for(size_t aj = 0; aj < atomIndex.size(); aj++)
                {
                    int    atomJ = atomIndex[aj];
                    for(int d = 0; d < DIM; d++)
                    {
                        int    row   = aj*DIM+d;
                        double value = -(mol->f_[atomJ][d]-fneg[row])/(2*stepSize);
                        if (false && debug)
                        {
                            fprintf(debug, "Setting H[%2d][%2d] = %10g\n", row, col, value); 
                        }
                        hessian->set(row, col, value);
                    }
                }
            }
        }
    }
    return epot0;
}

} // namespace alexandria
