/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
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
 * The GPU interface of gpu.h for ACT built without GPU support.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gpu.h"

namespace alexandria
{

namespace gpu
{

std::string implementation()
{
    return "";
}

std::vector<std::string> devices()
{
    return {};
}

std::vector<DevicePeak> devicePeaks()
{
    return {};
}

std::vector<std::string> noGpuReasons()
{
    return {};
}

void selectGpuOfProcess()
{
}

std::string gpuOfProcess()
{
    return "";
}

int gpuIndex()
{
    return -1;
}

bool useGpu()
{
    return false;
}

bool setUseGpu(bool)
{
    return false;
}

void turnOffGpu(const std::string &)
{
}

bool gpuTurnedOff()
{
    return false;
}

std::string gpuOffReason()
{
    return "";
}

bool relaxShells(const RelaxMolecule *, size_t, const RelaxSettings &, size_t *, bool *)
{
    return false;
}

class System
{
};

System *createSystem(const SystemBlock *, size_t)
{
    return nullptr;
}

void destroySystem(System *system)
{
    delete system;
}

void updateShells(System *, size_t, const MovingShell *, size_t)
{
}

void updateParameters(System *, const float *, const double *, double, double) {}

void createRowRecipes(System *, const SetRecipe *, size_t, const RowRecipe *, size_t, const double *, size_t, size_t) {}

void updateRecipeValues(System *, const double *) {}

void createShellVsiteRecipes(System *, const ShellRecipe *, size_t, const VsiteRecipe *, size_t) {}

void updateMovingShells(System *, const MovingShell *, size_t) {}

void updateVsites(System *, const VsiteEntry *, size_t) {}

void createAcm(System *, const AcmSystem &) {}

void updateAcm(System *, const double *) {}

void updateAcmCoordinates(System *, const double *) {}

void downloadCharges(System *, double *) {}

void downloadAcmCharges(System *, double *) {}

void constructVsites(System *)
{
}

void computeAll(System *, double *const *)
{
}

void computeBlocksInPlace(System *, const char *, const double **, bool)
{
}

void computeBlocksKeep(System *, const char *, bool, bool)
{
}

void computeBlocks(System *, const char *, double *const *, bool)
{
}

void resetShells(System *, const char *)
{
}

void keepShells(System *, const char *, bool)
{
}

void relaxBlocks(System *, const char *, const char *, const RelaxSettings &, size_t *nComputed, bool *notFinite,
                 size_t *nFailed)
{
    *nComputed = 0;
    *notFinite = false;
    *nFailed   = 0;
}

void copyCoordinates(System *, size_t, size_t, const uint32_t *, size_t, const uint32_t *)
{
}

void saveCoordinates(System *, size_t)
{
}

void restoreCoordinates(System *, size_t)
{
}

void downloadBlock(System *, size_t, int64_t *, float *)
{
}

void createDeviations(System *, const DeviationSystem &)
{
}

void keepEnergies(System *, const char *, bool)
{
}

void computeDeviations(System *, const char *, double *, double *)
{
}

bool computePairBlocks(const PairBatchBlock *, size_t, bool)
{
    return false;
}

} // namespace gpu

} // namespace alexandria
