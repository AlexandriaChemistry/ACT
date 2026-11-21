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
/*! \internal \file
 * \brief
 * An output environment for writing xvg files, for the parts of the tools
 * that a program without a command line uses.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_UTILITY_OUTPUT_ENVIRONMENT_H
#define ACT_UTILITY_OUTPUT_ENVIRONMENT_H

#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/programcontext.h"

namespace alexandria
{

/*! \brief The output environment of a command line, or one with its defaults
 *
 * The functions that write xvg files need the output environment that the
 * parser of the command line makes. A program without a command line gets
 * one with the defaults of the command line (-xvg xmgrace), which is
 * freed with this object.
 */
class OutputEnvironment
{
public:
    /*! \brief Use the given environment, or make one with the defaults
     * \param[in] given The environment of a command line, may be nullptr
     */
    explicit OutputEnvironment(gmx_output_env_t *given) : oenv_(given)
    {
        if (nullptr == oenv_)
        {
            output_env_init(&own_, gmx::getProgramContext(), time_ps, FALSE, exvgXMGRACE, 0);
            oenv_ = own_;
        }
    }
    //! Frees the environment that was made, if any
    ~OutputEnvironment()
    {
        if (own_)
        {
            output_env_done(own_);
        }
    }
    //! Not copied, since it may own the environment
    OutputEnvironment(const OutputEnvironment &)            = delete;
    //! Not copied, since it may own the environment
    OutputEnvironment &operator=(const OutputEnvironment &) = delete;
    //! \return the environment
    gmx_output_env_t *get() const { return oenv_; }
private:
    //! The environment that is used
    gmx_output_env_t *oenv_ = nullptr;
    //! The environment that was made here, if any
    gmx_output_env_t *own_  = nullptr;
};

} // namespace alexandria

#endif
