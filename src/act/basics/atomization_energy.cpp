/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "atomization_energy.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "act/utility/stringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
    
namespace alexandria
{

class AtomizationEnergyTerm
{
 private:
    //! Element / atom
    std::string elem_;
    //! Charge of the particle
    int         charge_;
    //! Data source
    std::string source_;
    //! The property
    std::string prop_;
    //! Temperature
    double      T_;
    //! The value of the property
    double      value_;
    //! The spin multiplicity
    int         mult_;
    //! The unit
    std::string unit_;
    //! The reference
    std::string ref_;
 public:
    //! Constructor, see descriptions above
    AtomizationEnergyTerm(const std::string &elem,
                          int                charge,
                          const std::string &source,
                          const std::string &prop,
                          double             T,
                          double             value,
                          int                mult,
                          const std::string &unit,
                          const std::string &ref) :
    elem_(elem), charge_(charge), source_(source), prop_(prop),
        T_(T), value_(value), mult_(mult), unit_(unit), ref_(ref) {}
    
    const std::string &elem() const { return elem_; }

    int charge() const { return charge_; }

    const std::string &source() const { return source_; }

    const std::string prop() const { return prop_; }

    double T() const { return T_; }

    double value() const { return value_; }

    int mult() const { return mult_; }

    const std::string unit() const { return unit_; }

    const std::string ref() const { return ref_; }
};

AtomizationEnergy::AtomizationEnergy()
{
    const char *acd     = "ACTDATA";
    auto        actdata = std::getenv(acd);
    if (nullptr == actdata || strlen(actdata) == 0)
    {
        fprintf(stderr, "Please set the environment variable %s to the correct directory.\n", acd);
        return;
    }
    std::vector<const char *> ahof = {
        "atomization-energies.csv", "atomization-energies-dft.csv"
    };
    for (const auto &aa : ahof)
    {
        std::string infile = gmx::formatString("%s/%s", actdata, aa);
        try
        {
            gmx::TextReader tr(infile);
            std::string     line;
            while (tr.readLine(&line))
            {
                if (line.find("#") == std::string::npos)
                {
                    auto words = gmx::splitDelimitedString(line, '|');
                    if (words.size() >= 8)
                    {
                        int    charge      = my_atoi(words[1].c_str(), "charge");
                        double temperature = my_atof(words[4].c_str(), "temperature");
                        double value       = my_atof(words[5].c_str(), "value");
                        int    mult        = my_atoi(words[1].c_str(), "multiplicity");
                        terms_.push_back(new AtomizationEnergyTerm(words[0], charge, words[2], words[3], temperature,
                                                                   value, mult, words[7], words[8]));
                    }
                    else if (debug)
                    {
                        fprintf(debug, "Funny line '%s' in '%s'\n", line.c_str(), infile.c_str());
                    }
                }
            }
            tr.close();
        }
        catch (gmx::FileIOError)
        {
            fprintf(stderr, "Cannot find %s\n", infile.c_str());
        }
    }
    if (debug)
    {
        fprintf(debug, "Read %zu terms for atomization energy\n", terms_.size());
    }
}

AtomizationEnergy::~AtomizationEnergy()
{
    for(auto &t : terms_)
    {
        delete t;
    }
}
        
double AtomizationEnergy::term(const std::string &elem,
                               int                charge,
                               const std::string &source,
                               const std::string &prop,
                               double             T,
                               std::string       *unit,
                               std::string       *ref) const
{
    //    AtomizationEnergyTerm aet(elem, charge, source, prop, T, 0, 1, "", "");
    auto term = std::find_if(terms_.begin(), terms_.end(),
                             [&](const AtomizationEnergyTerm *x)
                             { return (elem   == x->elem() &&
                                       charge == x->charge() &&
                                       source == x->source() &&
                                       prop   == x->prop() &&
                                       std::abs(T - x->T()) < 1e-3); });
    if (terms_.end() == term)
    {
        return 0;
    }
    else
    {
        if (nullptr != unit)
        {
            *unit = (*term)->unit();
        }
        if (nullptr != ref)
        {
            *ref = (*term)->ref();
        }
        return (*term)->value();
    }
}

}
