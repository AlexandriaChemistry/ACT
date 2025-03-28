/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef MOLPROP_UTIL_H
#define MOLPROP_UTIL_H

#include <string>
#include <vector>

#include "act/alexandria/molselect.h"
#include "act/molprop/molprop.h"
#include "gromacs/utility/arrayref.h"

/*! \brief
 * Enumerated type for MolPropSort function
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum MolPropSortAlgorithm {
    MPSA_MOLNAME     = 0,
    MPSA_FORMULA     = 1,
    MPSA_COMPOSITION = 2,
    MPSA_SELECTION   = 3,
    MPSA_NR          = 4
};

namespace alexandria
{
class QmCalc
{
    private:
        std::string method_, basis_, type_, lot_;
        int         count_;
    public:
        QmCalc(const std::string &method,
               const std::string &basis,
               const std::string &type) :
            method_(method), basis_(basis), type_(type), lot_(method), count_(1)
        {
            lot_.append("/");
            lot_.append(basis);
        };

        const std::string &method() const { return method_; }
        const std::string &basis() const { return basis_; }
        const std::string &type() const { return type_; }
        const std::string &lot() const { return lot_; }
        int count() const { return count_; }
        void increment() { count_++; }
        void decrement()
        {
            GMX_RELEASE_ASSERT(count_ > 0, "Trying to reduce count below zero");
            count_--;
        }
};

class QmCount
{
    private:
        std::vector<QmCalc>      qmc_;
        std::vector<std::string> conf_;

        std::vector<QmCalc>::const_iterator findCalc(const std::string &method,
                                                     const std::string &basis,
                                                     const std::string &type) const
        {
            return std::find_if(qmc_.begin(), qmc_.end(),
                                [method, basis, type](const QmCalc &qmc)
                                { return (qmc.method().compare(method) == 0 &&
                                          qmc.basis().compare(basis) == 0 &&
                                          qmc.type().compare(type) == 0); });
        }
        std::vector<QmCalc>::iterator findCalc(const std::string &method,
                                               const std::string &basis,
                                               const std::string &type)
        {
            return std::find_if(qmc_.begin(), qmc_.end(),
                                [method, basis, type](const QmCalc &qmc)
                                { return (qmc.method().compare(method) == 0 &&
                                          qmc.basis().compare(basis) == 0 &&
                                          qmc.type().compare(type) == 0); });
        }

    public:
        QmCount() {};

        std::vector<QmCalc>::iterator beginCalc() { return qmc_.begin(); }
        std::vector<QmCalc>::iterator endCalc() { return qmc_.end(); }
        std::vector<QmCalc>::const_iterator beginCalc() const { return qmc_.begin(); }
        std::vector<QmCalc>::const_iterator endCalc() const { return qmc_.end(); }
        size_t nCalc() const { return qmc_.size(); }
        int qmCalcCount(const std::string &method,
                        const std::string &basis,
                        const std::string &type) const;

        void addConf(const std::string &conformation);

        void addCalc(const std::string &method,
                     const std::string &basis,
                     const std::string &type);

};

/*! \brief Reset the molprop index number
 * \param[in] mp The molprops
 */
void generate_index(std::vector<MolProp> *mp);

/*! \brief Will merge multiple molprop files.
 * \param[in] infiles     MolProp filenames
 * \param[in] mp          The molprop vector
 * \return List of warning messages, if empty all is fine
 */
std::vector<std::string> merge_xml(const gmx::ArrayRef<const std::string> &infiles,
                                   std::vector<alexandria::MolProp>       *mp);

/*! \brief Check the available molprops to see what kind of calculations are stored.
 * \param[in]  mp       The molprops
 * \param[in]  mpo      The observable to check
 * \param[in]  fc_types Observable measurement types. Multiple strings may be
 *                      provided, separated by a colon (:). Ignored if nullptr.
 * \param[out] qmc      Overview of the calculations
 */
void find_calculations(const std::vector<alexandria::MolProp> &mp,
                       MolPropObservable                       mpo,
                       const char                             *fc_types,
                       QmCount                                *qmc);

/*! \brief
 * Sorts a vector of molprops
 *
 * Function that uses the std::sort routine and can apply different sorting
 * keys.
 *
 * \param[inout]  mp        The vector of MolProp
 * \param[in]     mpsa      The algorithm used for sorting
 * \param[in]     apt       Database of atom properties
 * \param[in]     gms       Optional structure containing selection criteria
 * \ingroup module_alexandria
 */
void MolPropSort(std::vector<MolProp> *mp,
                 MolPropSortAlgorithm  mpsa,
                 gmx_atomprop_t        apt,
                 const MolSelect      &gms);

} // namespace alexandria

/*! \brief Utility to split a user-provided lot
 *
 * \param[in]  lot    Level of theory
 * \param[out] method QM method
 * \param[out] basis  QM basis set
 */
void splitLot(const char  *lot,
              std::string *method,
              std::string *basis);

/*! \brief Utility to generate a lot
 *
 * \param[out] method QM method
 * \param[out] basis  QM basis set
 * \return  Level of theory
 */
std::string makeLot(const std::string &method,
                    const std::string &basis);
#endif
