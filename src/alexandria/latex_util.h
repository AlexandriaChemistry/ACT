/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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
 */
 
#ifndef LATEX_UTIL_H
#define LATEX_UTIL_H

//#include <math.h>
#include <stdio.h>
//#include <stdlib.h>
#include <string.h>

#include <string>
#include <vector>
//#include "categories.h"
//#include "molprop.h"
//#include "molprop_util.h"


namespace alexandria
{
/*! \brief Class for making LaTeX tables for publications
 */
class LongTable
{
private:
    //! File pointer to write to
    FILE                    *fp_ = nullptr;
    //! Whether or not the file was opened by us
    bool                     openedFile_ = false;
    //! Font type to use
    const char              *font_ = nullptr;
    //! Caption text
    std::string              caption_;
    //! Designator for the amount and type of columns
    std::string              columns_;
    //! Label for referincing the table
    std::string              label_;
    //! Header text
    std::vector<std::string> headLines_;
    //! Whether or not to print in landscape orientation
    bool                     bLandscape_ = false;
public:
    /*!\brief Constructor with a file pointer
     * \param[in] fp         The file pointer to write to
     * \param[in] bLandscape Whether or not to print in landscape orientation
     * \param[in] font       The font to use
     */
    LongTable(FILE *fp, bool bLandscape, const char *font);
    
    /*! Constructor with a file name
     * \param[in] fn         The file name to write to
     * \param[in] bLandscape Whether or not to print in landscape orientation
     */
    LongTable(const char *fn, bool bLandscape);
    
    /*! \brief Destructor
     * Will close the file pointer if it is open and is was
     * this structure that opened it. If not, the control is
     * left to the calling functions.
     */
    ~LongTable();

    /*! \brief Set the caption of the table
     * \param[in] caption The text for the caption
     */
    void setCaption(const char *caption) { caption_.assign(caption); }

    /*! \brief Set the label of the table
     * \param[in] label The string for the label
     */
    void setLabel(const char *label) { label_.assign(label); }

    /*! \brief Generate columns entry 
     * with first column left aligned and other center
     * \param[in] nColumns The number of columns
     */
    void setColumns(int nColumns);

    /*! \brief Generate columns entry 
     * according to the string as expected by LaTeX, e.g. 'lccrc'
     * \param[in] columns The column string
     */
    void setColumns(const char *columns) { columns_.assign(columns); }

    /*! \brief Add a header line
     * The header can consist of multiple lines, and this allows to
     * add one at a time.
     * \param[in] headline Another line
     */
    void addHeadLine(const std::string &headline) { headLines_.push_back(headline); }

    //! \brief Print the header text to the file
    void printHeader();
    
    //! \brief Print the footer text to the file
    void printFooter();
    
    /*! \brief Print a line to the table
     * \param[in] line The line to print, include LaTeX column markers
     */
    void printLine(const std::string &line);

    /*! \brief Print a line to the table based on columns
     * The number of columns should be the same as what was defined
     * using the setColumns function.
     * \param[in] columns The column entries
     */
    void printColumns(const std::vector<std::string> &columns);
    
    //! \brief Print a horizontal line to the table
    void printHLine();
};

}// namespace

#endif
