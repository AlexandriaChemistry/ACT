/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022,2025
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
 * Implements JSON utility
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef JSONTREE_H
#define JSONTREE_H

#include <cstdlib>

#include <string>
#include <vector>

namespace alexandria
{
/*! \brief Simple class to produce JSON or text outputs
 * The input is a tree structure that is filled in a program.
 * Only strings allowed for now.
 */
class JsonTree
{
private:
    // The key
    std::string           key_;
    // The value
    std::string           value_;
    // More objects;
    std::vector<JsonTree> objects_;
    
    
public:
    /*! \brief Constructor that sets the key for this one
     * \param[in] key   The key
     */
    JsonTree(const std::string &key)
    {
        key_ = key;
    }
    
    /*! \brief Constructor that sets the key and value for this one
     * \param[in] key   The key
     * \param[in] value The value
     */
    JsonTree(const std::string &key, const std::string &value)
    {
        key_ = key;
        value_ = value;
    }
    
    //! \return if there is nothing here
    bool empty() const { return objects_.empty() && value_.empty(); }
    
    //! \return whether there are complex trees
    bool branched() const { return !objects_.empty(); }

    /*! \brief Write everything in a single string
     * \param[in]    json   Whether or not to write JSON format
     * \param[inout] indent Indentation for printing
     * \return One big string
     */
    const std::string writeString(bool  json,
                                  int  *indent) const;
    /*! \brief Add a tree element
     * \param[in] object The object to add
     */
    void addObject(const JsonTree object)
    {
        objects_.push_back(object);
    }
    
    /*! \brief Add a tree element
     * \param[in] key   The key
     * \param[in] value The value
     */
    void addObject(const std::string &key, const std::string &value)
    {
        objects_.push_back(JsonTree(key, value));
    }
    
    /*! \brief Add a tree element with a subelement
     * \param[in] key   The key
     * \param[in] value The value
     * \param[in] unit  T
     */
    void addValueUnit(const std::string &key, const std::string &value, const std::string &unit)
    {
        JsonTree jt(key);
        jt.addObject("value", value);
        jt.addObject("unit", unit);
        objects_.push_back(jt);
    }
    
    /*! \brief Add a tree element
     * \param[in] key    The key
     * \param[in] object A complex object
     */
    void addObject(const std::string &key, const JsonTree object)
    {
        JsonTree jnew(key);
        jnew.addObject(object);
        objects_.push_back(jnew);
    }
    
    /*! \brief Write the tree to a file
     * \param[in] fileName The name of the file to write
     * \param[in] json     Whether or not to use JSON format
     */
    void write(const std::string &fileName,
               bool               json);
    
};

} // namespace alexandria

#endif
