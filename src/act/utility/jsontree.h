/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022,2025,2026
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
 * Implements JSON utility using nlohmann/json for serialization
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef JSONTREE_H
#define JSONTREE_H

#include <cstdlib>
#include <string>
#include <vector>

namespace alexandria
{
/*! \brief Text for a number that reads back as the same number
 * The shortest text is chosen, such that 0.1 is written as 0.1.
 * \param[in] value The number
 * \return the text, nan or inf for numbers that are not finite
 */
std::string jsonNumberText(double value);

/*! \brief Text for a single precision number that reads back as the same number
 * A float is not converted to double first, since that would add digits
 * that carry no information, e.g. 0.10000000149 for 0.1.
 * \param[in] value The number
 * \return the text, nan or inf for numbers that are not finite
 */
std::string jsonNumberText(float value);

/*! \brief Simple class to produce JSON or text outputs
 * The input is a tree structure that is filled in a program.
 * A value is a string, a number or a boolean. Numbers are stored as text as
 * well, such that the program decides on the number of digits, but they are
 * written as JSON numbers, without quotes. Numbers that are not finite
 * are written as null, since JSON has no way to write them.
 * Each JsonTree can have either a value or objects but not both.
 * If both are added a warning is printed. A tree without value and without
 * objects is written as null, a tree with an empty string as value as "".
 * The objects of a tree are written as the members of a JSON object, or, for
 * a tree made with array(), as the elements of a JSON array. In an array the
 * keys of the objects are not written.
 * JSON serialization is performed using the nlohmann/json library.
 */
class JsonTree
{
private:
    //! The key
    std::string           key_;
    //! The value
    std::string           value_;
    //! Whether a value was given, which may be an empty string
    bool                  hasValue_ = false;
    //! Whether the value is a number
    bool                  number_  = false;
    //! Whether the value is a boolean
    bool                  boolean_ = false;
    //! Whether the objects are the elements of an array
    bool                  array_   = false;
    //! Child objects
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
        key_      = key;
        value_    = value;
        hasValue_ = true;
    }

    /*! \brief Make an element that holds a number given as text
     * Use this when the program decides on the number of digits. Text
     * that is not a number is written as a string.
     * \param[in] key  The key
     * \param[in] text The number as text, spaces around it are ignored
     * \return the element
     */
    static JsonTree number(const std::string &key, const std::string &text)
    {
        JsonTree jt(key, text);
        jt.number_ = true;
        return jt;
    }

    /*! \brief Make an element that holds a boolean
     * \param[in] key   The key
     * \param[in] value The value, written as true or false
     * \return the element
     */
    static JsonTree boolean(const std::string &key, bool value)
    {
        JsonTree jt(key, value ? "true" : "false");
        jt.boolean_ = true;
        return jt;
    }

    /*! \brief Make an element of which the objects are written as an array
     * \param[in] key The key
     * \return the element, to which the elements of the array can be added
     * with addObject, addNumber, addBool or addArray. Their keys are not
     * written in JSON.
     */
    static JsonTree array(const std::string &key)
    {
        JsonTree jt(key);
        jt.array_ = true;
        return jt;
    }

    /*! \brief Constructor that sets the key and a numerical value
     * \param[in] key   The key
     * \param[in] value The value, written as a JSON number
     */
    JsonTree(const std::string &key, double value)
    {
        key_    = key;
        value_  = jsonNumberText(value);
        number_ = true;
    }

    /*! \brief Constructor that sets the key and a single precision value
     * \param[in] key   The key
     * \param[in] value The value, written as a JSON number
     */
    JsonTree(const std::string &key, float value)
    {
        key_    = key;
        value_  = jsonNumberText(value);
        number_ = true;
    }

    /*! \brief Constructor that sets the key and an integer value
     * \param[in] key   The key
     * \param[in] value The value, written as a JSON number
     */
    JsonTree(const std::string &key, int value)
    {
        key_    = key;
        value_  = std::to_string(value);
        number_ = true;
    }

    //! \return key
    const std::string &key() const { return key_; }

    //! \return value
    const std::string &value() const { return value_; }

    //! \return whether a value was given, which may be an empty string
    bool hasValue() const { return hasValue_; }

    //! \return whether the value is a number
    bool isNumber() const { return number_; }

    //! \return whether the value is a boolean
    bool isBool() const { return boolean_; }

    //! \return whether the objects are the elements of an array
    bool isArray() const { return array_; }

    //! \return objects
    const std::vector<JsonTree> &objects() const { return objects_; }
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
        objects_.emplace_back(key, value);
    }

    /*! \brief Add a tree element with a number
     * \param[in] key   The key
     * \param[in] value The value, written as a JSON number
     */
    void addObject(const std::string &key, double value)
    {
        objects_.emplace_back(key, value);
    }

    /*! \brief Add a tree element with a single precision number
     * \param[in] key   The key
     * \param[in] value The value, written as a JSON number
     */
    void addObject(const std::string &key, float value)
    {
        objects_.emplace_back(key, value);
    }

    /*! \brief Add a tree element with an integer
     * \param[in] key   The key
     * \param[in] value The value, written as a JSON number
     */
    void addObject(const std::string &key, int value)
    {
        objects_.emplace_back(key, value);
    }

    /*! \brief Add a tree element with a number that is given as text
     * Use this when the program decides on the number of digits. Text
     * that is not a number is written as a string.
     * \param[in] key  The key
     * \param[in] text The number as text, spaces around it are ignored
     */
    void addNumber(const std::string &key, const std::string &text)
    {
        objects_.push_back(number(key, text));
    }

    /*! \brief Add a tree element with a boolean
     * \param[in] key   The key
     * \param[in] value The value, written as true or false
     */
    void addBool(const std::string &key, bool value)
    {
        objects_.push_back(boolean(key, value));
    }

    /*! \brief Add an array of numbers, integers or strings
     * \tparam T double, float, int or std::string
     * \param[in] key    The key
     * \param[in] values The values, an empty vector gives an empty array
     */
    template<typename T>
    void addArray(const std::string &key, const std::vector<T> &values)
    {
        JsonTree jt = array(key);
        for (const auto &value : values)
        {
            jt.objects_.emplace_back("", value);
        }
        objects_.push_back(jt);
    }

    /*! \brief Add a tree element with a value and a unit
     * \param[in] key   The key
     * \param[in] value The value as text. If it is a number it is written
     *                  as a JSON number.
     * \param[in] unit  The unit
     */
    void addValueUnit(const std::string &key, const std::string &value, const std::string &unit);

    /*! \brief Add a tree element with a value and a unit
     * \param[in] key   The key
     * \param[in] value The value, written as a JSON number
     * \param[in] unit  The unit
     */
    void addValueUnit(const std::string &key, double value, const std::string &unit)
    {
        addValueUnit(key, jsonNumberText(value), unit);
    }

    /*! \brief Add a tree element with a value and a unit
     * \param[in] key   The key
     * \param[in] value The value, written as a JSON number
     * \param[in] unit  The unit
     */
    void addValueUnit(const std::string &key, float value, const std::string &unit)
    {
        addValueUnit(key, jsonNumberText(value), unit);
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
               bool               json,
               int                indent=2);
    
};

} // namespace alexandria

#endif
