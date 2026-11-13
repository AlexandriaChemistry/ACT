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
 * A small reader of JSON (RFC 8259) in UTF-8, for data files such as
 * alexandria.json, see src/act/python/mol_json.py. JsonTree writes JSON.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_UTILITY_JSONREADER_H
#define ACT_UTILITY_JSONREADER_H

#include <map>
#include <string>
#include <vector>

namespace alexandria
{

//! A value of a JSON document
class JsonValue
{
public:
    //! The kinds of values
    enum class Type
    {
        Null, Bool, Number, String, Array, Object
    };

    //! A null value
    JsonValue() = default;

    //! \return the kind of value
    Type type() const { return type_; }

    //! \return whether it is null
    bool isNull() const { return type_ == Type::Null; }

    //! \return the boolean, throws if it is not one
    bool asBool() const;

    //! \return the number, throws if it is not one
    double asNumber() const;

    //! \return the number as an integer, throws if it is not an integer
    long asInt() const;

    //! \return the string, throws if it is not one
    const std::string &asString() const;

    //! \return the elements of an array, throws if it is not one
    const std::vector<JsonValue> &asArray() const;

    //! \return the members of an object in their order, throws if it is not one
    const std::vector<std::pair<std::string, JsonValue>> &asObject() const;

    /*! \brief A member of an object
     * \param[in] key The name of the member
     * \return The value, nullptr if the object has no such member,
     *         throws if this is not an object
     */
    const JsonValue *find(const std::string &key) const;

    //! Make values, see parseJson
    static JsonValue makeBool(bool b);
    static JsonValue makeNumber(double d);
    static JsonValue makeString(std::string s);
    static JsonValue makeArray(std::vector<JsonValue> a);
    static JsonValue makeObject(std::vector<std::pair<std::string, JsonValue>> o);

private:
    Type                                           type_   = Type::Null;
    bool                                           bool_   = false;
    double                                         number_ = 0;
    std::string                                    string_;
    std::vector<JsonValue>                         array_;
    std::vector<std::pair<std::string, JsonValue>> object_;
};

/*! \brief Parse a JSON document
 * Strings are kept in UTF-8, escapes such as \\u00b1 are converted to it.
 * \param[in] text The document
 * \return The value of the document
 * \throws gmx::InvalidInputError with the position of an error
 */
JsonValue parseJson(const std::string &text);

/*! \brief Read and parse a JSON file in UTF-8
 * \param[in] fileName The file
 * \return The value of the document
 * \throws gmx::FileIOError if it cannot be read, gmx::InvalidInputError if it is not JSON
 */
JsonValue readJsonFile(const std::string &fileName);

} // namespace alexandria

#endif
