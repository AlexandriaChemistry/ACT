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
 * Implements the JSON reader, see jsonreader.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "jsonreader.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

namespace
{

//! The names of the types, for messages
const char *typeName(JsonValue::Type t)
{
    switch (t)
    {
    case JsonValue::Type::Null:   return "null";
    case JsonValue::Type::Bool:   return "a boolean";
    case JsonValue::Type::Number: return "a number";
    case JsonValue::Type::String: return "a string";
    case JsonValue::Type::Array:  return "an array";
    case JsonValue::Type::Object: return "an object";
    }
    return "unknown";
}

//! Throw if a value is not of a type
void expectType(const JsonValue &v, JsonValue::Type t)
{
    if (v.type() != t)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("JSON value is %s, not %s",
                                                           typeName(v.type()), typeName(t))));
    }
}

//! A recursive descent parser of RFC 8259
class Parser
{
public:
    //! \param[in] text The document
    explicit Parser(const std::string &text) : text_(text) {}

    //! \return The value of the whole document
    JsonValue document()
    {
        skipSpace();
        auto v = value(0);
        skipSpace();
        if (pos_ != text_.size())
        {
            fail("text after the end of the document");
        }
        return v;
    }

private:
    const std::string &text_;
    size_t             pos_ = 0;
    //! Deeper documents are refused, against a stack overflow
    static constexpr int maxDepth = 512;

    [[noreturn]] void fail(const std::string &what) const
    {
        // The line and column, for the message
        size_t line = 1, col = 1;
        for (size_t i = 0; i < pos_ && i < text_.size(); i++)
        {
            if (text_[i] == '\n')
            {
                line += 1;
                col   = 1;
            }
            else
            {
                col += 1;
            }
        }
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Invalid JSON at line %zu, column %zu: %s",
                                                           line, col, what.c_str())));
    }

    void skipSpace()
    {
        while (pos_ < text_.size() &&
               (text_[pos_] == ' ' || text_[pos_] == '\t' || text_[pos_] == '\n' || text_[pos_] == '\r'))
        {
            pos_++;
        }
    }

    bool consume(const char *word)
    {
        size_t n = std::char_traits<char>::length(word);
        if (text_.compare(pos_, n, word) == 0)
        {
            pos_ += n;
            return true;
        }
        return false;
    }

    JsonValue value(int depth)
    {
        if (depth > maxDepth)
        {
            fail("nested too deeply");
        }
        if (pos_ >= text_.size())
        {
            fail("unexpected end");
        }
        char c = text_[pos_];
        if (c == '{')
        {
            return object(depth);
        }
        if (c == '[')
        {
            return array(depth);
        }
        if (c == '"')
        {
            return JsonValue::makeString(string());
        }
        if (c == '-' || (c >= '0' && c <= '9'))
        {
            return JsonValue::makeNumber(number());
        }
        if (consume("true"))
        {
            return JsonValue::makeBool(true);
        }
        if (consume("false"))
        {
            return JsonValue::makeBool(false);
        }
        if (consume("null"))
        {
            return JsonValue();
        }
        fail(gmx::formatString("unexpected character '%c'", c));
    }

    JsonValue object(int depth)
    {
        pos_++; // {
        std::vector<std::pair<std::string, JsonValue>> members;
        skipSpace();
        if (consume("}"))
        {
            return JsonValue::makeObject(std::move(members));
        }
        while (true)
        {
            skipSpace();
            if (pos_ >= text_.size() || text_[pos_] != '"')
            {
                fail("expected the name of a member");
            }
            auto key = string();
            skipSpace();
            if (!consume(":"))
            {
                fail("expected ':'");
            }
            skipSpace();
            auto v = value(depth + 1);
            members.emplace_back(std::move(key), std::move(v));
            skipSpace();
            if (consume(","))
            {
                continue;
            }
            if (consume("}"))
            {
                return JsonValue::makeObject(std::move(members));
            }
            fail("expected ',' or '}'");
        }
    }

    JsonValue array(int depth)
    {
        pos_++; // [
        std::vector<JsonValue> elements;
        skipSpace();
        if (consume("]"))
        {
            return JsonValue::makeArray(std::move(elements));
        }
        while (true)
        {
            skipSpace();
            elements.push_back(value(depth + 1));
            skipSpace();
            if (consume(","))
            {
                continue;
            }
            if (consume("]"))
            {
                return JsonValue::makeArray(std::move(elements));
            }
            fail("expected ',' or ']'");
        }
    }

    //! Four hexadecimal digits of \\u
    unsigned hex4()
    {
        if (pos_ + 4 > text_.size())
        {
            fail("incomplete \\u escape");
        }
        unsigned v = 0;
        for (int i = 0; i < 4; i++)
        {
            char     c = text_[pos_++];
            unsigned d = 0;
            if (c >= '0' && c <= '9')
            {
                d = c - '0';
            }
            else if (c >= 'a' && c <= 'f')
            {
                d = c - 'a' + 10;
            }
            else if (c >= 'A' && c <= 'F')
            {
                d = c - 'A' + 10;
            }
            else
            {
                fail("invalid \\u escape");
            }
            v = 16*v + d;
        }
        return v;
    }

    //! Append a code point as UTF-8
    static void appendUtf8(unsigned cp, std::string *s)
    {
        if (cp < 0x80)
        {
            s->push_back(static_cast<char>(cp));
        }
        else if (cp < 0x800)
        {
            s->push_back(static_cast<char>(0xC0 | (cp >> 6)));
            s->push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        }
        else if (cp < 0x10000)
        {
            s->push_back(static_cast<char>(0xE0 | (cp >> 12)));
            s->push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
            s->push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        }
        else
        {
            s->push_back(static_cast<char>(0xF0 | (cp >> 18)));
            s->push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
            s->push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
            s->push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        }
    }

    std::string string()
    {
        pos_++; // "
        std::string s;
        while (true)
        {
            if (pos_ >= text_.size())
            {
                fail("unterminated string");
            }
            char c = text_[pos_++];
            if (c == '"')
            {
                return s;
            }
            if (static_cast<unsigned char>(c) < 0x20)
            {
                fail("control character in a string");
            }
            if (c != '\\')
            {
                // UTF-8 bytes are kept as they are
                s.push_back(c);
                continue;
            }
            if (pos_ >= text_.size())
            {
                fail("unterminated escape");
            }
            char e = text_[pos_++];
            switch (e)
            {
            case '"':  s.push_back('"');  break;
            case '\\': s.push_back('\\'); break;
            case '/':  s.push_back('/');  break;
            case 'b':  s.push_back('\b'); break;
            case 'f':  s.push_back('\f'); break;
            case 'n':  s.push_back('\n'); break;
            case 'r':  s.push_back('\r'); break;
            case 't':  s.push_back('\t'); break;
            case 'u':
                {
                    unsigned cp = hex4();
                    // A pair of surrogates for code points above 0xFFFF
                    if (cp >= 0xD800 && cp <= 0xDBFF)
                    {
                        if (!consume("\\u"))
                        {
                            fail("missing low surrogate");
                        }
                        unsigned lo = hex4();
                        if (lo < 0xDC00 || lo > 0xDFFF)
                        {
                            fail("invalid low surrogate");
                        }
                        cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                    }
                    else if (cp >= 0xDC00 && cp <= 0xDFFF)
                    {
                        fail("unexpected low surrogate");
                    }
                    appendUtf8(cp, &s);
                }
                break;
            default:
                fail(gmx::formatString("invalid escape '\\%c'", e));
            }
        }
    }

    double number()
    {
        size_t start = pos_;
        const auto digits = [this]()
        {
            size_t n = 0;
            while (pos_ < text_.size() && text_[pos_] >= '0' && text_[pos_] <= '9')
            {
                pos_++;
                n++;
            }
            return n;
        };
        if (text_[pos_] == '-')
        {
            pos_++;
        }
        if (pos_ < text_.size() && text_[pos_] == '0')
        {
            pos_++;
        }
        else if (digits() == 0)
        {
            fail("invalid number");
        }
        if (pos_ < text_.size() && text_[pos_] == '.')
        {
            pos_++;
            if (digits() == 0)
            {
                fail("invalid number");
            }
        }
        if (pos_ < text_.size() && (text_[pos_] == 'e' || text_[pos_] == 'E'))
        {
            pos_++;
            if (pos_ < text_.size() && (text_[pos_] == '+' || text_[pos_] == '-'))
            {
                pos_++;
            }
            if (digits() == 0)
            {
                fail("invalid number");
            }
        }
        // strtod on the exact characters, independent of the locale for these
        const std::string num = text_.substr(start, pos_ - start);
        std::istringstream iss(num);
        iss.imbue(std::locale::classic());
        double d = 0;
        iss >> d;
        return d;
    }
};

} // namespace

bool JsonValue::asBool() const
{
    expectType(*this, Type::Bool);
    return bool_;
}

double JsonValue::asNumber() const
{
    expectType(*this, Type::Number);
    return number_;
}

long JsonValue::asInt() const
{
    expectType(*this, Type::Number);
    if (number_ != std::floor(number_))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("JSON number %g is not an integer", number_)));
    }
    return static_cast<long>(number_);
}

const std::string &JsonValue::asString() const
{
    expectType(*this, Type::String);
    return string_;
}

const std::vector<JsonValue> &JsonValue::asArray() const
{
    expectType(*this, Type::Array);
    return array_;
}

const std::vector<std::pair<std::string, JsonValue>> &JsonValue::asObject() const
{
    expectType(*this, Type::Object);
    return object_;
}

const JsonValue *JsonValue::find(const std::string &key) const
{
    for (const auto &m : asObject())
    {
        if (m.first == key)
        {
            return &m.second;
        }
    }
    return nullptr;
}

JsonValue JsonValue::makeBool(bool b)
{
    JsonValue v;
    v.type_ = Type::Bool;
    v.bool_ = b;
    return v;
}

JsonValue JsonValue::makeNumber(double d)
{
    JsonValue v;
    v.type_   = Type::Number;
    v.number_ = d;
    return v;
}

JsonValue JsonValue::makeString(std::string s)
{
    JsonValue v;
    v.type_   = Type::String;
    v.string_ = std::move(s);
    return v;
}

JsonValue JsonValue::makeArray(std::vector<JsonValue> a)
{
    JsonValue v;
    v.type_  = Type::Array;
    v.array_ = std::move(a);
    return v;
}

JsonValue JsonValue::makeObject(std::vector<std::pair<std::string, JsonValue>> o)
{
    JsonValue v;
    v.type_   = Type::Object;
    v.object_ = std::move(o);
    return v;
}

JsonValue parseJson(const std::string &text)
{
    Parser p(text);
    return p.document();
}

JsonValue readJsonFile(const std::string &fileName)
{
    const std::ifstream in(fileName, std::ios::binary);
    if (!in)
    {
        GMX_THROW(gmx::FileIOError(gmx::formatString("Cannot read JSON file %s", fileName.c_str())));
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    std::string text = ss.str();
    // A byte order mark is allowed, and ignored
    if (text.compare(0, 3, "\xEF\xBB\xBF") == 0)
    {
        text.erase(0, 3);
    }
    try
    {
        return parseJson(text);
    }
    catch (gmx::InvalidInputError &e)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("%s in %s", e.what(), fileName.c_str())));
    }
}

} // namespace alexandria
