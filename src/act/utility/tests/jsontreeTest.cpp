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
 * Tests for the JsonTree utility class.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <iterator>
#include <limits>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "../jsontree.h"
#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

// ---- Accessor and state tests ----

TEST(JsonTreeTest, KeyOnlyConstructorIsEmpty)
{
    const JsonTree jt("myKey");
    EXPECT_EQ("myKey", jt.key());
    EXPECT_TRUE(jt.value().empty());
    EXPECT_TRUE(jt.objects().empty());
    EXPECT_TRUE(jt.empty());
    EXPECT_FALSE(jt.branched());
}

TEST(JsonTreeTest, StringValueConstructorAccessors)
{
    const JsonTree jt("molecule", std::string("water"));
    EXPECT_EQ("molecule", jt.key());
    EXPECT_EQ("water", jt.value());
    EXPECT_FALSE(jt.empty());
    EXPECT_FALSE(jt.branched());
    EXPECT_TRUE(jt.objects().empty());
}

TEST(JsonTreeTest, DoubleValueConstructorNotEmpty)
{
    const JsonTree jt("pi", 3.14);
    EXPECT_EQ("pi", jt.key());
    EXPECT_FALSE(jt.value().empty());
    EXPECT_FALSE(jt.empty());
    EXPECT_FALSE(jt.branched());
}

TEST(JsonTreeTest, IntValueConstructorStoresAsString)
{
    const JsonTree jt("count", 42);
    EXPECT_EQ("count", jt.key());
    EXPECT_EQ("42", jt.value());
    EXPECT_FALSE(jt.empty());
    EXPECT_FALSE(jt.branched());
}

// ---- addObject tests ----

TEST(JsonTreeTest, AddObjectStringMakesBranched)
{
    JsonTree jt("parent");
    EXPECT_FALSE(jt.branched());
    EXPECT_TRUE(jt.empty());

    jt.addObject("child", std::string("value1"));

    EXPECT_TRUE(jt.branched());
    EXPECT_FALSE(jt.empty());
    ASSERT_EQ(1u, jt.objects().size());
    EXPECT_EQ("child", jt.objects()[0].key());
    EXPECT_EQ("value1", jt.objects()[0].value());
}

TEST(JsonTreeTest, AddObjectDoubleMakesBranched)
{
    JsonTree jt("parent");
    jt.addObject("temperature", 298.15);
    EXPECT_TRUE(jt.branched());
    ASSERT_EQ(1u, jt.objects().size());
    EXPECT_EQ("temperature", jt.objects()[0].key());
    EXPECT_FALSE(jt.objects()[0].value().empty());
}

TEST(JsonTreeTest, AddObjectJsonTreeMakesBranched)
{
    JsonTree jt("parent");
    const JsonTree child("child", std::string("val"));
    jt.addObject(child);
    EXPECT_TRUE(jt.branched());
    ASSERT_EQ(1u, jt.objects().size());
    EXPECT_EQ("child", jt.objects()[0].key());
}

TEST(JsonTreeTest, AddObjectWithKeyWrapsInNewNode)
{
    JsonTree jt("root");
    JsonTree child("child");
    child.addObject("x", std::string("1"));
    jt.addObject("nested", child);

    ASSERT_EQ(1u, jt.objects().size());
    // The wrapper node has key "nested" and contains child as its only object
    const auto &wrapper = jt.objects()[0];
    EXPECT_EQ("nested", wrapper.key());
    ASSERT_EQ(1u, wrapper.objects().size());
    EXPECT_EQ("child", wrapper.objects()[0].key());
}

TEST(JsonTreeTest, MultipleAddObjectsAccumulate)
{
    JsonTree jt("root");
    jt.addObject("a", std::string("1"));
    jt.addObject("b", std::string("2"));
    jt.addObject("c", std::string("3"));
    EXPECT_EQ(3u, jt.objects().size());
}

// ---- addValueUnit tests ----

TEST(JsonTreeTest, AddValueUnitCreatesSubstructure)
{
    JsonTree jt("result");
    jt.addValueUnit("energy", "-285.8", "kJ/mol");

    EXPECT_TRUE(jt.branched());
    ASSERT_EQ(1u, jt.objects().size());
    const auto &sub = jt.objects()[0];
    EXPECT_EQ("energy", sub.key());
    ASSERT_EQ(2u, sub.objects().size());
    EXPECT_EQ("value", sub.objects()[0].key());
    EXPECT_EQ("-285.8", sub.objects()[0].value());
    EXPECT_EQ("unit", sub.objects()[1].key());
    EXPECT_EQ("kJ/mol", sub.objects()[1].value());
}

// ---- writeString JSON tests ----

TEST(JsonTreeTest, WriteStringJsonSimpleValue)
{
    const JsonTree jt("name", std::string("water"));
    int indent = -1;
    const std::string result = jt.writeString(true, &indent);
    EXPECT_EQ(R"({"name":"water"})", result);
    EXPECT_EQ(-1, indent);
}

TEST(JsonTreeTest, WriteStringJsonIntValue)
{
    const JsonTree jt("count", 7);
    int indent = -1;
    const std::string result = jt.writeString(true, &indent);
    EXPECT_EQ(R"({"count":7})", result);
    EXPECT_EQ(-1, indent);
}

TEST(JsonTreeTest, WriteStringJsonNestedLeaves)
{
    JsonTree jt("molecule");
    jt.addObject("name", std::string("water"));
    jt.addObject("formula", std::string("H2O"));
    int indent = -1;
    const std::string result = jt.writeString(true, &indent);
    // Proper JSON: children become members of an object
    const std::string expected =
              "{\"molecule\":{\"name\":\"water\",\"formula\":\"H2O\"}}";
    EXPECT_EQ(expected, result);
    EXPECT_EQ(-1, indent);
}

TEST(JsonTreeTest, WriteStringJsonValueUnitStructure)
{
    JsonTree jt("data");
    jt.addValueUnit("energy", "42.0", "kJ/mol");
    int indent = 0;
    const std::string result = jt.writeString(true, &indent);
    EXPECT_NE(std::string::npos, result.find("\"energy\""));
    EXPECT_NE(std::string::npos, result.find("\"value\""));
    // The value is a number, not a string
    EXPECT_EQ(std::string::npos, result.find("\"42.0\""));
    EXPECT_NE(std::string::npos, result.find("42.0"));
    EXPECT_NE(std::string::npos, result.find("\"unit\""));
    EXPECT_NE(std::string::npos, result.find("\"kJ/mol\""));
    EXPECT_EQ(0, indent);
}

// ---- numbers in JSON ----

/*! \brief Write a tree as JSON and parse the result
 * \param[in] jt The tree
 * \return the content below the key of the tree
 */
nlohmann::json parsed(const JsonTree &jt)
{
    int indent = -1;
    return nlohmann::json::parse(jt.writeString(true, &indent)).at(jt.key());
}

TEST(JsonTreeTest, NumbersAreJsonNumbers)
{
    JsonTree jt("n");
    jt.addObject("double", 3.25);
    jt.addObject("float", 0.1F);
    jt.addObject("int", 42);
    jt.addNumber("formatted", "1.5000");
    jt.addObject("text", std::string("12"));
    const auto js = parsed(jt);

    ASSERT_TRUE(js.at("double").is_number_float());
    EXPECT_EQ(3.25, js.at("double").get<double>());
    // A float is written with the digits of a float, not with the noise
    // that appears when it is converted to double.
    ASSERT_TRUE(js.at("float").is_number_float());
    EXPECT_EQ(0.1, js.at("float").get<double>());
    ASSERT_TRUE(js.at("int").is_number_integer());
    EXPECT_EQ(42, js.at("int").get<int>());
    ASSERT_TRUE(js.at("formatted").is_number());
    EXPECT_EQ(1.5, js.at("formatted").get<double>());
    // A string that looks like a number stays a string
    ASSERT_TRUE(js.at("text").is_string());
    EXPECT_EQ("12", js.at("text").get<std::string>());
}

TEST(JsonTreeTest, NumbersKeepTheirPrecision)
{
    const double small = 1.2345678901234e-9;
    const double large = -3553.1988123456;
    const float  fl    = 1234.5677F;
    JsonTree     jt("n");
    jt.addObject("small", small);
    jt.addObject("large", large);
    jt.addObject("float", fl);
    const auto js = parsed(jt);
    EXPECT_EQ(small, js.at("small").get<double>());
    EXPECT_EQ(large, js.at("large").get<double>());
    EXPECT_EQ(fl, js.at("float").get<float>());
}

TEST(JsonTreeTest, NumberConstructorsMarkNumbers)
{
    EXPECT_TRUE(JsonTree("a", 1.5).isNumber());
    EXPECT_TRUE(JsonTree("a", 1.5F).isNumber());
    EXPECT_TRUE(JsonTree("a", 3).isNumber());
    EXPECT_FALSE(JsonTree("a", std::string("3")).isNumber());
    EXPECT_FALSE(JsonTree("a").isNumber());
}

TEST(JsonTreeTest, ValueUnitWithNumberIsNumber)
{
    JsonTree jt("data");
    jt.addValueUnit("energy", "-285.8", "kJ/mol");
    jt.addValueUnit("entropy", 12.5, "J/mol K");
    const auto js = parsed(jt);
    ASSERT_TRUE(js.at("energy").at("value").is_number());
    EXPECT_EQ(-285.8, js.at("energy").at("value").get<double>());
    EXPECT_EQ("kJ/mol", js.at("energy").at("unit").get<std::string>());
    ASSERT_TRUE(js.at("entropy").at("value").is_number());
    EXPECT_EQ(12.5, js.at("entropy").at("value").get<double>());
}

TEST(JsonTreeTest, ValueUnitWithTextStaysText)
{
    JsonTree jt("data");
    jt.addValueUnit("phase", "liquid", "");
    const auto js = parsed(jt);
    ASSERT_TRUE(js.at("phase").at("value").is_string());
    EXPECT_EQ("liquid", js.at("phase").at("value").get<std::string>());
}

TEST(JsonTreeTest, NotANumberBecomesNull)
{
    JsonTree jt("n");
    jt.addObject("nan", std::nan(""));
    jt.addObject("inf", std::numeric_limits<double>::infinity());
    jt.addNumber("nantext", "nan");
    const auto js = parsed(jt);
    EXPECT_TRUE(js.at("nan").is_null());
    EXPECT_TRUE(js.at("inf").is_null());
    EXPECT_TRUE(js.at("nantext").is_null());
}

TEST(JsonTreeTest, AddNumberWithTextThatIsNoNumberWritesString)
{
    JsonTree jt("n");
    jt.addNumber("x", "12 apples");
    const auto js = parsed(jt);
    ASSERT_TRUE(js.at("x").is_string());
    EXPECT_EQ("12 apples", js.at("x").get<std::string>());
}

TEST(JsonTreeTest, AddNumberIgnoresSpaces)
{
    JsonTree jt("n");
    jt.addNumber("x", "   7.25 ");
    const auto js = parsed(jt);
    ASSERT_TRUE(js.at("x").is_number());
    EXPECT_EQ(7.25, js.at("x").get<double>());
}

TEST(JsonTreeTest, NumbersInTextOutput)
{
    JsonTree jt("n");
    jt.addObject("count", 7);
    jt.addObject("x", 2.5);
    int indent = 0;
    const std::string result = jt.writeString(false, &indent);
    EXPECT_NE(std::string::npos, result.find("count: 7"));
    EXPECT_NE(std::string::npos, result.find("x: 2.5"));
}

TEST(JsonTreeTest, EmptyStringIsWrittenAsString)
{
    JsonTree jt("n");
    jt.addObject("name", std::string(""));
    jt.addObject(JsonTree("nothing"));
    const auto js = parsed(jt);
    ASSERT_TRUE(js.at("name").is_string());
    EXPECT_EQ("", js.at("name").get<std::string>());
    // Without value and without objects there is nothing to write
    EXPECT_TRUE(js.at("nothing").is_null());
}

// ---- booleans and arrays in JSON ----

TEST(JsonTreeTest, BooleansAreJsonBooleans)
{
    JsonTree jt("b");
    jt.addBool("converged", true);
    jt.addBool("failed", false);
    jt.addObject("text", std::string("true"));
    const auto js = parsed(jt);
    ASSERT_TRUE(js.at("converged").is_boolean());
    EXPECT_TRUE(js.at("converged").get<bool>());
    ASSERT_TRUE(js.at("failed").is_boolean());
    EXPECT_FALSE(js.at("failed").get<bool>());
    // A string that looks like a boolean stays a string
    EXPECT_TRUE(js.at("text").is_string());
}

TEST(JsonTreeTest, BooleansInTextOutput)
{
    JsonTree jt("b");
    jt.addBool("converged", true);
    jt.addBool("failed", false);
    int indent = 0;
    const std::string result = jt.writeString(false, &indent);
    EXPECT_NE(std::string::npos, result.find("converged: true"));
    EXPECT_NE(std::string::npos, result.find("failed: false"));
}

TEST(JsonTreeTest, ArrayOfNumbers)
{
    JsonTree jt("root");
    jt.addArray("x", std::vector<double>({ 1.5, -2.25, 0.0 }));
    jt.addArray("n", std::vector<int>({ 3, 1, 2 }));
    jt.addArray("names", std::vector<std::string>({ "c2", "h1" }));
    const auto js = parsed(jt);
    ASSERT_TRUE(js.at("x").is_array());
    ASSERT_EQ(3U, js.at("x").size());
    EXPECT_EQ(1.5, js.at("x")[0].get<double>());
    EXPECT_EQ(-2.25, js.at("x")[1].get<double>());
    EXPECT_TRUE(js.at("x")[2].is_number());
    ASSERT_TRUE(js.at("n").is_array());
    EXPECT_TRUE(js.at("n")[0].is_number_integer());
    EXPECT_EQ(3, js.at("n")[0].get<int>());
    ASSERT_TRUE(js.at("names").is_array());
    EXPECT_EQ("h1", js.at("names")[1].get<std::string>());
}

TEST(JsonTreeTest, EmptyArrayIsWrittenAsArray)
{
    JsonTree jt("root");
    jt.addArray("x", std::vector<double>());
    jt.addObject("after", 1);
    const auto js = parsed(jt);
    ASSERT_TRUE(js.at("x").is_array());
    EXPECT_TRUE(js.at("x").empty());
}

TEST(JsonTreeTest, ArrayOfObjects)
{
    // The elements of an array keep their order, their keys are not written
    JsonTree atoms = JsonTree::array("atoms");
    for (int i = 0; i < 3; i++)
    {
        JsonTree atom("atom");
        atom.addObject("index", i);
        atom.addArray("x", std::vector<double>({ 0.1*i, 0.0, 1.0 }));
        atoms.addObject(atom);
    }
    EXPECT_TRUE(atoms.isArray());
    JsonTree jt("root");
    jt.addObject(atoms);
    const auto js = parsed(jt);
    ASSERT_TRUE(js.at("atoms").is_array());
    ASSERT_EQ(3U, js.at("atoms").size());
    EXPECT_EQ(2, js.at("atoms")[2].at("index").get<int>());
    ASSERT_TRUE(js.at("atoms")[1].at("x").is_array());
    EXPECT_EQ(1.0, js.at("atoms")[1].at("x")[2].get<double>());
}

TEST(JsonTreeTest, ArrayOfArrays)
{
    // For instance coordinates: one array of three numbers per atom
    JsonTree coords = JsonTree::array("coordinates");
    coords.addArray("", std::vector<double>({ 1, 2, 3 }));
    coords.addArray("", std::vector<double>({ 4, 5, 6 }));
    JsonTree jt("root");
    jt.addObject(coords);
    const auto js = parsed(jt);
    ASSERT_TRUE(js.at("coordinates").is_array());
    ASSERT_EQ(2U, js.at("coordinates").size());
    ASSERT_TRUE(js.at("coordinates")[1].is_array());
    EXPECT_EQ(6.0, js.at("coordinates")[1][2].get<double>());
}

TEST(JsonTreeTest, ArrayInTextOutput)
{
    JsonTree jt("root");
    jt.addArray("x", std::vector<double>({ 1.5, 2.5 }));
    int indent = 0;
    const std::string result = jt.writeString(false, &indent);
    EXPECT_NE(std::string::npos, result.find("x:"));
    EXPECT_NE(std::string::npos, result.find("1.5"));
    EXPECT_NE(std::string::npos, result.find("2.5"));
}

// ---- writeString text tests ----

TEST(JsonTreeTest, WriteStringTextSimpleValue)
{
    const JsonTree jt("name", std::string("water"));
    int indent = 0;
    const std::string result = jt.writeString(false, &indent);
    EXPECT_EQ("name: water", result);
    EXPECT_EQ(0, indent);
}

TEST(JsonTreeTest, WriteStringTextNestedLeaves)
{
    JsonTree jt("molecule");
    jt.addObject("name", std::string("water"));
    jt.addObject("formula", std::string("H2O"));
    int indent = 0;
    const std::string result = jt.writeString(false, &indent);
    // Leaf children are joined inline after the parent key
    EXPECT_NE(std::string::npos, result.find("molecule:"));
    EXPECT_NE(std::string::npos, result.find("name: water"));
    EXPECT_NE(std::string::npos, result.find("formula: H2O"));
    EXPECT_EQ(0, indent);
}

TEST(JsonTreeTest, WriteStringTextDeepTree)
{
    JsonTree jt("properties");
    jt.addValueUnit("energy", "-285.8", "kJ/mol");
    int indent = 0;
    const std::string result = jt.writeString(false, &indent);
    EXPECT_NE(std::string::npos, result.find("properties:"));
    EXPECT_NE(std::string::npos, result.find("energy:"));
    EXPECT_NE(std::string::npos, result.find("-285.8"));
    EXPECT_NE(std::string::npos, result.find("kJ/mol"));
    EXPECT_EQ(0, indent);
}

// ---- write() to file tests ----

TEST(JsonTreeTest, WriteJsonToFileProducesContent)
{
    gmx::test::TestFileManager files;
    const std::string tmpFile = files.getTemporaryFilePath(".json");

    JsonTree jt("molecule");
    jt.addObject("name", std::string("ethanol"));
    jt.addObject("formula", std::string("C2H5OH"));
    jt.write(tmpFile, true);

    std::ifstream ifs(tmpFile);
    ASSERT_TRUE(ifs.good()) << "File should be readable: " << tmpFile;
    const std::string content((std::istreambuf_iterator<char>(ifs)),
                               std::istreambuf_iterator<char>());
    EXPECT_FALSE(content.empty());
    EXPECT_NE(std::string::npos, content.find("\"molecule\""));
    EXPECT_NE(std::string::npos, content.find("\"ethanol\""));
}

TEST(JsonTreeTest, WriteTextToFileProducesContent)
{
    gmx::test::TestFileManager files;
    const std::string tmpFile = files.getTemporaryFilePath(".txt");

    JsonTree jt("molecule");
    jt.addObject("name", std::string("methane"));
    jt.write(tmpFile, false);

    std::ifstream ifs(tmpFile);
    ASSERT_TRUE(ifs.good()) << "File should be readable: " << tmpFile;
    const std::string content((std::istreambuf_iterator<char>(ifs)),
                               std::istreambuf_iterator<char>());
    EXPECT_FALSE(content.empty());
    EXPECT_NE(std::string::npos, content.find("molecule:"));
    EXPECT_NE(std::string::npos, content.find("methane"));
}

} // namespace

} // namespace alexandria
