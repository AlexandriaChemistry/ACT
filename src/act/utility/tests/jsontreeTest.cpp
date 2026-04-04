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

#include "../jsontree.h"

#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

// ---- Accessor and state tests ----

TEST(JsonTreeTest, KeyOnlyConstructorIsEmpty)
{
    JsonTree jt("myKey");
    EXPECT_EQ("myKey", jt.key());
    EXPECT_TRUE(jt.value().empty());
    EXPECT_TRUE(jt.objects().empty());
    EXPECT_TRUE(jt.empty());
    EXPECT_FALSE(jt.branched());
}

TEST(JsonTreeTest, StringValueConstructorAccessors)
{
    JsonTree jt("molecule", std::string("water"));
    EXPECT_EQ("molecule", jt.key());
    EXPECT_EQ("water", jt.value());
    EXPECT_FALSE(jt.empty());
    EXPECT_FALSE(jt.branched());
    EXPECT_TRUE(jt.objects().empty());
}

TEST(JsonTreeTest, DoubleValueConstructorNotEmpty)
{
    JsonTree jt("pi", 3.14);
    EXPECT_EQ("pi", jt.key());
    EXPECT_FALSE(jt.value().empty());
    EXPECT_FALSE(jt.empty());
    EXPECT_FALSE(jt.branched());
}

TEST(JsonTreeTest, IntValueConstructorStoresAsString)
{
    JsonTree jt("count", 42);
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
    JsonTree child("child", std::string("val"));
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
    JsonTree jt("name", std::string("water"));
    int indent = 0;
    std::string result = jt.writeString(true, &indent);
    EXPECT_EQ(R"({"name":"water"})", result);
    EXPECT_EQ(0, indent);
}

TEST(JsonTreeTest, WriteStringJsonIntValue)
{
    JsonTree jt("count", 7);
    int indent = 0;
    std::string result = jt.writeString(true, &indent);
    EXPECT_EQ(R"({"count":"7"})", result);
    EXPECT_EQ(0, indent);
}

TEST(JsonTreeTest, WriteStringJsonNestedLeaves)
{
    JsonTree jt("molecule");
    jt.addObject("name", std::string("water"));
    jt.addObject("formula", std::string("H2O"));
    int indent = 0;
    std::string result = jt.writeString(true, &indent);
    // Proper JSON: children become members of an object
    std::string expected =
        "{\"molecule\":{\"name\":\"water\",\"formula\":\"H2O\"}}";
    EXPECT_EQ(expected, result);
    EXPECT_EQ(0, indent);
}

TEST(JsonTreeTest, WriteStringJsonValueUnitStructure)
{
    JsonTree jt("data");
    jt.addValueUnit("energy", "42.0", "kJ/mol");
    int indent = 0;
    std::string result = jt.writeString(true, &indent);
    EXPECT_NE(std::string::npos, result.find("\"energy\""));
    EXPECT_NE(std::string::npos, result.find("\"value\""));
    EXPECT_NE(std::string::npos, result.find("\"42.0\""));
    EXPECT_NE(std::string::npos, result.find("\"unit\""));
    EXPECT_NE(std::string::npos, result.find("\"kJ/mol\""));
    EXPECT_EQ(0, indent);
}

// ---- writeString text tests ----

TEST(JsonTreeTest, WriteStringTextSimpleValue)
{
    JsonTree jt("name", std::string("water"));
    int indent = 0;
    std::string result = jt.writeString(false, &indent);
    EXPECT_EQ("name: water", result);
    EXPECT_EQ(0, indent);
}

TEST(JsonTreeTest, WriteStringTextNestedLeaves)
{
    JsonTree jt("molecule");
    jt.addObject("name", std::string("water"));
    jt.addObject("formula", std::string("H2O"));
    int indent = 0;
    std::string result = jt.writeString(false, &indent);
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
    std::string result = jt.writeString(false, &indent);
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
    std::string tmpFile = files.getTemporaryFilePath(".json");

    JsonTree jt("molecule");
    jt.addObject("name", std::string("ethanol"));
    jt.addObject("formula", std::string("C2H5OH"));
    jt.write(tmpFile, true);

    std::ifstream ifs(tmpFile);
    ASSERT_TRUE(ifs.good()) << "File should be readable: " << tmpFile;
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());
    EXPECT_FALSE(content.empty());
    EXPECT_NE(std::string::npos, content.find("\"molecule\""));
    EXPECT_NE(std::string::npos, content.find("\"ethanol\""));
}

TEST(JsonTreeTest, WriteTextToFileProducesContent)
{
    gmx::test::TestFileManager files;
    std::string tmpFile = files.getTemporaryFilePath(".txt");

    JsonTree jt("molecule");
    jt.addObject("name", std::string("methane"));
    jt.write(tmpFile, false);

    std::ifstream ifs(tmpFile);
    ASSERT_TRUE(ifs.good()) << "File should be readable: " << tmpFile;
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());
    EXPECT_FALSE(content.empty());
    EXPECT_NE(std::string::npos, content.find("molecule:"));
    EXPECT_NE(std::string::npos, content.find("methane"));
}

} // namespace

} // namespace alexandria
