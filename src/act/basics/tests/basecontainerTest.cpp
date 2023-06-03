/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2023
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
#include <gtest/gtest.h>

#include "../basecontainer.h"

#include <vector>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

class TopologyEntry
{
public:
    TopologyEntry() {};
    
    const TopologyEntry *self() const { return this; }
    
    int print() const
    {
        return 2;
    }
};

class TopologyInherited : public TopologyEntry
{
public:
    TopologyInherited() {}
    
    int koko() const
    {
        return 3;
    }
};


class BaseContainerTest : public gmx::test::CommandLineTestBase
{
protected:
    std::vector<BaseContainer<TopologyEntry>> entries_;
    BaseContainerTest ()
    {
        entries_.push_back(TopologyEntry{});
        entries_.push_back(TopologyInherited{});
    }
    
    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {
    }
    
    static void TearDownTestCase()
    {
    }
    
    void runTest(int i, int expected)
    {
        int result;
        if (i == 0)
        {
            result = entries_[0]->print();
        }
        else
        {
            auto f = static_cast<const TopologyInherited *>(entries_[1]->self());
            result = f->koko();
        }
        EXPECT_TRUE(result == expected);
    }
    
    void printThem(const std::vector<BaseContainer<TopologyEntry>> &entries)
    {
        size_t i = 0;
        for (const auto &e : entries)
        {
            EXPECT_TRUE(2 == e->print());
            if (i > 0)
            {
                auto f = static_cast<const TopologyInherited *>(e->self());
                EXPECT_TRUE(3 == f->koko());
            }
            i++;
        }
    }
    
    void doPrint()
    {
        printThem(entries_);
    }
};

TEST_F(BaseContainerTest, Base) 
{
    runTest(0, 2);
}

TEST_F(BaseContainerTest, Inherited) 
{
    runTest(1, 3);
}

TEST_F(BaseContainerTest, ConstVector) 
{
    doPrint();
}

} // namespace

} // namespace alexandria

