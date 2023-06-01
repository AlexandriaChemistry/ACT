/*
 * This source file is part of the Alexandria Chemistry Toolkit.
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

#ifndef ACT_BASICS_BASECONTAINER_H
#define ACT_BASICS_BASECONTAINER_H

#include <any>

namespace alexandria
{

// This code was inspired by:
// https://www.fluentcpp.com/2021/01/29/inheritance-without-pointers/ 
// It requires C++17.

template<typename ActBase>
struct BaseContainer
{
public:
    template<typename ConcreteType>
    BaseContainer(ConcreteType&& object)
        : storage{std::forward<ConcreteType>(object)}
        , getter{ [](std::any &storage) -> ActBase& { return std::any_cast<ConcreteType&>(storage); } }
    {}
    
    const ActBase *operator->() const { return &getter(storage); }
    
private:
    std::any storage;
    const ActBase& (*getter)(const std::any&);
};

} // namespace alexandria

#endif
