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
#include <type_traits>
#include <utility>

namespace alexandria
{

// This code was inspired by:
// https://www.fluentcpp.com/2021/01/29/inheritance-without-pointers/
// See also
// https://stackoverflow.com/questions/76379814/how-to-use-const-variables-using-c17-stdany-in-a-stdvector
// It requires C++17.
//
// The address of the stored object is kept, since std::any_cast checks the
// type on every access, which is expensive in inner loops. The address
// changes when the storage is copied or moved, therefore it is recomputed
// there.

/*! \brief Container holding an object of any class derived from ActBase by value
 * \tparam ActBase The base class that the stored objects derive from
 */
template<typename ActBase>
struct BaseContainer
{
public:
    /*! \brief Constructor storing a copy of, or moving, an object
     * \param[in] object The object of a class derived from ActBase
     */
    template<typename ConcreteType>
        requires (!std::is_same_v<std::remove_cvref_t<ConcreteType>, BaseContainer>)
    BaseContainer(ConcreteType&& object)
        : storage{std::forward<ConcreteType>(object)}
        , getter{ [](std::any &storage) -> ActBase& { return std::any_cast<std::remove_cvref_t<ConcreteType>&>(storage); } }
        , object_{ &getter(storage) }
    {}

    /*! \brief Copy constructor
     * \param[in] other The container to copy
     */
    BaseContainer(const BaseContainer &other)
        : storage{other.storage}, getter{other.getter}, object_{ &getter(storage) }
    {}

    /*! \brief Move constructor
     * \param[inout] other The container to move from
     */
    BaseContainer(BaseContainer &&other) noexcept
        : storage{std::move(other.storage)}, getter{other.getter}, object_{ &getter(storage) }
    {}

    /*! \brief Copy assignment
     * \param[in] other The container to copy
     * \return this container
     */
    BaseContainer &operator=(const BaseContainer &other)
    {
        if (this != &other)
        {
            storage = other.storage;
            getter  = other.getter;
            object_ = &getter(storage);
        }
        return *this;
    }

    /*! \brief Move assignment
     * \param[inout] other The container to move from
     * \return this container
     */
    BaseContainer &operator=(BaseContainer &&other) noexcept
    {
        if (this != &other)
        {
            storage = std::move(other.storage);
            getter  = other.getter;
            object_ = &getter(storage);
        }
        return *this;
    }

    //! \return a const pointer to the stored object
    const ActBase *operator->() const { return object_; }
    //! \return a pointer to the stored object
    ActBase *operator->() { return object_; }

private:
    //! Type-erased storage holding the object
    std::any storage;
    //! Function returning a reference to the object in storage, cast to ActBase
    ActBase& (*getter)(std::any&);
    //! The stored object
    ActBase *object_;
};

} // namespace alexandria

#endif
