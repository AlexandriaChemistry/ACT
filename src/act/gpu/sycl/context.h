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
 * The SYCL queue of the GPU and the device memory, kept between calls,
 * such that the memory is only allocated when more is needed.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_GPU_SYCL_CONTEXT_H
#define ACT_GPU_SYCL_CONTEXT_H

#include <cstddef>

#include <sycl/sycl.hpp>

namespace alexandria
{

namespace gpu
{

/*! \brief Device memory that grows when needed
 * The arrays are meant to be kept until the end of the program and are
 * not freed then, see gpuQueue.
 * \tparam T The type of the elements
 */
template <typename T>
class DeviceArray
{
public:
    DeviceArray() = default;
    DeviceArray(const DeviceArray &) = delete;
    DeviceArray &operator=(const DeviceArray &) = delete;

    /*! \brief Make room for at least n elements, the contents are lost when it grows
     * \param[in] queue The queue of the GPU
     * \param[in] n     Number of elements
     * \return the device memory
     */
    T *reserve(sycl::queue &queue, size_t n)
    {
        if (n > size_)
        {
            release();
            // Half more, such that slowly growing sizes allocate rarely
            size_   = n + n/2;
            queue_  = &queue;
            data_   = sycl::malloc_device<T>(size_, queue);
        }
        return data_;
    }

private:
    void release()
    {
        if (data_)
        {
            sycl::free(data_, *queue_);
        }
        data_ = nullptr;
        size_ = 0;
    }
    T           *data_  = nullptr;
    size_t       size_  = 0;
    sycl::queue *queue_ = nullptr;
};

/*! \brief The queue of the GPU used, nullptr without GPU
 * The first GPU found is used, with an in-order queue. It is never
 * deleted, since the order of destruction at the end of the program is
 * unknown relative to that of the SYCL runtime.
 * \return the queue
 */
sycl::queue *gpuQueue();

} // namespace gpu

} // namespace alexandria

#endif
