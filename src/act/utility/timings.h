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
 * Timing of the different parts of a program.
 *
 * Regions of the code are timed with ScopedTimer, or with Timings::start
 * and Timings::stop. Regions may be nested, the time of a region excludes
 * that of the regions it contains, such that the times of all regions add
 * up to the total time. Each region is identified by a name, which is
 * looked up once per call site:
 * \code
 *     static const size_t timer = Timings::instance().index("Charge generation");
 *     ScopedTimer         scoped(timer);
 * \endcode
 * The timings are per process and not thread safe, only the main thread
 * should use them.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_UTILITY_TIMINGS_H
#define ACT_UTILITY_TIMINGS_H

#include <chrono>
#include <string>
#include <vector>

namespace gmx
{
class TextWriter;
}

namespace alexandria
{

class CommunicationRecord;

//! \brief Exclusive times of named regions of the code in this process
class Timings
{
private:
    //! The clock used
    typedef std::chrono::steady_clock Clock;
    //! A region that has been started, but not stopped
    struct Active
    {
        //! Index of the region
        size_t            id;
        //! When it was started
        Clock::time_point start;
        //! Time spent in the regions it contains
        double            children;
    };
    //! Names of the regions
    std::vector<std::string> names_;
    //! Time of each region, excluding the regions it contains, in seconds
    std::vector<double>      seconds_;
    //! Number of times each region was started
    std::vector<double>      calls_;
    //! The regions that have been started but not stopped, innermost last
    std::vector<Active>      active_;

    //! Use instance() instead
    Timings() {}
public:
    //! \return the timings of this process
    static Timings &instance();

    /*! \brief Find a region
     * \param[in] name The name of the region, it is added if it is new
     * \return the index of the region
     */
    size_t index(const std::string &name);

    //! \brief Start timing region id
    void start(size_t id);

    //! \brief Stop timing region id, it must be the last one started
    void stop(size_t id);

    //! \return the names of the regions
    const std::vector<std::string> &names() const { return names_; }

    //! \return the time of each region in seconds, excluding the regions it contains
    const std::vector<double> &seconds() const { return seconds_; }

    //! \return how often each region was started
    const std::vector<double> &calls() const { return calls_; }
};

//! \brief Times a region from construction until destruction
class ScopedTimer
{
private:
    //! Index of the region
    size_t id_;
public:
    /*! \brief Start timing
     * \param[in] id Index of the region, from Timings::index
     */
    explicit ScopedTimer(size_t id) : id_(id) { Timings::instance().start(id_); }

    //! \brief Stop timing
    ~ScopedTimer() { Timings::instance().stop(id_); }

    //! Copying would stop the timer twice
    ScopedTimer(const ScopedTimer &) = delete;
    //! Copying would stop the timer twice
    ScopedTimer &operator=(const ScopedTimer &) = delete;
};

/*! \brief Collect the timings of all processes and print them on the master
 * All processes in the world communicator have to call this.
 * \param[in] cr       The communication record
 * \param[in] tw       Where to print, only used on the master and may be nullptr
 * \param[in] wallTime Wall time of this process in seconds
 */
void printTimings(const CommunicationRecord *cr,
                  gmx::TextWriter           *tw,
                  double                     wallTime);

} // namespace alexandria

#endif
