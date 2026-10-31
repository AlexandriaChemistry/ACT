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
class JsonTree;

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
    //! Floating point operations counted in each region, see addFlops
    std::vector<double>      flops_;
    //! Floating point operations counted outside of all regions
    double                   untimedFlops_ = 0;
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

    /*! \brief Start timing a region
     * \param[in] id The index of the region
     */
    void start(size_t id);

    /*! \brief Stop timing a region
     * \param[in] id The index of the region, it must be the last one started
     */
    void stop(size_t id);

    /*! \brief Set the times and numbers of calls of all regions to zero
     * The regions are kept, since their indices are stored where they
     * are used. No region may be active.
     */
    void reset();

    /*! \brief Count floating point operations
     * They are attributed to the innermost region that is active, such
     * that the rate in a region follows from its time. See printTimings
     * for what is counted.
     * \param[in] flops The number of floating point operations
     */
    void addFlops(double flops)
    {
        if (active_.empty())
        {
            untimedFlops_ += flops;
        }
        else
        {
            flops_[active_.back().id] += flops;
        }
    }

    //! Returned by parent if no region is active
    static constexpr size_t none = static_cast<size_t>(-1);

    //! \return the innermost region that is active, or none
    size_t parent() const { return active_.empty() ? none : active_.back().id; }

    //! \return the names of the regions
    const std::vector<std::string> &names() const { return names_; }

    //! \return the time of each region in seconds, excluding the regions it contains
    const std::vector<double> &seconds() const { return seconds_; }

    //! \return how often each region was started
    const std::vector<double> &calls() const { return calls_; }

    //! \return the floating point operations counted in each region
    const std::vector<double> &flops() const { return flops_; }

    //! \return the floating point operations counted outside of all regions
    double untimedFlops() const { return untimedFlops_; }
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

/*! \brief Times the consecutive phases of a program and its wall time
 * \code
 *     PhaseTimer phase("Startup");
 *     ...
 *     phase.next("Output");
 *     ...
 *     phase.stop();
 *     printTimings(cr, tw, phase.wallTime());
 * \endcode
 */
class PhaseTimer
{
private:
    //! When the program started
    std::chrono::steady_clock::time_point wallStart_;
    //! Index of the current phase
    size_t                                id_;
    //! Whether a phase is being timed
    bool                                  running_ = true;
public:
    /*! \brief Start the wall time and the first phase
     * If no region is active, the timings are reset first, such that they
     * are those of this program only, even if it is called more than once
     * in the same process.
     * \param[in] name Name of the first phase
     */
    explicit PhaseTimer(const std::string &name);

    //! \brief Stop the current phase if it is still running, e.g. after an early return
    ~PhaseTimer();

    //! Copying would stop the phase twice
    PhaseTimer(const PhaseTimer &) = delete;
    //! Copying would stop the phase twice
    PhaseTimer &operator=(const PhaseTimer &) = delete;

    /*! \brief Stop the current phase and start the next one
     * \param[in] name Name of the next phase, the current phase
     *                 continues if it has this name
     */
    void next(const std::string &name);

    //! \brief Stop the current phase, all regions it contains must have been stopped
    void stop();

    //! \return the wall time since construction in seconds
    double wallTime() const;
};

//! \brief Theoretical floating point performance of one core
struct PeakPerformance
{
    //! Double precision GFlop/s, zero if unknown
    double      gflops = 0;
    //! How it was determined
    std::string basis;
};

/*! \brief Collect the timings of all processes and print them on the master
 * All processes in the world communicator have to call this. If floating
 * point operations were counted, see Timings::addFlops, the floating point
 * performance is printed as well, per core on average over the wall time
 * of all processes and in the parts in which they were counted, and as a
 * fraction of the theoretical performance if that is given. The peak is
 * for double precision, which all calculations use if GMX_DOUBLE is set,
 * otherwise only the potentials with explicit SIMD, as is printed.
 * \param[in]  cr       The communication record
 * \param[in]  tw       Where to print, only used on the master and may be nullptr
 * \param[in]  wallTime Wall time of this process in seconds
 * \param[out] jtree    If not nullptr, the same table is added to it on the master
 * \param[in]  peak     Theoretical performance of one core, may be nullptr
 */
void printTimings(const CommunicationRecord *cr,
                  gmx::TextWriter           *tw,
                  double                     wallTime,
                  JsonTree                  *jtree = nullptr,
                  const PeakPerformance     *peak  = nullptr);

} // namespace alexandria

#endif
