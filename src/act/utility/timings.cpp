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
 * Implements the timing of the different parts of a program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include "timings.h"

#include <algorithm>
#include <map>
#include <numeric>

#include "act/utility/communicationrecord.h"
#include "act/utility/jsontree.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace alexandria
{

Timings &Timings::instance()
{
    static Timings timings;
    return timings;
}

size_t Timings::index(const std::string &name)
{
    auto it = std::find(names_.begin(), names_.end(), name);
    if (it != names_.end())
    {
        return it - names_.begin();
    }
    names_.push_back(name);
    seconds_.push_back(0);
    calls_.push_back(0);
    return names_.size() - 1;
}

void Timings::start(size_t id)
{
    active_.push_back({ id, Clock::now(), 0.0 });
}

void Timings::stop(size_t id)
{
    if (active_.empty() || active_.back().id != id)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Stopping timer %s, which is not the last one started",
                                                       id < names_.size() ? names_[id].c_str() : "unknown")));
    }
    auto   &a       = active_.back();
    double  elapsed = std::chrono::duration<double>(Clock::now() - a.start).count();
    seconds_[id]   += elapsed - a.children;
    calls_[id]     += 1;
    active_.pop_back();
    if (!active_.empty())
    {
        active_.back().children += elapsed;
    }
}

void Timings::reset()
{
    if (!active_.empty())
    {
        GMX_THROW(gmx::InternalError("Cannot reset the timings while a region is active"));
    }
    std::fill(seconds_.begin(), seconds_.end(), 0.0);
    std::fill(calls_.begin(), calls_.end(), 0.0);
}

PhaseTimer::PhaseTimer(const std::string &name)
    : wallStart_(std::chrono::steady_clock::now())
{
    auto &timings = Timings::instance();
    if (Timings::none == timings.parent())
    {
        timings.reset();
    }
    id_ = timings.index(name);
    timings.start(id_);
}

PhaseTimer::~PhaseTimer()
{
    // Only stop if that is possible, a destructor should not throw
    if (running_ && Timings::instance().parent() == id_)
    {
        Timings::instance().stop(id_);
    }
}

void PhaseTimer::next(const std::string &name)
{
    auto id = Timings::instance().index(name);
    if (running_ && id == id_)
    {
        return;
    }
    stop();
    id_      = id;
    Timings::instance().start(id_);
    running_ = true;
}

void PhaseTimer::stop()
{
    if (running_)
    {
        Timings::instance().stop(id_);
        running_ = false;
    }
}

double PhaseTimer::wallTime() const
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - wallStart_).count();
}

namespace
{

//! The timings of one process
struct ProcessTimings
{
    //! Role of the process
    NodeType                 role;
    //! Wall time of the process
    double                   wallTime;
    //! Names of the regions
    std::vector<std::string> names;
    //! Time of each region
    std::vector<double>      seconds;
    //! Number of calls of each region
    std::vector<double>      calls;
};

//! Separates the names of the regions when they are sent
const char nameSeparator = '\n';

} // namespace

void printTimings(const CommunicationRecord *cr,
                  gmx::TextWriter           *tw,
                  double                     wallTime,
                  JsonTree                  *jtree)
{
    auto          &timings = Timings::instance();
    ProcessTimings mine    = { cr->nodeType(), wallTime, timings.names(),
                               timings.seconds(), timings.calls() };
    std::vector<ProcessTimings> all;
    if (!cr->isParallel())
    {
        all.push_back(mine);
    }
    else
    {
        // Collective operations are used, since point to point communication
        // is turned off when the master has finished, see CommunicationRecord::done.
        MPI_Comm comm  = cr->comm_world();
        int      nproc = cr->size();
        bool     root  = cr->isMaster();
        // Role and wall time
        double              info[2] = { static_cast<double>(mine.role), mine.wallTime };
        std::vector<double> infos(root ? 2*nproc : 0);
        MPI_Gather(info, 2, MPI_DOUBLE, infos.data(), 2, MPI_DOUBLE, 0, comm);
        // Gather a vector of variable length per process
        auto gather = [comm, nproc, root](const auto &mydata, MPI_Datatype type,
                                          std::vector<int> *counts, auto *alldata)
        {
            int count = mydata.size();
            counts->resize(root ? nproc : 0);
            MPI_Gather(&count, 1, MPI_INT, counts->data(), 1, MPI_INT, 0, comm);
            std::vector<int> displs(counts->size(), 0);
            for (size_t i = 1; i < counts->size(); i++)
            {
                displs[i] = displs[i-1] + (*counts)[i-1];
            }
            if (root)
            {
                alldata->resize(displs.empty() ? 0 : displs.back() + counts->back());
            }
            MPI_Gatherv(mydata.data(), count, type, alldata->data(), counts->data(),
                        displs.data(), type, 0, comm);
        };
        std::string names;
        for (const auto &name : mine.names)
        {
            names += name + nameSeparator;
        }
        std::vector<int>    nameCounts, secondCounts, callCounts;
        std::string         allNames;
        std::vector<double> allSeconds, allCalls;
        gather(names, MPI_CHAR, &nameCounts, &allNames);
        gather(mine.seconds, MPI_DOUBLE, &secondCounts, &allSeconds);
        gather(mine.calls, MPI_DOUBLE, &callCounts, &allCalls);
        if (!root)
        {
            return;
        }
        size_t nameOffset = 0, valueOffset = 0;
        for (int p = 0; p < nproc; p++)
        {
            ProcessTimings pt;
            pt.role     = static_cast<NodeType>(static_cast<int>(infos[2*p]));
            pt.wallTime = infos[2*p+1];
            std::string pnames = allNames.substr(nameOffset, nameCounts[p]);
            size_t      begin  = 0;
            for (size_t end = pnames.find(nameSeparator); end != std::string::npos;
                 end = pnames.find(nameSeparator, begin))
            {
                pt.names.push_back(pnames.substr(begin, end - begin));
                begin = end + 1;
            }
            pt.seconds.assign(allSeconds.begin() + valueOffset, allSeconds.begin() + valueOffset + secondCounts[p]);
            pt.calls.assign(allCalls.begin() + valueOffset, allCalls.begin() + valueOffset + callCounts[p]);
            nameOffset  += nameCounts[p];
            valueOffset += secondCounts[p];
            all.push_back(std::move(pt));
        }
    }
    if (nullptr == tw && nullptr == jtree)
    {
        return;
    }
    // Sum over the processes, in total and per role
    const std::vector<NodeType>    roles     = { NodeType::Master, NodeType::MiddleMan, NodeType::Helper };
    const std::vector<std::string> roleNames = { "Master", "Middlemen", "Helpers" };
    std::map<NodeType, int>        nproc;
    struct Sum
    {
        double                     seconds = 0;
        double                     calls   = 0;
        std::map<NodeType, double> perRole;
    };
    std::map<std::string, Sum> sums;
    const std::string          notTimed  = "Not timed";
    double                     totalWall = 0;
    for (const auto &pt : all)
    {
        nproc[pt.role] += 1;
        totalWall      += pt.wallTime;
        double timed    = 0;
        for (size_t i = 0; i < pt.names.size() && i < pt.seconds.size() && i < pt.calls.size(); i++)
        {
            if (0 == pt.calls[i])
            {
                // Not used since the timings were reset
                continue;
            }
            auto &s = sums[pt.names[i]];
            s.seconds           += pt.seconds[i];
            s.calls             += pt.calls[i];
            s.perRole[pt.role]  += pt.seconds[i];
            timed               += pt.seconds[i];
        }
        auto &s = sums[notTimed];
        s.seconds          += pt.wallTime - timed;
        s.perRole[pt.role] += pt.wallTime - timed;
    }
    std::vector<std::pair<std::string, Sum>> rows(sums.begin(), sums.end());
    std::sort(rows.begin(), rows.end(),
              [](const auto &a, const auto &b) { return a.second.seconds > b.second.seconds; });

    std::string counts;
    for (size_t r = 0; r < roles.size(); r++)
    {
        if (nproc[roles[r]] > 0)
        {
            counts += gmx::formatString("%s%d %s", counts.empty() ? "" : ", ",
                                        nproc[roles[r]], roleNames[r].c_str());
        }
    }
    if (tw)
    {
        tw->writeStringFormatted("\nTime spent in the different parts of the program, summed over %zu process(es)",
                                 all.size());
        tw->writeStringFormatted(" (%s).\n", counts.c_str());
        tw->writeStringFormatted("The time of each part excludes the parts that it calls. MPI times are given per\n");
        tw->writeStringFormatted("part of the program they are called from. 'MPI waiting' is waiting for a message\n");
        tw->writeStringFormatted("to arrive, 'MPI collectives' includes waiting for the other processes. The\n");
        tw->writeStringFormatted("columns per role are averages per process in seconds.\n");
        tw->writeStringFormatted("Wall time of the master %.2f s.\n\n", mine.wallTime);
        std::string header = gmx::formatString("%-48s %10s %6s", "Part", "Seconds", "%");
        for (size_t r = 0; r < roles.size(); r++)
        {
            if (nproc[roles[r]] > 0)
            {
                header += gmx::formatString(" %10s", roleNames[r].c_str());
            }
        }
        header += gmx::formatString(" %12s", "Calls");
        tw->writeStringFormatted("%s\n", header.c_str());
    }
    // The same table in JSON, with the per role columns as averages
    // per process in seconds
    JsonTree jparts("Parts");
    auto printRow = [&](const std::string &name, const Sum &s)
    {
        double      percent = totalWall > 0 ? 100*s.seconds/totalWall : 0.0;
        std::string line    = gmx::formatString("%-48s %10.2f %6.2f", name.c_str(), s.seconds, percent);
        JsonTree    jrow(name);
        jrow.addObject("Seconds", gmx::formatString("%.3f", s.seconds));
        jrow.addObject("Percent", gmx::formatString("%.2f", percent));
        for (size_t r = 0; r < roles.size(); r++)
        {
            if (nproc[roles[r]] > 0)
            {
                auto   pr      = s.perRole.find(roles[r]);
                double average = (pr == s.perRole.end() ? 0.0 : pr->second)/nproc[roles[r]];
                line += gmx::formatString(" %10.3f", average);
                jrow.addObject(roleNames[r], gmx::formatString("%.3f", average));
            }
        }
        if (s.calls > 0)
        {
            line += gmx::formatString(" %12.0f", s.calls);
            jrow.addObject("Calls", gmx::formatString("%.0f", s.calls));
        }
        if (tw)
        {
            tw->writeStringFormatted("%s\n", line.c_str());
        }
        jparts.addObject(jrow);
    };
    for (const auto &row : rows)
    {
        printRow(row.first, row.second);
    }
    // Subtotals of the potentials and of MPI
    auto subtotal = [&rows, &printRow](const std::string &prefix, const std::string &name)
    {
        Sum sum;
        for (const auto &row : rows)
        {
            if (row.first.compare(0, prefix.size(), prefix) == 0)
            {
                sum.seconds += row.second.seconds;
                sum.calls   += row.second.calls;
                for (const auto &pr : row.second.perRole)
                {
                    sum.perRole[pr.first] += pr.second;
                }
            }
        }
        if (sum.seconds > 0)
        {
            printRow(name, sum);
        }
    };
    subtotal("Potential ", "All potentials");
    subtotal("MPI ", "All MPI");
    Sum total;
    total.seconds = totalWall;
    for (const auto &pt : all)
    {
        total.perRole[pt.role] += pt.wallTime;
    }
    printRow("Total", total);
    if (tw)
    {
        tw->writeStringFormatted("\n");
    }
    if (jtree)
    {
        JsonTree jtimings("Timings");
        jtimings.addObject("Processes", counts);
        jtimings.addValueUnit("Wall time of the master", gmx::formatString("%.3f", mine.wallTime), "s");
        jtimings.addObject(jparts);
        jtree->addObject(jtimings);
    }
}

} // namespace alexandria
