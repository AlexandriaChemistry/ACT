/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021
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
 
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

#include <sys/time.h>
#include <sys/resource.h>

//
// Code taken from
// https://www.tutorialspoint.com/how-to-get-memory-usage-under-linux-in-cplusplus
//
#include "memory_check.h"

using namespace std;

static void mem_usage(double *vm_usage, double *resident_set) 
{
   *vm_usage = 0.0;
   *resident_set = 0.0;
   //get info from proc directory
   ifstream stat_stream("/proc/self/stat",ios_base::in);
   // Works on linux only, so check whether we have a stream
   if (stat_stream)
   {
       //create some variables to get info
       string pid, comm, state, ppid, pgrp, session, tty_nr;
       string tpgid, flags, minflt, cminflt, majflt, cmajflt;
       string utime, stime, cutime, cstime, priority, nice;
       string O, itrealvalue, starttime;
       unsigned long vsize;
       long rss;
       stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                   >> utime >> stime >> cutime >> cstime >> priority >> nice
                   >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
       stat_stream.close();
       long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // for x86-64 is configured to use 2MB pages
       *vm_usage = vsize / 1024.0;
       *resident_set = rss * page_size_kb;
   }
   else
   {
       // Apple
       // https://stackoverflow.com/questions/1543157/how-can-i-find-out-how-much-memory-my-c-app-is-using-on-the-mac
       struct rusage usage;
       if (0 == getrusage(RUSAGE_SELF, &usage))
       {
           *resident_set = usage.ru_maxrss; // bytes
       }
   }
}

void print_memory_usage_low(FILE *fp, const char *file, int line)
{
   if (fp)
   {
       double vm, rss;
       mem_usage(&vm, &rss);
       if (rss > 0)
       {
           fprintf(fp, "%s %d VMEM %g Resident %g\n", file, line, vm, rss);
           fflush(fp);
       }
   }
}
