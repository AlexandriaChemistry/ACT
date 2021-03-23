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
