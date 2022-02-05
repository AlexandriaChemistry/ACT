#ifndef _memory_check_h
#define _memory_check_h

#include <cstdio>

extern void print_memory_usage_low(FILE *fp, const char *file, int line);

#define print_memory_usage(fp) print_memory_usage_low(fp, __FILE__, __LINE__)

#endif

