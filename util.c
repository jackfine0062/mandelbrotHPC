/**
 * Several utility functions for displaying results.
 */

#if defined(linux)
#define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "util.h"

/**
 * Prints a positive number with the given number of sigfigs and a unit. The
 * value is scaled to the correct unit (which are mult apart - 1000 for SI and
 * 1024 for digit prefixes).
 */
void print_with_unit(double val, int sigfigs, int mult,
                     const char** units, size_t n_units) {
    size_t i_unit = 0;
    while (i_unit < n_units && val >= mult) { val /= mult; i_unit++; }
    if (i_unit == 0) { sigfigs = 0; }
    else if (val < 10) { sigfigs -= 1; }
    else if (val < 100) { sigfigs -= 2; }
    else { sigfigs -= 3; }
    printf("%.*f %s", sigfigs, val, units[i_unit]);
}

/**
 * Prints a number of bytes after converting to a nicer unit.
 */
void print_bytes(size_t n) {
    static const char* units[4] = {"bytes", "KiB", "MiB", "GiB"};
    print_with_unit(n, 3, 1024, units, 4);
}

/**
 * Print the time (in seconds) with the right units and 3 significant digits.
 */
void print_time(double seconds) {
    static const char* units[4] = {"ns", "us", "ms", "s"};
    print_with_unit(seconds * 1000000000.0, 3, 1000, units, 4);
}

/**
 * Get the difference between two times.
 */
double get_time_diff(struct timespec* start, struct timespec* end) {
    double diff = end->tv_sec - start->tv_sec;
    diff += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return diff;
}

// get_num_physical_cores() and get_num_logical_cores() have to be specialized
// for each OS.
#if defined(__APPLE__)
#include <sys/sysctl.h>
size_t __get_sysctl_size_t(const char* name) {
    size_t var = 0, sizeof_var = sizeof(size_t);
    sysctlbyname(name, &var, &sizeof_var, 0, 0);
    return var;
}
size_t get_num_physical_cores() { return __get_sysctl_size_t("hw.physicalcpu"); }
size_t get_num_logical_cores() { return __get_sysctl_size_t("hw.logicalcpu"); }
size_t get_num_cores_affinity() { return get_num_logical_cores(); } // macOS doesn't really support affinity
#elif defined(linux)
#include <unistd.h>
#include <sched.h>
size_t get_num_physical_cores() { return (sysconf(_SC_NPROCESSORS_ONLN) + 1) / 2; } // TODO: this assumes processor has 2 threads per core
size_t get_num_logical_cores() { return sysconf(_SC_NPROCESSORS_ONLN); }
size_t get_num_cores_affinity() { cpu_set_t cs; CPU_ZERO(&cs); sched_getaffinity(0, sizeof(cs), &cs); return CPU_COUNT(&cs); }
#else
#error Unrecognized OS
#endif
