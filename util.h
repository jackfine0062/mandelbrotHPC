/**
 * Several utility functions for displaying results.
 */

#pragma once

#include <stdlib.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Prints a positive number with the given number of sigfigs and a unit. The
 * value is scaled to the correct unit (which are mult apart - 1000 for SI and
 * 1024 for digit prefixes).
 */
void print_with_unit(double val, int sigfigs, int mult,
                     const char** units, size_t n_units);

/**
 * Prints a number of bytes after converting to a nicer unit.
 */
void print_bytes(size_t n);

/**
 * Print the time (in seconds) with the right units and 3 significant digits.
 */
void print_time(double seconds);

/**
 * Get the difference between two times.
 */
double get_time_diff(struct timespec* start, struct timespec* end);

/**
 * Get the number of physical cores on the machine.
 */
size_t get_num_physical_cores();

/**
 * Get the number of logical cores on the machine.
 * This includes hardware threads (i.e. hyperthreads).
 */
size_t get_num_logical_cores();

/**
 * Get the number of cores dedicted to this process.
 */
size_t get_num_cores_affinity();

#ifdef __cplusplus
}
#endif