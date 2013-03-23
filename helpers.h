/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef HELPERS_H
#define HELPERS_H

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <time.h>

// M_PI isn't defined on windows
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
#define popen _popen
#define pclose _pclose
#define strdup _strdup

#if defined(_USE_32BIT_TIME_T)
#define _gmtime_s _gmtime32_s
#else
#define _gmtime_s _gmtime64_s
#endif

#endif

#define error_jump(label, ret, ...) do { ret = error(__VA_ARGS__); goto label; } while(0)

typedef struct
{
    time_t time;
    uint16_t ms;
} ts_time;

void free_2d_array(char **array, size_t len);
float *cast_double_array_to_float(double *d_ptr, size_t count);
size_t get_matching_files(const char *pattern, char ***files);
char *get_first_matching_file(char *pattern);
int compare_double(const void *a, const void *b);
int compare_float(const void *a, const void *b);
int error(const char * format, ...);
void millisleep(int ms);
int ts_exec_read(const char *cmd, char **output);
int ts_exec_write(const char *cmd, const void *restrict data, size_t size);
int init_ds9();

void calculate_amplitude_spectrum(double *time, double *mma, size_t count,
                                  double freq_min, double freq_max,
                                  double *freq, double *ampl, size_t length);

char *canonicalize_path(const char *path);

ts_time parse_time(const char *string);
ts_time parse_time_ccdops(const char *string);
ts_time parse_date_time(const char *date, const char *time);
void serialize_time(ts_time t, char buf[24]);
double ts_time_to_utc_hour(ts_time t);
double ts_time_to_bjd(ts_time t, double ra, double dec, double epoch);
double ts_difftime(ts_time a, ts_time b);
char *prompt_user_input(char *message, char *fallback);

bool region_contains(uint16_t r[4], size_t x, size_t y);
double region_mean(uint16_t r[4], double *data, size_t stride);

double mean_exclude_sigma(double *data, size_t n, double sigma);

#endif
