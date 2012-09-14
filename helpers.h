/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#ifndef HELPERS_H
#define HELPERS_H

#include <stdarg.h>
#include <time.h>

// M_PI isn't defined on windows
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define error_jump(label, ret, ...) do { ret = error(__VA_ARGS__); goto label; } while(0)

void free_2d_array(char **array, int len);
float *cast_double_array_to_float(double *d_ptr, size_t count);
int get_matching_files(const char *pattern, char ***files);
char *get_first_matching_file(char *pattern);
int compare_double(const void *a, const void *b);
int error(const char * format, ...);
void die(const char * format, ...) __attribute__ ((noreturn));
int ts_exec_read(const char *cmd, char **output);
int ts_exec_write(const char *cmd, const void *restrict data, size_t size);
int init_ds9();

void calculate_amplitude_spectrum(double fmin, double fmax, double *t, double *data, int numData, double *outFreq, double *outAmpl, int numOut);
void calculate_amplitude_spectrum_float(float fmin, float fmax, float *time, float *data, int numData, float *outFreq, float *outAmpl, int numOut);

#if (defined _WIN32)
int vasprintf (char **resultp, const char *format, va_list args);
int asprintf(char **resultp, const char *format, ...);
char *strndup(const char *s, size_t max);
#endif

char *canonicalize_path(const char *path);
time_t ts_timegm(struct tm *t);
void ts_gmtime(time_t in, struct tm *out);
time_t parse_time_t(const char *string);
struct tm parse_date_time_tm(const char *date, const char *time);
void serialize_time_t(time_t t, char buf[20]);

#endif