/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#include <stdarg.h>

#ifndef HELPERS_H
#define HELPERS_H

#define error_jump(label, ret, ...) do { ret = error(__VA_ARGS__); goto label; } while(0)

void free_2d_array(char **array, int len);
float *cast_double_array_to_float(double *d_ptr, size_t count);
int get_matching_files(const char *pattern, char ***files);
int get_first_matching_file(char *pattern, char *filenamebuf, int buflen);
int compare_double(const void *a, const void *b);
int error(const char * format, ...);
void die(const char * format, ...) __attribute__ ((noreturn));
int init_ds9(char *);
int tell_ds9(char *title, char *command, void *data, int dataSize);
int ask_ds9(char *title, char *command, char **outbuf);

void calculate_amplitude_spectrum(double fmin, double fmax, double *t, double *data, int numData, double *outFreq, double *outAmpl, int numOut);
void calculate_amplitude_spectrum_float(float fmin, float fmax, float *time, float *data, int numData, float *outFreq, float *outAmpl, int numOut);
#endif