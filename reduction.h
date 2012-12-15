/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef REDUCTION_H
#define REDUCTION_H

#include "datafile.h"

int create_flat(const char *pattern, size_t minmax, const char *masterdark, const char *outname);
int create_dark(const char *pattern, size_t minmax, const char *outname);
int display_frame(char *data_path, char *frame_name);
int print_frame_metadata(char *frame_path);

int update_reduction(char *dataPath);
int create_reduction_file(char *filePath);
int update_preview(char *preview_filename, char *ds9_title);
int calculate_bjd(char *date, char *time, char *ra_string, char *dec_string, double epoch);
int create_ts(char *reference_date, char *reference_time, char **filenames, size_t num_datafiles, char *ts_filename);

int display_tracer(char *dataPath);
int calculate_profile(char *dataPath, int obsIndex, int targetIndex);
int detect_repeats(char *dataPath);
int reduce_aperture_range(char *base_name, double min, double max, double step, char *prefix);

#endif
