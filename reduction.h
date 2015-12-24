/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef REDUCTION_H
#define REDUCTION_H

#include "datafile.h"

int create_flat(const char *pattern, size_t minmax, const char *masterbias, const char *masterdark, const char *outname);
int create_bias(const char *pattern, size_t minmax, double bias_fudge, const char *outname);
int create_dark(const char *pattern, size_t minmax, const char *masterbias, const char *outname);
int display_frame(char *data_path, char *frame_name);
int print_frame_metadata(char *frame_path);

int update_reduction(char *dataPath);
int create_reduction_file(char *filePath);
int update_preview(char *preview_filename, char *ds9_title, char *autoguide_output);
int calculate_bjd(char *date, char *time, char *ra_string, char *dec_string);
int create_ts(char *reference_date, char *reference_time, char **filenames, size_t num_datafiles, char *ts_filename, bool use_ratio, float comparison_magnitude);
int create_jdratio(char *reference_date, char *reference_time, char **filenames, size_t num_datafiles, char *ratio_filename);

int display_tracer(char *dataPath);
int frame_translation(const char *frame, const char *bias, const char *reference, const char *dark_path, const char *flat_path);
int reduce_aperture_range(char *base_name, double min, double max, double step, char *prefix);

#endif
