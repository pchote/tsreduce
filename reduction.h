/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include "datafile.h"

#ifndef REDUCTION_H
#define REDUCTION_H

int generate_photometry_dft_data(datafile *data,
                                 double **raw_time, double **raw, double **mean_sky, size_t *num_raw,
                                 double **time, double **ratio, double **polyfit, double **mma, double **fwhm, size_t *num_filtered,
                                 double **ratio_noise, double **mma_noise,
                                 double *ratio_mean, double *ratio_std, double *mma_mean, double *mma_std,
                                 double **freq, double **ampl, size_t *num_dft);

int create_flat(const char *pattern, int minmax, const char *masterdark, const char *outname);
int create_dark(const char *pattern, int minmax, const char *outname);

int reduce_single_frame(char *framePath, char *darkPath, char *flatPath, char *outPath);
int update_reduction(char *dataPath);
int create_reduction_file(char *filePath);
int update_preview(char *preview_filename, char *ds9_title);
int calculate_bjd(char *date, char *time, char *ra_string, char *dec_string, double epoch);
int create_ts(char *reference_date, char *reference_time, char **filenames, int num_datafiles, char *ts_filename);

#endif
