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
                                 double **raw_time, double **raw, size_t *num_raw,
                                 double **time, double **ratio, double **polyfit, double **mmi, size_t *num_filtered,
                                 double **ratio_noise, double **mmi_noise,
                                 double *ratio_mean, double *ratio_std, double *mmi_mean, double *mmi_std,
                                 double **freq, double **ampl, size_t *num_dft);

int create_flat(const char *pattern, int minmax, const char *masterdark, const char *outname);
int create_dark(const char *pattern, int minmax, const char *outname);

int ccd_readnoise(const char *framePattern);
int reduce_single_frame(char *framePath, char *darkPath, char *flatPath, char *outPath);
int update_reduction(char *dataPath);
int create_reduction_file(char *framePath, char *framePattern, char *darkTemplate, char *flatTemplate, char *filePath);
int create_mmi(char *dataPath);
int evaluate_aperture_snr(char *dataPath, double minAperture, double maxAperture, int numApertures);

#endif
