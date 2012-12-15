/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef DATAFILE_H
#define DATAFILE_H

#include <stdio.h>
#include "aperture.h"
#include "hashmap.h"
#include "helpers.h"

struct observation
{
    struct observation *next;
    struct observation *prev;

    char *filename;
    double time;
    double ratio;
    double ratio_noise;

    // Pointers inside data array
    double *star;
    double *noise;
    double *sky;
    double2 *pos;
    double *fwhm;

    // Dynamically sized by allocating
    // extra space with malloc
    char data[1];
};

struct target_data
{
    aperture aperture;
    char *label;
    double scale;
};

typedef struct
{
    int version;
    char *frame_dir;
    char *frame_pattern;
    char *dark_template;
    char *flat_template;
    struct target_data *targets;
    size_t target_count;
    ts_time reference_time;

    struct observation *obs_start;
    struct observation *obs_end;
    size_t obs_count;
    map_t filename_map;

    uint8_t ratio_fit_degree;
    double plot_max_raw;
    double plot_max_dft;
    int plot_num_uhz;
    double plot_min_uhz;
    double plot_max_uhz;
    bool plot_error_bars;

    double ccd_gain;
    double ccd_readnoise;
    double ccd_platescale;

    char *coord_ra;
    char *coord_dec;
    double coord_epoch;

    double2 *blocked_ranges;
    int num_blocked_ranges;
} datafile;

struct photometry_data
{
    double scaled_target_max;
    double ratio_mean;
    double ratio_std;
    double mma_mean;
    double mma_std;
    double fwhm_mean;
    double fwhm_std;

    // Raw (unfiltered) photometry
    double *target_time;
    double *target_intensity;
    double *target_noise;

    size_t *target_count;
    double *target_snr;

    // Data averaged over all raw photometry
    double *raw_time;
    double *sky;
    double *fwhm;
    size_t raw_count;

    // Processed photometry
    // Filtered by blocked ranges and NaN ratio
    double *time;
    double *ratio;
    double *ratio_noise;
    double *ratio_fit;
    double *mma;
    double *mma_noise;
    size_t filtered_count;

    // Polynomial fit coefficients to ratio
    double *fit_coeffs;
    size_t fit_coeffs_count;

    // Time scale
    double time_offset;
    double time_min;
    double time_max;
    int time_exponent;
    double time_scale;
};

struct dft_data
{
    double *freq;
    double *ampl;
    size_t count;

    double max_ampl;
    double mean_ampl;

    double min_freq;
    double max_freq;
};

datafile *datafile_alloc();
datafile *datafile_load(char *dataFile);
void datafile_free(datafile *data);
struct observation *datafile_new_observation(datafile *data);
void datafile_append_observation(datafile *data, struct observation *obs);
void datafile_discard_observations(datafile *data);
int datafile_save(datafile *data, char *filename);

struct photometry_data *datafile_generate_photometry(datafile *data);
void datafile_free_photometry(struct photometry_data *data);

struct dft_data *datafile_generate_dft(datafile *data, struct photometry_data *pd);
struct dft_data *datafile_generate_window(datafile *data, struct photometry_data *pd, double freq, double range, size_t count);
void datafile_free_dft(struct dft_data *data);

#endif
