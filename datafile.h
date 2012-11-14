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
    double fwhm;

    // Pointers inside data array
    double *star;
    double *sky;
    double2 *pos;

    // Dynamically sized by allocating
    // extra space with malloc
    char data[1];
};

typedef struct
{
    int version;
    char *frame_dir;
    char *frame_pattern;
    char *dark_template;
    char *flat_template;
    target *targets;
    int num_targets;
    ts_time reference_time;

    struct observation *obs_start;
    struct observation *obs_end;
    size_t obs_count;
    map_t filename_map;

    int plot_fit_degree;
    double plot_max_raw;
    int plot_num_uhz;
    double plot_min_uhz;
    double plot_max_uhz;

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
    bool has_noise;
    bool has_fwhm;

    double ratio_mean;
    double ratio_std;
    double ratio_snr;
    double mma_mean;
    double mma_std;

    // Raw (unfiltered) photometry
    double *raw_time;
    double *raw;
    double *sky;
    size_t raw_count;

    // Processed photometry
    // Filtered by blocked ranges and NaN ratio
    double *time;
    double *ratio;
    double *ratio_noise;
    double *ratio_fit;
    double *mma;
    double *mma_noise;
    double *fwhm;
    size_t filtered_count;

    // Polynomial fit coefficients to ratio
    double *fit_coeffs;
    size_t fit_coeffs_count;
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
void datafile_free_dft(struct dft_data *data);

#endif
