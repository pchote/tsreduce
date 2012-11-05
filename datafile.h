/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include "aperture.h"
#include "hashmap.h"
#include "helpers.h"

#ifndef DATAFILE_H
#define DATAFILE_H

// Maximum number of frames to load when creating a flat or dark frame
#define MAX_FRAMES 100

// Maximum number of targets to track
#define MAX_TARGETS 10

// Maximum number of ranges to block
#define MAX_BLOCKED_RANGES 10

struct observation
{
    double star[MAX_TARGETS];
    double sky[MAX_TARGETS];
    double2 pos[MAX_TARGETS];
    double time;
    double ratio;
    double ratio_noise;
    double fwhm;
    char filename[64];

    struct observation *next;
    struct observation *prev;
};

typedef struct
{
    FILE *file;
    int version;
    char *frame_dir;
    char *frame_pattern;
    char *dark_template;
    char *flat_template;
    target targets[MAX_TARGETS];
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

    double2 blocked_ranges[MAX_BLOCKED_RANGES];
    int num_blocked_ranges;
} datafile;

datafile *datafile_alloc();
datafile *datafile_load(char *dataFile);
void datafile_free(datafile *data);
void datafile_append_observation(datafile *data, struct observation *obs);
int datafile_save_header(datafile *data, char *filename);

#endif
