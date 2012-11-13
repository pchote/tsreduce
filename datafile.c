/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "datafile.h"
#include "helpers.h"

#define PLOT_FIT_DEGREE_DEFAULT 2
#define PLOT_MAX_RAW_DEFAULT 0
#define PLOT_NUM_UHZ_DEFAULT 1000
#define PLOT_MIN_UHZ_DEFAULT 0
#define PLOT_MAX_UHZ_DEFAULT 10000
#define CCD_GAIN_DEFAULT 0
#define CCD_READNOISE_DEFAULT 0
#define CCD_PLATESCALE_DEFAULT 1.0

/*
 * Allocate a datafile on the heap and set default values
 */
datafile *datafile_alloc()
{
    datafile *dp = calloc(1, sizeof(datafile));
    dp->version = 6;
    dp->plot_fit_degree = PLOT_FIT_DEGREE_DEFAULT;
    dp->plot_max_raw = PLOT_MAX_RAW_DEFAULT;
    dp->plot_num_uhz = PLOT_NUM_UHZ_DEFAULT;
    dp->plot_min_uhz = PLOT_MIN_UHZ_DEFAULT;
    dp->plot_max_uhz = PLOT_MAX_UHZ_DEFAULT;
    dp->ccd_gain = CCD_GAIN_DEFAULT;
    dp->ccd_readnoise = CCD_READNOISE_DEFAULT;
    dp->ccd_platescale = CCD_PLATESCALE_DEFAULT;

    dp->filename_map = hashmap_new();
    return dp;
}

/*
 * Allocate a new datafile and load with data from the specified file
 */
datafile* datafile_load(char *filename)
{
    datafile *dp = datafile_alloc();

    // Open the data file (created with `tsreduce init`)
    dp->file = fopen(filename, "r+");
    if (dp->file == NULL)
        return NULL;

    char linebuf[1024];
    char stringbuf[1024];
    while (fgets(linebuf, 1024, dp->file) != NULL)
    {
        //
        // Header
        //
        if (!strncmp(linebuf,"# FramePattern:", 15))
        {
            sscanf(linebuf, "# FramePattern: %1024s\n", stringbuf);
            dp->frame_pattern = strdup(stringbuf);
        }
        else if (!strncmp(linebuf,"# FrameDir:", 11))
        {
            sscanf(linebuf, "# FrameDir: %1024s\n", stringbuf);
            dp->frame_dir = strdup(stringbuf);
        }
        else if (!strncmp(linebuf,"# DarkTemplate:", 15))
        {
            sscanf(linebuf, "# DarkTemplate: %1024s\n", stringbuf);
            dp->dark_template = strdup(stringbuf);
        }
        else if (!strncmp(linebuf,"# FlatTemplate:", 15))
        {
            sscanf(linebuf, "# FlatTemplate: %1024s\n", stringbuf);
            dp->flat_template = strdup(stringbuf);
        }
        else if (!strncmp(linebuf,"# Target:", 9))
        {
            if (dp->version == 3)
            {
                sscanf(linebuf, "# Target: (%lf, %lf, %lf, %lf, %lf)\n",
                       &dp->targets[dp->num_targets].x,
                       &dp->targets[dp->num_targets].y,
                       &dp->targets[dp->num_targets].r,
                       &dp->targets[dp->num_targets].s1,
                       &dp->targets[dp->num_targets].s2);
                dp->targets[dp->num_targets].plot_scale = 1;
                dp->num_targets++;
            }
            else
            {
                sscanf(linebuf, "# Target: (%lf, %lf, %lf, %lf, %lf) [%lf]\n",
                       &dp->targets[dp->num_targets].x,
                       &dp->targets[dp->num_targets].y,
                       &dp->targets[dp->num_targets].r,
                       &dp->targets[dp->num_targets].s1,
                       &dp->targets[dp->num_targets].s2,
                       &dp->targets[dp->num_targets].plot_scale);
                dp->num_targets++;
            }
        }
        else if (!strncmp(linebuf,"# ReferenceTime:", 16))
            dp->reference_time = parse_time(linebuf+17);

        else if (!strncmp(linebuf,"# Version:", 10))
            sscanf(linebuf, "# Version: %d\n", &dp->version);
        else if (!strncmp(linebuf,"# PlotFitDegree:", 16))
            sscanf(linebuf, "# PlotFitDegree: %d\n", &dp->plot_fit_degree);
        else if (!strncmp(linebuf,"# PlotMaxRaw:", 13))
            sscanf(linebuf, "# PlotMaxRaw: %lf\n", &dp->plot_max_raw);
        else if (!strncmp(linebuf,"# PlotMinUhz:", 13))
            sscanf(linebuf, "# PlotMinUhz: %lf\n", &dp->plot_min_uhz);
        else if (!strncmp(linebuf,"# PlotMaxUhz:", 13))
            sscanf(linebuf, "# PlotMaxUhz: %lf\n", &dp->plot_max_uhz);
        else if (!strncmp(linebuf,"# PlotNumUhz:", 13))
            sscanf(linebuf, "# PlotNumUhz: %d\n", &dp->plot_num_uhz);
        else if (!strncmp(linebuf,"# CCDGain:", 10))
            sscanf(linebuf, "# CCDGain: %lf\n", &dp->ccd_gain);
        else if (!strncmp(linebuf,"# CCDReadNoise:", 15))
            sscanf(linebuf, "# CCDReadNoise: %lf\n", &dp->ccd_readnoise);
        else if (!strncmp(linebuf,"# CCDPlateScale:", 16))
            sscanf(linebuf, "# CCDPlateScale: %lf\n", &dp->ccd_platescale);
        else if (!strncmp(linebuf,"# RA:", 5))
        {
            sscanf(linebuf, "# RA: %1024s\n", stringbuf);
            dp->coord_ra = strdup(stringbuf);
        }
        else if (!strncmp(linebuf,"# DEC:", 6))
        {
            sscanf(linebuf, "# DEC: %1024s\n", stringbuf);
            dp->coord_dec = strdup(stringbuf);
        }
        else if (!strncmp(linebuf,"# Epoch:", 8))
            sscanf(linebuf, "# Epoch: %lf\n", &dp->coord_epoch);
        else if (!strncmp(linebuf,"# BlockRange:", 13))
        {
            sscanf(linebuf, "# BlockRange: (%lf, %lf)\n",
                   &dp->blocked_ranges[dp->num_blocked_ranges].x,
                   &dp->blocked_ranges[dp->num_blocked_ranges].y);
            dp->num_blocked_ranges++;
        }

        // Skip header / comment lines
        if (linebuf[0] == '#')
            continue;

        // Skip empty lines
        if (linebuf[0] == '\n')
            continue;

        //
        // Observations
        //

        struct observation *obs = malloc(sizeof(struct observation));
        if (!obs)
            return NULL;

        // TODO: Make this safe against malformed lines

        // Time
        // Note: strtok_r is unsupported under windows
        char *token;

        if (!(token = strtok(linebuf, " ")))
            goto parse_error;
        obs->time = atof(token);

        // Target intensity / sky / aperture x / aperture y
        for (int i = 0; i < dp->num_targets; i++)
        {
            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->star[i] = atof(token);

            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->sky[i] = atof(token);

            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->pos[i].x = atof(token);

            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->pos[i].y = atof(token);
        }

        if (!(token = strtok(NULL, " ")))
            goto parse_error;
        obs->ratio = atof(token);

        if (dp->version >= 5)
        {
            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->ratio_noise = atof(token);
        }

        if (dp->version >= 6)
        {
            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->fwhm = atof(token);
        }

        // Filename
        if (!(token = strtok(NULL, " ")))
            goto parse_error;
        strncpy(obs->filename, token, sizeof(obs->filename));

        // Strip newline
        obs->filename[strlen(obs->filename)-1] = '\0';

        datafile_append_observation(dp, obs);

        continue;
    parse_error:
        error("Skipping malformed observations line");
        free(obs);
    }
    return dp;
}

void datafile_free(datafile *data)
{
    if (data->file)
        fclose(data->file);

    free(data->frame_dir);
    free(data->frame_pattern);
    free(data->dark_template);
    free(data->flat_template);
    free(data->coord_ra);
    free(data->coord_dec);

    struct observation *obs, *next;
    for (obs = data->obs_start; obs; obs = next)
    {
        next = obs->next;
        free(obs);
    }

    hashmap_free(data->filename_map);

    free(data);
}

void datafile_append_observation(datafile *data, struct observation *obs)
{
    obs->prev = data->obs_end;
    obs->next = NULL;
    if (data->obs_end)
        data->obs_end->next = obs;
    else
        data->obs_start = obs;

    data->obs_end = obs;
    data->obs_count++;

    hashmap_put(data->filename_map, obs->filename, obs);
}

/*
 * Save the header (but not data) from a struct datafile
 */
int datafile_save_header(datafile *data, char *filename)
{
    // data->file may point at a different file, so
    // create a new output
    FILE *out = fopen(filename, "w+");
    if (!out)
        return error("Error opening file %s", filename);

    // Write the file
    fprintf(out, "# Puoko-nui Online reduction output\n");
    if (data->version)
        fprintf(out, "# Version: %d\n", data->version);
    if (data->frame_dir)
        fprintf(out, "# FrameDir: %s\n", data->frame_dir);
    if (data->frame_pattern)
        fprintf(out, "# FramePattern: %s\n", data->frame_pattern);
    if (data->dark_template)
        fprintf(out, "# DarkTemplate: %s\n", data->dark_template);
    if (data->flat_template)
        fprintf(out, "# FlatTemplate: %s\n", data->flat_template);
    if (data->plot_fit_degree != PLOT_FIT_DEGREE_DEFAULT)
        fprintf(out, "# PlotFitDegree: %d\n", data->plot_fit_degree);
    if (data->plot_max_raw != PLOT_MAX_RAW_DEFAULT)
        fprintf(out, "# PlotMaxRaw: %f\n", data->plot_max_raw);
    if (data->plot_min_uhz != PLOT_MIN_UHZ_DEFAULT)
        fprintf(out, "# PlotMinUhz: %f\n", data->plot_min_uhz);
    if (data->plot_max_uhz != PLOT_MAX_UHZ_DEFAULT)
        fprintf(out, "# PlotMaxUhz: %f\n", data->plot_max_uhz);
    if (data->ccd_gain != CCD_GAIN_DEFAULT)
        fprintf(out, "# CCDGain: %f\n", data->ccd_gain);
    if (data->ccd_readnoise != CCD_READNOISE_DEFAULT)
        fprintf(out, "# CCDReadNoise: %f\n", data->ccd_readnoise);
    if (data->ccd_platescale != CCD_PLATESCALE_DEFAULT)
        fprintf(out, "# CCDPlateScale: %f\n", data->ccd_platescale);

    if (data->reference_time.time)
    {
        char datetimebuf[24];
        serialize_time(data->reference_time, datetimebuf);
        fprintf(out, "# ReferenceTime: %s\n", datetimebuf);
    }

    for (int i = 0; i < data->num_targets; i++)
        fprintf(out, "# Target: (%f, %f, %f, %f, %f) [1.0]\n",
                data->targets[i].x, data->targets[i].y,
                data->targets[i].r, data->targets[i].s1, data->targets[i].s2);

    for (int i = 0; i < data->num_blocked_ranges; i++)
        fprintf(out, "# BlockRange: (%f, %f)\n",
                data->blocked_ranges[i].x, data->blocked_ranges[i].y);

    fclose(out);
    return 0;
}
