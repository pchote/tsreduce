/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "datafile.h"
#include "fit.h"
#include "helpers.h"

#define PLOT_FIT_DEGREE_DEFAULT 2
#define PLOT_MAX_RAW_DEFAULT 0
#define PLOT_NUM_UHZ_DEFAULT 1000
#define PLOT_MIN_UHZ_DEFAULT 0
#define PLOT_MAX_UHZ_DEFAULT 10000
#define CCD_GAIN_DEFAULT 0
#define CCD_READNOISE_DEFAULT 0
#define CCD_PLATESCALE_DEFAULT 1.0

extern int verbosity;

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
    FILE *input = fopen(filename, "r+");
    if (!input)
        return NULL;

    char linebuf[1024];
    char stringbuf[1024];

    // Count number of targets and blocked ranges
    size_t target_count = 0;
    size_t block_count = 0;
    while (fgets(linebuf, 1024, input) != NULL)
    {
        if (!strncmp(linebuf,"# Target:", 9))
            target_count++;
        else if (!strncmp(linebuf,"# BlockRange:", 13))
            block_count++;
    }

    dp->targets = malloc(target_count*sizeof(target));
    dp->blocked_ranges = malloc(block_count*sizeof(double2));
    if (!dp->targets || !dp->blocked_ranges)
    {
        fclose(input);
        return NULL;
    }

    rewind(input);
    while (fgets(linebuf, 1024, input) != NULL)
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
        else if (!strncmp(linebuf,"# Target:", 9) &&
                 dp->num_targets < target_count)
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
        else if (!strncmp(linebuf,"# BlockRange:", 13) &&
                 dp->num_blocked_ranges < block_count)
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

        struct observation *obs = datafile_new_observation(dp);
        if (!obs)
            return NULL;

        // Time
        // Note: strtok_r is unsupported under windows
        char *token;

        if (!(token = strtok(linebuf, " ")))
            goto parse_error;
        obs->time = atof(token);

        // Target intensity / sky / aperture x / aperture y
        for (size_t i = 0; i < dp->num_targets; i++)
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
        obs->filename = strdup(token);

        size_t len = strlen(obs->filename);
        if (len < 1)
            goto parse_error;

        // Strip newline
        obs->filename[len-1] = '\0';

        datafile_append_observation(dp, obs);

        continue;
    parse_error:
        error("Skipping malformed observations line");
        free(obs);
    }

    fclose(input);
    return dp;
}

void datafile_free(datafile *data)
{
    free(data->frame_dir);
    free(data->frame_pattern);
    free(data->dark_template);
    free(data->flat_template);
    free(data->coord_ra);
    free(data->coord_dec);

    datafile_discard_observations(data);

    hashmap_free(data->filename_map);

    free(data);
}

struct observation *datafile_new_observation(datafile *data)
{
    size_t star_size = data->num_targets*sizeof(double);
    size_t sky_size = data->num_targets*sizeof(double);
    size_t pos_size = data->num_targets*sizeof(double2);
    size_t s = sizeof(struct observation) - 1 + star_size + sky_size + pos_size;

    struct observation *obs = calloc(1, s);
    if (!obs)
        return NULL;

    obs->star = (double *)obs->data;
    obs->sky = (double *)(obs->star + data->num_targets);
    obs->pos = (double2 *)(obs->sky + data->num_targets);

    return obs;
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

void datafile_discard_observations(datafile *data)
{
    struct observation *obs, *next;
    for (obs = data->obs_start; obs; obs = next)
    {
        next = obs->next;
        free(obs->filename);
        free(obs);
    }
    data->obs_start = NULL;
    data->obs_end = NULL;
}

int datafile_save(datafile *data, char *filename)
{
    // data->file may point at a different file, so
    // create a new output
    FILE *out = fopen(filename, "w");
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

    for (size_t i = 0; i < data->num_targets; i++)
        fprintf(out, "# Target: (%f, %f, %f, %f, %f) [1.0]\n",
                data->targets[i].x, data->targets[i].y,
                data->targets[i].r, data->targets[i].s1, data->targets[i].s2);

    for (size_t i = 0; i < data->num_blocked_ranges; i++)
        fprintf(out, "# BlockRange: (%f, %f)\n",
                data->blocked_ranges[i].x, data->blocked_ranges[i].y);

    struct observation *obs;
    for (obs = data->obs_start; obs; obs = obs->next)
    {
        fprintf(out, "%.1f ", obs->time);
        for (size_t i = 0; i < data->num_targets; i++)
        {
            fprintf(out, "%.2f ", obs->star[i]);
            fprintf(out, "%.2f ", obs->sky[i]);
            fprintf(out, "%.2f %.2f ", obs->pos[i].x, obs->pos[i].y);
        }

        fprintf(out, "%.5e ", obs->ratio);
        if (data->version >= 5)
            fprintf(out, "%.5e ", obs->ratio_noise);

        if (data->version >= 6)
            fprintf(out, "%.3f ", obs->fwhm*data->ccd_platescale);

        fprintf(out, "%s\n", obs->filename);
    }

    fclose(out);
    return 0;
}

struct photometry_data *datafile_generate_photometry(datafile *data)
{
    // TODO: Expose this as a parameter
    double mma_filter_sigma = 3;

    if (!data->obs_start)
        return NULL;

    struct photometry_data *p = malloc(sizeof(struct photometry_data));
    if (!p)
        return NULL;

    p->raw_time = calloc(data->obs_count, sizeof(double));
    p->raw = calloc(data->obs_count*data->num_targets, sizeof(double));
    p->sky = calloc(data->obs_count, sizeof(double));

    p->time = calloc(data->obs_count, sizeof(double));
    p->ratio = calloc(data->obs_count, sizeof(double));
    p->ratio_noise = calloc(data->obs_count, sizeof(double));
    p->ratio_fit = calloc(data->obs_count, sizeof(double));
    p->mma = calloc(data->obs_count, sizeof(double));
    p->mma_noise = calloc(data->obs_count, sizeof(double));

    p->fwhm = calloc(data->obs_count, sizeof(double));

    p->fit_coeffs_count = data->plot_fit_degree + 1;
    p->fit_coeffs = calloc(p->fit_coeffs_count, sizeof(double));

    if (!p->raw_time || !p->raw || !p->sky ||
        !p->time || !p->ratio || !p->ratio_noise ||
        !p->ratio_fit || !p->mma || !p->mma_noise ||
        !p->fwhm || !p->fit_coeffs)
    {
        datafile_free_photometry(p);
        error("Allocation error");
        return NULL;
    }

    //
    // Parse raw data
    //

    p->has_noise = (data->version >= 5 && data->dark_template);
    p->has_fwhm = (data->version >= 6);
    p->raw_count = data->obs_count;
    p->filtered_count = 0;
    p->ratio_mean = 0;

    // External code may modify obs_count to restrict data processing,
    // so both checks are required
    struct observation *obs = data->obs_start;
    for (size_t i = 0; obs && i < data->obs_count; obs = obs->next, i++)
    {
        p->raw_time[i] = obs->time;
        p->sky[i] = 0;

        for (size_t j = 0; j < data->num_targets; j++)
        {
            p->raw[j*data->obs_count + i] = obs->star[j];
            p->sky[i] += obs->sky[j]/data->num_targets;
        }

        // Filter bad observations
        bool skip = false;
        for (size_t j = 0; j < data->num_blocked_ranges; j++)
            if (p->raw_time[i] >= data->blocked_ranges[j].x && p->raw_time[i] <= data->blocked_ranges[j].y)
            {
                skip = true;
                break;
            }

        // Invalid observations have noise = nan
        if (skip || isnan(obs->ratio_noise))
            continue;

        p->time[p->filtered_count] = obs->time;
        p->ratio[p->filtered_count] = obs->ratio;
        p->ratio_mean += obs->ratio;

        // Read noise and fwhm from data file if available
        if (p->has_noise)
            p->ratio_noise[p->filtered_count] = obs->ratio_noise;

        if (p->has_fwhm)
            p->fwhm[p->filtered_count] = obs->fwhm;

        p->filtered_count++;
    }

    p->ratio_mean /= p->filtered_count;

    // Ratio standard deviation
    p->ratio_std = 0;
    for (size_t i = 0; i < p->filtered_count; i++)
        p->ratio_std += (p->ratio[i] - p->ratio_mean)*(p->ratio[i] - p->ratio_mean);
    p->ratio_std = sqrt(p->ratio_std/p->filtered_count);

    //
    // Calculate polynomial fit
    //
    if (p->filtered_count < p->fit_coeffs_count)
    {
        datafile_free_photometry(p);
        error("Insufficient data for polynomial fit");
        return NULL;
    }

    if (fit_polynomial(p->time, p->ratio, p->has_noise ? p->ratio_noise : NULL, p->filtered_count, p->fit_coeffs, data->plot_fit_degree))
    {
        datafile_free_photometry(p);
        error("Polynomial fit failed");
        return NULL;
    }

    //
    // Calculate mma
    //

    p->mma_mean = 0;
    for (size_t i = 0; i < p->filtered_count; i++)
    {
        // Subtract polynomial fit and convert to mma
        p->ratio_fit[i] = 0;
        double pow = 1;
        for (size_t j = 0; j < p->fit_coeffs_count; j++)
        {
            p->ratio_fit[i] += pow*p->fit_coeffs[j];
            pow *= p->time[i];
        }
        p->mma[i] = 1000*(p->ratio[i] - p->ratio_fit[i])/p->ratio_fit[i];

        if (p->has_noise)
        {
            double numer_error = fabs(p->ratio_noise[i]/(p->ratio[i] - p->ratio_fit[i]));
            double denom_error = fabs(p->ratio_noise[i]/p->ratio[i]);
            p->mma_noise[i] = (numer_error + denom_error)*fabs(p->mma[i]);
        }

        p->mma_mean += p->mma[i];
    }
    p->mma_mean /= p->filtered_count;

    // mma standard deviation
    p->mma_std = 0;
    for (size_t i = 0; i < p->filtered_count; i++)
        p->mma_std += (p->mma[i] - p->mma_mean)*(p->mma[i] - p->mma_mean);
    p->mma_std = sqrt(p->mma_std/p->filtered_count);

    double mma_corrected_mean = 0;
    size_t mma_corrected_count = 0;

    // Discard outliers and recalculate mean
    for (size_t i = 0; i < p->filtered_count; i++)
    {
        if (fabs(p->mma[i] - p->mma_mean) > mma_filter_sigma*p->mma_std)
        {
            if (verbosity >= 1)
                error("%f is an outlier, setting to 0", p->time[i]);
            p->mma[i] = 0;
        }
        else
        {
            mma_corrected_mean += p->mma[i];
            mma_corrected_count++;
        }
    }

    mma_corrected_mean /= mma_corrected_count;
    p->mma_mean = mma_corrected_mean;

    return p;
}

void datafile_free_photometry(struct photometry_data *data)
{
    free(data->raw_time);
    free(data->raw);
    free(data->sky);
    free(data->time);
    free(data->ratio);
    free(data->ratio_noise);
    free(data->ratio_fit);
    free(data->mma);
    free(data->mma_noise);
    free(data->fwhm);
    free(data->fit_coeffs);
    free(data);
}