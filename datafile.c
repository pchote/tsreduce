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

#define CUR_DATAFILE_VERSION 7
#define MIN_DATAFILE_VERSION 7

#define RATIO_FIT_DEGREE_DEFAULT 2
#define PLOT_MAX_RAW_DEFAULT 0
#define PLOT_NUM_UHZ_DEFAULT 1000
#define PLOT_MIN_UHZ_DEFAULT 0
#define PLOT_MAX_UHZ_DEFAULT 10000
#define PLOT_ERROR_BARS_DEFAULT false
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
    dp->version = CUR_DATAFILE_VERSION;
    dp->ratio_fit_degree = RATIO_FIT_DEGREE_DEFAULT;
    dp->plot_max_raw = PLOT_MAX_RAW_DEFAULT;
    dp->plot_num_uhz = PLOT_NUM_UHZ_DEFAULT;
    dp->plot_min_uhz = PLOT_MIN_UHZ_DEFAULT;
    dp->plot_max_uhz = PLOT_MAX_UHZ_DEFAULT;
    dp->plot_error_bars = PLOT_ERROR_BARS_DEFAULT;
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

    dp->targets = malloc(target_count*sizeof(struct target_data));
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
                 dp->target_count < target_count)
        {
            aperture *t = &dp->targets[dp->target_count].aperture;
            sscanf(linebuf, "# Target: (%lf, %lf, %lf, %lf, %lf) [%lf]\n",
                   &t->x, &t->y, &t->r, &t->s1, &t->s2,
                   &dp->targets[dp->target_count].scale);
            dp->target_count++;
        }
        else if (!strncmp(linebuf,"# ReferenceTime:", 16))
            dp->reference_time = parse_time(linebuf+17);

        else if (!strncmp(linebuf,"# Version:", 10))
        {
            sscanf(linebuf, "# Version: %d\n", &dp->version);
            if (dp->version < MIN_DATAFILE_VERSION)
            {
                fprintf(stderr, "Obsolete file format (Version %d).\n", dp->version);
                fprintf(stderr, "Re-reduce data in Version %d or later to continue.\n", MIN_DATAFILE_VERSION);

                // TODO Cleanup resources properly
                return NULL;
            }
        }
        else if (!strncmp(linebuf,"# RatioFitDegree:", 17))
            sscanf(linebuf, "# RatioFitDegree: %hhu\n", &dp->ratio_fit_degree);
        else if (!strncmp(linebuf,"# PlotMaxRaw:", 13))
            sscanf(linebuf, "# PlotMaxRaw: %lf\n", &dp->plot_max_raw);
        else if (!strncmp(linebuf,"# PlotMinUhz:", 13))
            sscanf(linebuf, "# PlotMinUhz: %lf\n", &dp->plot_min_uhz);
        else if (!strncmp(linebuf,"# PlotErrorBars:", 16))
        {
            uint8_t temp;
            sscanf(linebuf, "# PlotErrorBars: %hhu\n", &temp);
            dp->plot_error_bars = temp;
        }
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

        // Note: strtok_r is unsupported under windows
        char *token;

        // Filename
        if (!(token = strtok(linebuf, " ")))
            goto parse_error;
        obs->filename = strdup(token);

        // Time
        if (!(token = strtok(NULL, " ")))
            goto parse_error;
        obs->time = atof(token);

        // Target intensity / sky / aperture x / aperture y
        for (size_t i = 0; i < dp->target_count; i++)
        {
            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->star[i] = atof(token);

            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->noise[i] = atof(token);

            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->sky[i] = atof(token);

            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->pos[i].x = atof(token);

            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->pos[i].y = atof(token);

            if (!(token = strtok(NULL, " ")))
                goto parse_error;
            obs->fwhm[i] = atof(token);
        }

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
    size_t star_size = data->target_count*sizeof(double);
    size_t noise_size = data->target_count*sizeof(double);
    size_t sky_size = data->target_count*sizeof(double);
    size_t pos_size = data->target_count*sizeof(double2);
    size_t fwhm_size = data->target_count*sizeof(double);
    size_t s = sizeof(struct observation) - 1 + star_size + noise_size + sky_size + pos_size + fwhm_size;

    struct observation *obs = calloc(1, s);
    if (!obs)
        return NULL;

    obs->star = (double *)obs->data;
    obs->noise = (double *)(obs->star + data->target_count);
    obs->sky = (double *)(obs->noise + data->target_count);
    obs->pos = (double2 *)(obs->sky + data->target_count);
    obs->fwhm = (double *)(obs->pos + data->target_count);

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
    fprintf(out, "### tsreduce reduction data\n");
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
    if (data->ratio_fit_degree != RATIO_FIT_DEGREE_DEFAULT)
        fprintf(out, "# RatioFitDegree: %d\n", data->ratio_fit_degree);
    if (data->plot_max_raw != PLOT_MAX_RAW_DEFAULT)
        fprintf(out, "# PlotMaxRaw: %g\n", data->plot_max_raw);
    if (data->plot_min_uhz != PLOT_MIN_UHZ_DEFAULT)
        fprintf(out, "# PlotMinUhz: %g\n", data->plot_min_uhz);
    if (data->plot_max_uhz != PLOT_MAX_UHZ_DEFAULT)
        fprintf(out, "# PlotMaxUhz: %g\n", data->plot_max_uhz);
    if (data->ccd_gain != CCD_GAIN_DEFAULT)
        fprintf(out, "# CCDGain: %g\n", data->ccd_gain);
    if (data->ccd_readnoise != CCD_READNOISE_DEFAULT)
        fprintf(out, "# CCDReadNoise: %g\n", data->ccd_readnoise);
    if (data->ccd_platescale != CCD_PLATESCALE_DEFAULT)
        fprintf(out, "# CCDPlateScale: %g\n", data->ccd_platescale);
    if (data->plot_error_bars != PLOT_ERROR_BARS_DEFAULT)
        fprintf(out, "# PlotErrorBars: %d\n", data->plot_error_bars);

    if (data->reference_time.time)
    {
        char datetimebuf[24];
        serialize_time(data->reference_time, datetimebuf);
        fprintf(out, "# ReferenceTime: %s\n", datetimebuf);
    }

    fprintf(out, "### (x, y, Radius, Inner Sky Radius, Outer Sky Radius) [plot scale]\n");
    for (size_t i = 0; i < data->target_count; i++)
    {
        aperture *t = &data->targets[i].aperture;
        fprintf(out, "# Target: (%6.2f, %6.2f, %6.2f, %6.2f, %6.2f) [%.2f]\n",
                t->x, t->y, t->r, t->s1, t->s2,
                data->targets[i].scale);
    }
    if (data->num_blocked_ranges > 0)
        fprintf(out, "### (Start (s), End (s))\n");
    for (size_t i = 0; i < data->num_blocked_ranges; i++)
        fprintf(out, "# BlockRange: (%g, %g)\n",
                data->blocked_ranges[i].x, data->blocked_ranges[i].y);

    fprintf(out, "### Filename,          Time  ");
    for (size_t i = 0; i < data->target_count; i++)
        fprintf(out, " |  Star    Noise     Sky     x      y     FWHM");
    fprintf(out, "\n");

    fprintf(out, "###                     (s)  ");
    for (size_t i = 0; i < data->target_count; i++)
        fprintf(out, " | (ADU/s) (ADU/s)  (ADU/s)  (px)   (px)   (px)");
    fprintf(out, "\n");

    struct observation *obs;
    for (obs = data->obs_start; obs; obs = obs->next)
    {
        fprintf(out, "%s ", obs->filename);
        fprintf(out, "%9.3f ", obs->time);
        for (size_t i = 0; i < data->target_count; i++)
        {
            fprintf(out, "%10.2f ", obs->star[i]);
            fprintf(out, "%6.2f ", obs->noise[i]);
            fprintf(out, "%9.2f ", obs->sky[i]);
            fprintf(out, "%6.2f %6.2f ", obs->pos[i].x, obs->pos[i].y);
            fprintf(out, "%5.2f", obs->fwhm[i]);
        }
        fprintf(out, "\n");

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

    struct photometry_data *p = calloc(1, sizeof(struct photometry_data));
    if (!p)
        return NULL;

    p->target_time = calloc(data->obs_count*data->target_count, sizeof(double));
    p->target_intensity = calloc(data->obs_count*data->target_count, sizeof(double));
    p->target_noise = calloc(data->obs_count*data->target_count, sizeof(double));

    p->target_count = calloc(data->target_count, sizeof(size_t));
    p->target_snr = calloc(data->target_count, sizeof(double));

    p->raw_time = calloc(data->obs_count, sizeof(double));
    p->sky = calloc(data->obs_count, sizeof(double));
    p->fwhm = calloc(data->obs_count, sizeof(double));

    p->time = calloc(data->obs_count, sizeof(double));
    p->ratio = calloc(data->obs_count, sizeof(double));
    p->ratio_noise = calloc(data->obs_count, sizeof(double));
    p->ratio_fit = calloc(data->obs_count, sizeof(double));
    p->mma = calloc(data->obs_count, sizeof(double));
    p->mma_noise = calloc(data->obs_count, sizeof(double));

    p->fit_coeffs_count = data->ratio_fit_degree + 1;
    p->fit_coeffs = calloc(p->fit_coeffs_count, sizeof(double));

    if (!p->target_time || !p->target_intensity || !p->target_noise ||
        !p->target_count || !p->target_snr ||
        !p->raw_time || !p->sky || !p->fwhm ||
        !p->time || !p->ratio || !p->ratio_noise ||
        !p->ratio_fit || !p->mma || !p->mma_noise ||
        !p->fit_coeffs)
    {
        datafile_free_photometry(p);
        error("Allocation error");
        return NULL;
    }

    //
    // Parse raw data
    //

    for (size_t i = 0; i < data->target_count; i++)
    {
        p->target_count[i] = 0;
        p->target_snr[i] = 0;
    }

    p->scaled_target_max = 0;
    p->filtered_count = 0;
    p->ratio_mean = 0;
    p->fwhm_mean = 0;

    // External code may modify obs_count to restrict data processing,
    // so both checks are required
    struct observation *obs = data->obs_start;
    for (size_t i = 0; obs && i < data->obs_count; obs = obs->next, i++)
    {
        p->raw_time[p->raw_count] = obs->time;
        p->sky[p->raw_count] = 0;
        p->fwhm[p->raw_count] = 0;
        size_t target_count = 0;

        double comparison_intensity = 0;
        double comparison_noise = 0;

        for (size_t j = 0; j < data->target_count; j++)
        {
            if (isnan(obs->star[j]) || isnan(obs->noise[j]) || isnan(obs->fwhm[j]))
                continue;

            size_t k = j*data->obs_count + p->target_count[j];
            p->target_time[k] = obs->time;
            p->target_intensity[k] = obs->star[j];
            p->target_noise[k] = obs->noise[j];

            p->target_count[j]++;
            p->target_snr[j] += obs->star[j] / obs->noise[j];

            p->sky[p->raw_count] += obs->sky[j];
            p->fwhm[p->raw_count] += obs->fwhm[j];
            target_count++;


            if (j > 0)
            {
                comparison_intensity += obs->star[j];
                comparison_noise += obs->noise[j];
            }

            double r = obs->star[j]*data->targets[j].scale;
            if (r > p->scaled_target_max)
                p->scaled_target_max = r;
        }

        if (target_count > 0)
        {
            p->sky[p->raw_count] /= target_count;
            p->fwhm[p->raw_count] /= target_count;
            p->fwhm_mean += p->fwhm[p->raw_count];
            p->raw_count++;
        }

        // Cannot calculate ratio if we've lost one or more targets
        // (each contributes a factor proportional to its relative intensity)
        if (target_count != data->target_count)
            continue;

        bool skip = false;
        for (size_t j = 0; j < data->num_blocked_ranges; j++)
            if (obs->time >= data->blocked_ranges[j].x && obs->time <= data->blocked_ranges[j].y)
            {
                skip = true;
                break;
            }

        if (skip)
            continue;

        p->time[p->filtered_count] = obs->time;
        p->ratio[p->filtered_count] = obs->star[0] / comparison_intensity;
        p->ratio_mean += p->ratio[p->filtered_count];

        p->ratio_noise[p->filtered_count] = (obs->noise[0]/obs->star[0] + comparison_noise/comparison_intensity)*p->ratio[p->filtered_count];

        p->filtered_count++;
    }

    for (size_t i = 0; i < data->target_count; i++)
        p->target_snr[i] /= p->target_count[i];

    p->ratio_mean /= p->filtered_count;
    p->fwhm_mean /= p->raw_count;

    // Ratio and fwhm standard deviation
    p->ratio_std = 0;
    for (size_t i = 0; i < p->filtered_count; i++)
    {
        p->ratio_std += (p->ratio[i] - p->ratio_mean)*(p->ratio[i] - p->ratio_mean);
        p->fwhm_std += (p->fwhm[i] - p->fwhm_mean)*(p->fwhm[i] - p->fwhm_mean);
    }
    p->ratio_std = sqrt(p->ratio_std/p->filtered_count);
    p->fwhm_std = sqrt(p->fwhm_std/p->filtered_count);

    //
    // Calculate polynomial fit
    //
    if (p->filtered_count < p->fit_coeffs_count)
    {
        datafile_free_photometry(p);
        error("Insufficient data for polynomial fit");
        return NULL;
    }

    if (fit_polynomial(p->time, p->ratio, p->ratio_noise, p->filtered_count, p->fit_coeffs, data->ratio_fit_degree))
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

        double numer_error = fabs(p->ratio_noise[i]/(p->ratio[i] - p->ratio_fit[i]));
        double denom_error = fabs(p->ratio_noise[i]/p->ratio[i]);
        p->mma_noise[i] = (numer_error + denom_error)*fabs(p->mma[i]);

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

void datafile_free_photometry(struct photometry_data *p)
{
    free(p->target_time);
    free(p->target_intensity);
    free(p->target_noise);

    free(p->target_count);
    free(p->target_snr);

    free(p->raw_time);
    free(p->sky);
    free(p->fwhm);

    free(p->time);
    free(p->ratio);
    free(p->ratio_noise);
    free(p->ratio_fit);
    free(p->mma);
    free(p->mma_noise);
    free(p->fit_coeffs);

    free(p);
}

struct dft_data *datafile_generate_dft(datafile *data, struct photometry_data *pd)
{
    struct dft_data *d = calloc(1, sizeof(struct dft_data));
    if (!d)
        return NULL;

    d->min_freq = data->plot_min_uhz*1e-6;
    d->max_freq = data->plot_max_uhz*1e-6;
    d->count = data->plot_num_uhz;
    d->freq = calloc(d->count, sizeof(double));
    d->ampl = calloc(d->count, sizeof(double));

    if (!d->freq || !d->count)
    {
        datafile_free_dft(d);
        error("Allocation error");
        return NULL;
    }

    calculate_amplitude_spectrum(pd->time, pd->mma, pd->filtered_count,
                                 d->min_freq, d->max_freq,
                                 d->freq, d->ampl, d->count);

    d->max_ampl = 0;
    d->mean_ampl = 0;
    for (size_t i = 0; i < d->count; i++)
    {
        d->max_ampl = fmax(d->max_ampl, d->ampl[i]);
        d->mean_ampl += d->ampl[i];
    }
    d->mean_ampl /= d->count;

    return d;
}

void datafile_free_dft(struct dft_data *data)
{
    free(data->freq);
    free(data->ampl);
    free(data);
}