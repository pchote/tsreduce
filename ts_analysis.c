/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <cpgplot.h>

#include "reduction.h"
#include "datafile.h"
#include "helpers.h"
#include "fit.h"

int display_tracer(char *dataPath)
{
    int ret;

    // Read file header
    datafile *data = datafile_load(dataPath);
    if (data == NULL)
        return error("Error opening data file");

    if (chdir(data->frame_dir))
        error_jump(setup_error, ret, "Invalid frame path: %s", data->frame_dir);

    if (init_ds9())
        error_jump(setup_error, ret, "Unable to launch ds9");

    if (!data->obs_start)
        error_jump(setup_error, ret, "No observations to display");

    char command[1024];
    snprintf(command, 1024, "xpaset tsreduce file %s/%s", data->frame_dir, data->obs_end->filename);
    ts_exec_write(command, NULL, 0);

    // Set scaling mode
    ts_exec_write("xpaset tsreduce scale mode zscale", NULL, 0);

    // Flip X axis
    ts_exec_write("xpaset tsreduce orient x", NULL, 0);

    // Display results in ds9 - errors are non-fatal
    ts_exec_write("xpaset tsreduce regions delete all", NULL, 0);

    if (data->obs_start->next)
        for (struct observation *obs = data->obs_start->next; obs; obs = obs->next)
            for (int i = 0; i < data->num_targets; i++)
            {
                snprintf(command, 1024, "xpaset tsreduce regions command '{line %f %f %f %f # line= 0 0 color=red select=0}'",
                         obs->prev->pos[i].x + 1, obs->prev->pos[i].y + 1,
                         obs->pos[i].x + 1, obs->pos[i].y + 1);
                ts_exec_write(command, NULL, 0);
            }

    // Draw apertures
    for (int i = 0; i < data->num_targets; i++)
    {
        double2 xy = data->obs_end->pos[i];
        target t = data->targets[i];
        snprintf(command, 1024, "xpaset tsreduce regions command '{circle %f %f %f #color=red select=0}'", xy.x + 1, xy.y + 1, t.r);
        ts_exec_write(command, NULL, 0);

        snprintf(command, 1024, "xpaset tsreduce regions command '{annulus %f %f %f %f #select=0}'",
                 xy.x + 1, xy.y + 1, t.s1, t.s2);
        ts_exec_write(command, NULL, 0);
    }

    ts_exec_write(command, NULL, 0);
    ts_exec_write("xpaset -p tsreduce update now", NULL, 0);

setup_error:
    datafile_free(data);
    return ret;
}

// Output radial profile information for the given targetIndex, obsIndex in
// the reduction file at dataPath
int calculate_profile(char *dataPath, int obsIndex, int targetIndex)
{
    int ret = 0;

    // Read file header
    datafile *data = datafile_load(dataPath);
    if (data == NULL)
        return error("Error opening data file");

    if (obsIndex >= data->obs_count)
        error_jump(setup_error, ret, "Requested observation is out of range: max is %d", data->obs_count-1);

    if (chdir(data->frame_dir))
        error_jump(setup_error, ret, "Invalid frame path: %s", data->frame_dir);

    char *filename;
    if (obsIndex >= 0)
    {
        struct observation *obs = data->obs_start;
        for (size_t i = 0; obs && i < obsIndex; i++, obs = obs->next);
            // Do nothing

        filename = strdup(obs->filename);
    }
    else
        filename = get_first_matching_file(data->frame_pattern);

    if (!filename)
        error_jump(setup_error, ret, "No matching files found");

    framedata *frame = framedata_load(filename);
    if (!frame)
        error_jump(setup_error, ret, "Error loading frame %s", filename);
    subtract_bias(frame);

    framedata *dark = framedata_load(data->dark_template);
    if (!dark)
        error_jump(dark_error, ret, "Error loading frame %s", data->dark_template);
    framedata_subtract(frame, dark);

    framedata *flat = framedata_load(data->flat_template);
    if (!flat)
        error_jump(flat_error, ret, "Error loading frame %s", data->flat_template);
    framedata_divide(frame, flat);

    if (targetIndex < 0 || targetIndex >= data->num_targets)
        error_jump(process_error, ret, "Invalid target `%d' selected", targetIndex);

    target t = data->targets[targetIndex];
    double2 xy;
    if (center_aperture(t, frame, &xy))
        error_jump(process_error, ret, "Aperture centering failed");
    t.x = xy.x; t.y = xy.y;

    double sky_intensity, sky_std_dev;
    if (calculate_background(t, frame, &sky_intensity, &sky_std_dev))
        error_jump(process_error, ret, "Background calculation failed");

    // Calculation lives in its own scope to ensure jumping to the error handling is safe
    // TODO: this is a mess - either tidy this up or remove the function completely
    {
        const int numIntensity = 46;
        double intensity[numIntensity];
        double noise[numIntensity];
        double radii[numIntensity];
        double profile[numIntensity];

        double readnoise, gain;
        if (framedata_get_header_dbl(flat, "CCD-READ", &readnoise))
            readnoise = data->ccd_readnoise;

        if (readnoise <= 0)
            error_jump(process_error, ret, "CCD Read noise unknown. Define CCDReadNoise in %s.", dataPath);

        if (framedata_get_header_dbl(flat, "CCD-GAIN", &gain))
            gain = data->ccd_gain;

        if (gain <= 0)
            error_jump(process_error, ret, "CCD Gain unknown. Define CCDGain in %s.", dataPath);

        printf("# Read noise: %f\n", readnoise);
        printf("# Gain: %f\n", gain);

        // Calculate the remaining integrated intensities
        for (int i = 1; i < numIntensity; i++)
        {
            radii[i] = i/5.0 + 1;
            integrate_aperture_and_noise(xy, radii[i], frame, dark, readnoise, gain, &intensity[i], &noise[i]);
            intensity[i] -= sky_intensity*M_PI*radii[i]*radii[i];
        }

        // Normalize integrated count by area to give an intensity profile
        // r = 0 value is sampled from the central pixel directly
        radii[0] = 0;
        noise[0] = 0;
        intensity[0] = 0;
        profile[0] = frame->data[frame->cols*((int)xy.y) + (int)xy.x];

        // Central integrated value is a disk
        profile[1] = intensity[1]/(M_PI*radii[1]*radii[1]);

        // Remaining areas are annuli
        for (int i = 2; i < numIntensity; i++)
        {
            double area = M_PI*(radii[i]*radii[i] - radii[i-1]*radii[i-1]);
            profile[i] = (intensity[i] - intensity[i-1])/area;
        }

        // Print sky value
        printf("# Sky background: %f\n", sky_intensity);
        printf("# Sky stddev: %f\n", sky_std_dev);

        // Estimate FWHM by linear interpolation between points
        for (int i = 1; i < numIntensity; i++)
            if (profile[i] < profile[0]/2)
            {
                double fwhm = 2*(radii[i - 1] + (radii[i] - radii[i-1])*(profile[0]/2 - profile[i-1])/(profile[i] - profile[i-1]));
                printf("# Estimated FWHM: %f px (%f arcsec)\n", fwhm, fwhm*0.68083798727887213);
                break;
            }

        // Estimate radius that encloses 85%, 90%, 95% intensity
        for (int i = 0; i < numIntensity - 1; i++)
            if (intensity[i + 1] > 0.85*intensity[numIntensity-1])
            {
                printf("# Estimated 85%%: %f\n", i + (0.85*intensity[numIntensity-1] - intensity[i]) / (intensity[i+1] - intensity[i]));
                break;
            }
        for (int i = 0; i < numIntensity - 1; i++)
            if (intensity[i + 1] > 0.90*intensity[numIntensity-1])
            {
                printf("# Estimated 90%%: %f\n", i + (0.90*intensity[numIntensity-1] - intensity[i]) / (intensity[i+1] - intensity[i]));
                break;
            }
        for (int i = 0; i < numIntensity - 1; i++)
            if (intensity[i + 1] > 0.95*intensity[numIntensity-1])
            {
                printf("# Estimated 95%%: %f\n", i + (0.95*intensity[numIntensity-1] - intensity[i]) / (intensity[i+1] - intensity[i]));
                break;
            }

        // Estimate radius where signal reaches 5x,10x sky sigma
        for (int i = 1; i < numIntensity; i++)
            if (profile[i] < 5*sky_std_dev)
            {
                printf("# Estimated 5 sky sigma: %f\n", i - 1 + (5*sky_std_dev - profile[i-1])/(profile[i] - profile[i-1]));
                break;
            }
        for (int i = 1; i < numIntensity; i++)
            if (profile[i] < 10*sky_std_dev)
            {
                printf("# Estimated 10 sky sigma: %f\n", i - 1 + (10*sky_std_dev - profile[i-1])/(profile[i] - profile[i-1]));
                break;
            }

        // Print profile values
        for (int i = 0; i < numIntensity; i++)
            printf("%f %f %f %f\n", radii[i], profile[i], intensity[i], intensity[i]/noise[i]);
    }

process_error:
    framedata_free(flat);
flat_error:
    framedata_free(dark);
dark_error:
    framedata_free(frame);
    free(filename);
setup_error:
    datafile_free(data);
    return ret;
}


// List the timestamps/filenames that are corrupted by continous downloads
int detect_repeats(char *dataPath)
{
    // Read file header
    datafile *data = datafile_load(dataPath);
    if (data == NULL)
        return error("Error opening data file");

    // No data
    if (!data->obs_start)
    {
        datafile_free(data);
        return error("File specifies no observations");
    }

    if (!data->obs_start->next)
    {
        datafile_free(data);
        return error("File requires at least two observations");
    }

    bool bad_range = false;
    for (struct observation *obs = data->obs_start->next; obs; obs = obs->next)
    {
        double time = obs->time;
        if (obs->time - obs->prev->time < 0.1f)
            bad_range = true;

        if (bad_range)
        {
            printf("%s @ %.1f", obs->filename, time);
            for (int j = 0; j < data->num_blocked_ranges; j++)
                if (obs->time >= data->blocked_ranges[j].x && obs->time <= data->blocked_ranges[j].y)
                {
                    printf(" blocked");
                    break;
                }
            printf("\n");
        }

        if (bad_range && obs->time - obs->prev->time >= 0.1f)
            bad_range = false;
    }

    datafile_free(data);
    return 0;
}

int plot_fits_internal(datafile *data, char *tsDevice, double tsSize, char *dftDevice, double dftSize)
{
    int plot_colors_max = 8;
    int plot_colors[] = {4,2,8,3,5,6,7,9};

    struct photometry_data *pd = datafile_generate_photometry(data);
    if (!pd)
        return error("Photometry calculation failed");

    struct dft_data *dd = datafile_generate_dft(data, pd);
    if (!dd)
    {
        datafile_free_photometry(pd);
        return error("DFT calculation failed");
    }

    float min_seconds = pd->raw_time[0];
    float max_seconds = pd->raw_time[pd->raw_count - 1];
    int secexp = (int)(log10(max_seconds) / 3)*3;

    // Time in hours
    float min_time = (float)ts_time_to_utc_hour(data->reference_time) + min_seconds / 3600;
    float max_time = min_time + (max_seconds - min_seconds)/3600;

    double snr_ratio = 0;
    if (data->version >= 5)
    {
        double total_ratio = 0;
        double total_ratio_noise = 0;
        for (int i = 0; i < pd->filtered_count; i++)
        {
            total_ratio += pd->ratio[i];
            total_ratio_noise += pd->ratio_noise[i];
        }
        snr_ratio = total_ratio/total_ratio_noise;
    }

    // Cast the double arrays to float, does not allocate any new memory.
    float *raw_time = cast_double_array_to_float(pd->raw_time, pd->raw_count);
    float *raw = cast_double_array_to_float(pd->raw, pd->raw_count*data->num_targets);
    float *mean_sky = cast_double_array_to_float(pd->sky, pd->raw_count);

    float min_mma = pd->mma_mean - 5*pd->mma_std;
    float max_mma = pd->mma_mean + 5*pd->mma_std;
    float min_ratio = pd->ratio_mean - 5*pd->ratio_std;
    float max_ratio = pd->ratio_mean + 5*pd->ratio_std;

    if (cpgopen(tsDevice) <= 0)
        return error("Unable to open PGPLOT window");
    cpgpap(tsSize, 0.6);

    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(1);

    // Plot blocked ranges
    cpgsvp(0.1, 0.9, 0.075, 0.93);
    cpgswin(min_seconds, max_seconds, 0, 1);
    cpgsci(14);
    for (int j = 0; j < data->num_blocked_ranges; j++)
    {
        cpgmove(data->blocked_ranges[j].x, 0);
        cpgdraw(data->blocked_ranges[j].x, 1);
        cpgmove(data->blocked_ranges[j].y, 0);
        cpgdraw(data->blocked_ranges[j].y, 1);
    }
    cpgsci(1);

    double mean_fwhm = 0;
    if (pd->filtered_count > 0)
    {
        float *time = cast_double_array_to_float(pd->time, pd->filtered_count);
        float *ratio = cast_double_array_to_float(pd->ratio, pd->filtered_count);
        float *ratio_noise = cast_double_array_to_float(pd->ratio_noise, pd->filtered_count);
        float *fwhm = cast_double_array_to_float(pd->fwhm, pd->filtered_count);
        float *polyfit = cast_double_array_to_float(pd->ratio_fit, pd->filtered_count);
        float *mma = cast_double_array_to_float(pd->mma, pd->filtered_count);
        float *mma_noise = cast_double_array_to_float(pd->mma_noise, pd->filtered_count);

        // Fitted MMA
        cpgsvp(0.1, 0.9, 0.8, 0.93);
        cpgsch(1.0);

        // Top axis in UTC Hour
        cpgswin(min_time, max_time, min_mma, max_mma);
        cpgsch(0.7);
        cpgbox("cstm", 1, 4, "bcstnv", 0, 0);

        // Bottom axis in seconds
        cpgswin(min_seconds, max_seconds, min_mma, max_mma);
        cpgbox("bst", 0, 0, "0", 0, 0);
        cpgsch(1.0);

        cpgmtxt("t", 1.75, 0.5, 0.5, "UTC Hour");
        cpgmtxt("l", 2.75, 0.5, 0.5, "mma");

        if (data->version >= 5)
            cpgerrb(6, pd->filtered_count, time, mma, mma_noise, 0.0);
        else
            cpgpt(pd->filtered_count, time, mma, 229);

        // Ratio
        cpgsvp(0.1, 0.9, 0.67, 0.8);
        cpgmtxt("l", 2.75, 0.5, 0.5, "Ratio");

        // Top axis in UTC Hour
        cpgswin(min_time, max_time, min_ratio, max_ratio);
        cpgsch(0.7);
        cpgbox("cst", 1, 4, "bcstnv", 0, 0);

        // Bottom axis in seconds
        cpgswin(min_seconds, max_seconds, min_ratio, max_ratio);
        cpgbox("bst", 0, 0, "0", 0, 0);
        cpgsch(1.0);

        // Plot error bars if ratio_noise is available (data->version >= 5)
        if (pd->has_noise)
            cpgerrb(6, pd->filtered_count, time, ratio, ratio_noise, 0.0);
        else
            cpgpt(pd->filtered_count, time, ratio, 229);

        // Plot the polynomial fit
        cpgsci(2);
        cpgline(pd->filtered_count, time, polyfit);
        cpgsci(1);

        // FWHM
        if (pd->has_fwhm)
        {
            float min_fwhm = FLT_MAX;
            float max_fwhm = -FLT_MAX;
            for (int i = 0; i < pd->filtered_count; i++)
            {
                min_fwhm = fmin(min_fwhm, fwhm[i]);
                max_fwhm = fmax(max_fwhm, fwhm[i]);
                mean_fwhm += fwhm[i]/pd->filtered_count;
            }
            float mid_fwhm = (max_fwhm + min_fwhm)/2;
            float fwhm_range = (max_fwhm - min_fwhm)/2;
            min_fwhm = mid_fwhm - 1.2*fwhm_range;
            max_fwhm = mid_fwhm + 1.2*fwhm_range;

            cpgsvp(0.1, 0.9, 0.54, 0.67);
            cpgmtxt("l", 2.75, 0.5, 0.5, "FWHM (\")");

            // Top axis in UTC Hour
            cpgswin(min_time, max_time, min_fwhm, max_fwhm);
            cpgsch(0.7);
            cpgbox("cst", 1, 4, "bcstnv", 0, 0);

            // Bottom axis in seconds
            cpgswin(min_seconds, max_seconds, min_fwhm, max_fwhm);
            cpgbox("bst", 0, 0, "0", 0, 0);
            cpgsch(1.0);
            cpgpt(pd->filtered_count, time, fwhm, 229);
        }
    }

    // Raw Data
    double max_raw = data->plot_max_raw;

    // Maximum not specified - calculate from data
    if (max_raw == 0)
    {
        for (int i = 0; i < pd->raw_count; i++)
            for (int j = 0; j < data->num_targets; j++)
            {
                double r = raw[j*pd->raw_count + i]*data->targets[j].plot_scale;
                if (r > max_raw)
                    max_raw = r;
            }
        max_raw *= 1.2;
    }

    // Calculate base 10 exponent to reduce label length
    int rawexp = (int)log10(max_raw);

    cpgsvp(0.1, 0.9, 0.075, 0.54);

    // Top axis in UTC Hour
    cpgswin(min_time, max_time, 0, max_raw/pow(10, rawexp));
    cpgsch(0.7);

    cpgbox("cst", 1, 4, "bcstnv", 0, 0);

    // Bottom axis in seconds
    cpgswin(min_seconds/pow(10, secexp), max_seconds/pow(10, secexp), 0, max_raw);
    cpgbox("bstn", 0, 0, "0", 0, 0);
    cpgsch(1.0);

    char label[64];
    snprintf(label, 64, "Count Rate (10\\u%d\\d ADU/s)", rawexp);
    cpgmtxt("l", 2.75, 0.5, 0.5, label);

    if (secexp == 0)
        strcpy(label, "Run Time (s)");
    else
        snprintf(label, 64, "Run Time (10\\u%d\\d s)", secexp);

    cpgmtxt("b", 2.5, 0.5, 0.5, label);

    for (int j = 0; j < data->num_targets; j++)
    {
        cpgswin(min_seconds, max_seconds, 0, max_raw/data->targets[j].plot_scale);
        cpgsci(plot_colors[j%plot_colors_max]);
        cpgpt(pd->raw_count, raw_time, &raw[j*pd->raw_count], 229);
    }

    cpgsci(15);
    cpgswin(min_seconds, max_seconds, 0, max_raw);
    cpgpt(pd->raw_count, raw_time, mean_sky, 229);

    if (pd->filtered_count > 0)
    {
        cpgsvp(0.1, 0.9, 0.015, 0.05);
        cpgswin(0, 1, 0, 1);
        cpgsci(1);
        cpgsch(0.9);
        snprintf(label, 32, "Ratio SNR: %.2f", snr_ratio);
        cpgptxt(0.95, 0, 0, 1.0, label);

        snprintf(label, 32, "Mean FWHM: %.2f\": (%.2fpx)", mean_fwhm, mean_fwhm/data->ccd_platescale);
        cpgptxt(0.05, 0, 0, 0.0, label);
        cpgsch(1.0);
    }

    // Type labels
    cpgsvp(0.1, 0.9, 0.45, 0.55);
    int num_labels = data->num_targets + 1;
    cpgswin(0, num_labels, 0, 1);
    for (int j = 0; j < num_labels; j++)
    {
        cpgsci(plot_colors[j%plot_colors_max]);

        char label[30];
        if (j == num_labels - 1)
        {
            cpgsci(15);
            strcpy(label, "Mean Sky");
        }
        else if (j == 0)
        {
            if (data->targets[j].plot_scale == 1.0)
                strcpy(label, "Target");
            else
                snprintf(label, 30, "%g \\x Target", data->targets[j].plot_scale);
        }
        else
        {
            if (data->targets[j].plot_scale == 1.0)
                snprintf(label, 30, "Comparison %d", j);
            else
                snprintf(label, 30, "%g \\x Comparison %d", data->targets[j].plot_scale, j);
        }
        cpgptxt(j+0.5, 0.5, 0, 0.5, label);
    }
    cpgend();

    if (dd->count > 0)
    {
        // Calculate baseline scale
        char uhzlabel[20];
        double scale = 1e6;
        char *unit = "\\gm";

        if (dd->max_freq > 1.5e4)
        {
            scale = 1.0e-3;
            unit = "k";
        }
        else if (dd->max_freq > 1.5e1)
        {
            scale = 1.0;
            unit = "";
        }
        else if (dd->max_freq > 1.5e-2)
        {
            scale = 1.0e3;
            unit = "m";
        }
        snprintf(uhzlabel, 20, "Frequency (%sHz)", unit);

        char ampl_label[32];
        snprintf(ampl_label, 32, "Mean amplitude: %.2f mma", dd->mean_ampl);

        // Convert DFT data to float arrays for plotting
        cast_double_array_to_float(dd->freq, dd->count);
        cast_double_array_to_float(dd->ampl, dd->count);

        // Draw plot
        if (cpgopen(dftDevice) <= 0)
            return error("Unable to open PGPLOT window");

        cpgpap(dftSize, 0.6);
        cpgask(0);
        cpgslw(1);
        cpgsfs(2);
        cpgscf(1);

        // DFT
        cpgsvp(0.1, 0.9, 0.075, 0.93);
        cpgswin(scale*dd->min_freq, scale*dd->max_freq, 0, 1);
        cpgbox("bcstn", 0, 0, "0", 0, 0);
        cpgswin(dd->min_freq, dd->max_freq, 0, 1.1*dd->max_ampl);
        cpgbox("0", 0, 0, "bcnst", 0, 0);

        cpgsci(12);
        cpgline(dd->count, (float *)dd->freq, (float *)dd->ampl);
        cpgsci(1);

        cpgswin(0, 1, 0, 1);
        cpgptxt(0.97, 0.93, 0, 1.0, ampl_label);
        cpgmtxt("b", 2.5, 0.5, 0.5, uhzlabel);
        cpgmtxt("l", 2, 0.5, 0.5, "Amplitude (mma)");
        cpgend();
    }

    datafile_free_dft(dd);
    datafile_free_photometry(pd);

    return 0;
}

int plot_fits(char *dataPath, char *tsDevice, double tsSize, char *dftDevice, double dftSize)
{
    datafile *data = datafile_load(dataPath);
    if (!data)
        return error("Error opening data file %s", dataPath);

    int ret = plot_fits_internal(data, tsDevice, tsSize, dftDevice, dftSize);
    datafile_free(data);
    return ret;
}

int playback_reduction(char *dataPath, int delay, int step, char *tsDevice, double tsSize, char *dftDevice, double dftSize)
{
    int ret = 0;

    // Count the number of observations
    datafile *data = datafile_load(dataPath);
    if (!data)
        return error("Error opening data file %s", dataPath);

    size_t limit = data->obs_count;
    for (size_t i = 2; i < limit; i+= step)
    {
        // Limit the data to the first N observations
        data->obs_count = i;

        clock_t start = clock();
        if (plot_fits_internal(data, tsDevice, tsSize, dftDevice, dftSize))
            error_jump(plot_error, ret, "Plotting error");

        // Attempt to compensate for calculation time
        int ms = (clock() - start)*1000/CLOCKS_PER_SEC;
        if (ms < delay)
            millisleep(delay - ms);
    }

    // Finish by plotting with all data
    data->obs_count = limit;
    if (plot_fits_internal(data, tsDevice, tsSize, dftDevice, dftSize))
        error_jump(plot_error, ret, "Plotting error");

plot_error:
    datafile_free(data);
    return ret;
}

int reduce_aperture_range(char *base_name, double min, double max, double step, char *prefix)
{
    int ret = 0;
    datafile *data = datafile_load(base_name);
    if (data == NULL)
        return error("Error opening data file");

    datafile_discard_observations(data);

    // Force at least version 5 to include noise calculation
    if (data->version < 5)
        data->version = 5;

    char *dir = getcwd(NULL, 0);
    // Create and update a datafile for each aperture
    double radius = min;
    do
    {
        for (int i = 0; i < data->num_targets; i++)
            data->targets[i].r = radius;

        size_t filename_len = strlen(prefix) + 11;
        char *filename = malloc(filename_len*sizeof(char));
        snprintf(filename, filename_len, "%s-%0.2f.dat", prefix, radius);

        chdir(dir);
        // Errors are non-fatal -> proceeed to the next file
        if (!datafile_save(data, filename))
            update_reduction(filename);

        free(filename);
        radius += step;
    } while (radius < max);

    free(dir);
    datafile_free(data);
    return ret;
}

int plot_range(char *datafile_pattern)
{
    char **datafile_names;
    size_t num_files = get_matching_files(datafile_pattern, &datafile_names);
    if (num_files < 0)
        return error("Error matching pattern %s", datafile_pattern);

    if (num_files == 0)
    {
        free(datafile_names);
        return error("No datafiles found matching pattern %s", datafile_pattern);
    }
    int i = 0;
    float x,y;
    char c;
    while (true)
    {
        plot_fits(datafile_names[i], "5/xs", 9.41, "6/xs", 9.41);

        if (cpgopen("9/xs") <= 0)
        {
            error("Unable to open PGPLOT window");
            break;
        }

        cpgask(0);
        cpgpap(3, 0.1);
        cpgsch(20);
        cpgswin(0, 1, 0, 1);
        cpgsci(1);
        cpgptxt(0.5, 0.7, 0, 0.5, datafile_names[i]);
        cpgptxt(0.5, 0.0, 0, 0.5, "Next file: [n] Prev file: [p] Quit: [q]");

        cpgcurs(&x,&y,&c);
        cpgend();

        switch (c)
		{
			case 'n':
				if (i < num_files - 1) i++;
				break;
			case 'p':
				if (i > 0) i--;
				break;
			case 'q':
				goto endLoop;
			default:
				break;
		}
    }

endLoop:
    free_2d_array(datafile_names, num_files);
    return 0;
}

/*
 * Report the start time, data length, and number of observations in a data file
 */
int report_time(char *dataPath)
{
    // Read file header
    datafile *data = datafile_load(dataPath);
    if (data == NULL)
        return error("Error opening data file %s", dataPath);

    double *raw_time, *raw, *mean_sky, *time, *ratio, *fwhm, *polyfit, *mma, *ratio_noise, *mma_noise;
    double ratio_mean, ratio_std, mma_mean, mma_std;
    size_t num_raw, num_filtered;
    if (generate_photometry_dft_data(data,
                                     &raw_time, &raw, &mean_sky, &num_raw,
                                     &time, &ratio, &polyfit, &mma, &fwhm, &num_filtered,
                                     &ratio_noise, &mma_noise,
                                     &ratio_mean, &ratio_std, &mma_mean, &mma_std,
                                     NULL, NULL, NULL))
    {
        datafile_free(data);
        return error("Error generating data");
    }

    char datetimebuf[24];
    serialize_time(data->reference_time, datetimebuf);
    printf("%s %.2f %zu\n", datetimebuf, (time[num_filtered-1] - time[0])/3600, num_filtered);

    free(raw_time);
    free(raw);
    free(mean_sky);
    free(time);
    free(ratio);
    free(polyfit);
    free(mma);
    free(fwhm);
    datafile_free(data);
    return 0;
}

