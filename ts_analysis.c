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
    int ret = 0;
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
            for (size_t i = 0; i < data->target_count; i++)
            {
                snprintf(command, 1024, "xpaset tsreduce regions command '{line %f %f %f %f # line= 0 0 color=red select=0}'",
                         obs->prev->pos[i].x + 1, obs->prev->pos[i].y + 1,
                         obs->pos[i].x + 1, obs->pos[i].y + 1);
                ts_exec_write(command, NULL, 0);
            }

    // Draw apertures
    for (size_t i = 0; i < data->target_count; i++)
    {
        double2 xy = data->obs_end->pos[i];
        aperture *t = &data->targets[i].aperture;
        snprintf(command, 1024, "xpaset tsreduce regions command '{circle %f %f %f #color=red select=0}'", xy.x + 1, xy.y + 1, t->r);
        ts_exec_write(command, NULL, 0);

        snprintf(command, 1024, "xpaset tsreduce regions command '{annulus %f %f %f %f #select=0}'",
                 xy.x + 1, xy.y + 1, t->s1, t->s2);
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

    if (targetIndex < 0 || targetIndex >= data->target_count)
        error_jump(process_error, ret, "Invalid target `%d' selected", targetIndex);

    aperture a = data->targets[targetIndex].aperture;
    double2 xy;
    if (center_aperture(a, frame, &xy))
        error_jump(process_error, ret, "Aperture centering failed");

    a.x = xy.x;
    a.y = xy.y;

    double sky_intensity, sky_std_dev;
    if (calculate_background(a, frame, &sky_intensity, &sky_std_dev))
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

static int plot_internal(datafile *data, const char *tsDevice, double tsSize, const char *dftDevice, double dftSize)
{
    size_t plot_colors_max = 8;
    uint8_t plot_colors[] = {4,2,8,3,5,6,7,9};

    struct photometry_data *pd = datafile_generate_photometry(data);
    if (!pd)
        return error("Photometry calculation failed");

    struct dft_data *dd = datafile_generate_dft(data, pd);
    if (!dd)
    {
        datafile_free_photometry(pd);
        return error("DFT calculation failed");
    }

    double min_seconds = pd->raw_time[0];
    double max_seconds = pd->raw_time[pd->raw_count - 1];
    int secexp = (int)(log10(max_seconds) / 3)*3;
    double sec_scale = 1.0/pow(10, secexp);

    // Time in hours
    double min_time = ts_time_to_utc_hour(data->reference_time) + min_seconds / 3600;
    double max_time = min_time + (max_seconds - min_seconds)/3600;

    if (cpgopen(tsDevice) <= 0)
        return error("Unable to open PGPLOT window");

    cpgpap(tsSize, 0.6);
    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(1);

    //
    // Plot blocked ranges
    //
    cpgsvp(0.1, 0.9, 0.075, 0.93);
    cpgswin(min_seconds, max_seconds, 0, 1);
    cpgsci(14);
    for (size_t j = 0; j < data->num_blocked_ranges; j++)
    {
        cpgmove(data->blocked_ranges[j].x, 0);
        cpgdraw(data->blocked_ranges[j].x, 1);
        cpgmove(data->blocked_ranges[j].y, 0);
        cpgdraw(data->blocked_ranges[j].y, 1);
    }
    cpgsci(1);

    // Generic label buffer for passing data to pgplot
    char label[64];
    const size_t label_len = 64;

    //
    // Plot raw data
    //
    cast_double_array_to_float(pd->target_time, data->obs_count*data->target_count);
    cast_double_array_to_float(pd->target_intensity, data->obs_count*data->target_count);
    cast_double_array_to_float(pd->target_noise, data->obs_count*data->target_count);

    cast_double_array_to_float(pd->raw_time, pd->raw_count);
    cast_double_array_to_float(pd->sky, pd->raw_count);
    cast_double_array_to_float(pd->fwhm, pd->raw_count);
    {
        double min_raw = 0;
        double max_raw = data->plot_max_raw ? data->plot_max_raw : 1.25*pd->scaled_target_max;
        int raw_exp = (int)log10(max_raw);
        double raw_scale = 1.0/pow(10, raw_exp);

        cpgsvp(0.1, 0.9, 0.075, 0.54);

        // Plot top axis markers in UTC hour, bottom axis markers in seconds
        cpgsch(0.7);
        cpgswin(min_time, max_time, raw_scale*min_raw, raw_scale*max_raw);
        cpgbox("cst", 1, 4, "bcstnv", 0, 0);
        cpgswin(sec_scale*min_seconds, sec_scale*max_seconds, min_raw, max_raw);
        cpgbox("bstn", 0, 0, "0", 0, 0);
        cpgsch(1.0);

        snprintf(label, label_len, "Count Rate (10\\u%d\\d ADU/s)", raw_exp);
        cpgmtxt("l", 2.75, 0.5, 0.5, label);

        if (secexp == 0)
            strncpy(label, "Run Time (s)", label_len);
        else
            snprintf(label, label_len, "Run Time (10\\u%d\\d s)", secexp);
        cpgmtxt("b", 2.5, 0.5, 0.5, label);

        // Raw intensities
        for (size_t j = 0; j < data->target_count; j++)
        {
            cpgswin(min_seconds, max_seconds, min_raw, max_raw/data->targets[j].scale);
            cpgsci(plot_colors[j%plot_colors_max]);

            size_t k = j*data->obs_count;
            float *time = &((float *)pd->target_time)[k];
            float *intensity =  &((float *)pd->target_intensity)[k];
            float *noise =  &((float *)pd->target_noise)[k];

            if (data->plot_error_bars)
                cpgerrb(6, pd->target_count[j], time, intensity, noise, 0.0);
            else
                cpgpt(pd->target_count[j], time, intensity, 229);
        }

        // Mean sky intensity
        cpgswin(min_seconds, max_seconds, min_raw, max_raw);
        cpgsci(15);
        cpgpt(pd->raw_count, (float *)pd->raw_time, (float *)pd->sky, 229);

        // Labels
        cpgsvp(0.1, 0.9, 0.45, 0.55);
        cpgswin(0, data->target_count + 1, 0, 1);
        for (size_t j = 0; j <= data->target_count; j++)
        {
            cpgsci(plot_colors[j%plot_colors_max]);

            if (j == data->target_count)
            {
                cpgsci(15);
                strncpy(label, "Mean Sky", label_len);
            }
            else
            {
                if (data->targets[j].scale == 1.0)
                    strncpy(label, data->targets[j].label, label_len);
                else
                    snprintf(label, label_len, "%g \\x %s", data->targets[j].scale, data->targets[j].label);
            }

            cpgptxt(j+0.5, 0.5, 0, 0.5, label);
            if (j < data->target_count)
            {
                snprintf(label, label_len, "SNR: %.0f", pd->target_snr[j]);
                cpgsch(0.8);
                cpgptxt(j+0.5, 0.25, 0, 0.5, label);
                cpgsch(1.0);
            }
        }
        cpgsci(1);
    }

    //
    // Plot FWHM
    //
    {
        double min_fwhm = pd->fwhm_mean - 5*pd->fwhm_std;
        double max_fwhm = pd->fwhm_mean + 5*pd->fwhm_std;

        cpgsvp(0.1, 0.9, 0.54, 0.67);
        cpgmtxt("l", 2.75, 0.5, 0.5, "FWHM (\")");

        // Plot top axis markers in UTC hour, bottom axis markers in seconds
        cpgsch(0.7);
        cpgswin(min_time, max_time, min_fwhm*data->ccd_platescale, max_fwhm*data->ccd_platescale);
        cpgbox("cst", 1, 4, "bcstnv", 0, 0);
        cpgswin(sec_scale*min_seconds, sec_scale*max_seconds, min_fwhm, max_fwhm);
        cpgbox("bst", 0, 0, "0", 0, 0);
        cpgsch(1.0);

        cpgswin(min_seconds, max_seconds, min_fwhm, max_fwhm);
        cpgpt(pd->raw_count, (float *)pd->raw_time, (float *)pd->fwhm, 229);
    }

    //
    // Plot Ratio
    //
    cast_double_array_to_float(pd->time, pd->filtered_count);
    {
        double min_ratio = pd->ratio_mean - 5*pd->ratio_std;
        double max_ratio = pd->ratio_mean + 5*pd->ratio_std;
        cast_double_array_to_float(pd->ratio, pd->filtered_count);
        cast_double_array_to_float(pd->ratio_noise, pd->filtered_count);
        cast_double_array_to_float(pd->ratio_fit, pd->filtered_count);

        cpgsvp(0.1, 0.9, 0.67, 0.8);
        cpgmtxt("l", 2.75, 0.5, 0.5, "Ratio");

        // Plot top axis markers in UTC hour, bottom axis markers in seconds
        cpgsch(0.7);
        cpgswin(min_time, max_time, min_ratio, max_ratio);
        cpgbox("cst", 1, 4, "bcstnv", 0, 0);
        cpgswin(sec_scale*min_seconds, sec_scale*max_seconds, min_ratio, max_ratio);
        cpgbox("bst", 0, 0, "0", 0, 0);
        cpgsch(1.0);

        cpgswin(min_seconds, max_seconds, min_ratio, max_ratio);
        if (data->plot_error_bars)
            cpgerrb(6, pd->filtered_count, (float *)pd->time, (float *)pd->ratio, (float *)pd->ratio_noise, 0.0);
        else
            cpgpt(pd->filtered_count, (float *)pd->time, (float *)pd->ratio, 229);

        // Plot the polynomial fit
        cpgsci(2);

        // Don't plot the fit through blocked regions
        float *t = (float *)pd->time;
        float *t_end = &t[pd->filtered_count-1];
        float *f = (float *)pd->ratio_fit;
        do
        {
            double min_time = t[0];
            double max_time = max_seconds;

            // Find the first blocked range that affects the data
            for (size_t j = 0; j < data->num_blocked_ranges; j++)
            {
                double2 r = data->blocked_ranges[j];
                if (min_time > r.x && min_time < r.y)
                    min_time = r.y;

                if (max_time > r.x && min_time < r.y)
                    max_time = r.x;
            }

            // Find the first point to draw
            while (t < t_end && *t < min_time)
            {
                t++;
                f++;
            }

            if (t == t_end)
                break;

            // Find the last point to draw
            size_t count = 0;
            while (t + count < t_end && *(t + count) < max_time)
                count++;

            cpgline(count, t, f);

            // Start next search from end of this section
            t += count + 1;
            f += count + 1;
        } while (t + 1 < t_end);

        cpgsci(1);
    }

    //
    // Plot mma
    //
    {
        double min_mma = pd->mma_mean - 5*pd->mma_std;
        double max_mma = pd->mma_mean + 5*pd->mma_std;
        cast_double_array_to_float(pd->mma, pd->filtered_count);
        cast_double_array_to_float(pd->mma_noise, pd->filtered_count);

        cpgsvp(0.1, 0.9, 0.8, 0.93);
        cpgsch(1.0);

        cpgmtxt("t", 1.75, 0.5, 0.5, "UTC Hour");
        cpgmtxt("l", 2.75, 0.5, 0.5, "mma");

        // Plot top axis markers in UTC hour, bottom axis markers in seconds
        cpgsch(0.7);
        cpgswin(min_time, max_time, min_mma, max_mma);
        cpgbox("cstm", 1, 4, "bcstnv", 0, 0);
        cpgswin(sec_scale*min_seconds, sec_scale*max_seconds, min_mma, max_mma);
        cpgbox("bst", 0, 0, "0", 0, 0);
        cpgsch(1.0);
        
        cpgswin(min_seconds, max_seconds, min_mma, max_mma);
        if (data->plot_error_bars)
            cpgerrb(6, pd->filtered_count, (float *)pd->time, (float *)pd->mma, (float *)pd->mma_noise, 0.0);
        else
            cpgpt(pd->filtered_count, (float *)pd->time, (float *)pd->mma, 229);
    }

    //
    // Plot mean FWHM label
    //
    {
        cpgsvp(0.1, 0.9, 0.015, 0.05);
        cpgswin(0, 1, 0, 1);
        cpgsci(1);
        cpgsch(0.9);

        snprintf(label, 32, "Mean FWHM: %.2f\": (%.2fpx)", pd->fwhm_mean*data->ccd_platescale, pd->fwhm_mean);
        cpgptxt(0.05, 0, 0, 0.0, label);

        cpgsch(1.0);
    }
    cpgend();

    //
    // Plot DFT
    //
    if (cpgopen(dftDevice) <= 0)
        return error("Unable to open PGPLOT window");

    cpgpap(dftSize, 0.6);
    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(1);
    // Calculate baseline scale
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

    // DFT
    cpgsvp(0.1, 0.9, 0.075, 0.93);
    cpgswin(scale*dd->min_freq, scale*dd->max_freq, 0, 1);
    cpgbox("bcstn", 0, 0, "0", 0, 0);
    cpgswin(dd->min_freq, dd->max_freq, 0, 1.1*dd->max_ampl);
    cpgbox("0", 0, 0, "bcnst", 0, 0);

    cpgsci(12);
    cast_double_array_to_float(dd->freq, dd->count);
    cast_double_array_to_float(dd->ampl, dd->count);
    cpgline(dd->count, (float *)dd->freq, (float *)dd->ampl);
    cpgsci(1);

    // Calculate median intensity
    qsort(dd->ampl, dd->count, sizeof(float), compare_float);

    snprintf(label, label_len, "Frequency (%sHz)", unit);
    cpgmtxt("b", 2.5, 0.5, 0.5, label);

    cpgmtxt("l", 2, 0.5, 0.5, "Amplitude (mma)");

    cpgswin(0, 1, 0, 1);
    snprintf(label, label_len, "Mean: %.2f mma", dd->mean_ampl);
    cpgptxt(0.97, 0.94, 0, 1.0, label);
    snprintf(label, label_len, "Median: %.2f mma", ((float *)dd->ampl)[dd->count/2]);
    cpgptxt(0.97, 0.90, 0, 1.0, label);

    cpgend();

    datafile_free_dft(dd);
    datafile_free_photometry(pd);

    return 0;
}

int online_plot(char *dataPath, char *tsDevice, double tsSize, char *dftDevice, double dftSize)
{
    datafile *data = datafile_load(dataPath);
    if (!data)
        return error("Error opening data file %s", dataPath);

    int ret = plot_internal(data, tsDevice, tsSize, dftDevice, dftSize);
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
    for (size_t i = data->ratio_fit_degree + 1; i < limit; i += step)
    {
        // Limit the data to the first N observations
        data->obs_count = i;

        clock_t start = clock();
        if (plot_internal(data, tsDevice, tsSize, dftDevice, dftSize))
            error_jump(plot_error, ret, "Plotting error");

        // Attempt to compensate for calculation time
        int ms = (clock() - start)*1000/CLOCKS_PER_SEC;
        if (ms < delay)
            millisleep(delay - ms);
    }

    // Finish by plotting with all data
    data->obs_count = limit;
    if (plot_internal(data, tsDevice, tsSize, dftDevice, dftSize))
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

    char *dir = getcwd(NULL, 0);
    // Create and update a datafile for each aperture
    double radius = min;
    do
    {
        for (size_t i = 0; i < data->target_count; i++)
            data->targets[i].aperture.r = radius;

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
        online_plot(datafile_names[i], "5/xs", 9.41, "6/xs", 9.41);

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

