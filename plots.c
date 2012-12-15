/*
 * Copyright 2010, 2011, 2012 Paul Chote
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

#include "plots.h"
#include "datafile.h"
#include "helpers.h"
#include "fit.h"

static size_t plot_colors_max = 8;
static uint8_t plot_colors[] = {2,8,3,12,5,6,7,9};

static void plot_time_axes(float x1, float x2, float y1, float y2,
                           datafile *data, struct photometry_data *pd)
{
    cpgsvp(x1, x2, y1, y2);
    cpgsch(0.9);
    cpgswin(pd->time_offset + pd->time_min, pd->time_offset + pd->time_max, 0, 1);
    cpgtbox("cstmXYZH", 0, 0, "0", 0, 0);
    cpgswin(pd->time_scale*pd->time_min, pd->time_scale*pd->time_max, 0, 1);
    cpgbox("bstn", 0, 0, "0", 0, 0);
    cpgsch(1.0);

    cpgmtxt("t", 2.0, 0.5, 0.5, "UTC Time");

    char label[64];
    if (pd->time_exponent == 0)
        strncpy(label, "Run Time (s)", 64);
    else
        snprintf(label, 64, "Run Time (10\\u%d\\d s)", pd->time_exponent);

    cpgmtxt("b", 2.5, 0.5, 0.5, label);
}

static int plot_raw_panel(float x1, float x2, float y1, float y2,
                          datafile *data, struct photometry_data *pd)
{
    int ret = 0;
    char label[64];
    const size_t label_len = 64;

    // PGPLOT requires float arrays
    float *time = malloc(data->obs_count*data->target_count*sizeof(float));
    float *intensity = malloc(data->obs_count*data->target_count*sizeof(float));
    float *noise = malloc(data->obs_count*data->target_count*sizeof(float));
    if (!time || !intensity || !noise)
        error_jump(target_allocation_error, ret, "Allocation error");

    for (size_t i = 0; i < data->obs_count*data->target_count; i++)
    {
        time[i] = pd->target_time[i];
        intensity[i] = pd->target_intensity[i];
        noise[i] = pd->target_noise[i];
    }

    float *sky_time = malloc(pd->raw_count*sizeof(float));
    float *sky = malloc(pd->raw_count*sizeof(float));
    if (!sky_time || !sky)
        error_jump(sky_allocation_error, ret, "Allocation error");

    for (size_t i = 0; i < pd->raw_count; i++)
    {
        sky_time[i] = pd->raw_time[i];
        sky[i] = pd->sky[i];
    }

    double min_raw = 0;
    double max_raw = data->plot_max_raw ? data->plot_max_raw : 1.25*pd->scaled_target_max;
    int raw_exp = (int)log10(max_raw);
    double raw_scale = 1.0/pow(10, raw_exp);

    cpgsvp(x1, x2, y1, y2);

    // Plot top axis markers in UTC hour, bottom axis markers in seconds
    cpgsch(0.9);
    cpgswin(pd->time_offset + pd->time_min, pd->time_offset + pd->time_max, raw_scale*min_raw, raw_scale*max_raw);
    cpgtbox("cstZ", 0, 0, "bcstnv", 0, 0);
    cpgswin(pd->time_scale*pd->time_min, pd->time_scale*pd->time_max, min_raw, max_raw);
    cpgbox("bst", 0, 0, "0", 0, 0);
    cpgsch(1.0);

    snprintf(label, label_len, "Count Rate (10\\u%d\\d ADU/s)", raw_exp);
    cpgmtxt("l", 2.75, 0.5, 0.5, label);

    // Raw intensities
    for (size_t j = 0; j < data->target_count; j++)
    {
        cpgswin(pd->time_min, pd->time_max, min_raw, max_raw/data->targets[j].scale);
        cpgsci(plot_colors[j%plot_colors_max]);

        size_t k = j*data->obs_count;
        if (data->plot_error_bars)
            cpgerrb(6, pd->target_count[j], &time[k], &intensity[k], &noise[k], 0.0);
        else
            cpgpt(pd->target_count[j], &time[k], &intensity[k], 229);
    }

    // Mean sky intensity
    cpgswin(pd->time_min, pd->time_max, min_raw, max_raw);
    cpgsci(15);
    cpgpt(pd->raw_count, sky_time, sky, 229);

    // Labels
    cpgsvp(x1, x2, 0.82*y2, y2);
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
            cpgsch(0.9);
            cpgptxt(j+0.5, 0.25, 0, 0.5, label);
            cpgsch(1.0);
        }
    }
    cpgsci(1);

sky_allocation_error:
    free(sky_time);
    free(sky);

target_allocation_error:
    free(time);
    free(intensity);
    free(noise);

    return ret;
}

static int plot_fwhm_panel(float x1, float x2, float y1, float y2,
                           datafile *data, struct photometry_data *pd)
{
    int ret = 0;

    // PGPLOT requires float arrays
    float *time = malloc(pd->raw_count*sizeof(float));
    float *fwhm = malloc(pd->raw_count*sizeof(float));
    if (!time || !fwhm)
        error_jump(allocation_error, ret, "Allocation error");

    for (size_t i = 0; i < pd->raw_count; i++)
    {
        time[i] = pd->raw_time[i];
        fwhm[i] = pd->fwhm[i];
    }

    double min_fwhm = pd->fwhm_mean - 5*pd->fwhm_std;
    double max_fwhm = pd->fwhm_mean + 7.5*pd->fwhm_std;

    cpgsvp(x1, x2, y1, y2);

    // Plot top axis markers in UTC hour, bottom axis markers in seconds
    // Reserve the top 25% of the plot for the mean FWHM label
    cpgmtxt("l", 2.75, 0.5, 0.5, "FWHM (\")");

    cpgsch(0.9);
    cpgswin(pd->time_offset + pd->time_min, pd->time_offset + pd->time_max, min_fwhm*data->ccd_platescale, max_fwhm*data->ccd_platescale);
    cpgtbox("cstZ", 0, 0, "bcstnv", 0, 0);
    cpgswin(pd->time_scale*pd->time_min, pd->time_scale*pd->time_max, min_fwhm, max_fwhm);
    cpgbox("bst", 0, 0, "0", 0, 0);
    cpgsch(1.0);

    cpgswin(pd->time_min, pd->time_max, min_fwhm, max_fwhm);
    cpgpt(pd->raw_count, time, fwhm, 229);
    cpgsci(1);

allocation_error:
    free(time);
    free(fwhm);

    return ret;
}

static int plot_ratio_panel(float x1, float x2, float y1, float y2,
                            datafile *data, struct photometry_data *pd)
{
    int ret = 0;

    // PGPLOT requires float arrays
    float *time = malloc(pd->filtered_count*sizeof(float));
    float *ratio = malloc(pd->filtered_count*sizeof(float));
    float *noise = malloc(pd->filtered_count*sizeof(float));
    float *fit = malloc(pd->filtered_count*sizeof(float));
    if (!time || !ratio || !noise || !fit)
        error_jump(allocation_error, ret, "Allocation error");

    for (size_t i = 0; i < pd->filtered_count; i++)
    {
        time[i] = pd->time[i];
        ratio[i] = pd->ratio[i];
        noise[i] = pd->ratio_noise[i];
        fit[i] = pd->ratio_fit[i];
    }

    double min_ratio = pd->ratio_mean - 5*pd->ratio_std;
    double max_ratio = pd->ratio_mean + 5*pd->ratio_std;

    cpgsvp(x1, x2, y1, y2);
    cpgmtxt("l", 2.75, 0.5, 0.5, "Ratio");

    // Plot top axis markers in UTC hour, bottom axis markers in seconds
    cpgsch(0.9);
    cpgswin(pd->time_offset + pd->time_min, pd->time_offset + pd->time_max, min_ratio, max_ratio);
    cpgtbox("cstZ", 0, 0, "bc", 0, 0);
    cpgswin(pd->time_scale*pd->time_min, pd->time_scale*pd->time_max, min_ratio, max_ratio);
    cpgbox("bst", 0, 0, "0", 0, 0);
    cpgsch(1.0);

    cpgswin(pd->time_min, pd->time_max, min_ratio, max_ratio);
    if (data->plot_error_bars)
        cpgerrb(6, pd->filtered_count, time, ratio, noise, 0.0);
    else
        cpgpt(pd->filtered_count, time, ratio, 229);

    // Plot the polynomial fit
    cpgsci(2);

    // Don't plot the fit through blocked regions
    float *t = time;
    float *t_end = &t[pd->filtered_count];
    float *f = fit;
    do
    {
        double min_time = t[0];
        double max_time = pd->time_max;

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
        while (t + count < t_end && *(t + count) <= max_time)
            count++;

        cpgline(count, t, f);

        // Start next search from end of this section
        t += count + 1;
        f += count + 1;
    } while (t + 1 < t_end);

    cpgsci(1);

allocation_error:
    free(time);
    free(ratio);
    free(noise);
    free(fit);

    return ret;
}

static int plot_mma_panel(float x1, float x2, float y1, float y2,
                          datafile *data, struct photometry_data *pd)
{
    int ret = 0;

    // PGPLOT requires float arrays
    float *time = malloc(pd->filtered_count*sizeof(float));
    float *mma = malloc(pd->filtered_count*sizeof(float));
    float *mma_noise = malloc(pd->filtered_count*sizeof(float));
    if (!time || !mma || !mma_noise)
        error_jump(allocation_error, ret, "Allocation error");

    for (size_t i = 0; i < pd->filtered_count; i++)
    {
        time[i] = pd->time[i];
        mma[i] = pd->mma[i];
        mma_noise[i] = pd->mma_noise[i];
    }

    double min_mma = pd->mma_mean - 5*pd->mma_std;
    double max_mma = pd->mma_mean + 5*pd->mma_std;

    cpgsvp(x1, x2, y1, y2);

    cpgsch(1.0);
    cpgmtxt("l", 2.75, 0.5, 0.5, "mma");

    // Plot top axis markers in UTC hour, bottom axis markers in seconds
    cpgsch(0.9);
    cpgswin(pd->time_offset + pd->time_min, pd->time_offset + pd->time_max, min_mma, max_mma);
    cpgtbox("cstZ", 0, 0, "bcstnv", 0, 0);
    cpgswin(pd->time_scale*pd->time_min, pd->time_scale*pd->time_max, min_mma, max_mma);
    cpgbox("bst", 0, 0, "0", 0, 0);
    cpgsch(1.0);

    cpgswin(pd->time_min, pd->time_max, min_mma, max_mma);
    if (data->plot_error_bars)
        cpgerrb(6, pd->filtered_count, time, mma, mma_noise, 0.0);
    else
        cpgpt(pd->filtered_count, time, mma, 229);

allocation_error:
    free(time);
    free(mma);
    free(mma_noise);

    return ret;
}

int online_focus_plot(char *data_path, const char *device, double size)
{
    int ret = 0;

    datafile *data = datafile_load(data_path);
    if (!data)
        return error("Error opening data file %s", data_path);

    struct photometry_data *pd = datafile_generate_photometry(data);
    if (!pd)
        error_jump(photometry_error, ret, "Photometry calculation failed");

    if (cpgopen(device) <= 0)
        error_jump(setup_error, ret, "Unable to open PGPLOT window");

    cpgpap(size, 0.6);
    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(1);

    if (plot_raw_panel(0.065, 0.98, 0.075, 0.55, data, pd))
        error_jump(plot_error, ret, "Error plotting raw panel");

    if (plot_fwhm_panel(0.065, 0.98, 0.55, 0.91, data, pd))
        error_jump(plot_error, ret, "Error plotting fwhm panel");

    plot_time_axes(0.065, 0.98, 0.075, 0.91, data, pd);

plot_error:
    cpgend();
setup_error:
    datafile_free_photometry(pd);
photometry_error:
    datafile_free(data);
    return ret;
}

static int plot_internal(datafile *data, const char *ts_device, const char *dft_device, double size)
{
    int ret = 0;

    struct photometry_data *pd = datafile_generate_photometry(data);
    if (!pd)
        error_jump(photometry_error, ret, "Photometry calculation failed");

    struct dft_data *dd = datafile_generate_dft(data, pd);
    if (!dd)
        error_jump(dft_error, ret, "DFT calculation failed");

    double window_range = (data->plot_max_uhz - data->plot_min_uhz)/16;
    double window_freq = (data->plot_max_uhz + data->plot_min_uhz)/2;
    size_t window_count = data->plot_num_uhz/5;

    struct dft_data *wd = datafile_generate_window(data, pd, window_freq, window_range, window_count);
    if (!wd)
        error_jump(window_error, ret, "DFT window calculation failed");

    if (cpgopen(ts_device) <= 0)
        error_jump(setup_error, ret, "Unable to open PGPLOT window");

    cpgpap(size, 0.6);
    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(1);

    //
    // Plot blocked ranges
    //
    cpgsvp(0.075, 0.975, 0.075, 0.93);
    cpgswin(pd->time_min, pd->time_max, 0, 1);
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

    if (plot_raw_panel(0.065, 0.98, 0.075, 0.55, data, pd))
        error_jump(plot_error, ret, "Error plotting raw panel");

    if (plot_fwhm_panel(0.065, 0.98, 0.55, 0.67, data, pd))
        error_jump(plot_error, ret, "Error plotting fwhm panel");

    // Mean FWHM label
    snprintf(label, 32, "Mean: %.2f\" (%.2fpx)", pd->fwhm_mean*data->ccd_platescale, pd->fwhm_mean);
    cpgsch(0.9);
    cpgmtxt("t", -1.2, 1.0, 1.05, label);
    cpgsch(1.0);

    if (plot_ratio_panel(0.065, 0.98, 0.67, 0.79, data, pd))
        error_jump(plot_error, ret, "Error plotting ratio panel");

    if (plot_mma_panel(0.065, 0.98, 0.79, 0.91, data, pd))
        error_jump(plot_error, ret, "Error plotting fwhm panel");

    plot_time_axes(0.065, 0.98, 0.075, 0.91, data, pd);

    cpgend();

    //
    // Plot DFT
    //
    if (cpgopen(dft_device) <= 0)
        error_jump(setup_error, ret, "Unable to open PGPLOT window");

    cpgpap(size, 0.6);
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
    double max_dft = data->plot_max_dft ? data->plot_max_dft : 1.1*dd->max_ampl;
    cpgsvp(0.065, 0.98, 0.07, 0.97);
    cpgswin(scale*dd->min_freq, scale*dd->max_freq, 0, 1);
    cpgsch(0.9);
    cpgbox("bcstn", 0, 0, "0", 0, 0);
    cpgswin(dd->min_freq, dd->max_freq, 0, max_dft);
    cpgbox("0", 0, 0, "bcstnv", 0, 0);
    cpgsch(1.0);

    cpgsci(2);
    cast_double_array_to_float(dd->freq, dd->count);
    cast_double_array_to_float(dd->ampl, dd->count);
    cpgline(dd->count, (float *)dd->freq, (float *)dd->ampl);
    cpgsci(1);

    // Calculate median intensity
    qsort(dd->ampl, dd->count, sizeof(float), compare_float);

    snprintf(label, label_len, "Frequency (%sHz)", unit);
    cpgmtxt("b", 2.5, 0.5, 0.5, label);

    cpgmtxt("l", 2.75, 0.5, 0.5, "Amplitude (mma)");

    cpgswin(0, 1, 0, 1);
    snprintf(label, label_len, "Mean: %.2f mma", dd->mean_ampl);
    cpgptxt(0.97, 0.94, 0, 1.0, label);
    snprintf(label, label_len, "Median: %.2f mma", ((float *)dd->ampl)[dd->count/2]);
    cpgptxt(0.97, 0.90, 0, 1.0, label);

    // DFT Window
    cpgsvp(0.09, 0.2025, 0.75, 0.94);
    cpgswin(wd->min_freq, wd->max_freq, 0, 1.1);
    cpgsci(2);
    cast_double_array_to_float(wd->freq, wd->count);
    cast_double_array_to_float(wd->ampl, wd->count);
    cpgline(wd->count, (float *)wd->freq, (float *)wd->ampl);
    cpgsci(1);
    cpgbox("bc", 0, 0, "bc", 0, 0);
    cpgmtxt("b", 1.25, 0.5, 0.5, "Window");


plot_error:
    cpgend();
setup_error:
    datafile_free_dft(wd);
window_error:
    datafile_free_dft(dd);
dft_error:
    datafile_free_photometry(pd);
photometry_error:
    return ret;
}

int online_plot(char *data_path, char *ts_device, char *dft_device, double size)
{
    datafile *data = datafile_load(data_path);
    if (!data)
        return error("Error opening data file %s", data_path);

    int ret = plot_internal(data, ts_device, dft_device, size);
    datafile_free(data);
    return ret;
}

int playback_reduction(char *data_path, int delay, int step, char *ts_device, char *dft_device, double size)
{
    int ret = 0;

    // Count the number of observations
    datafile *data = datafile_load(data_path);
    if (!data)
        return error("Error opening data file %s", data_path);

    size_t limit = data->obs_count;
    for (size_t i = data->ratio_fit_degree + 1; i < limit; i += step)
    {
        // Limit the data to the first N observations
        data->obs_count = i;

        clock_t start = clock();
        if (plot_internal(data, ts_device, dft_device, size))
            error_jump(plot_error, ret, "Plotting error");

        // Attempt to compensate for calculation time
        int ms = (clock() - start)*1000/CLOCKS_PER_SEC;
        if (ms < delay)
            millisleep(delay - ms);
    }

    // Finish by plotting with all data
    data->obs_count = limit;
    if (plot_internal(data, ts_device, dft_device, size))
        error_jump(plot_error, ret, "Plotting error");

plot_error:
    datafile_free(data);
    return ret;
}

int plot_range(char *datafile_pattern)
{
    char **datafile_names;
    size_t num_files = get_matching_files(datafile_pattern, &datafile_names);

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
        online_plot(datafile_names[i], "5/xs", "6/xs", 9.41);

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

