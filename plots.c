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

#define PL_DOUBLE
#include <plplot/plplot.h>

#include "plots.h"
#include "datafile.h"
#include "helpers.h"
#include "fit.h"

// Set the color table to match PGPLOT
static void set_color_table()
{
    PLINT r[] = {0,248,248,  0, 20,  0,249,247,248,125,  0, 15,130,249, 80,168};
    PLINT g[] = {0,248,  0,248, 19,248,  8,247,127,248,248,130, 17,  0, 80,168};
    PLINT b[] = {0,248, 25, 79,237,247,239, 84, 46, 80,143,240,238,125, 80,168};
    plscmap0(r,g,b, 16);
}

static size_t plot_colors_max = 8;
static uint8_t plot_colors[] = {2,8,9,11,5,6,7,12};

static int plot_error_bars(double *x, double *y, double *dy, size_t n)
{
    int ret = 0;
    double *y_min = malloc(n*sizeof(double));
    double *y_max = malloc(n*sizeof(double));
    if (!y_min || !y_max)
        error_jump(allocation_error, ret, "Allocation error");

    for (size_t i = 0; i < n; i++)
    {
        y_min[i] = y[i] - dy[i];
        y_max[i] = y[i] + dy[i];
    }

    plsmin(0, 0);
    plerry(n, x, y_min, y_max);
    plsmin(0, 1);
allocation_error:
    free(y_min);
    free(y_max);
    return ret;
}

static void plot_time_axes(float x1, float x2, float y1, float y2,
                           datafile *data, struct photometry_data *pd)
{
    plvpor(x1, x2, y1, y2);
    plschr(0, 0.9);
    plwind(pd->time_offset + pd->time_min, pd->time_offset + pd->time_max, 0, 1);
    plbox("cstmd", 0, 0, "0", 0, 0);
    plwind(pd->time_scale*pd->time_min, pd->time_scale*pd->time_max, 0, 1);
    plbox("bstn", 0, 0, "0", 0, 0);
    plschr(0, 1.0);

    plmtex("t", 2.0, 0.5, 0.5, "UTC Time");

    char label[64];
    if (pd->time_exponent == 0)
        strncpy(label, "Run Time (s)", 64);
    else
        snprintf(label, 64, "Run Time (10\\u%d\\d s)", pd->time_exponent);

    plmtex("b", 2.5, 0.5, 0.5, label);

}

static int plot_raw_panel(float x1, float x2, float y1, float y2,
                          datafile *data, struct photometry_data *pd)
{
    int ret = 0;
    char label[64];
    const size_t label_len = 64;

    double min_raw = 0;
    double max_raw = data->plot_max_raw ? data->plot_max_raw : 1.25*pd->scaled_target_max;
    int raw_exp = (int)log10(max_raw);
    double raw_scale = 1.0/pow(10, raw_exp);

    plvpor(x1, x2, y1, y2);

    // Plot top axis markers in UTC hour, bottom axis markers in seconds
    plschr(0, 0.9);
    plwind(pd->time_offset + pd->time_min, pd->time_offset + pd->time_max, raw_scale*min_raw, raw_scale*max_raw);
    plbox("cstd", 0, 0, "bcstnv", 0, 0);

    plwind(pd->time_scale*pd->time_min, pd->time_scale*pd->time_max, min_raw, max_raw);
    plbox("bst", 0, 0, "0", 0, 0);
    plschr(0, 1.0);

    snprintf(label, label_len, "Count Rate (10\\u%d\\d ADU/s)", raw_exp);
    plmtex("l", 2.75, 0.5, 0.5, label);

    // Raw intensities
    for (size_t j = 0; j < data->target_count; j++)
    {
        plwind(pd->time_min, pd->time_max, min_raw, max_raw/data->targets[j].scale);
        plcol0(plot_colors[j%plot_colors_max]);

        size_t k = j*data->obs_count;
        if (data->plot_error_bars)
            plot_error_bars(&pd->target_time[k], &pd->target_intensity[k], &pd->target_noise[k], data->obs_count);
        else
            plpoin(pd->target_count[j], &pd->target_time[k], &pd->target_intensity[k], 20);
    }

    // Mean sky intensity
    plwind(pd->time_min, pd->time_max, min_raw, max_raw);
    plcol0(15);
    plpoin(pd->raw_count, pd->raw_time, pd->sky, 20);

    // Labels
    plvpor(x1, x2, 0.82*y2, y2);
    plwind(0, data->target_count + 1, 0, 1);
    for (size_t j = 0; j <= data->target_count; j++)
    {
        plcol0(plot_colors[j%plot_colors_max]);

        if (j == data->target_count)
        {
            plcol0(15);
            strncpy(label, "Mean Sky", label_len);
        }
        else
        {
            if (data->targets[j].scale == 1.0)
                strncpy(label, data->targets[j].label, label_len);
            else
                snprintf(label, label_len, "%g \\x %s", data->targets[j].scale, data->targets[j].label);
        }

        plptex(j+0.5, 0.5, 0, 0, 0.5, label);
        if (j < data->target_count)
        {
            snprintf(label, label_len, "SNR: %.0f", pd->target_snr[j]);
            plschr(0, 0.9);
            plptex(j+0.5, 0.25, 0, 0, 0.5, label);
            plschr(0, 1.0);
        }

    }
    plcol0(1);
    return ret;
}

static int plot_fwhm_panel(float x1, float x2, float y1, float y2,
                           datafile *data, struct photometry_data *pd)
{
    int ret = 0;

    // Smoothed FWHM estimate
    double *smooth = malloc(pd->raw_count*sizeof(double));
    if (!smooth)
        error_jump(allocation_error, ret, "Allocation error");

    for (size_t i = 0; i < pd->raw_count; i++)
        smooth[i] = pd->fwhm[i];

    double min_fwhm = pd->fwhm_mean - 5*pd->fwhm_std;
    double max_fwhm = pd->fwhm_mean + 8*pd->fwhm_std;

    plvpor(x1, x2, y1, y2);

    // Plot top axis markers in UTC hour, bottom axis markers in seconds
    // Reserve the top 25% of the plot for the mean FWHM label
    plmtex("l", 2.75, 0.5, 0.5, "FWHM (\")");

    plschr(0, 0.9);
    plwind(pd->time_offset + pd->time_min, pd->time_offset + pd->time_max, min_fwhm*data->ccd_platescale, max_fwhm*data->ccd_platescale);
    plbox("cstd", 0, 0, "bcstnv", 0, 0);
    plwind(pd->time_scale*pd->time_min, pd->time_scale*pd->time_max, min_fwhm, max_fwhm);
    plbox("bst", 0, 0, "0", 0, 0);
    plschr(0, 1.0);

    plwind(pd->time_min, pd->time_max, min_fwhm, max_fwhm);
    plpoin(pd->raw_count, pd->raw_time, pd->fwhm, 20);

    // Calculate running mean
    for (size_t i = 0; i < pd->raw_count; i++)
    {
        // Handle start and end of data
        size_t run = data->plot_fwhm_smooth;
        size_t start = i - run/2;
        if (run/2 > i)
        {
            start = 0;
            run -= run/2 - i;
        }

        if (start + run >= pd->raw_count)
            run = pd->raw_count - start - 1;

        smooth[i] = mean_exclude_sigma(&smooth[start], run, 3);
    }

    plcol0(2);
    plline(pd->raw_count, pd->raw_time, smooth);
    plcol0(1);

    // Mean FWHM label
    plschr(0, 0.9);
    char label[32];
    snprintf(label, 32, "Mean: %.2f\" (%.2fpx)", pd->fwhm_mean*data->ccd_platescale, pd->fwhm_mean);
    plmtex("t", -1.25, 0.0, -0.1, label);

    snprintf(label, 32, "Current: %.2f\" (%.2fpx)", smooth[pd->raw_count - 1]*data->ccd_platescale, smooth[pd->raw_count - 1]);
    plmtex("t", -1.25, 1.0, 1.1, label);
    plschr(0, 1.0);

    free(smooth);
allocation_error:
    return ret;
}

static int plot_ratio_panel(float x1, float x2, float y1, float y2,
                            datafile *data, struct photometry_data *pd)
{
    int ret = 0;

    double min_ratio = pd->ratio_mean - 5*pd->ratio_std;
    double max_ratio = pd->ratio_mean + 5*pd->ratio_std;

    plvpor(x1, x2, y1, y2);
    plmtex("l", 2.75, 0.5, 0.5, "Ratio");

    // Plot top axis markers in UTC hour, bottom axis markers in seconds
    plschr(0, 0.9);
    plwind(pd->time_offset + pd->time_min, pd->time_offset + pd->time_max, min_ratio, max_ratio);
    plbox("cstd", 0, 0, "bc", 0, 0);
    plwind(pd->time_scale*pd->time_min, pd->time_scale*pd->time_max, min_ratio, max_ratio);
    plbox("bst", 0, 0, "0", 0, 0);
    plschr(0, 1.0);

    plwind(pd->time_min, pd->time_max, min_ratio, max_ratio);
    if (data->plot_error_bars)
        plot_error_bars(pd->time, pd->ratio, pd->ratio_noise, pd->filtered_count);
    else
        plpoin(pd->filtered_count, pd->time, pd->ratio, 20);

    // Plot the polynomial fit
    plcol0(2);

    // Don't plot the fit through blocked regions
    double *t = pd->time;
    double *t_end = &t[pd->filtered_count];
    double *f = pd->ratio_fit;
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

        plline(count, t, f);

        // Start next search from end of this section
        t += count + 1;
        f += count + 1;
    } while (t + 1 < t_end);

    plcol0(1);
    return ret;
}

static int plot_mma_panel(float x1, float x2, float y1, float y2,
                          datafile *data, struct photometry_data *pd)
{
    int ret = 0;

    double min_mma = pd->mma_mean - 5*pd->mma_std;
    double max_mma = pd->mma_mean + 5*pd->mma_std;

    plvpor(x1, x2, y1, y2);

    plschr(0, 1.0);
    plmtex("l", 2.75, 0.5, 0.5, "mma");

    // Plot top axis markers in UTC hour, bottom axis markers in seconds
    plschr(0, 0.9);
    plwind(pd->time_offset + pd->time_min, pd->time_offset + pd->time_max, min_mma, max_mma);
    plbox("cstd", 0, 0, "bcstnv", 0, 0);
    plwind(pd->time_scale*pd->time_min, pd->time_scale*pd->time_max, min_mma, max_mma);
    plbox("bst", 0, 0, "0", 0, 0);
    plschr(0, 1.0);

    plwind(pd->time_min, pd->time_max, min_mma, max_mma);
    if (data->plot_error_bars)
        plot_error_bars(pd->time, pd->mma, pd->mma_noise, pd->filtered_count);
    else
        plpoin(pd->filtered_count, pd->time, pd->mma, 20);

    return ret;
}

static int plot_dft_panel(float x1, float x2, float y1, float y2,
                          datafile *data, struct dft_data *dft, struct dft_data *window)
{
    int ret = 0;
    char label[64];
    const size_t label_len = 64;

    // Calculate baseline scale
    double scale = 1e6;
    char *unit = "\\gm";

    if (dft->max_freq > 1.5e4)
    {
        scale = 1.0e-3;
        unit = "k";
    }
    else if (dft->max_freq > 1.5e1)
    {
        scale = 1.0;
        unit = "";
    }
    else if (dft->max_freq > 1.5e-2)
    {
        scale = 1.0e3;
        unit = "m";
    }

    double max_dft = data->plot_max_dft ? data->plot_max_dft : 1.1*dft->max_ampl;

    // Check for overlap with the DFT panel
    double max = 0;
    double threshold_ampl = 0.65*max_dft;
    double threshold_freq = 0.16*dft->max_freq;
    for (size_t i = 0; i < dft->count && dft->freq[i] < threshold_freq; i++)
        max = fmax(max, dft->ampl[i]);
    if (max > threshold_ampl)
        max_dft *= max/threshold_ampl;

    plvpor(x1, x2, y1, y2);
    plwind(scale*dft->min_freq, scale*dft->max_freq, 0, 1);
    plschr(0, 0.9);
    plbox("bcstn", 0, 0, "0", 0, 0);
    plwind(dft->min_freq, dft->max_freq, 0, max_dft);
    plbox("0", 0, 0, "bcstnv", 0, 0);
    plschr(0, 1.0);

    plcol0(2);
    plline(dft->count, dft->freq, dft->ampl);
    plcol0(1);

    snprintf(label, label_len, "Frequency (%sHz)", unit);
    plmtex("b", 2.5, 0.5, 0.5, label);

    plmtex("l", 2.75, 0.5, 0.5, "Amplitude (mma)");

    // DFT Window
    float wx1 = x1 + 0.025;
    float wx2 = wx1 + (x2 - x1)/8;
    float wy2 = y2 - 0.04;
    float wy1 = y2 - (y2 - y1)/3.5;

    plvpor(wx1, wx2, wy1, wy2);
    plwind(window->min_freq, window->max_freq, 0, 1.1);
    plcol0(2);
    plline(window->count, window->freq, window->ampl);
    plcol0(1);
    plbox("bc", 0, 0, "bc", 0, 0);
    plmtex("b", 1.25, 0.5, 0.5, "Window");

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
/*
    if (cpgopen(device) <= 0)
        error_jump(setup_error, ret, "Unable to open PGPLOT window");

    cpgpap(size, 0.6);
    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(1);

    if (plot_raw_panel(0.065, 0.98, 0.075, 0.55, data, pd))
        error_jump(plot_error, ret, "Error plotting raw panel");

    if (plot_fwhm_panel(0.065, 0.98, 0.55, 0.93, data, pd))
        error_jump(plot_error, ret, "Error plotting fwhm panel");

    plot_time_axes(0.065, 0.98, 0.075, 0.93, data, pd);

plot_error:
    cpgend();
 */
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

    plinit();
    pladv(0);
/*
    cpgpap(size, 0.6);
    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(1);
*/
    pltimefmt("%H:%M");
    set_color_table();
    plssym(0, 0.3);

    //
    // Plot blocked ranges
    //
    plvpor(0.065, 0.98, 0.075, 0.93);
    plwind(pd->time_min, pd->time_max, 0, 1);
    plcol0(14);
    for (size_t j = 0; j < data->num_blocked_ranges; j++)
    {
        pljoin(data->blocked_ranges[j].x, 0, data->blocked_ranges[j].x, 1);
        pljoin(data->blocked_ranges[j].y, 0, data->blocked_ranges[j].y, 1);
    }
    plcol0(1);

    // Generic label buffer for passing data to pgplot
    char label[64];
    const size_t label_len = 64;

    if (plot_raw_panel(0.065, 0.98, 0.075, 0.55, data, pd))
        error_jump(plot_error, ret, "Error plotting raw panel");

    if (plot_fwhm_panel(0.065, 0.98, 0.55, 0.67, data, pd))
        error_jump(plot_error, ret, "Error plotting fwhm panel");

    if (plot_ratio_panel(0.065, 0.98, 0.67, 0.79, data, pd))
        error_jump(plot_error, ret, "Error plotting ratio panel");

    if (plot_mma_panel(0.065, 0.98, 0.79, 0.93, data, pd))
        error_jump(plot_error, ret, "Error plotting fwhm panel");

    plot_time_axes(0.065, 0.98, 0.075, 0.93, data, pd);
    pladv(0);
/*
    cpgpap(size, 0.6);
    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(1);
*/
    // DFT
    if (plot_dft_panel(0.065, 0.98, 0.07, 0.97, data, dd, wd))
        error_jump(plot_error, ret, "Error plotting dft panel");

    // Calculate median intensity
    qsort(dd->ampl, dd->count, sizeof(double), compare_double);

    plvpor(0.065, 0.98, 0.07, 0.97);
    plwind(0, 1, 0, 1);
    snprintf(label, label_len, "Mean: %.2f mma", dd->mean_ampl);
    plptex(0.97, 0.94, 0, 0, 1.0, label);
    snprintf(label, label_len, "Median: %.2f mma", dd->ampl[dd->count/2]);
    plptex(0.97, 0.90, 0, 0, 1.0, label);

plot_error:
    plend();
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
/*
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
        plschr(0, 20);
        plwind(0, 1, 0, 1);
        plcol0(1);
        plptex(0.5, 0.7, 0, 0, 0.5, datafile_names[i]);
        plptex(0.5, 0.0, 0, 0, 0.5, "Next file: [n] Prev file: [p] Quit: [q]");

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
*/
endLoop:
    free_2d_array(datafile_names, num_files);
    return 0;
}

