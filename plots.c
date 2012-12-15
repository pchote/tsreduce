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

#include "plots.h"
#include "datafile.h"
#include "helpers.h"
#include "fit.h"

int online_focus_plot(char *data_path, const char *device, double size)
{
    size_t plot_colors_max = 8;
    uint8_t plot_colors[] = {2,8,3,12,5,6,7,9};

    datafile *data = datafile_load(data_path);
    if (!data)
        return error("Error opening data file %s", data_path);

    struct photometry_data *pd = datafile_generate_photometry(data);
    if (!pd)
    {
        datafile_free(data);
        return error("Photometry calculation failed");
    }

    double base_seconds = 3600*ts_time_to_utc_hour(data->reference_time);
    double min_seconds = pd->raw_time[0];
    double max_seconds = pd->raw_time[pd->raw_count - 1];
    int secexp = (int)(log10(max_seconds) / 3)*3;
    double sec_scale = 1.0/pow(10, secexp);

    if (cpgopen(device) <= 0)
    {
        datafile_free(data);
        return error("Unable to open PGPLOT window");
    }

    cpgpap(size, 0.6);
    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(1);

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

        cpgsvp(0.065, 0.98, 0.075, 0.55);

        // Plot top axis markers in UTC hour, bottom axis markers in seconds
        cpgsch(0.9);
        cpgswin(base_seconds + min_seconds, base_seconds + max_seconds, raw_scale*min_raw, raw_scale*max_raw);
        cpgtbox("cstZ", 0, 0, "bcstnv", 0, 0);
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
        cpgsvp(0.065, 0.98, 0.45, 0.55);
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
    }

    //
    // Plot FWHM
    //
    {
        double min_fwhm = pd->fwhm_mean - 5*pd->fwhm_std;
        double max_fwhm = pd->fwhm_mean + 7.5*pd->fwhm_std;

        cpgsvp(0.065, 0.98, 0.55, 0.91);

        // Plot top axis markers in UTC hour, bottom axis markers in seconds
        // Reserve the top 25% of the plot for the mean FWHM label
        cpgmtxt("t", 2.0, 0.5, 0.5, "UTC Time");
        cpgmtxt("l", 2.75, 0.5, 0.5, "FWHM (\")");

        cpgsch(0.9);
        cpgswin(base_seconds + min_seconds, base_seconds + max_seconds, min_fwhm*data->ccd_platescale, max_fwhm*data->ccd_platescale);
        cpgtbox("cstmXYZH", 0, 0, "bcstnv", 0, 0);
        cpgswin(sec_scale*min_seconds, sec_scale*max_seconds, min_fwhm, max_fwhm);
        cpgbox("bst", 0, 0, "0", 0, 0);
        cpgsch(1.0);

        cpgswin(min_seconds, max_seconds, min_fwhm, max_fwhm);
        cpgpt(pd->raw_count, (float *)pd->raw_time, (float *)pd->fwhm, 229);
    }

    cpgend();

    datafile_free_photometry(pd);
    datafile_free(data);

    return 0;
}

static int plot_internal(datafile *data, const char *tsDevice, double tsSize, const char *dftDevice, double dftSize)
{
    size_t plot_colors_max = 8;
    uint8_t plot_colors[] = {2,8,3,12,5,6,7,9};

    struct photometry_data *pd = datafile_generate_photometry(data);
    if (!pd)
        return error("Photometry calculation failed");

    struct dft_data *dd = datafile_generate_dft(data, pd);
    if (!dd)
    {
        datafile_free_photometry(pd);
        return error("DFT calculation failed");
    }

    double window_range = (data->plot_max_uhz - data->plot_min_uhz)/16;
    double window_freq = (data->plot_max_uhz + data->plot_min_uhz)/2;
    size_t window_count = data->plot_num_uhz/5;

    struct dft_data *wd = datafile_generate_window(data, pd, window_freq, window_range, window_count);
    if (!wd)
    {
        datafile_free_photometry(pd);
        datafile_free_dft(dd);
        return error("DFT window calculation failed");
    }

    double base_seconds = 3600*ts_time_to_utc_hour(data->reference_time);
    double min_seconds = pd->raw_time[0];
    double max_seconds = pd->raw_time[pd->raw_count - 1];
    int secexp = (int)(log10(max_seconds) / 3)*3;
    double sec_scale = 1.0/pow(10, secexp);

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
    cpgsvp(0.075, 0.975, 0.075, 0.93);
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

        cpgsvp(0.065, 0.98, 0.075, 0.55);

        // Plot top axis markers in UTC hour, bottom axis markers in seconds
        cpgsch(0.9);
        cpgswin(base_seconds + min_seconds, base_seconds + max_seconds, raw_scale*min_raw, raw_scale*max_raw);
        cpgtbox("cstZ", 0, 0, "bcstnv", 0, 0);
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
        cpgsvp(0.065, 0.98, 0.45, 0.55);
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
    }

    //
    // Plot FWHM
    //
    {
        double min_fwhm = pd->fwhm_mean - 5*pd->fwhm_std;
        double max_fwhm = pd->fwhm_mean + 7.5*pd->fwhm_std;

        cpgsvp(0.065, 0.98, 0.55, 0.67);
        cpgmtxt("l", 2.75, 0.5, 0.5, "FWHM (\")");

        // Plot top axis markers in UTC hour, bottom axis markers in seconds
        // Reserve the top 25% of the plot for the mean FWHM label
        cpgsch(0.9);
        cpgswin(base_seconds + min_seconds, base_seconds + max_seconds, min_fwhm*data->ccd_platescale, max_fwhm*data->ccd_platescale);
        cpgtbox("cstZ", 0, 0, "bcstnv", 0, 0);
        cpgswin(sec_scale*min_seconds, sec_scale*max_seconds, min_fwhm, max_fwhm);
        cpgbox("bst", 0, 0, "0", 0, 0);
        cpgsch(1.0);

        cpgswin(min_seconds, max_seconds, min_fwhm, max_fwhm);
        cpgpt(pd->raw_count, (float *)pd->raw_time, (float *)pd->fwhm, 229);

        snprintf(label, 32, "Mean: %.2f\" (%.2fpx)", pd->fwhm_mean*data->ccd_platescale, pd->fwhm_mean);
        cpgsch(0.9);
        cpgptxt(max_seconds, max_fwhm - 3*pd->fwhm_std, 0, 1.05, label);
        cpgsch(1.0);
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

        cpgsvp(0.065, 0.98, 0.67, 0.79);
        cpgmtxt("l", 2.75, 0.5, 0.5, "Ratio");

        // Plot top axis markers in UTC hour, bottom axis markers in seconds
        cpgsch(0.9);
        cpgswin(base_seconds + min_seconds, base_seconds + max_seconds, min_ratio, max_ratio);
        cpgtbox("cstZ", 0, 0, "bc", 0, 0);
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
        float *t_end = &t[pd->filtered_count];
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
            while (t + count < t_end && *(t + count) <= max_time)
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

        cpgsvp(0.065, 0.98, 0.79, 0.91);
        cpgsch(1.0);

        cpgmtxt("t", 2.0, 0.5, 0.5, "UTC Time");
        cpgmtxt("l", 2.75, 0.5, 0.5, "mma");

        // Plot top axis markers in UTC hour, bottom axis markers in seconds
        cpgsch(0.9);
        cpgswin(base_seconds + min_seconds, base_seconds + max_seconds, min_mma, max_mma);
        cpgtbox("cstmXYZH", 0, 0, "bcstnv", 0, 0);
        cpgswin(sec_scale*min_seconds, sec_scale*max_seconds, min_mma, max_mma);
        cpgbox("bst", 0, 0, "0", 0, 0);
        cpgsch(1.0);
        
        cpgswin(min_seconds, max_seconds, min_mma, max_mma);
        if (data->plot_error_bars)
            cpgerrb(6, pd->filtered_count, (float *)pd->time, (float *)pd->mma, (float *)pd->mma_noise, 0.0);
        else
            cpgpt(pd->filtered_count, (float *)pd->time, (float *)pd->mma, 229);
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

    cpgend();

    datafile_free_dft(wd);
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

