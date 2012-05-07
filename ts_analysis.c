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
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <cpgplot.h>

#include "reduction.h"
#include "datafile.h"
#include "helpers.h"
#include "fit.h"

// Display the apertures specified by the reduction file at dataPath for the
// observation index obsIndex over an image frame in ds9
int display_targets(char *dataPath, int obsIndex)
{
    if (!init_ds9("tsreduce"))
        return error("Unable to launch ds9");

    // Read file header
    datafile *data = datafile_load(dataPath);
    if (data == NULL)
        return error("Error opening data file");

    if (obsIndex >= data->num_obs)
    {
        datafile_free(data);
        return error("Requested observation is out of range: max is %d", data->num_obs-1);
    }

    if (chdir(data->frame_dir))
    {
        datafile_free(data);
        return error("Invalid frame path: %s", data->frame_dir);
    }

    char command[128];
    char filenamebuf[NAME_MAX];
    filenamebuf[0] = '\0';
    if (obsIndex >= 0)
        strncpy(filenamebuf, data->obs[obsIndex].filename, NAME_MAX);
    else if (!get_first_matching_file(data->frame_pattern, filenamebuf, NAME_MAX))
    {
        datafile_free(data);
        return error("No matching files found");
    }

    // DS9 errors are nonfatal
    snprintf(command, 128, "file %s/%s", data->frame_dir, filenamebuf);
    if (tell_ds9("tsreduce", command, NULL, 0))
        error("ds9 command failed: %s", command);

    // Set scaling mode
    if (tell_ds9("tsreduce", "scale mode 99.5", NULL, 0))
        error("ds9 command failed: scale mode 99.5");

    // Flip X axis
    if (tell_ds9("tsreduce", "orient x", NULL, 0))
        error("ds9 command failed: orient x");

    for (int i = 0; i < data->num_targets; i++)
    {
        double x, y;
        // Use the specified target aperture
        // DS9 starts labelling pixels at 1
        if (obsIndex < 0)
        {
            x = data->targets[i].x + 1;
            y = data->targets[i].y + 1;
        }
        // Use the centered observation
        else
        {
            x = data->obs[obsIndex].pos[i].x + 1;
            y = data->obs[obsIndex].pos[i].y + 1;
        }

        snprintf(command, 128, "regions command {circle %f %f %f #color=red}", x, y, data->targets[i].r);
        if (tell_ds9("tsreduce", command, NULL, 0))
            error("ds9 command failed: %s\n", command);

        snprintf(command, 128, "regions command {circle %f %f %f #background}", x, y, data->targets[i].s1);
        if (tell_ds9("tsreduce", command, NULL, 0))
            error("ds9 command failed: %s\n", command);

        snprintf(command, 128, "regions command {circle %f %f %f #background}", x, y, data->targets[i].s2);
        if (tell_ds9("tsreduce", command, NULL, 0))
            error("ds9 command failed: %s\n", command);
    }

    datafile_free(data);
    return 0;
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

    if (obsIndex >= data->num_obs)
        error_jump(setup_error, ret, "Requested observation is out of range: max is %d", data->num_obs-1);

    if (chdir(data->frame_dir))
        error_jump(setup_error, ret, "Invalid frame path: %s", data->frame_dir);

    char filenamebuf[NAME_MAX];
    if (obsIndex >= 0)
        strncpy(filenamebuf, data->obs[obsIndex].filename, NAME_MAX);
    else if (!get_first_matching_file(data->frame_pattern, filenamebuf, NAME_MAX))
        error_jump(setup_error, ret, "No matching files found");

    framedata *frame = framedata_load(filenamebuf);
    if (!frame)
        error_jump(setup_error, ret, "Error loading frame %s", filenamebuf);
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
        const int numIntensity = 21;
        double intensity[numIntensity];
        double noise[numIntensity];
        double radii[numIntensity];
        double profile[numIntensity];

        double readnoise = framedata_get_header_dbl(flat, "CCD-READ");
        double gain = framedata_get_header_dbl(flat, "CCD-GAIN");

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
    if (data->num_obs <= 0)
    {
        datafile_free(data);
        return error("File specifies no observations");
    }

    int last_bad = 0;
    double last_time = data->obs[0].time;
    for (int i = 1; i < data->num_obs; i++)
    {
        double time = data->obs[i].time;
        if (time - last_time < 0.1f)
            last_bad = 1;

        if (last_bad)
            printf("%s @ %.1f\n", data->obs[i].filename, time);

        if (last_bad && time - last_time >= 0.1f)
            last_bad = 0;

        last_time = time;
    }

    datafile_free(data);
    return 0;
}

/*
 * Cast an array of doubles to an array of floats reusing the same memory
 */
float *cast_double_array_to_float(double *d_ptr, size_t count)
{
    float *f_ptr = (float *)d_ptr;
    for (size_t i = 0; i < count; i++)
        f_ptr[i] = d_ptr[i];
    return f_ptr;
}

int plot_fits(char *dataPath, char *tsDevice, char *dftDevice)
{
    int plot_colors_max = 8;
    int plot_colors[] = {4,2,8,3,5,6,7,9};

    // Read file header
    datafile *data = datafile_load(dataPath);
    if (data == NULL)
        return error("Error opening data file %s", dataPath);

    // Start time in hours
    struct tm starttime;
    gmtime_r(&data->reference_time, &starttime);
    float min_time = starttime.tm_hour + starttime.tm_min / 60.0 + starttime.tm_sec / 3600.0;

    double *raw_time_d, *raw_d, *time_d, *ratio_d, *polyfit_d, *mmi_d, *ratio_noise_d, *mmi_noise_d, *freq_d, *ampl_d;
    double ratio_mean, ratio_std, mmi_mean, mmi_std;
    size_t num_raw, num_filtered, num_dft;
    if (generate_photometry_dft_data(data,
                                     &raw_time_d, &raw_d, &num_raw,
                                     &time_d, &ratio_d, &polyfit_d, &mmi_d, &num_filtered,
                                     &ratio_noise_d, &mmi_noise_d,
                                     &ratio_mean, &ratio_std, &mmi_mean, &mmi_std,
                                     &freq_d, &ampl_d, &num_dft))
    {
        datafile_free(data);
        return error("Error generating data");
    }

    float max_time = min_time + raw_time_d[num_raw-1]/3600;

    double snr_ratio = 0;
    if (ratio_noise_d)
    {
        double total_ratio = 0;
        double total_ratio_noise = 0;
        for (int i = 0; i < num_filtered; i++)
        {
            total_ratio += ratio_d[i];
            total_ratio_noise += ratio_noise_d[i];
        }
        snr_ratio = total_ratio/total_ratio_noise;
    }

    // Cast the double arrays to float, does not allocate any new memory.
    float *raw_time = cast_double_array_to_float(raw_time_d, num_raw);
    float *raw = cast_double_array_to_float(raw_d, num_raw*data->num_targets);
    float *time = cast_double_array_to_float(time_d, num_filtered);
    float *ratio = cast_double_array_to_float(ratio_d, num_filtered);
    float *polyfit = cast_double_array_to_float(polyfit_d, num_filtered);
    float *mmi = cast_double_array_to_float(mmi_d, num_filtered);
    float *freq = cast_double_array_to_float(freq_d, num_dft);
    float *ampl = cast_double_array_to_float(ampl_d, num_dft);

    float *ratio_noise = NULL;
    if (ratio_noise_d)
        ratio_noise = cast_double_array_to_float(ratio_noise_d, num_filtered);

    float *mmi_noise = NULL;
    if (mmi_noise_d)
        mmi_noise = cast_double_array_to_float(mmi_noise_d, num_filtered);

    float min_mmi = mmi_mean - 5*mmi_std;
    float max_mmi = mmi_mean + 5*mmi_std;
    float min_ratio = ratio_mean - 5*ratio_std;
    float max_ratio = ratio_mean + 5*ratio_std;

    // Determine max dft ampl
    float max_dft_ampl = 0;
    for (int i = 0; i < num_dft; i++)
        max_dft_ampl = fmax(max_dft_ampl, ampl[i]);
    max_dft_ampl *= 1.1;

    if (cpgopen(tsDevice ? tsDevice : "5/xs") <= 0)
        return error("Unable to open PGPLOT window");

    // 800 x 480
    if (!tsDevice)
        cpgpap(9.41, 0.6);

    cpgask(0);
    cpgslw(3);
    cpgsfs(2);
    cpgscf(2);

    // Fitted MMI
    cpgsvp(0.1, 0.9, 0.75, 0.9);
    cpgsch(1.25);
    cpgmtxt("t", 2, 0.5, 0.5, "Time Series Data");
    cpgsch(1.0);

    cpgmtxt("l", 2.5, 0.5, 0.5, "mma");
    cpgswin(min_time, max_time, min_mmi, max_mmi);
    cpgbox("bcstm", 1, 4, "bcstn", 0, 0);

    cpgswin(0, raw_time[num_raw-1], min_mmi, max_mmi);

    if (mmi_noise)
        cpgerrb(6, num_filtered, time, mmi, mmi_noise, 0.0);
    else
        cpgpt(num_filtered, time, mmi, 20);

    // Ratio
    cpgsvp(0.1, 0.9, 0.55, 0.75);
    cpgmtxt("l", 2.5, 0.5, 0.5, "Ratio");
    cpgswin(min_time, max_time, min_ratio, max_ratio);
    cpgbox("bcst", 1, 4, "bcstn", 0, 0);
    cpgswin(0, raw_time[num_raw-1], min_ratio, max_ratio);

    // Plot error bars if ratio_noise is available (data->version >= 5)
    if (ratio_noise)
        cpgerrb(6, num_filtered, time, ratio, ratio_noise, 0.0);
    else
        cpgpt(num_filtered, time, ratio, 20);

    // Plot the polynomial fit
    cpgsci(2);
    cpgline(num_filtered, time, polyfit);
    cpgsci(1);

    // Raw Data
    cpgsvp(0.1, 0.9, 0.075, 0.55);
    cpgmtxt("l", 2.5, 0.5, 0.5, "Counts Per Second");
    cpgmtxt("b", 2.5, 0.5, 0.5, "UTC Hour");
    cpgswin(min_time, max_time, 0, data->plot_max_raw);
    cpgbox("bcstn", 1, 4, "bcstn", 0, 0);

    for (int j = 0; j < data->num_targets; j++)
    {
        cpgswin(0, raw_time[num_raw-1], 0, data->plot_max_raw/data->targets[j].plot_scale);
        cpgsci(plot_colors[j%plot_colors_max]);
        cpgpt(num_raw, raw_time, &raw[j*num_raw], 20);
    }

    // SNR label
    cpgswin(0, 1, 0, 1);
    cpgsci(1);
    char *snr_label;
    asprintf(&snr_label, "Ratio SNR: %.2f", snr_ratio);
    cpgptxt(0.97, 0.9, 0, 1.0, snr_label);
    free(snr_label);
    cpgend();

    if (cpgopen(dftDevice ? dftDevice : "6/xs") <= 0)
        return error("Unable to open PGPLOT window");

    // 800 x 480
    if (!dftDevice)
        cpgpap(9.41, 0.6);

    cpgask(0);
    cpgslw(3);
    cpgsfs(2);
    cpgscf(2);

    // DFT
    cpgsvp(0.1, 0.9, 0.075, 0.87);
    cpgswin(data->plot_min_uhz, data->plot_max_uhz, 0, 1);
    cpgbox("bstn", 0, 0, "bcnst", 5, 5);

    // Plot period on top axis
    cpgsci(1);
    cpgmove(data->plot_min_uhz, 1);
    cpgdraw(data->plot_max_uhz, 1);
    int default_periods[] = {100, 125, 150, 200, 300, 500, 2000};
    int default_period_count = 7;
    char buf[20];
    for (int i = 0; i < default_period_count; i++)
    {
        double freq = 1e6/default_periods[i];
        if (freq < data->plot_min_uhz || freq > data->plot_max_uhz)
            continue;

        cpgmove(freq, 1);
        cpgdraw(freq, 0.98);
        snprintf(buf, 20, "%d", default_periods[i]);
        cpgptxt(freq, 1.02, 0, 0.5, buf);
    }

    cpgswin(data->plot_min_uhz*1e-6, data->plot_max_uhz*1e-6, 0, max_dft_ampl);
    cpgsci(12);
    cpgline(num_dft, freq, ampl);
    cpgsci(1);

    cpgmtxt("b", 2.5, 0.5, 0.5, "Frequency (\\gmHz)");
    cpgmtxt("l", 2, 0.5, 0.5, "Amplitude (mma)");
    cpgmtxt("t", 2, 0.5, 0.5, "Period (s)");

    cpgsch(1.25);
    cpgmtxt("t", 3.2, 0.5, 0.5, "Fourier Amplitude Spectrum");
    cpgend();

    free(raw_time);
    free(raw);
    free(time);
    free(ratio);
    free(polyfit);
    free(mmi);
    free(freq);
    free(ampl);
    datafile_free(data);
    return 0;
}

int amplitude_spectrum(char *dataPath)
{
    datafile *data = datafile_load(dataPath);
    if (data == NULL)
        return error("Error opening data file");

    double *raw_time, *raw, *time, *ratio, *polyfit, *mmi, *freq, *ampl;
    size_t num_raw, num_filtered, num_dft;
    if (generate_photometry_dft_data(data,
                                     &raw_time, &raw, &num_raw,
                                     &time, &ratio, &polyfit, &mmi, &num_filtered,
                                     NULL, NULL,
                                     NULL, NULL, NULL, NULL,
                                     &freq, &ampl, &num_dft))
    {
        datafile_free(data);
        return error("Error generating data");
    }

    for (int i = 0; i < num_dft; i++)
        printf("%f %f\n", freq[i], ampl[i]);

    free(raw_time);
    free(raw);
    free(time);
    free(ratio);
    free(polyfit);
    free(mmi);
    free(freq);
    free(ampl);
    datafile_free(data);
    return 0;
}


int reduce_aperture_range(char *base_name, double min, double max, double step, char *prefix)
{
    int ret = 0;
    datafile *data = datafile_load(base_name);
    if (data == NULL)
        return error("Error opening data file");

    // Create and update a datafile for each aperture
    double radius = min;
    do
    {
        for (int i = 0; i < data->num_targets; i++)
            data->targets[i].r = radius;

        char *filename;
        asprintf(&filename, "%s-%0.2f.dat", prefix, radius);

        // Errors are non-fatal -> proceeed to the next file
        if (!datafile_save_header(data, filename))
            update_reduction(filename);

        free(filename);
        radius += step;
    } while (radius < max);

    datafile_free(data);
    return ret;
}

int plot_range(char *datafile_pattern)
{
    char **datafile_names;
    int num_files = get_matching_files(datafile_pattern, &datafile_names);
    if (num_files == 0)
        return error("No datafiles found matching pattern %s", datafile_pattern);

    int i = 0;
    float x,y;
    char c;
    while (true)
    {
        plot_fits(datafile_names[i], NULL, NULL);

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
 * Apply a constant time offset to an mmi file
 */
int offset_time(char *mmiFile, double offset)
{
    char linebuf[1024];
    FILE *file = fopen(mmiFile, "r+");
    if (file == NULL)
        return error("Unable to open file: %s", mmiFile);

    double time, mmi;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
        {
            printf("%s", linebuf);
            continue;
        }

        sscanf(linebuf, "%lf %lf\n", &time, &mmi);
        printf("%14.6f %14.2f\n", time+offset, mmi);
    }
    return 0;
}