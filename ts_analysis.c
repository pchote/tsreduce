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

#include "tsreduce.h"
#include "helpers.h"
#include "fit.h"

// Display the apertures specified by the reduction file at dataPath for the
// observation index obsIndex over an image frame in ds9
int display_targets(char *dataPath, int obsIndex)
{
    // Read file header
    datafile data = read_data_header(dataPath);
    if (data.file == NULL)
        return error("Error opening data file");

    if (obsIndex >= data.num_obs)
    {
        fclose(data.file);
        return error("Requested observation is out of range: max is %d", data.num_obs-1);
    }

    chdir(data.frame_dir);

    if (!init_ds9("tsreduce"))
    {
        fclose(data.file);
        return error("Unable to launch ds9");
    }

    char command[128];
    char filenamebuf[NAME_MAX];
    filenamebuf[0] = '\0';

    // Observation index not specified: open the first image that matches
    if (obsIndex >= 0)
        strncpy(filenamebuf, data.obs[obsIndex].filename, NAME_MAX);
    else if (!get_first_matching_file(data.frame_pattern, filenamebuf, NAME_MAX))
    {
        fclose(data.file);
        return error("No matching files found");
    }

    snprintf(command, 128, "file %s/%s", data.frame_dir, filenamebuf);
    if (tell_ds9("tsreduce", command, NULL, 0))
    {
        fclose(data.file);
        return error("ds9 command failed: %s", command);
    }

    // Set scaling mode
    if (tell_ds9("tsreduce", "scale mode 99.5", NULL, 0))
    {
        fclose(data.file);
        return error("ds9 command failed: scale mode 99.5");
    }

    // Flip X axis
    if (tell_ds9("tsreduce", "orient x", NULL, 0))
    {
        fclose(data.file);
        return error("ds9 command failed: orient x");
    }

    for (int i = 0; i < data.num_targets; i++)
    {
        double x, y;
        // Use the specified target aperture
        // DS9 starts labelling pixels at 1
        if (obsIndex < 0)
        {
            x = data.targets[i].x + 1;
            y = data.targets[i].y + 1;
        }
        // Use the centered observation
        else
        {
            x = data.obs[obsIndex].pos[i].x + 1;
            y = data.obs[obsIndex].pos[i].y + 1;
        }

        snprintf(command, 128, "regions command {circle %f %f %f #color=red}", x, y, data.targets[i].r);
        if (tell_ds9("tsreduce", command, NULL, 0))
            fprintf(stderr, "ds9 command failed: %s\n", command);
        snprintf(command, 128, "regions command {circle %f %f %f #background}", x, y, data.targets[i].s1);
        if (tell_ds9("tsreduce", command, NULL, 0))
            fprintf(stderr, "ds9 command failed: %s\n", command);
        snprintf(command, 128, "regions command {circle %f %f %f #background}", x, y, data.targets[i].s2);
        if (tell_ds9("tsreduce", command, NULL, 0))
            fprintf(stderr, "ds9 command failed: %s\n", command);
    }

    fclose(data.file);
    return 0;
}

// Output radial profile information for the given targetIndex, obsIndex in
// the reduction file at dataPath
int calculate_profile(char *dataPath, int obsIndex, int targetIndex)
{
    // Read file header
    datafile data = read_data_header(dataPath);
    if (data.file == NULL)
        return error("Error opening data file");

    if (obsIndex >= data.num_obs)
    {
        fclose(data.file);
        return error("Requested observation is out of range: max is %d", data.num_obs-1);
    }

    chdir(data.frame_dir);

    char filenamebuf[NAME_MAX];
    if (obsIndex >= 0)
        strncpy(filenamebuf, data.obs[obsIndex].filename, NAME_MAX);
    else if (!get_first_matching_file(data.frame_pattern, filenamebuf, NAME_MAX))
    {
        fclose(data.file);
        return error("No matching files found");
    }

    framedata frame = framedata_new(filenamebuf, FRAMEDATA_DBL);
    if (data.dark_template != NULL)
    {
        framedata dark = framedata_new(data.dark_template, FRAMEDATA_DBL);
        framedata_subtract(&frame, &dark);
        framedata_free(dark);
    }

    if (data.flat_template != NULL)
    {
        framedata flat = framedata_new(data.flat_template, FRAMEDATA_DBL);
        framedata_divide(&frame, &flat);
        framedata_free(flat);
    }

    if (targetIndex < 0 || targetIndex >= data.num_targets)
    {
        fclose(data.file);
        return error("Invalid target `%d' selected", targetIndex);
    }

    target t = data.targets[targetIndex];
    double2 xy = converge_aperture(t, &frame);
    t.x = xy.x; t.y = xy.y;
    double2 bg = calculate_background(t, &frame);

    const int numIntensity = 21;
    double intensity[numIntensity];
    double radii[numIntensity];
    double profile[numIntensity];

    // Calculate the remaining integrated intensities
    for (int i = 1; i < numIntensity; i++)
    {
        radii[i] = i;
        intensity[i] = integrate_aperture(xy, radii[i], &frame) - bg.x*M_PI*radii[i]*radii[i];
    }

    // Normalize integrated count by area to give an intensity profile
    // r = 0 value is sampled from the central pixel directly
    radii[0] = 0;
    intensity[0] = 0;
    profile[0] = frame.dbl_data[frame.cols*((int)xy.y) + (int)xy.x];

    // Central integrated value is a disk
    profile[1] = intensity[1]/(M_PI*radii[1]*radii[1]);

    // Remaining areas are annuli
    for (int i = 2; i < numIntensity; i++)
    {
        double area = M_PI*(radii[i]*radii[i] - radii[i-1]*radii[i-1]);
        profile[i] = (intensity[i] - intensity[i-1])/area;
    }

    // Print sky value
    printf("# Sky background: %f\n", bg.x);
    printf("# Sky stddev: %f\n", bg.y);

    // Estimate FWHM by linear interpolation between points
    for (int i = 1; i < numIntensity; i++)
        if (profile[i] < profile[0]/2)
        {
            double fwhm = i - 1 + (profile[0]/2 - profile[i-1])/(profile[i] - profile[i-1]);
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
        if (profile[i] < 5*bg.y)
        {
            printf("# Estimated 5 sky sigma: %f\n", i - 1 + (5*bg.y - profile[i-1])/(profile[i] - profile[i-1]));
            break;
        }
    for (int i = 1; i < numIntensity; i++)
        if (profile[i] < 10*bg.y)
        {
            printf("# Estimated 10 sky sigma: %f\n", i - 1 + (10*bg.y - profile[i-1])/(profile[i] - profile[i-1]));
            break;
        }

    // Print profile values
    for (int i = 0; i < numIntensity; i++)
        printf("%f %f %f\n", radii[i], profile[i], intensity[i]);

    return 0;
}


// List the timestamps/filenames that are corrupted by continous downloads
int detect_repeats(char *dataPath)
{
    // Read file header
    datafile data = read_data_header(dataPath);
    if (data.file == NULL)
        return error("Error opening data file");

    // No data
    if (data.num_obs <= 0)
    {
        fclose(data.file);
        return error("File specifies no observations");
    }

    int last_bad = 0;
    double last_time = data.obs[0].time;
    for (int i = 1; i < data.num_obs; i++)
    {
        double time = data.obs[i].time;
        if (time - last_time < 0.1f)
            last_bad = 1;

        if (last_bad)
            printf("%s @ %.1f\n", data.obs[i].filename, time);

        if (last_bad && time - last_time >= 0.1f)
            last_bad = 0;

        last_time = time;
    }

    fclose(data.file);
    return 0;
}

int plot_fits(char *dataPath)
{
    int plot_colors_max = 8;
    int plot_colors[] = {4,2,8,3,5,6,7,9};

    // Read file header
    datafile data = read_data_header(dataPath);
    if (data.file == NULL)
        return error("Error opening data file");

    chdir(data.frame_dir);

    // No data
    if (data.num_obs <= 0)
    {
        fclose(data.file);
        return error("File specifies no observations");
    }

    // Start time in hours
    struct tm starttime;
    gmtime_r(&data.reference_time, &starttime);
    float min_time = starttime.tm_hour + starttime.tm_min / 60.0 + starttime.tm_sec / 3600.0;

    // Time Series data
    float *time = (float *)malloc(data.num_obs*sizeof(float));
    if (time == NULL)
        return error("malloc failed");

    float *raw = (float *)malloc(data.num_obs*data.num_targets*sizeof(float));
    if (raw == NULL)
        return error("malloc failed");

    float *ratio = (float *)malloc(data.num_obs*sizeof(float));
    if (ratio == NULL)
        return error("malloc failed");

    float *polyfit = (float *)malloc(data.num_obs*sizeof(float));
    if (polyfit == NULL)
        return error("malloc failed");

    // Calculate polynomial fit to the ratio
    double *coeffs = (double *)malloc((data.plot_fit_degree+1)*sizeof(double));
    if (coeffs == NULL)
        return error("malloc failed");

    double ratio_mean = 0;
    for (int i = 0; i < data.num_obs; i++)
    {
        time[i] = data.obs[i].time;
        ratio[i] = data.obs[i].ratio;
        ratio_mean += ratio[i];

        for (int j = 0; j < data.num_targets; j++)
            raw[j*data.num_obs + i] = data.obs[i].star[j];
    }
    float max_time = min_time + time[data.num_obs-1]/3600;

    ratio_mean /= data.num_obs;

    // Calculate standard deviation
    double ratio_std = 0;
    for (int i = 0; i < data.num_obs; i++)
        ratio_std += (ratio[i] - ratio_mean)*(ratio[i] - ratio_mean);
    ratio_std = sqrt(ratio_std/data.num_obs);
    float min_ratio = (ratio_mean - 5*ratio_std);
    float max_ratio = (ratio_mean + 5*ratio_std);


    if (fit_polynomial(time, ratio, data.num_obs, coeffs, data.plot_fit_degree))
    {
        free(coeffs);
        free(ratio);
        free(raw);
        free(time);
        fclose(data.file);
        return error("Fit failed");
    }

    float *mmi = (float *)malloc(data.num_obs*sizeof(float));
    if (mmi == NULL)
        return error("malloc failed");

    double mmi_mean = 0;
    for (int i = 0; i < data.num_obs; i++)
    {
        // Subtract polynomial fit and convert to mmi
        polyfit[i] = 0;
        double pow = 1;
        for (int j = 0; j <= data.plot_fit_degree; j++)
        {
            polyfit[i] += pow*coeffs[j];
            pow *= time[i];
        }
        mmi[i] = ratio[i] > 0 ? 1000*(ratio[i] - polyfit[i])/ratio[i] : 0;
        mmi_mean += mmi[i];
    }
    mmi_mean /= data.num_obs;

    // Calculate standard deviation
    double mmi_std = 0;
    for (int i = 0; i < data.num_obs; i++)
        mmi_std += (mmi[i] - mmi_mean)*(mmi[i] - mmi_mean);
    mmi_std = sqrt(mmi_std/data.num_obs);

    double mmi_corrected_mean = 0;
    int mmi_corrected_count = 0;

    // Discard outliers and recalculate mean
    for (int i = 0; i < data.num_obs; i++)
    {
        if (fabs(mmi[i] - mmi_mean) > 3*mmi_std)
            mmi[i] = 0;
        else
        {
            mmi_corrected_mean += mmi[i];
            mmi_corrected_count++;
        }
    }
    mmi_corrected_mean /= mmi_corrected_count;

    float min_mmi = (mmi_corrected_mean - 5*mmi_std);
    float max_mmi = (mmi_corrected_mean + 5*mmi_std);

    if (cpgopen("5/xs") <= 0)
        return error("Unable to open PGPLOT window");

    // 800 x 480
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

    cpgmtxt("l", 2.5, 0.5, 0.5, "MMI");
    cpgswin(min_time, max_time, min_mmi, max_mmi);
    cpgbox("bcstm", 1, 4, "bcstn", 0, 0);

    cpgswin(0, time[data.num_obs-1], min_mmi, max_mmi);
    cpgpt(data.num_obs, time, mmi, 2);

    // Ratio
    cpgsvp(0.1, 0.9, 0.55, 0.75);
    cpgmtxt("l", 2.5, 0.5, 0.5, "Ratio");
    cpgswin(min_time, max_time, min_ratio, max_ratio);
    cpgbox("bcst", 1, 4, "bcstn", 0, 0);

    cpgswin(0, time[data.num_obs-1], min_ratio, max_ratio);
    cpgpt(data.num_obs, time, ratio, 2);

    // Plot the polynomial fit
    cpgsci(2);
    cpgline(data.num_obs, time, polyfit);
    cpgsci(1);

    // Raw Data
    cpgsvp(0.1, 0.9, 0.075, 0.55);
    cpgmtxt("l", 2.5, 0.5, 0.5, "Counts Per Second");
    cpgmtxt("b", 2.5, 0.5, 0.5, "UTC Hour");
    cpgswin(min_time, max_time, 0, data.plot_max_raw);
    cpgbox("bcstn", 1, 4, "bcstn", 0, 0);

    for (int j = 0; j < data.num_targets; j++)
    {
        cpgswin(0,time[data.num_obs-1], 0, data.plot_max_raw/data.targets[j].plot_scale);
        cpgsci(plot_colors[j%plot_colors_max]);
        cpgpt(data.num_obs, time, &raw[j*data.num_obs], 2);
    }
    cpgend();

    // DFT data
    float *freq = (float *)malloc(data.plot_num_uhz*sizeof(float));
    if (freq == NULL)
        return error("malloc failed");

    float *ampl = (float *)malloc(data.plot_num_uhz*sizeof(float));
    if (ampl == NULL)
        return error("malloc failed");

    calculate_amplitude_spectrum_float(data.plot_min_uhz*1e-6, data.plot_max_uhz*1e-6, time, mmi, data.num_obs, freq, ampl, data.plot_num_uhz);

    // Determine max dft ampl
    float max_dft_ampl = 0;
    for (int i = 0; i < data.plot_num_uhz; i++)
        max_dft_ampl = fmax(max_dft_ampl, ampl[i]);
    max_dft_ampl *= 1.1;

    if (cpgopen("6/xs") <= 0)
        return error("Unable to open PGPLOT window");

    // 800 x 480
    cpgpap(9.41, 0.6);
    cpgask(0);
    cpgslw(3);
    cpgsfs(2);
    cpgscf(2);

    // DFT
    cpgsvp(0.1, 0.9, 0.075, 0.9);
    cpgswin(data.plot_min_uhz, data.plot_max_uhz, 0, max_dft_ampl);
    cpgbox("bcstn", 0, 0, "bcnst", 10, 2);

    cpgswin(data.plot_min_uhz*1e-6, data.plot_max_uhz*1e-6, 0, max_dft_ampl);
    cpgsci(12);
    cpgline(data.plot_num_uhz, freq, ampl);
    cpgsci(1);

    cpgmtxt("b", 2.5, 0.5, 0.5, "Frequency (\\gmHz)");
    cpgmtxt("l", 2, 0.5, 0.5, "Amplitude (mma)");

    cpgsch(1.25);
    cpgmtxt("t", 1.5, 0.5, 0.5, "Amplitude Spectrum");

    free(coeffs);
    free(ampl);
    free(freq);
    free(mmi);
    free(ratio);
    free(raw);
    free(time);
    fclose(data.file);

    return 0;
}

int amplitude_spectrum(char *dataPath)
{
    // Read file header
    datafile data = read_data_header(dataPath);
    if (data.file == NULL)
        return error("Error opening data file");

    chdir(data.frame_dir);

    // No data
    if (data.num_obs <= 0)
    {
        fclose(data.file);
        return error("File specifies no observations");
    }

    // Time Series data
    float *time = (float *)malloc(data.num_obs*sizeof(float));
    if (time == NULL)
        return error("malloc failed");

    float *raw = (float *)malloc(data.num_obs*data.num_targets*sizeof(float));
    if (raw == NULL)
        return error("malloc failed");

    float *ratio = (float *)malloc(data.num_obs*sizeof(float));
    if (ratio == NULL)
        return error("malloc failed");

    float *polyfit = (float *)malloc(data.num_obs*sizeof(float));
    if (polyfit == NULL)
        return error("malloc failed");

    // Calculate polynomial fit to the ratio
    double *coeffs = (double *)malloc((data.plot_fit_degree+1)*sizeof(double));
    if (coeffs == NULL)
        return error("malloc failed");

    double ratio_mean = 0;
    for (int i = 0; i < data.num_obs; i++)
    {
        time[i] = data.obs[i].time;
        ratio[i] = data.obs[i].ratio;
        ratio_mean += ratio[i];

        for (int j = 0; j < data.num_targets; j++)
            raw[j*data.num_obs + i] = data.obs[i].star[j];
    }

    ratio_mean /= data.num_obs;

    // Calculate standard deviation
    double ratio_std = 0;
    for (int i = 0; i < data.num_obs; i++)
        ratio_std += (ratio[i] - ratio_mean)*(ratio[i] - ratio_mean);
    ratio_std = sqrt(ratio_std/data.num_obs);

    if (fit_polynomial(time, ratio, data.num_obs, coeffs, data.plot_fit_degree))
    {
        free(coeffs);
        free(ratio);
        free(raw);
        free(time);
        fclose(data.file);
        return error("Fit failed");
    }

    float *mmi = (float *)malloc(data.num_obs*sizeof(float));
    if (mmi == NULL)
        return error("malloc failed");

    double mmi_mean = 0;
    for (int i = 0; i < data.num_obs; i++)
    {
        // Subtract polynomial fit and convert to mmi
        polyfit[i] = 0;
        double pow = 1;
        for (int j = 0; j <= data.plot_fit_degree; j++)
        {
            polyfit[i] += pow*coeffs[j];
            pow *= time[i];
        }
        mmi[i] = 1000*(ratio[i] - polyfit[i])/ratio[i];
        mmi_mean += mmi[i];
    }
    mmi_mean /= data.num_obs;

    // Calculate standard deviation
    double mmi_std = 0;
    for (int i = 0; i < data.num_obs; i++)
        mmi_std += (mmi[i] - mmi_mean)*(mmi[i] - mmi_mean);
    mmi_std = sqrt(mmi_std/data.num_obs);

    double mmi_corrected_mean = 0;
    int mmi_corrected_count = 0;

    // Discard outliers and recalculate mean
    for (int i = 0; i < data.num_obs; i++)
    {
        if (fabs(mmi[i] - mmi_mean) > 3*mmi_std)
            mmi[i] = 0;
        else
        {
            mmi_corrected_mean += mmi[i];
            mmi_corrected_count++;
        }
    }
    mmi_corrected_mean /= mmi_corrected_count;

    // DFT data
    float *freq = (float *)malloc(data.plot_num_uhz*sizeof(float));
    if (freq == NULL)
        return error("malloc failed");

    float *ampl = (float *)malloc(data.plot_num_uhz*sizeof(float));
    if (ampl == NULL)
        return error("malloc failed");

    calculate_amplitude_spectrum_float(data.plot_min_uhz*1e-6, data.plot_max_uhz*1e-6, time, mmi, data.num_obs, freq, ampl, data.plot_num_uhz);

    for (int i = 0; i < data.plot_num_uhz; i++)
        printf("%f %f\n", freq[i], ampl[i]);

    free(coeffs);
    free(ampl);
    free(freq);
    free(mmi);
    free(ratio);
    free(raw);
    free(time);
    fclose(data.file);

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
