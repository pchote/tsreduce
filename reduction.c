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
#include <dirent.h>
#include <regex.h>
#include <stdbool.h>

#include "datafile.h"
#include "framedata.h"
#include "helpers.h"
#include "aperture.h"
#include "fit.h"
#include "astro_convert.h"

extern int verbosity;

int generate_photometry_dft_data(datafile *data,
    // Raw data, unfiltered by blocked ranges
    double **raw_time, double **raw, double **mean_sky, size_t *num_raw,
    // Filtered data: polynomial fit, ratio, mma
    double **time, double **ratio, double **polyfit, double **mma, double **fwhm, size_t *num_filtered,
    double **ratio_noise, double **mma_noise,
    double *ratio_mean_out, double *ratio_std_out, double *mma_mean_out, double *mma_std_out,
    // Fourier transform; set to NULL to skip calculation
    double **freq, double **ampl, size_t *num_dft)
{
    int ret = 0;

    // No data
    if (data->num_obs <= 0)
        return error("File specifies no observations");

    // Time Series data
    *raw_time = (double *)malloc(data->num_obs*sizeof(double));
    if (*raw_time == NULL)
        error_jump(raw_time_alloc_error, ret, "raw_time malloc failed");

    *raw = (double *)malloc(data->num_obs*data->num_targets*sizeof(double));
    if (*raw == NULL)
        error_jump(raw_alloc_error, ret, "raw malloc failed");

    *mean_sky = (double *)malloc(data->num_obs*sizeof(double));
    if (*mean_sky == NULL)
        error_jump(mean_sky_alloc_error, ret, "mean_sky malloc failed");

    *time = (double *)malloc(data->num_obs*sizeof(double));
    if (*time == NULL)
        error_jump(time_alloc_error, ret, "time malloc failed");

    *ratio = (double *)malloc(data->num_obs*sizeof(double));
    if (*ratio == NULL)
        error_jump(ratio_alloc_error, ret, "ratio malloc failed");

    *ratio_noise = (double *)malloc(data->num_obs*sizeof(double));
    if (*ratio == NULL)
        error_jump(ratio_noise_alloc_error, ret, "ratio malloc failed");

    *fwhm = (double *)malloc(data->num_obs*sizeof(double));
    if (*fwhm == NULL)
        error_jump(fwhm_alloc_error, ret, "ratio malloc failed");

    *polyfit = (double *)malloc(data->num_obs*sizeof(double));
    if (*polyfit == NULL)
        error_jump(polyfit_alloc_error, ret, "polyfit malloc failed");

    // Calculate polynomial fit to the ratio
    double *coeffs = (double *)calloc(data->plot_fit_degree + 1, sizeof(double));
    if (coeffs == NULL)
        error_jump(coeffs_alloc_error, ret, "coeffs malloc failed");

    double ratio_mean = 0;
    *num_filtered = 0;
    *num_raw = data->num_obs;
    for (int i = 0; i < data->num_obs; i++)
    {
        (*raw_time)[i] = data->obs[i].time;

        for (int j = 0; j < data->num_targets; j++)
        {
            (*raw)[j*data->num_obs + i] = data->obs[i].star[j];
            (*mean_sky)[i] += data->obs[i].sky[j]/data->num_targets;
        }

        // Filter bad observations
        bool skip = false;
        for (int j = 0; j < data->num_blocked_ranges; j++)
            if ((*raw_time)[i] >= data->blocked_ranges[j].x && (*raw_time)[i] <= data->blocked_ranges[j].y)
            {
                skip = true;
                break;
            }

        // Invalid observations have noise = nan
        if (isnan(data->obs[i].ratio_noise))
        {
            printf("Skipping observation at %f\n", data->obs[i].time);
            skip = true;
        }
        if (!skip)
        {
            // Calculate ratio from raw data, ignoring the value in the data file
            double target = data->obs[i].star[0];
            double comparison = 0;
            for (int j = 1; j < data->num_targets; j++)
                comparison += data->obs[i].star[j];

            // Skip invalid observations
            if (target == 0 || comparison == 0)
            {
                if (verbosity >= 1)
                    error("Ignoring bad observation at %f", data->obs[i].time);
                continue;
            }

            (*time)[*num_filtered] = data->obs[i].time;
            (*ratio)[*num_filtered] = comparison > 0 ? target/comparison : 0;

            // Read noise and fwhm from data file if available
            (*ratio_noise)[*num_filtered] = (data->version >= 5 && data->dark_template) ? data->obs[i].ratio_noise : 0;
            (*fwhm)[*num_filtered] = (data->version >= 6) ? data->obs[i].fwhm : 0;

            ratio_mean += (*ratio)[*num_filtered];
            (*num_filtered)++;
        }
    }

    ratio_mean /= *num_filtered;
    if (ratio_mean_out)
        *ratio_mean_out = ratio_mean;

    // Calculate standard deviation
    double ratio_std = 0;
    for (int i = 0; i < *num_filtered; i++)
        ratio_std += ((*ratio)[i] - ratio_mean)*((*ratio)[i] - ratio_mean);
    ratio_std = sqrt(ratio_std/(*num_filtered));
    if (ratio_std_out)
        *ratio_std_out = ratio_std;

    // Only attempt fit if it makes sense
    if (*num_filtered <= data->plot_fit_degree)
    {
        *num_filtered = 0;
        *mma_mean_out = 0;
        *mma = NULL;
        *mma_noise = NULL;
        *freq = NULL;
        *ampl = NULL;
        *num_dft = 0;
        return 0;
    }

    double *fit_noise = (data->version >= 5 && data->dark_template) ? *ratio_noise : NULL;
    if (fit_polynomial(*time, *ratio, fit_noise, *num_filtered, coeffs, data->plot_fit_degree))
        error_jump(poly_fit_error, ret, "Polynomial fit failed");

    *mma = (double *)malloc(*num_filtered*sizeof(double));
    if (*mma == NULL)
        error_jump(mma_alloc_error, ret, "mma malloc failed");

    *mma_noise = (double *)malloc(*num_filtered*sizeof(double));
    if (*mma_noise == NULL)
        error_jump(mma_noise_alloc_error, ret, "ratio malloc failed");

    double mma_mean = 0;
    for (int i = 0; i < *num_filtered; i++)
    {
        // Subtract polynomial fit and convert to mma
        (*polyfit)[i] = 0;
        double pow = 1;
        for (int j = 0; j <= data->plot_fit_degree; j++)
        {
            (*polyfit)[i] += pow*coeffs[j];
            pow *= (*time)[i];
        }
        (*mma)[i] = 1000*((*ratio)[i] - (*polyfit)[i])/(*polyfit)[i];

        if (data->version >= 5 && data->dark_template)
        {
            double numer_error = fabs((*ratio_noise)[i]/((*ratio)[i] - (*polyfit)[i]));
            double denom_error = fabs((*ratio_noise)[i]/(*ratio)[i]);
            (*mma_noise)[i] = (numer_error + denom_error)*fabs((*mma)[i]);
        }
        else
            (*mma_noise)[i] = 0;

        mma_mean += (*mma)[i];
    }
    mma_mean /= *num_filtered;

    // Calculate standard deviation
    double mma_std = 0;
    for (int i = 0; i < *num_filtered; i++)
        mma_std += ((*mma)[i] - mma_mean)*((*mma)[i] - mma_mean);
    mma_std = sqrt(mma_std/(*num_filtered));

    double mma_corrected_mean = 0;
    int mma_corrected_count = 0;

    // Discard outliers and recalculate mean
    for (int i = 0; i < *num_filtered; i++)
    {
        if (fabs((*mma)[i] - mma_mean) > 3*mma_std)
        {
            if (verbosity >= 1)
                error("%f is an outlier, setting to 0", (*time)[i]);
            (*mma)[i] = 0;
        }
        else
        {
            mma_corrected_mean += (*mma)[i];
            mma_corrected_count++;
        }
    }
    mma_corrected_mean /= mma_corrected_count;
    if (mma_mean_out)
        *mma_mean_out = mma_corrected_mean;

    if (mma_std_out)
        *mma_std_out = mma_std;

    // Calculate DFT
    if (freq && ampl && num_dft)
    {
        *freq = (double *)malloc(data->plot_num_uhz*sizeof(double));
        if (*freq == NULL)
            error_jump(freq_alloc_error, ret, "freq malloc failed");

        *ampl = (double *)malloc(data->plot_num_uhz*sizeof(double));
        if (*ampl == NULL)
            error_jump(ampl_alloc_error, ret, "ampl malloc failed");
        *num_dft = data->plot_num_uhz;
        calculate_amplitude_spectrum(data->plot_min_uhz*1e-6, data->plot_max_uhz*1e-6,
                                           *time, *mma, *num_filtered,
                                           *freq, *ampl, *num_dft);
    }

    // Cleanup temporary memory
    free(coeffs);
    return 0;

    // Free allocated memory on error
ampl_alloc_error:
    free(*freq);
    *freq = NULL;
freq_alloc_error:
    free(*mma_noise);
    *mma_noise = NULL;
mma_noise_alloc_error:
    free(*mma);
    *mma = NULL;
mma_alloc_error:
poly_fit_error:
    free(coeffs);
coeffs_alloc_error:
    free(*polyfit);
    *polyfit = NULL;
polyfit_alloc_error:
    free(*fwhm);
    *fwhm = NULL;
fwhm_alloc_error:
    free(*ratio_noise);
    *ratio_noise = NULL;
ratio_noise_alloc_error:
    free(*ratio);
    *ratio = NULL;
ratio_alloc_error:
    free(*time);
    *time = NULL;
time_alloc_error:
    free(*mean_sky);
    *mean_sky = NULL;
mean_sky_alloc_error:
    free(*raw);
    *raw = NULL;
raw_alloc_error:
    free(*raw_time);
    *raw_time = NULL;
raw_time_alloc_error:
    return ret;
}

static time_t get_frame_time(framedata *frame)
{
    // Fits strings have a max length of FLEN_VALUE = 71
    char datebuf[128], timebuf[128], datetimebuf[257];
    if (framedata_has_header_string(frame, "UTC-BEG"))
    {
        framedata_get_header_string(frame, "UTC-DATE", datebuf);
        framedata_get_header_string(frame, "UTC-BEG", timebuf);
        snprintf(datetimebuf, 257, "%s %s", datebuf, timebuf);
    }
    else if (framedata_has_header_string(frame, "GPSTIME"))
        framedata_get_header_string(frame, "GPSTIME", datetimebuf);

    return parse_time_t(datetimebuf);
}

// Prepare a raw flat frame for combining into the master-flat
// Frame is dark (and bias if overscanned) subtracted, and normalized so the mean count is 1
//
// dark is expected to be bias subtracted if overscan is present
// This ensures that the bias is correctly removed, as the scaling of the dark intensity would otherwise prevent this
//
// Returns the mean intensity before normalization
double prepare_flat(framedata *flat, framedata *dark, double *mean_out)
{
    long flatexp, darkexp;
    if (framedata_get_header_long(flat, "EXPTIME", &flatexp))
        return error("EXPTIME is undefined in flat frame");

    if (framedata_get_header_long(dark, "EXPTIME", &darkexp))
        return error("EXPTIME is undefined in dark frame");

    // Subtract bias
    subtract_bias(flat);

    // Subtract dark, normalized to the flat exposure time
    double exp_ratio = flatexp*1.0/darkexp;
    for (int i = 0; i < flat->rows*flat->cols; i++)
        flat->data[i] -= exp_ratio*dark->data[i];

    // Calculate mean in image area only
    int *ir = flat->regions.image_region;
    double mean = mean_in_region(flat, ir);

    // Calculate standard deviation
    double std = 0;
    for (int j = ir[2]; j < ir[3]; j++)
        for (int i = ir[0]; i < ir[1]; i++)
        {
            double temp = flat->data[flat->cols*j + i] - mean;
            std += temp*temp;
        }
    std = sqrt(std/flat->regions.image_px);

    // Recalculate the mean, excluding outliers at 3 sigma
    double mean_new = 0;
    int count = 0;
    for (int j = ir[2]; j < ir[3]; j++)
        for (int i = ir[0]; i < ir[1]; i++)
            if (fabs(flat->data[flat->cols*j + i] - mean) < 3*std)
            {
                mean_new += flat->data[flat->cols*j + i];
                count++;
            }
    mean_new /= count;

    // Normalize flat so the image region mean is 1
    for (int i = 0; i < flat->rows*flat->cols; i++)
        flat->data[i] /= mean_new;

    // Return original mean level
    *mean_out = mean_new;
    return 0;
}

// Create a flat field frame from the frames listed by the command `flatcmd',
// rejecting `minmax' highest and lowest pixel values after subtracting
// the dark frame `masterdark'.
// Save the resulting image to the file `outname'
//
// Also calculates the readnoise and gain if overscan is available
int create_flat(const char *pattern, int minmax, const char *masterdark, const char *outname)
{
    int ret = 0;

    // Find the filenames that match the specified pattern
    char **frame_paths;
    int num_frames = get_matching_files(pattern, &frame_paths);
    if (num_frames < 0)
        error_jump(match_error, ret, "Error matching files");

    // Ensure there are enough frames to discard the requested number of pixels
    if (num_frames <= 2*minmax)
        error_jump(insufficient_frames, ret,
            "Insufficient frames. %d found, %d will be discarded", num_frames, 2*minmax);

    // Load the master dark frame
    // Frame geometry for all subsequent frames is assumed to match the master-dark
    framedata *dark = framedata_load(masterdark);
    if (!dark)
        error_jump(insufficient_frames, ret, "Error loading frame %s", masterdark);

    // Data cube for processing the flat data
    //        data[0] = frame[0][0,0], data[1] = frame[1][0,0] ... data[num_frames] = frame[0][1,0] etc
    double *data_cube = (double *)malloc(num_frames*dark->cols*dark->rows*sizeof(double));
    if (data_cube == NULL)
        error_jump(datacube_failed, ret, "data_cube alloc failed");

    // Mean intensity of each flat frame
    double *mean_flat = (double *)malloc(num_frames*sizeof(double));
    if (mean_flat == NULL)
        error_jump(meanflat_failed, ret, "mean_flat alloc failed");

    // We need to iterate over the individual frames multiple times to calculate gain
    framedata **frames = (framedata **)malloc(num_frames*sizeof(framedata*));
    if (frames == NULL)
        error_jump(frames_failed, ret, "frames alloc failed");

    // Load flat field frames into a data cube for the image data, and a cube for the bias data
    // data[0] = flat[0][0,0], data[1] = flat[1][0,0] ... data[num_frames] = flat[0][0,1] etc
    for(int i = 0; i < num_frames; i++)
    {
        if (verbosity >= 1)
            printf("loading `%s`\n", frame_paths[i]);
        frames[i] = framedata_load(frame_paths[i]);
        if (!frames[i])
        {
            // Cleanup the frames that we have allocated
            for (int j = 0; j < i; j++)
                framedata_free(frames[j]);

            error_jump(loadflat_failed, ret, "Error loading frame %s", frame_paths[i]);
        }

        if (frames[i]->rows != dark->rows || frames[i]->cols != dark->cols)
        {
            // Cleanup the frames that we have allocated
            for (int j = 0; j <= i; j++)
                framedata_free(frames[j]);

            error_jump(loadflat_failed, ret,
                "Frame %s dimensions mismatch. Expected (%d,%d), was (%d, %d)",
                frames[i], dark->rows, dark->cols, frames[i]->rows, frames[i]->cols);
        }

        // Dark-subtract and normalize the frame, recording the mean level for calculating the gain
        if (prepare_flat(frames[i], dark, &mean_flat[i]))
        {
            // Cleanup the frames that we have allocated
            for (int j = 0; j <= i; j++)
                framedata_free(frames[j]);

            error_jump(loadflat_failed, ret, "prepare_flat failed on frame %s", frames[i]);
        }

        // Store normalized data in data cube
        for (int j = 0; j < dark->rows*dark->cols; j++)
            data_cube[num_frames*j+i] = frames[i]->data[j];
    }

    // Calculate read noise from overscan region
    // Subtracting master-dark removed the mean level, so just add
    double readnoise = 0;
    if (dark->regions.has_overscan)
    {
        int *br = dark->regions.bias_region;
        double var = 0;
        for (int j = br[2]; j < br[3]; j++)
            for (int i = br[0]; i < br[1]; i++)
                for (int k = 0; k < num_frames; k++)
                {
                    // Need to multiply by the mean level for the frame to give the variance in ADU
                    double temp = mean_flat[k]*data_cube[num_frames*(dark->cols*j + i) + k];
                    var += temp*temp;
                }
        readnoise = sqrt(var/(dark->regions.bias_px*num_frames));
    }

    // Calculate median image for the master-flat
    // Loop over the pixels, sorting the values from each image into increasing order
    double *median_flat = (double *)malloc(dark->rows*dark->cols*sizeof(double));
    if (median_flat == NULL)
        error_jump(median_failed, ret, "median_flat alloc failed");

    for (int j = 0; j < dark->rows*dark->cols; j++)
    {
        qsort(data_cube + num_frames*j, num_frames, sizeof(double), compare_double);

        // then average the non-rejected pixels into the output array
        median_flat[j] = 0;
        for (int i = minmax; i < num_frames - minmax; i++)
            median_flat[j] += data_cube[num_frames*j + i];
        median_flat[j] /= (num_frames - 2*minmax);
    }

    // Calculate gain from each image
    double *gain = (double *)malloc(num_frames*sizeof(double));
    if (gain == NULL)
        error_jump(gain_failed, ret, "gain alloc failed");

    if (dark->regions.has_overscan)
    {
        // Calculate mean dark level for gain calculation
        int *ir = dark->regions.image_region;
        double mean_dark = mean_in_region(dark, ir);

        for(int k = 0; k < num_frames; k++)
        {
            // Calculate the variance by taking the difference between
            // the (normalized) frame and the (normalized) master-flat
            double var = 0;
            for (int j = ir[2]; j < ir[3]; j++)
                for (int i = ir[0]; i < ir[1]; i++)
                {
                    double temp = mean_flat[k]*(frames[k]->data[dark->cols*j + i] - median_flat[dark->cols*j + i]);
                    var += temp*temp;
                }
            var /= dark->regions.image_px;
            gain[k] = (mean_flat[k] + mean_dark) / (var - readnoise*readnoise);

            if (verbosity >= 1)
                printf("%d var: %f mean: %f dark: %f gain: %f\n", k, var, mean_flat[k], mean_dark, gain[k]);
        }
    }

    // Find median value
    qsort(gain, num_frames, sizeof(double), compare_double);
    double median_gain = gain[num_frames/2];

    // Find mean value
    double mean_gain = 0;
    for(int k = 0; k < num_frames; k++)
        mean_gain += gain[k];
    mean_gain /= num_frames;

    // Replace values outside the image region with 1, so overscan survives flatfielding
    if (dark->regions.has_overscan)
    {
        int *br = dark->regions.image_region;
        for (int k = 0; k < dark->rows*dark->cols; k++)
        {
            int x = k % dark->cols;
            int y = k / dark->cols;
            bool in_image = x >= br[0] && x < br[1] && y >= br[2] && y < br[3];
            if (!in_image)
                median_flat[k] = 1;
        }
    }

    // Create a new fits file
    fitsfile *out;
    int status = 0;
    size_t filename_len = strlen(outname) + 2;
    char *filename = malloc(filename_len*sizeof(char));
    snprintf(filename, filename_len, "!%s", outname);
    fits_create_file(&out, filename, &status);
    free(filename);

    // Create the primary array image (16-bit short integer pixels
    fits_create_img(out, DOUBLE_IMG, 2, (long []){dark->cols, dark->rows}, &status);

    // Set header keys for readout noise and gain
    if (dark->regions.has_overscan)
    {
        fits_update_key(out, TDOUBLE, "CCD-READ", &readnoise, "Estimated read noise (ADU)", &status);
        fits_update_key(out, TDOUBLE, "CCD-GAIN", &median_gain, "Estimated gain (electrons/ADU)", &status);

        printf("Readnoise: %f\n", readnoise);
        printf("Gain: %f\n", median_gain);
    }

    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, dark->rows*dark->cols, median_flat, &status))
    {
        // Warn, but continue
        error("fits_write_img failed with status %d", status);
    }

    fits_close_file(out, &status);

    // Cleanup
    free(gain);
gain_failed:
    free(median_flat);
median_failed:
    for (int j = 0; j < num_frames; j++)
        framedata_free(frames[j]);
loadflat_failed:
    free(frames);
frames_failed:
    free(mean_flat);
meanflat_failed:
    free(data_cube);
datacube_failed:
    framedata_free(dark);
insufficient_frames:
    free_2d_array(frame_paths, num_frames);
match_error:
    return ret;
}

// Create a darkframe from the frames listed by the command `darkcmd',
// rejecting `minmax' highest and lowest pixel values.
// Save the resulting image to the file `outname'
int create_dark(const char *pattern, int minmax, const char *outname)
{
    int ret = 0;

    char **frame_paths;
    int num_frames = get_matching_files(pattern, &frame_paths);
    if (num_frames < 0)
        error_jump(match_error, ret, "Error matching files");

    if (num_frames < 2*minmax)
        error_jump(insufficient_frames, ret,
            "Insufficient frames. %d found, %d will be discarded", num_frames, 2*minmax);

    framedata *base = framedata_load(frame_paths[0]);
    if (!base)
        error_jump(insufficient_frames, ret, "Error loading frame %s", frame_paths[0]);

    long exptime;
    if (framedata_get_header_long(base, "EXPTIME", &exptime))
        error_jump(dark_failed, ret, "EXPTIME undefined in %s", frame_paths[0]);

    double *median_dark = (double *)malloc(base->rows*base->cols*sizeof(double));
    if (median_dark == NULL)
        error_jump(dark_failed, ret, "median_dark alloc failed");

    // Data cube for processing the flat data
    //        data[0] = frame[0][0,0], data[1] = frame[1][0,0] ... data[num_frames] = frame[0][1,0] etc
    double *data_cube = (double *)malloc(num_frames*base->cols*base->rows*sizeof(double));
    if (data_cube == NULL)
        error_jump(datacube_failed, ret, "data_cube alloc failed");

    for (int i = 0; i < num_frames; i++)
    {
        if (verbosity >= 1)
            printf("loading `%s`\n", frame_paths[i]);
        framedata *f = framedata_load(frame_paths[i]);
        if (!f)
            error_jump(loaddark_failed, ret, "Error loading frame %s", frame_paths[i]);

        if (f->rows != base->rows || f->cols != base->cols)
        {
            framedata_free(f);
            error_jump(loaddark_failed, ret,
                "Frame %s dimensions mismatch. Expected (%d,%d), was (%d, %d)",
                frame_paths[i], base->rows, base->cols, f->rows, f->cols);
        }

        subtract_bias(f);
        for (int j = 0; j < base->rows*base->cols; j++)
            data_cube[num_frames*j+i] = f->data[j];

        framedata_free(f);
    }
    
    // Loop over the pixels, sorting the values from each image into increasing order
    for (int j = 0; j < base->rows*base->cols; j++)
    {
        qsort(data_cube + num_frames*j, num_frames, sizeof(double), compare_double);

        // then average the non-rejected pixels into the output array
        median_dark[j] = 0;
        for (int i = minmax; i < num_frames - minmax; i++)
            median_dark[j] += data_cube[num_frames*j + i];
        median_dark[j] /= (num_frames - 2*minmax);
    }
    free(data_cube);

    // Create a new fits file
    fitsfile *out;
    int status = 0;
    
    size_t filename_len = strlen(outname) + 2;
    char *filename = malloc(filename_len*sizeof(char));
    snprintf(filename, filename_len, "!%s", outname);
    fits_create_file(&out, filename, &status);
    free(filename);
    
    // Create the primary array image (16-bit short integer pixels
    fits_create_img(out, DOUBLE_IMG, 2, (long []){base->cols, base->rows}, &status);
    fits_update_key(out, TLONG, "EXPTIME", &exptime, "Actual integration time (sec)", &status);
    
    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, base->rows*base->cols, median_dark, &status))
    {
        // Warn, but continue
        error("fits_write_img failed with status %d", status);
    }

    fits_close_file(out, &status);
loaddark_failed:
datacube_failed:
    free(median_dark);
dark_failed:
    framedata_free(base);
insufficient_frames:
    free_2d_array(frame_paths, num_frames);
match_error:
    return ret;
}

// Dark-subtract and flatfield the frame at `framePath' using the dark
// at `darkPath' and the flatfield at `flatPath'.
// Save the resulting image as a (double) floating-point fits file `outname'
int reduce_single_frame(char *framePath, char *darkPath, char *flatPath, char *outPath)
{
    int ret = 0;

    framedata *base = framedata_load(framePath);
    if (!base)
        error_jump(load_error, ret, "Error loading frame %s", framePath);

    framedata *dark = framedata_load(darkPath);
    if (!dark)
        error_jump(load_error, ret, "Error loading frame %s", darkPath);

    framedata *flat = framedata_load(flatPath);
    if (!flat)
        error_jump(load_error, ret, "Error loading frame %s", flatPath);

    // Process frame
    subtract_bias(base);
    framedata_subtract(base, dark);
    framedata_divide(base, flat);

    // Create a new fits file
    fitsfile *out;
    int status = 0;

    size_t filename_len = strlen(outPath) + 2;
    char *filename = malloc(filename_len*sizeof(char));
    snprintf(filename, filename_len, "!%s", outPath);
    fits_create_file(&out, filename, &status);
    free(filename);

    // Create the primary array image
    fits_create_img(out, DOUBLE_IMG, 2, (long []){base->cols, base->rows}, &status);

    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, base->rows*base->cols, base->data, &status))
        error_jump(write_error, ret, "fits_write_img failed with status %d", status);

write_error:
    fits_close_file(out, &status);

load_error:
    framedata_free(flat);
    framedata_free(dark);
    framedata_free(base);
    return ret;
}

// Load the reduction file at dataPath and reduce any new data
int update_reduction(char *dataPath)
{
    int ret = 0;

    // Read file header
    datafile *data = datafile_load(dataPath);
    if (data == NULL)
        return error("Error opening data file");

    if (data->version < 3)
        error_jump(data_error, ret, "Invalid data file version `%d'. Requires version >= 3", data->version);

    if (chdir(data->frame_dir))
        error_jump(data_error, ret, "Invalid frame path: %s", data->frame_dir);

    double readnoise = 0, gain = 1;
    framedata *flat = NULL;

    if (data->flat_template)
    {
        flat = framedata_load(data->flat_template);
        if (framedata_get_header_dbl(flat, "CCD-READ", &readnoise))
            readnoise = data->ccd_readnoise;

        if (readnoise <= 0)
            error_jump(flat_error, ret, "CCD Read noise unknown. Define CCDReadNoise in %s.", dataPath);

        if (framedata_get_header_dbl(flat, "CCD-GAIN", &gain))
            gain = data->ccd_gain;

        if (gain <= 0)
            error_jump(flat_error, ret, "CCD Gain unknown. Define CCDGain in %s.", dataPath);
    }

    framedata *dark = data->dark_template ? framedata_load(data->dark_template) : NULL;

    // Iterate through the files in the directory
    char **frame_paths;
    int start_obs = data->num_obs;
    int num_frames = get_matching_files(data->frame_pattern, &frame_paths);
    for (int i = 0; i < num_frames; i++)
    {
        // Check whether the frame has been processed
        int processed = FALSE;
        for (int j = 0; j < data->num_obs; j++)
            if (strcmp(frame_paths[i], data->obs[j].filename) == 0)
            {
                processed = TRUE;
                break;
            }

        if (processed)
            continue;

        printf("Reducing %s\n", frame_paths[i]);

        framedata *frame = framedata_load(frame_paths[i]);
        if (!frame)
        {
            framedata_free(frame);
            error_jump(process_error, ret, "Error loading frame %s", frame_paths[i]);
        }

        long exptime;
        if (framedata_get_header_long(frame, "EXPTIME", &exptime))
        {
            framedata_free(frame);
            error_jump(process_error, ret, "EXPTIME undefined in %s", frame_paths[i]);
        }

        // Calculate time at the start of the exposure relative to ReferenceTime
        time_t frame_time = get_frame_time(frame);
        double starttime = difftime(frame_time, data->reference_time);

        // Process frame
        subtract_bias(frame);
        if (dark)
            framedata_subtract(frame, dark);
        if (flat)
            framedata_divide(frame, flat);

        // Observation start time
        fprintf(data->file, "%.1f ", starttime);
        data->obs[data->num_obs].time = starttime;

        // Process frame
        double comparisonIntensity = 0;
        double targetIntensity = 0;
        double comparisonNoise = 0;
        double targetNoise = 0;
        double mean_fwhm = 0;

        bool failed = false;
        for (int i = 0; i < data->num_targets; i++)
        {
            // Use the aperture position from the previous frame
            // as a starting point if it is valid
            target t = data->targets[i];
            if (data->num_obs > 0)
            {
                double2 last = data->obs[data->num_obs-1].pos[i];
                if (last.x > 0 && last.x < frame->cols)
                    t.x = last.x;
                if (last.y > 0 && last.y < frame->rows)
                    t.y = last.y;
            }

            double sky = 0;
            double intensity = 0;
            double noise = 0;
            double2 xy = {0,0};

            if (!center_aperture(t, frame, &xy))
            {
                double bg = 0;
                if (calculate_background(t, frame, &bg, NULL))
                    bg = 0;

                // Integrate sky over the aperture and normalize per unit time
                sky = bg*M_PI*t.r*t.r / exptime;

                if (dark)
                    integrate_aperture_and_noise(xy, t.r, frame, dark, readnoise, gain, &intensity, &noise);
                else
                {
                    intensity = integrate_aperture(xy, t.r, frame);
                    noise = 0;
                }
                intensity = intensity/exptime - sky;
                noise /= exptime;

                if (data->version >= 6)
                {
                    // Calculate FWHM
                    double centerProfile = frame->data[frame->cols*((int)xy.y) + (int)xy.x] - bg;
                    double lastIntensity = 0;
                    double lastProfile = centerProfile;
                    int maxRadius = (int)t.s2 + 1;

                    double fwhm = 0;
                    for (int radius = 1; radius <= maxRadius; radius++)
                    {
                        double intensity = integrate_aperture(xy, radius, frame) - bg*M_PI*radius*radius;
                        double profile = (intensity - lastIntensity) / (M_PI*(2*radius-1));

                        if (profile < centerProfile/2)
                        {
                            double lastRadius = radius - 1;
                            fwhm = 2*(lastRadius + (radius - lastRadius)*(centerProfile/2 - lastProfile)/(profile - lastProfile));
                            break;
                        }

                        lastIntensity = intensity;
                        lastProfile = profile;
                    }

                    if (fwhm == 0)
                        failed = true;

                    mean_fwhm += fwhm / data->num_targets;
                }
            }
            else
                failed = true;

            fprintf(data->file, "%.2f ", intensity); // intensity (ADU/s)
            fprintf(data->file, "%.2f ", sky); // sky intensity (ADU/s)
            fprintf(data->file, "%.2f %.2f ", xy.x, xy.y); // Aperture center

            data->obs[data->num_obs].star[i] = intensity;
            data->obs[data->num_obs].sky[i] = sky;
            data->obs[data->num_obs].pos[i] = xy;

            if (i == 0)
            {
                targetIntensity = intensity;
                targetNoise = noise;
            }
            else
            {
                comparisonIntensity += intensity;
                comparisonNoise += noise;
            }
        }

        // Ratio
        double ratio = comparisonIntensity > 0 ? targetIntensity / comparisonIntensity : 0;
        double ratioNoise = failed ? sqrt(-1) : (targetNoise/targetIntensity + comparisonNoise/comparisonIntensity)*ratio;
        data->obs[data->num_obs].ratio = ratio;

        fprintf(data->file, "%.3e ", ratio);
        if (data->version >= 5)
            fprintf(data->file, "%.3e ", ratioNoise);

        if (data->version >= 6)
            fprintf(data->file, "%.3f ", mean_fwhm*data->ccd_platescale);

        // Filename
        fprintf(data->file, "%s\n", frame_paths[i]);
        strncpy(data->obs[data->num_obs].filename, frame_paths[i], 64);
        data->num_obs++;

        framedata_free(frame);
    }
    
    printf("Reduced %d observations\n", data->num_obs - start_obs);

process_error:
    free_2d_array(frame_paths, num_frames);
    if (dark)
        framedata_free(dark);
flat_error:
    if (flat)
        framedata_free(flat);
data_error:
    datafile_free(data);
    return ret;
}

int create_reduction_file(char *outname)
{
    int ret = 0;

    FILE *fileTest = fopen(outname, "w");
    if (fileTest == NULL)
        return error("Unable to create data file: %s. Does it already exist?", outname);
    fclose(fileTest);

    datafile *data = datafile_alloc();

    // Store the current directory so we can return before saving the data file
    char *datadir = getcwd(NULL, 0);

    while (true)
    {
        char *ret = prompt_user_input("Enter frame path", ".");
        data->frame_dir = canonicalize_path(ret);
        if (!chdir(data->frame_dir))
        {
            free(ret);
            break;
        }
        printf("Invalid frame path: %s\n", ret);
        free(data->frame_dir);
        free(ret);
    }

    char *input = prompt_user_input("Use calibration frames", "y");
    bool use_calibration = !strcmp(input, "y");
    free(input);

    if (use_calibration)
    {
        data->dark_template = prompt_user_input("Enter output master dark filename", "master-dark.fits.gz");

        // Create master-dark if necessary
        if (access(data->dark_template, F_OK) != -1)
            printf("Skipping master dark creation - file already exists\n");
        else
        {
            char dark_pattern[1039];
            int num_darks;
            while (true)
            {
                char *ret = prompt_user_input("Enter dark prefix", "dark");
                snprintf(dark_pattern, 1039, "%s-[0-9]+.fits.gz", ret);
                free(ret);

                char **dark_filenames;
                num_darks = get_matching_files(dark_pattern, &dark_filenames);
                if (num_darks > 0)
                {
                    free_2d_array(dark_filenames, num_darks);
                    break;
                }

                printf("No files found matching pattern: %s/%s\n", data->frame_dir, dark_pattern);
            }

            int minmax = 0;
            while (true)
            {
                char fallback[32];
                snprintf(fallback, 32, "%d", num_darks/2);
                char *ret = prompt_user_input("Enter number of darks around median to average", fallback);
                int count = atoi(ret);
                free(ret);

                if (count > 0 && count <= num_darks)
                {
                    minmax = (num_darks - count) / 2;
                    break;
                }
                printf("Number must be between 0 and %d\n", num_darks);
            }

            int failed = create_dark(dark_pattern, minmax, data->dark_template);
            if (failed)
                error_jump(create_dark_error, ret, "master dark generation failed");
        }

        data->flat_template = prompt_user_input("Enter output master flat filename", "master-flat.fits.gz");
    
        // Create master-flat if necessary
        if (access(data->flat_template, F_OK) != -1)
            printf("Skipping master flat creation - file already exists\n");
        else
        {
            char flat_pattern[1039];
            int num_flats;
            while (true)
            {
                char *ret = prompt_user_input("Enter flat prefix", "flat");
                snprintf(flat_pattern, 32, "%s-[0-9]+.fits.gz", ret);
                free(ret);

                char **flat_filenames;
                num_flats = get_matching_files(flat_pattern, &flat_filenames);
                if (num_flats > 0)
                {
                    free_2d_array(flat_filenames, num_flats);
                    break;
                }

                printf("No files found matching pattern: %s/%s\n", data->frame_dir, flat_pattern);
            }

            int minmax = 0;
            while (true)
            {
                char fallback[32];
                snprintf(fallback, 32, "%d", num_flats/2);
                char *ret = prompt_user_input("Enter number of flats around median to average", fallback);
                int count = atoi(ret);
                free(ret);

                if (count > 0 && count <= num_flats)
                {
                    minmax = (num_flats - count) / 2;
                    break;
                }
                printf("Number must be between 0 and %d\n", num_flats);
            }

            int failed = create_flat(flat_pattern, minmax, data->dark_template, data->flat_template);
            if (failed)
                error_jump(create_flat_error, ret, "master flat generation failed");
        }
    }
    else
    {
        data->dark_template = NULL;
        data->flat_template = NULL;
    }

    char *preview_filename;
    while (true)
    {
        char namebuf[1039];
        char *ret = prompt_user_input("Enter target prefix", NULL);
        snprintf(namebuf, 1039, "%s-[0-9]+.fits.gz", ret);
        free(ret);
        data->frame_pattern = strdup(namebuf);

        preview_filename = get_first_matching_file(data->frame_pattern);
        if (preview_filename)
            break;

        printf("No files found matching pattern: %s/%s\n", data->frame_dir, data->frame_pattern);
        free(preview_filename);
        free(data->frame_pattern);
    }

    // Open the file to find the reference time
    framedata *frame = framedata_load(preview_filename);
    if (!frame)
        error_jump(frameload_error, ret, "Error loading frame %s", preview_filename);

    subtract_bias(frame);
    data->reference_time = get_frame_time(frame);

    if (data->dark_template)
    {
        framedata *dark = framedata_load(data->dark_template);
        if (!dark)
            error_jump(frameload_error, ret, "Error loading frame %s", data->dark_template);

        framedata_subtract(frame, dark);
        framedata_free(dark);
    }

    if (data->flat_template)
    {
        framedata *flat = framedata_load(data->flat_template);
        if (!flat)
            error_jump(frameload_error, ret, "Error loading frame %s", data->flat_template);

        if (framedata_get_header_dbl(flat, "CCD-READ", (double[]){0}))
        {
            while (true)
            {
                char *ret = prompt_user_input("Enter CCD Readnoise (ADU):", "3.32");
                data->ccd_readnoise = atof(ret);
                free(ret);
                if (data->ccd_readnoise > 0)
                    break;

                printf("Number must be greater than 0\n");
            }
        }

        if (framedata_get_header_dbl(flat, "CCD-GAIN", (double[]){0}))
        {
            while (true)
            {
                char *ret = prompt_user_input("Enter CCD Gain (ADU):", "2.00");
                data->ccd_gain = atof(ret);
                free(ret);

                if (data->ccd_gain > 0)
                    break;

                printf("Number must be greater than 0\n");
            }
        }

        framedata_divide(frame, flat);
        framedata_free(flat);
    }

    while (true)
    {
        if (init_ds9())
            return error("Unable to launch ds9");

        {
            char command[128];
            snprintf(command, 128, "xpaset tsreduce array [xdim=%d,ydim=%d,bitpix=-64]", frame->cols, frame->rows);
            if (ts_exec_write(command, frame->data, frame->rows*frame->cols*sizeof(double)))
                error_jump(frameload_error, ret, "ds9 command failed: %s", command);
        }

        // Set scaling mode
        if (ts_exec_write("xpaset tsreduce scale mode zscale", NULL, 0))
            error_jump(frameload_error, ret, "ds9 command failed");

        // Flip X axis
        if (ts_exec_write("xpaset tsreduce orient x", NULL, 0))
            error_jump(frameload_error, ret, "ds9 command failed");

        // Set region mode to annulus
        if (ts_exec_write("xpaset tsreduce regions shape annulus", NULL, 0))
            error_jump(frameload_error, ret, "ds9 command failed");

        printf("Circle the target stars and surrounding sky in ds9\nPress enter in this terminal to continue...");
        getchar();

        // Read zoom from ds9
        char *ds9_zoom;
        if (ts_exec_read("xpaget tsreduce zoom", &ds9_zoom))
            error_jump(frameload_error, ret, "ds9 request zoom failed");
        float zoom = atof(ds9_zoom);
        free(ds9_zoom);

        char *ds9buf;
        if (ts_exec_read("xpaget tsreduce regions", &ds9buf))
            error_jump(frameload_error, ret, "ds9 request regions failed");

        // Parse the region definitions
        data->num_targets = 0;
        char *cur = ds9buf;
        double largest_aperture = 0;
        double sky[MAX_TARGETS];
        for (; (cur = strstr(cur, "annulus")) != NULL; cur++)
        {
            if (data->num_targets == MAX_TARGETS)
            {
                printf("Limit of %d targets reached. Remaining targets have been ignored", MAX_TARGETS);
                break;
            }

            // Read aperture coords
            target t;
            sscanf(cur, "annulus(%lf,%lf,%lf,%lf)", &t.x, &t.y, &t.s1, &t.s2);

            // ds9 denotes the bottom-left pixel (1,1), tsreduce uses (0,0)
            t.x -= 1;
            t.y -= 1;

            // Estimate initial aperture size as inner sky radius
            t.r = t.s1;

            if (verbosity >= 1)
                printf("Initial aperture xy: (%f,%f) r: %f s:(%f,%f)\n", t.x, t.y, t.r, t.s1, t.s2);

            double2 xy;
            if (center_aperture(t, frame, &xy))
            {
                printf("Centering failed to converge. Removing aperture.\n");
                continue;
            }

            t.x = xy.x;
            t.y = xy.y;

            double sky_intensity, sky_std_dev;
            if (calculate_background(t, frame, &sky_intensity, &sky_std_dev))
            {
                printf("Background calculation failed. Removing aperture.\n");
                continue;
            }
            sky[data->num_targets] = sky_intensity;

            // Estimate the radius where the star flux falls to 5 times the std. dev. of the background
            double lastIntensity = 0;
            double lastProfile = frame->data[frame->cols*((int)xy.y) + (int)xy.x] - sky_intensity;
            int maxRadius = (int)t.s2 + 1;
            for (int radius = 1; radius <= maxRadius; radius++)
            {
                double intensity = integrate_aperture(xy, radius, frame) - sky_intensity*M_PI*radius*radius;
                double profile = (intensity - lastIntensity) / (M_PI*(2*radius-1));
                if (profile < 5*sky_std_dev)
                {
                    // Linear interpolate radii to estimate the radius that gives 5*stddev
                    t.r = radius - 1 + (5*sky_std_dev - lastProfile) / (profile - lastProfile);
                    largest_aperture = fmax(largest_aperture, t.r);
                    break;
                }
                lastIntensity = intensity;
                lastProfile = profile;
            }

            // Set target parameters
            data->targets[data->num_targets++] = t;
        }
        free(ds9buf);

        // Set aperture radii to the same size, equal to the largest calculated above
        for (int i = 0; i < data->num_targets; i++)
            data->targets[i].r = largest_aperture;

        // Display results in ds9 - errors are non-fatal
        ts_exec_write("xpaset tsreduce regions delete all", NULL, 0);

        for (int i = 0; i < data->num_targets; i++)
        {
            double x = data->targets[i].x + 1;
            double y = data->targets[i].y + 1;
            double r = data->targets[i].r;
            double s1 = data->targets[i].s1;
            double s2 = data->targets[i].s2;

            char command[1024];
            snprintf(command, 1024, "xpaset tsreduce regions command '{circle %f %f %f #color=red select=0}'", x, y, r);
            ts_exec_write(command, NULL, 0);
            snprintf(command, 1024, "xpaset tsreduce regions command '{annulus %f %f %f %f #background select=0}'", x, y, s1, s2);
            ts_exec_write(command, NULL, 0);

            char msg[16];
            if (i == 0)
                strcpy(msg, "Target");
            else
                snprintf(msg, 16, "Comparison %d", i);

            double intensity = frame->data[frame->cols*((size_t)data->targets[i].y) + (size_t)data->targets[i].x] - sky[i];

            snprintf(command, 1024, "xpaset -p tsreduce regions command '{text %f %f #color=green select=0 text=\"%s\"}'", x, y - s2 - 10/zoom, msg);
            ts_exec_write(command, NULL, 0);
            snprintf(command, 1024, "xpaset -p tsreduce regions command '{text %f %f #color=green select=0 text=\"%.0f ADU\"}'", x, y - s2 - 25/zoom, intensity);
            ts_exec_write(command, NULL, 0);

            snprintf(command, 1024, "xpaset -p tsreduce regions command '{box %f %f %f %f #color=black select=0 width=%d}'",
                     x, y - s2 - 17.5/zoom, 115/zoom, 17.5/zoom, 18);
            ts_exec_write(command, NULL, 0);
        }
        ts_exec_write("xpaset -p tsreduce update now", NULL, 0);

        char *ret = prompt_user_input("Are the displayed apertures correct?:", "y");
        bool done = !strcmp(ret, "y");
        free(ret);
        if (done)
            break;
    }
    printf("Set %d targets\n", data->num_targets);

    // Save to disk
    if (chdir(datadir))
        error_jump(frameload_error, ret, "Invalid data path: %s", datadir);

    datafile_save_header(data, outname);
    printf("Saved to %s\n", outname);
    printf("Run `tsreduce update %s` to reduce existing files\n", outname);
    printf("Run `tsreduce plot %s` to preview reduced data\n", outname);
frameload_error:
    framedata_free(frame);
    free(preview_filename);
create_flat_error:
create_dark_error:
    datafile_free(data);
    return ret;
}

// Update the ds9 preview with a frame and FWHM calculations
int update_preview(char *preview_filename, char *ds9_title)
{
    int ret = 0;
    char ds9_command_buf[1024];

    framedata *frame = framedata_load(preview_filename);
    if (!frame)
        error_jump(frame_error, ret, "Error loading frame %s", preview_filename);
    subtract_bias(frame);

    double plate_scale = 1;
    if (framedata_get_header_dbl(frame, "IM-SCALE", &plate_scale))
        printf("IM-SCALE header key not found. Assuming 1px = 1arcsec.\n");

    // Read regions from ds9
    snprintf(ds9_command_buf, 1024, "xpaget %s regions", ds9_title);
    char *ds9_regions;
    if (ts_exec_read(ds9_command_buf, &ds9_regions))
        error_jump(region_error, ret, "ds9 request regions failed");

    // Read zoom from ds9
    snprintf(ds9_command_buf, 1024, "xpaget %s zoom", ds9_title);
    char *ds9_zoom;
    if (ts_exec_read(ds9_command_buf, &ds9_zoom))
        error_jump(region_error, ret, "ds9 request zoom failed");
    float zoom = atof(ds9_zoom);
    free(ds9_zoom);

    // Read binning from frame
    long binning = 1;
    framedata_get_header_long(frame, "CCD-BIN", &binning);

    // Set new frame
    snprintf(ds9_command_buf, 1024, "xpaset -p %s file %s", ds9_title, preview_filename);
    ts_exec_write(ds9_command_buf, NULL, 0);

    // Parse the region definitions
    char *cur = ds9_regions;
    for (; (cur = strstr(cur, "annulus")) != NULL; cur++)
    {
        // Read aperture coords
        target t;
        int select = 1;
        sscanf(cur, "annulus(%lf,%lf,%lf,%lf) # color=red select=%d", &t.x, &t.y, &t.s1, &t.s2, &select);
        if (!select)
            continue;

        // ds9 denotes the bottom-left pixel (1,1), tsreduce uses (0,0)
        t.x -= 1;
        t.y -= 1;

        // Estimate initial aperture size as inner sky radius
        t.r = t.s1;

        double2 xy;
        if (center_aperture(t, frame, &xy))
        {
            printf("Aperture converge failed. Removing target\n");
            continue;
        }
        t.x = xy.x;
        t.y = xy.y;

        double sky_intensity, sky_std_dev;
        if (calculate_background(t, frame, &sky_intensity, &sky_std_dev))
        {
            printf("Sky calculation failed. Removing target\n");
            continue;
        }

        double centerProfile = frame->data[frame->cols*((int)xy.y) + (int)xy.x] - sky_intensity;
        double lastIntensity = 0;
        double lastProfile = centerProfile;
        int maxRadius = (int)t.s2 + 1;

        double fwhm = 0;
        for (int radius = 1; radius <= maxRadius; radius++)
        {
            double intensity = integrate_aperture(xy, radius, frame) - sky_intensity*M_PI*radius*radius;
            double profile = (intensity - lastIntensity) / (M_PI*(2*radius-1));

            if (profile < centerProfile/2)
            {
                double lastRadius = radius - 1;
                fwhm = 2*(lastRadius + (radius - lastRadius)*(centerProfile/2 - lastProfile)/(profile - lastProfile));
                break;
            }

            lastIntensity = intensity;
            lastProfile = profile;
        }

        if (fwhm == 0)
        {
            printf("Invalid fwhm. Removing target\n");
            continue;
        }
        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{annulus %f %f %f %f}'", ds9_title, t.x + 1, t.y + 1, t.s1, t.s2);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{annulus %f %f %f %f #color=red select=0 background}'", ds9_title, t.x + 1, t.y + 1, t.s1, t.s2);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{circle %f %f %f #color=red select=0}'", ds9_title, t.x + 1, t.y + 1, fwhm);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{text %f %f #color=green select=0 text=\"FWHM: %.2f arcsec\"}'",
                 ds9_title, t.x + 1, t.y + 1 - t.s2 - 10/zoom, fwhm*binning*plate_scale);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{text %f %f #color=green select=0 text=\"Peak: %.0f ADU/px\"}'",
                 ds9_title, t.x + 1, t.y + 1 - t.s2 - 25/zoom, centerProfile);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{text %f %f #color=green select=0 text=\"BG: %.0f ADU/px\"}'",
                 ds9_title, t.x + 1, t.y + 1 - t.s2 - 40/zoom, sky_intensity);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{box %f %f %f %f #color=black select=0 width=%d}'",
                 ds9_title, t.x + 1, t.y + 1 - t.s2 - 25/zoom, 115/zoom, 25/zoom, 25);
        ts_exec_write(ds9_command_buf, NULL, 0);
    }

    // Display frame time
    char frame_end[128], frame_date[128], frame_object[128];
    long frame_exp;
    framedata_get_header_string(frame, "UTC-END", frame_end);
    framedata_get_header_string(frame, "UTC-DATE", frame_date);
    framedata_get_header_string(frame, "OBJECT", frame_object);
    framedata_get_header_long(frame, "EXPTIME", &frame_exp);

    snprintf(ds9_command_buf, 1024,
             "xpaset -p %s regions command '{text %f %f #color=green select=0 font=\"helvetica 12 bold roman\" text=\"%s @ %lds\"}'",
             ds9_title, frame->cols/2.0, frame->rows + 30/zoom, frame_object, frame_exp);
    ts_exec_write(ds9_command_buf, NULL, 0);

    snprintf(ds9_command_buf, 1024,
             "xpaset -p %s regions command '{text %f %f #color=green select=0 font=\"helvetica 12 bold roman\" text=\"Ending: %s %s\"}'",
             ds9_title, frame->cols/2.0, frame->rows + 10/zoom, frame_date, frame_end);
    ts_exec_write(ds9_command_buf, NULL, 0);

    // Force the display to update
    snprintf(ds9_command_buf, 1024, "xpaset -p %s update now", ds9_title);
    ts_exec_write(ds9_command_buf, NULL, 0);

region_error:
    framedata_free(frame);
frame_error:
    return ret;
}

/*
 * Calculate the BJD time for a given UTC timestamp and observation coordinates
 */
int calculate_bjd(char *date, char *time, char *ra_string, char *dec_string, double epoch)
{
    // Convert ra from HH:MM:SS to radians
    double a,b,c;
    sscanf(ra_string, "%lf:%lf:%lf", &a, &b, &c);
    double ra = (a + b/60 + c/3600)*M_PI/12;

    // Convert dec from DD:'':"" to radians
    sscanf(dec_string, "%lf:%lf:%lf", &a, &b, &c);
    double dec = (a + b/60 + c/3600)*M_PI/180;

    // Generate a struct tm for our reference time
    struct tm t = parse_date_time_tm(date, time);

    // Convert from UT to TT
    t.tm_sec += utcttoffset(&t);

    // Calculate BJD
    double2 reference_coords = precess((double2){ra, dec}, epoch, tmtoyear(&t));
    double reference_bjd = jdtobjd(tmtojd(&t), reference_coords);
    printf("%f\n", reference_bjd);

    return 0;
}

/*
 * Create a timeseries file from a list of datafiles.
 * Times are converted to BJD after the specified reference
 */
int create_ts(char *reference_date, char *reference_time, char **filenames, int num_datafiles, char *ts_filename)
{
    int ret = 0;
    if (num_datafiles < 1)
        return error("No datafiles specified");

    datafile **datafiles = (datafile **)malloc(num_datafiles*sizeof(datafile *));
    if (!datafiles)
        error_jump(datafile_error, ret, "Error allocating datafiles");

    for (int i = 0; i < num_datafiles; i++)
    {
        datafiles[i] = datafile_load(filenames[i]);
        if (!datafiles[i])
        {
            for (int j = 0; j < i; j++)
                datafile_free(datafiles[j]);
            free(datafiles);
            error_jump(datafile_error, ret, "Error loading datafile %s", filenames[i]);
        }
    }
    printf("Loaded %d datafiles\n", num_datafiles);

    if (!datafiles[0]->coord_ra || !datafiles[0]->coord_dec || datafiles[0]->coord_epoch == 0)
        error_jump(coord_error, ret, "Datafile %s doesn't specify star coordinates", filenames[0]);

    // Convert ra from HH:MM:SS to radians
    double a,b,c;
    sscanf(datafiles[0]->coord_ra, "%lf:%lf:%lf", &a, &b, &c);
    double ra = (a + b/60 + c/3600)*M_PI/12;

    // Convert dec from DD:'':"" to radians
    sscanf(datafiles[0]->coord_dec, "%lf:%lf:%lf", &a, &b, &c);
    double dec = (a + b/60 + c/3600)*M_PI/180;
    double2 coords = {ra, dec};
    double epoch = datafiles[0]->coord_epoch;

    // Generate a struct tm for our reference time
    struct tm t = parse_date_time_tm(reference_date, reference_time);

    // Convert from UT to TT
    t.tm_sec += utcttoffset(&t);

    // Calculate BJD
    double2 reference_coords = precess(coords, epoch, tmtoyear(&t));
    double reference_bjd = jdtobjd(tmtojd(&t), reference_coords);

    printf("Reference BJD: %f\n", reference_bjd);

    FILE *out = fopen(ts_filename, "w+");
    if (!out)
        error_jump(output_error, ret, "Error opening file %s", ts_filename);

    // Print file header
    fprintf(out, "# tsreduce create-ts output file\n");
    fprintf(out, "# Reference time: %s %s UTC; %f BJD\n", reference_date, reference_time, reference_bjd);
    fprintf(out, "# Files:\n");
    for (int i = 0; i < num_datafiles; i++)
        fprintf(out, "#   %s\n", filenames[i]);

    // Convert data to BJD relative to the reference time
    int num_saved = 0;
    for (int i = 0; i < num_datafiles; i++)
    {
        // Calculate time offset relative to the reference
        struct tm st;
        ts_gmtime(datafiles[i]->reference_time, &st);

        // Convert from UT to TT
        st.tm_sec += utcttoffset(&t);
        double start_jd = tmtojd(&st);

        // Calculate precessed RA and DEC at the start of each night
        // This is already far more accurate than we need
        double2 start_coords = precess(coords, epoch, tmtoyear(&st));
        printf("%s start BJD: %f\n", filenames[i], jdtobjd(start_jd, start_coords));

        double *raw_time, *raw, *mean_sky, *time, *ratio, *ratio_noise, *fwhm, *polyfit, *mma, *mma_noise;
        size_t num_raw, num_filtered;
        if (generate_photometry_dft_data(datafiles[i],
                                         &raw_time, &raw, &mean_sky, &num_raw,
                                         &time, &ratio, &polyfit, &mma, &fwhm, &num_filtered,
                                         &ratio_noise, &mma_noise,
                                         NULL, NULL, NULL, NULL,
                                         NULL, NULL, NULL))
        {
            error_jump(processing_error, ret, "Error generating photometry data for data %s", filenames[i]);
        }

        for (int j = 0; j < num_filtered; j++)
        {
            double bjd = jdtobjd(start_jd + time[j]/86400, start_coords);
            fprintf(out,"%f %f %f\n", bjd - reference_bjd, mma[j], mma_noise[j]);
            num_saved++;
        }

        free(raw_time);
        free(raw);
        free(mean_sky);
        free(time);
        free(ratio);
        free(ratio_noise);
        free(polyfit);
        free(mma);
        free(mma_noise);
        free(fwhm);
    }
    printf("Converted %d observations\n", num_saved);

processing_error:
    fclose(out);
output_error:
coord_error:
    for (int i = 0; i < num_datafiles; i++)
        datafile_free(datafiles[i]);
    free(datafiles);
datafile_error:
    return ret;
}
