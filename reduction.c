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

int generate_photometry_dft_data(datafile *data,
    // Raw data, unfiltered by blocked ranges
    double **raw_time, double **raw, size_t *num_raw,
    // Filtered data: polynomial fit, ratio, MMI
    double **time, double **ratio, double **polyfit, double **mma, size_t *num_filtered,
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
    if (raw == NULL)
        error_jump(raw_alloc_error, ret, "raw malloc failed");

    *time = (double *)malloc(data->num_obs*sizeof(double));
    if (*time == NULL)
        error_jump(time_alloc_error, ret, "time malloc failed");

    *ratio = (double *)malloc(data->num_obs*sizeof(double));
    if (*ratio == NULL)
        error_jump(ratio_alloc_error, ret, "ratio malloc failed");

    *ratio_noise = (double *)malloc(data->num_obs*sizeof(double));
    if (*ratio == NULL)
        error_jump(ratio_noise_alloc_error, ret, "ratio malloc failed");

    *polyfit = (double *)malloc(data->num_obs*sizeof(double));
    if (*polyfit == NULL)
        error_jump(polyfit_alloc_error, ret, "polyfit malloc failed");

    // Calculate polynomial fit to the ratio
    double *coeffs = (double *)malloc((data->plot_fit_degree+1)*sizeof(double));
    if (coeffs == NULL)
        error_jump(coeffs_alloc_error, ret, "coeffs malloc failed");

    double ratio_mean = 0;
    *num_filtered = 0;
    *num_raw = data->num_obs;
    for (int i = 0; i < data->num_obs; i++)
    {
        (*raw_time)[i] = data->obs[i].time;

        for (int j = 0; j < data->num_targets; j++)
            (*raw)[j*data->num_obs + i] = data->obs[i].star[j];

        // Filter bad observations
        bool skip = false;
        for (int j = 0; j < data->num_blocked_ranges; j++)
            if ((*raw_time)[i] >= data->blocked_ranges[j].x && (*raw_time)[i] <= data->blocked_ranges[j].y)
            {
                skip = true;
                break;
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
                error("Ignoring bad observation at %f", data->obs[i].time);
                continue;
            }

            (*time)[*num_filtered] = data->obs[i].time;
            (*ratio)[*num_filtered] = comparison > 0 ? target/comparison : 0;

            // Read noise from data file if available
            (*ratio_noise)[*num_filtered] = (data->version >= 5) ? data->obs[i].ratio_noise : 1;

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

    if (fit_polynomial(*time, *ratio, *ratio_noise, *num_filtered, coeffs, data->plot_fit_degree))
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

        if (data->version >= 5)
        {
            double numer_error = fabs((*ratio_noise)[i]/((*ratio)[i] - (*polyfit)[i]));
            double denom_error = fabs((*ratio_noise)[i]/(*ratio)[i]);
            (*mma_noise)[i] = (numer_error + denom_error)*fabs((*mma)[i]);
        }
        else
            (*mma_noise)[i] = 1;

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
    free(*ratio_noise);
ratio_noise_alloc_error:
    free(*ratio);
    *ratio = NULL;
ratio_alloc_error:
    free(*time);
    *time = NULL;
time_alloc_error:
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
    char datebuf[128], timebuf[128], datetimebuf[257];
    if (framedata_has_header_string(frame, "UTC-BEG"))
    {
        framedata_get_header_string(frame, "UTC-DATE", datebuf);
        framedata_get_header_string(frame, "UTC-BEG", timebuf);
        sprintf(datetimebuf, "%s %s", datebuf, timebuf);
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
double prepare_flat(framedata *flat, framedata *dark)
{
    int flatexp = framedata_get_header_int(flat, "EXPTIME");
    int darkexp = framedata_get_header_int(dark, "EXPTIME");

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
    return mean_new;
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
        mean_flat[i] = prepare_flat(frames[i], dark);

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
    char *outbuf;
    asprintf(&outbuf, "!%s", outname);
    fits_create_file(&out, outbuf, &status);
    free(outbuf);

    // Create the primary array image (16-bit short integer pixels
    fits_create_img(out, DOUBLE_IMG, 2, (long []){dark->cols, dark->rows}, &status);

    // Set header keys for readout noise and gain
    if (dark->regions.has_overscan)
    {
        fits_update_key(out, TDOUBLE, "CCD-READ", &readnoise, "Estimated read noise (ADU)", &status);
        fits_update_key(out, TDOUBLE, "CCD-GAIN", &median_gain, "Estimated gain (electrons/ADU)", &status);

        printf("Readnoise: %f\n", readnoise);
        printf("Gain mean: %f median: %f\n", mean_gain, median_gain);
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

    int exptime = framedata_get_header_int(base, "EXPTIME");
    double *median_dark = (double *)malloc(base->rows*base->cols*sizeof(double));
    if (median_dark == NULL)
        error_jump(dark_failed, ret, "median_dark alloc failed");

    // Data cube for processing the flat data
    //        data[0] = frame[0][0,0], data[1] = frame[1][0,0] ... data[num_frames] = frame[0][1,0] etc
    double *data_cube = (double *)malloc(num_frames*base->cols*base->rows*sizeof(double));
    if (data_cube == NULL)
        error_jump(datacube_failed, ret, "data_cube alloc failed");

    for( int i = 0; i < num_frames; i++)
    {
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
    char outbuf[2048];
    sprintf(outbuf, "!%s", outname);
    
    fits_create_file(&out, outbuf, &status);
    
    // Create the primary array image (16-bit short integer pixels
    fits_create_img(out, DOUBLE_IMG, 2, (long []){base->cols, base->rows}, &status);
    fits_update_key(out, TINT, "EXPTIME", &exptime, "Actual integration time (sec)", &status);
    
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
    char outbuf[2048];
    sprintf(outbuf, "!%s", outPath);
    fits_create_file(&out, outbuf, &status);

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

    framedata *flat = framedata_load(data->flat_template);
    if (!flat)
        error_jump(frame_error, ret, "Error loading frame %s", data->flat_template);

    double readnoise = data->ccd_readnoise > 0 ? data->ccd_readnoise : framedata_get_header_dbl(flat, "CCD-READ");
    double gain = data->ccd_gain > 0 ? data->ccd_gain : framedata_get_header_dbl(flat, "CCD-GAIN");

    framedata *dark = framedata_load(data->dark_template);
    if (!dark)
        error_jump(frame_error, ret, "Error loading frame %s", data->flat_template);

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
        int exptime = framedata_get_header_int(frame, "EXPTIME");

        // Calculate time at the start of the exposure relative to ReferenceTime
        time_t frame_time = get_frame_time(frame);
        double starttime = difftime(frame_time, data->reference_time);

        // Process frame
        subtract_bias(frame);
        framedata_subtract(frame, dark);
        framedata_divide(frame, flat);

        // Observation start time
        fprintf(data->file, "%.1f ", starttime);
        data->obs[data->num_obs].time = starttime;

        // Process frame
        double comparisonIntensity = 0;
        double targetIntensity = 0;
        double comparisonNoise = 0;
        double targetNoise = 0;
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
                if (calculate_background(t, frame, &sky, NULL))
                    sky = 0;

                // Integrate sky over the aperture and normalize per unit time
                sky *= M_PI*t.r*t.r / exptime;

                integrate_aperture_and_noise(xy, t.r, frame, dark, readnoise, gain, &intensity, &noise);
                intensity = intensity/exptime - sky;
                noise /= exptime;
            }

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
        double ratioNoise = ratio ? (targetNoise/targetIntensity + comparisonNoise/comparisonIntensity)*ratio : 1;
        data->obs[data->num_obs].ratio = ratio;

        fprintf(data->file, "%.3e ", ratio);
        if (data->version >= 5)
            fprintf(data->file, "%.3e ", ratioNoise);

        // Filename
        fprintf(data->file, "%s\n", frame_paths[i]);
        strncpy(data->obs[data->num_obs].filename, frame_paths[i], 64);
        data->num_obs++;

        framedata_free(frame);
    }
    
    printf("Reduced %d observations\n", data->num_obs - start_obs);

process_error:
    free_2d_array(frame_paths, num_frames);

frame_error:
    framedata_free(flat);
    framedata_free(dark);

data_error:
    datafile_free(data);
    return ret;
}


// Create a reduction file at filePath, with frames from the dir framePath
// matching the regex framePattern, with dark and flat frames given by
// darkTemplate and flatTemplate
int create_reduction_file(char *framePath, char *framePattern, char *darkTemplate, char *flatTemplate, char *outname)
{
    int ret = 0;

#if (defined _WIN32 || defined _WIN64)
    // Windows doesn't support the "x" mode
    const char *mode = "w";
#else
    // Non-rigorous test that we won't overwrite an existing file
    const char *mode = "wx";
#endif
    FILE *fileTest = fopen(outname, mode);
    if (fileTest == NULL)
        return error("Unable to create data file: %s. Does it already exist?", outname);
    fclose(fileTest);

    if (init_ds9())
        return error("Unable to launch ds9");

    char pathBuf[PATH_MAX];
    realpath(framePath, pathBuf);
    datafile *data = datafile_alloc();

    // Store the current directory so we can return before saving the data file
    char *datadir = getcwd(NULL, 0);

    if (chdir(pathBuf))
        error_jump(setup_error, ret, "Invalid frame path: %s", pathBuf);
    data->frame_dir = strdup(pathBuf);
    
    char *filename = get_first_matching_file(framePattern);
    if (!filename)
        error_jump(setup_error, ret, "No matching files found");
    data->frame_pattern = strdup(framePattern);
    
    // Open the file to find the reference time
    framedata *frame = framedata_load(filename);
    if (!frame)
        error_jump(frameload_error, ret, "Error loading frame %s", filename);

    subtract_bias(frame);
    data->reference_time = get_frame_time(frame);

    framedata *dark = framedata_load(darkTemplate);
    if (!dark)
        error_jump(frameload_error, ret, "Error loading frame %s", darkTemplate);

    framedata_subtract(frame, dark);
    framedata_free(dark);
    data->dark_template = strdup(darkTemplate);

    framedata *flat = framedata_load(flatTemplate);
    if (!flat)
        error_jump(frameload_error, ret, "Error loading frame %s", flatTemplate);

    framedata_divide(frame, flat);
    framedata_free(flat);
    data->flat_template = strdup(flatTemplate);

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

    printf("Circle the target stars and surrounding sky in ds9 then press enter to continue...\n");
    getchar();

    char *ds9buf;
    if (ts_exec_read("xpaget tsreduce regions", &ds9buf))
        error_jump(frameload_error, ret, "ds9 request regions failed");

    // Parse the region definitions
    char *cur = ds9buf;
    double largest_aperture = 0;
    while ((cur = strstr(cur, "circle")) != NULL)
    {
        if (data->num_targets == MAX_TARGETS)
        {
            printf("Limit of %d targets reached. Remaining targets have been ignored", MAX_TARGETS);
            break;
        }
        
        // Read aperture coords
        target t;
        sscanf(cur, "circle(%lf,%lf,%lf)", &t.x, &t.y, &t.s2);
        
        // ds9 denotes the bottom-left pixel (1,1), tsreduce uses (0,0)
        t.x -= 1;
        t.y -= 1;
        
        //
        // Calculate the optimum aperture positioning
        //
        
        // Set outer sky radius to selected radius, inner to outer - 5
        t.s1 = t.s2 - 5;
        
        // Estimate initial aperture size as inner sky radius
        t.r = t.s1;
        
        printf("Initial aperture xy: (%f,%f) r: %f s:(%f,%f)\n", t.x, t.y, t.r, t.s1, t.s2);
        
        double2 xy;
        if (center_aperture(t, frame, &xy))
            error_jump(aperture_converge_error, ret, "Aperture did not converge");
        t.x = xy.x;
        t.y = xy.y;

        double sky_intensity, sky_std_dev;
        if (calculate_background(t, frame, &sky_intensity, &sky_std_dev))
            error_jump(aperture_converge_error, ret, "Background calculation failed");

        // Estimate the radius where the star flux falls to 5 times the std. dev. of the background
        double lastIntensity = 0;
        double lastProfile = frame->data[frame->cols*((int)xy.y) + (int)xy.x];
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
        
        // Increment by one char so strstr will find the next instance
        cur++;
    }

    // Set aperture radii to the same size, equal to the largest calculated above
    for (int i = 0; i < data->num_targets; i++)
        data->targets[i].r = largest_aperture;
    
    printf("Founds %d targets\n", data->num_targets);
    
    // Display results in ds9 - errors are non-fatal
    ts_exec_write("xpaset tsreduce regions delete all", NULL, 0);

    for (int i = 0; i < data->num_targets; i++)
    {
        double x = data->targets[i].x + 1;
        double y = data->targets[i].y + 1;
        char command[128];
        snprintf(command, 128, "xpaset tsreduce regions command '{circle %f %f %f #color=red}'", x, y, data->targets[i].r);
        ts_exec_write(command, NULL, 0);
        snprintf(command, 128, "xpaset tsreduce regions command '{circle %f %f %f #background}'", x, y, data->targets[i].s1);
        ts_exec_write(command, NULL, 0);
        snprintf(command, 128, "xpaset tsreduce regions command '{circle %f %f %f #background}'", x, y, data->targets[i].s2);
        ts_exec_write(command, NULL, 0);
    }

    // Save to disk
    if (chdir(datadir))
        error_jump(aperture_converge_error, ret, "Invalid data path: %s", datadir);

    datafile_save_header(data, outname);
aperture_converge_error:
    free(ds9buf);
frameload_error:
    framedata_free(frame);
    free(filename);
setup_error:
    datafile_free(data);
    return ret;
}

/*
 * Convert a .dat file to .mma for analysis with ts3add etc
 */
int create_mma(char *dataPath)
{
    // Read file header
    datafile *data = datafile_load(dataPath);
    if (data == NULL)
        return error("Error opening data file");

    double *raw_time, *raw, *time, *ratio, *ratio_noise, *polyfit, *mma, *mma_noise;
    size_t num_raw, num_filtered;
    if (generate_photometry_dft_data(data,
                                     &raw_time, &raw, &num_raw,
                                     &time, &ratio, &polyfit, &mma, &num_filtered,
                                     &ratio_noise, &mma_noise,
                                     NULL, NULL, NULL, NULL,
                                     NULL, NULL, NULL))
    {
        datafile_free(data);
        return error("Error generating MMI data");
    }

    printf("# points = %d; dt = %f\n", (unsigned int)num_filtered, (time[num_filtered-1] - time[0])/3600.0);
    printf("# tgt = %s\n", data->frame_pattern);

    char buf[20];
    serialize_time_t(data->reference_time, buf);
    printf("# UT start = %s\n", buf);
    printf("# Time MMI error\n");

    for (int i = 0; i < num_filtered; i++)
        printf("%f %f %f\n", time[i]/3600.0, mma[i], mma_noise[i]);

    free(raw_time);
    free(raw);
    free(time);
    free(ratio);
    free(ratio_noise);
    free(polyfit);
    free(mma);
    free(mma_noise);
    datafile_free(data);
    return 0;
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

        double *raw_time, *raw, *time, *ratio, *ratio_noise, *polyfit, *mma, *mma_noise;
        size_t num_raw, num_filtered;
        if (generate_photometry_dft_data(datafiles[i],
                                         &raw_time, &raw, &num_raw,
                                         &time, &ratio, &polyfit, &mma, &num_filtered,
                                         &ratio_noise, &mma_noise,
                                         NULL, NULL, NULL, NULL,
                                         NULL, NULL, NULL))
        {
            error_jump(processing_error, ret, "Error generating MMI data for data %s", filenames[i]);
        }

        for (int j = 0; j < num_filtered; j++)
        {
            double bjd = jdtobjd(start_jd + time[j]/86400, start_coords);
            fprintf(out,"%f %f %f\n", bjd - reference_bjd, mma[j], mma_noise[j]);
            num_saved++;
        }

        free(raw_time);
        free(raw);
        free(time);
        free(ratio);
        free(ratio_noise);
        free(polyfit);
        free(mma);
        free(mma_noise);
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
