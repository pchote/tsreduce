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
#include <dirent.h>
#include <regex.h>
#include <stdbool.h>
#include <stdint.h>
#include <fitsio.h>

#include "datafile.h"
#include "framedata.h"
#include "helpers.h"
#include "aperture.h"
#include "fit.h"
#include "hashmap.h"

extern int verbosity;

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
int create_flat(const char *pattern, size_t minmax, const char *masterdark, const char *outname)
{
    int ret = 0;

    // Find the filenames that match the specified pattern
    char **frame_paths;
    size_t num_frames = get_matching_files(pattern, &frame_paths);
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
    for(size_t i = 0; i < num_frames; i++)
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
        for (size_t j = 0; j < dark->rows*dark->cols; j++)
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

    for (size_t j = 0; j < dark->rows*dark->cols; j++)
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

        for(size_t k = 0; k < num_frames; k++)
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
                printf("%zu var: %f mean: %f dark: %f gain: %f\n", k, var, mean_flat[k], mean_dark, gain[k]);
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
        for (size_t k = 0; k < dark->rows*dark->cols; k++)
        {
            size_t x = k % dark->cols;
            size_t y = k / dark->cols;
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
    for (size_t j = 0; j < num_frames; j++)
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
int create_dark(const char *pattern, size_t minmax, const char *outname)
{
    int ret = 0;

    char **frame_paths;
    size_t num_frames = get_matching_files(pattern, &frame_paths);
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

    for (size_t i = 0; i < num_frames; i++)
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
    for (size_t j = 0; j < base->rows*base->cols; j++)
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

    char *datadir = getcwd(NULL, 0);
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
    size_t start_obs = data->obs_count;
    size_t num_frames = get_matching_files(data->frame_pattern, &frame_paths);
    for (size_t i = 0; i < num_frames; i++)
    {
        // Check if file has been processed
        struct observation *unused;
        if (hashmap_get(data->filename_map, frame_paths[i], (void**)(&unused)) != MAP_MISSING)
            continue;

        printf("Reducing %s\n", frame_paths[i]);

        framedata *frame = framedata_load(frame_paths[i]);
        if (!frame)
        {
            framedata_free(frame);
            error_jump(process_error, ret, "Error loading frame %s", frame_paths[i]);
        }

        double exptime;
        if (framedata_get_header_dbl(frame, "EXPTIME", &exptime))
        {
            framedata_free(frame);
            error_jump(process_error, ret, "EXPTIME undefined in %s", frame_paths[i]);
        }

        // Calculate time at the start of the exposure relative to ReferenceTime
        ts_time frame_time = framedata_start_time(frame);
        double starttime = ts_difftime(frame_time, data->reference_time);

        // Process frame
        subtract_bias(frame);
        if (dark)
            framedata_subtract(frame, dark);
        if (flat)
            framedata_divide(frame, flat);

        struct observation *obs = datafile_new_observation(data);
        if (!obs)
            error_jump(process_error, ret, "Allocation error");

        // Observation start time
        obs->time = starttime;

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
            if (data->obs_end)
            {
                double2 last = data->obs_end->pos[i];
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
                    double fwhm = estimate_fwhm(frame, xy, bg, t.s1);
                    if (fwhm < 1)
                        failed = true;

                    mean_fwhm += fwhm / data->num_targets;
                }
            }
            else
                failed = true;

            obs->star[i] = intensity;
            obs->sky[i] = sky;
            obs->pos[i] = xy;

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

        obs->fwhm = mean_fwhm;
        obs->ratio = comparisonIntensity > 0 ? targetIntensity / comparisonIntensity : 0;
        obs->ratio_noise = failed ? sqrt(-1) : (targetNoise/targetIntensity + comparisonNoise/comparisonIntensity)*obs->ratio;
        obs->filename = strdup(frame_paths[i]);

        datafile_append_observation(data, obs);
        framedata_free(frame);
    }
    
    if (chdir(datadir))
        error_jump(process_error, ret, "Invalid data path: %s", datadir);

    datafile_save(data, dataPath);

    printf("Reduced %zu observations\n", data->obs_count - start_obs);

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
                snprintf(dark_pattern, 1039, "^%s-[0-9]+.fits.gz", ret);
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
                snprintf(flat_pattern, 32, "^%s-[0-9]+.fits.gz", ret);
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
        snprintf(namebuf, 1039, "^%s-[0-9]+.fits.gz", ret);
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
    data->reference_time = framedata_start_time(frame);

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

        if (framedata_get_header_dbl(flat, "IM-SCALE", &data->ccd_platescale))
        {
            char *ret = prompt_user_input("Enter CCD platescale (arcsec/px):", "0.66");
            data->ccd_platescale = atof(ret);
            free(ret);
        }

        framedata_divide(frame, flat);
        framedata_free(flat);
    }

    while (true)
    {
        uint8_t aperture_type = 0;
        double aperture_size = 0;
        while (true)
        {
            printf("Aperture types are  0: 5sigma,  1: 3FWHM,  2: Manual\n");
            char *ret = prompt_user_input("Select aperture type:", "0");
            aperture_type = atoi(ret);
            free(ret);
            if (aperture_type <= 2)
                break;

            printf("Number must be between 0 and 2\n");
        }

        if (aperture_type == 2)
        {
            while (true)
            {
                char *ret = prompt_user_input("Enter aperture radius (px):", "8");
                aperture_size = atof(ret);
                free(ret);
                if (aperture_size > 0)
                    break;

                printf("Number must be greater than 0\n");
            }
        }

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
        // Temporarily hardcode target limit to 10
        // TODO: Fix this
        data->targets = malloc(10*sizeof(target));
        double sky[10];

        for (; (cur = strstr(cur, "annulus")) != NULL; cur++)
        {
            // TODO: Fix this
            if (data->num_targets == 10)
            {
                printf("Limit of %d targets reached. Remaining targets have been ignored", 10);
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

            switch (aperture_type)
            {
                case 0: // 5*sky sigma
                {
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
                    break;
                }
                case 1: // 3*FWHM
                {
                    double fwhm = estimate_fwhm(frame, xy, sky_intensity, t.s1);
                    t.r = 3*fwhm/2;
                    largest_aperture = fmax(largest_aperture, t.r);

                    break;
                }
                case 2:
                default:
                    t.r = aperture_size;
                    largest_aperture = aperture_size;
                    break;
            }

            // Set target parameters
            data->targets[data->num_targets++] = t;
        }
        free(ds9buf);

        // Set aperture radii to the same size, equal to the largest calculated above
        for (int i = 0; i < data->num_targets; i++)
            data->targets[i].r = largest_aperture;

        printf("Aperture radius: %fpx\n", largest_aperture);

        // Display results in ds9 - errors are non-fatal
        ts_exec_write("xpaset tsreduce regions delete all", NULL, 0);

        for (size_t i = 0; i < data->num_targets; i++)
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
                snprintf(msg, 16, "Comparison %zu", i);

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

    datafile_save(data, outname);
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

        double fwhm = estimate_fwhm(frame, xy, sky_intensity, t.s1);
        if (fwhm < 1)
        {
            printf("Invalid fwhm. Removing target\n");
            continue;
        }

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{annulus %f %f %f %f}'", ds9_title, t.x + 1, t.y + 1, t.s1, t.s2);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{annulus %f %f %f %f #color=red select=0 background}'", ds9_title, t.x + 1, t.y + 1, t.s1, t.s2);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{circle %f %f %f #color=red select=0}'", ds9_title, t.x + 1, t.y + 1, fwhm/2);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{text %f %f #color=green select=0 text=\"FWHM: %.2f arcsec\"}'",
                 ds9_title, t.x + 1, t.y + 1 - t.s2 - 10/zoom, fwhm*binning*plate_scale);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{text %f %f #color=green select=0 text=\"Peak: %.0f ADU/px\"}'",
                 ds9_title, t.x + 1, t.y + 1 - t.s2 - 25/zoom, frame->data[frame->cols*((int)xy.y) + (int)xy.x] - sky_intensity);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{text %f %f #color=green select=0 text=\"BG: %.0f ADU/px\"}'",
                 ds9_title, t.x + 1, t.y + 1 - t.s2 - 40/zoom, sky_intensity);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{box %f %f %f %f #color=black select=0 width=%d}'",
                 ds9_title, t.x + 1, t.y + 1 - t.s2 - 25/zoom, 115/zoom, 25/zoom, 25);
        ts_exec_write(ds9_command_buf, NULL, 0);
    }

    // Display frame time
    double frame_exp;
    char *frame_end = framedata_get_header_string(frame, "UTC-END");
    char *frame_date = framedata_get_header_string(frame, "UTC-DATE");
    char *frame_object = framedata_get_header_string(frame, "OBJECT");
    framedata_get_header_dbl(frame, "EXPTIME", &frame_exp);

    if (frame_object)
    {
        snprintf(ds9_command_buf, 1024,
                 "xpaset -p %s regions command '{text %f %f #color=green select=0 font=\"helvetica 12 bold roman\" text=\"%s @ %gs\"}'",
                 ds9_title, frame->cols/2.0, frame->rows + 30/zoom, frame_object, frame_exp);
        ts_exec_write(ds9_command_buf, NULL, 0);
    }

    if (frame_date && frame_end)
    {
        snprintf(ds9_command_buf, 1024,
                 "xpaset -p %s regions command '{text %f %f #color=green select=0 font=\"helvetica 12 bold roman\" text=\"Ending: %s %s\"}'",
                 ds9_title, frame->cols/2.0, frame->rows + 10/zoom, frame_date, frame_end);
        ts_exec_write(ds9_command_buf, NULL, 0);
    }

    free(frame_object);
    free(frame_date);
    free(frame_end);

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

    printf("%f\n", ts_time_to_bjd(parse_date_time(date, time), ra, dec, epoch));
    return 0;
}

/*
 * Create a timeseries file from a list of datafiles.
 * Times are converted to BJD after the specified reference
 */
int create_ts(char *reference_date, char *reference_time, char **filenames, size_t num_datafiles, char *ts_filename)
{
    int ret = 0;
    if (num_datafiles < 1)
        return error("No datafiles specified");

    datafile **datafiles = (datafile **)malloc(num_datafiles*sizeof(datafile *));
    if (!datafiles)
        error_jump(datafile_error, ret, "Error allocating datafiles");

    for (size_t i = 0; i < num_datafiles; i++)
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
    printf("Loaded %zu datafiles\n", num_datafiles);

    if (!datafiles[0]->coord_ra || !datafiles[0]->coord_dec || datafiles[0]->coord_epoch == 0)
        error_jump(coord_error, ret, "Datafile %s doesn't specify star coordinates", filenames[0]);

    // Convert ra from HH:MM:SS to radians
    double a,b,c;
    sscanf(datafiles[0]->coord_ra, "%lf:%lf:%lf", &a, &b, &c);
    double ra = (a + b/60 + c/3600)*M_PI/12;

    // Convert dec from DD:'':"" to radians
    sscanf(datafiles[0]->coord_dec, "%lf:%lf:%lf", &a, &b, &c);
    double dec = (a + b/60 + c/3600)*M_PI/180;
    double epoch = datafiles[0]->coord_epoch;

    double reference_bjd = ts_time_to_bjd(parse_date_time(reference_date, reference_time), ra, dec, epoch);
    printf("Reference BJD: %f\n", reference_bjd);

    FILE *out = fopen(ts_filename, "w+");
    if (!out)
        error_jump(output_error, ret, "Error opening file %s", ts_filename);

    // Print file header
    fprintf(out, "# tsreduce create-ts output file\n");
    fprintf(out, "# Reference time: %s %s UTC; %f BJD\n", reference_date, reference_time, reference_bjd);
    fprintf(out, "# Files:\n");
    for (size_t i = 0; i < num_datafiles; i++)
        fprintf(out, "#   %s\n", filenames[i]);

    // Convert data to BJD relative to the reference time
    size_t num_saved = 0;
    for (size_t i = 0; i < num_datafiles; i++)
    {
        double start_bjd = ts_time_to_bjd(datafiles[i]->reference_time, ra, dec, epoch);

        // Calculate precessed RA and DEC at the start of each night
        // This is already far more accurate than we need
        printf("%s start BJD: %f\n", filenames[i], start_bjd);

        struct photometry_data *pd = datafile_generate_photometry(datafiles[i]);
        if (!pd)
            error_jump(processing_error, ret, "Error generating photometry data for data %s", filenames[i]);

        for (size_t j = 0; j < pd->filtered_count; j++)
        {
            ts_time obstime = datafiles[i]->reference_time;
            obstime.time += (time_t)(pd->time[j]);
            obstime.ms += round(1000*fmod(pd->time[j], 1));
            fprintf(out,"%f %f %f\n", ts_time_to_bjd(obstime, ra, dec, epoch) - reference_bjd, pd->mma[j], pd->mma_noise[j]);
            num_saved++;
        }

        datafile_free_photometry(pd);
    }
    printf("Converted %zu observations\n", num_saved);

processing_error:
    fclose(out);
output_error:
coord_error:
    for (size_t i = 0; i < num_datafiles; i++)
        datafile_free(datafiles[i]);
    free(datafiles);
datafile_error:
    return ret;
}
