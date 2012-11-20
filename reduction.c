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
#include <dirent.h>
#include <regex.h>
#include <stdbool.h>
#include <stdint.h>

#include "datafile.h"
#include "framedata.h"
#include "helpers.h"
#include "aperture.h"
#include "fit.h"
#include "hashmap.h"

extern int verbosity;

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

    // Ensure there are enough frames to discard the requested number of pixels
    if (num_frames <= 2*minmax)
        error_jump(insufficient_frames_error, ret,
            "Insufficient frames. %d found, %d will be discarded", num_frames, 2*minmax);

    // Load the master dark frame
    // Frame geometry for all subsequent frames is assumed to match the master-dark
    framedata *dark = framedata_load(masterdark);
    if (!dark)
        error_jump(setup_error, ret, "Error loading frame %s", masterdark);

    framedata *base = framedata_load(frame_paths[0]);
    if (!base)
        error_jump(setup_error, ret, "Error loading frame %s", frame_paths[0]);

    if (base->rows != dark->rows || base->cols != dark->cols)
        error_jump(setup_error, ret, "Dark and flat frame sizes don't match");


    uint16_t image_region[4] = {0, base->cols, 0, base->rows};
    uint16_t bias_region[4] = {0, 0, 0, 0};

    {
        char *str;
        if (framedata_get_metadata(base, "IMAG-RGN", FRAME_METADATA_STRING, &str) == FRAME_METADATA_OK)
            sscanf(str, "[%hu, %hu, %hu, %hu]", &image_region[0], &image_region[1],
                                                &image_region[2], &image_region[3]);

        if (framedata_get_metadata(base, "BIAS-RGN", FRAME_METADATA_STRING, &str) == FRAME_METADATA_OK)
            sscanf(str, "[%hu, %hu, %hu, %hu]", &bias_region[0], &bias_region[1],
                                                &bias_region[2], &bias_region[3]);
    }
    size_t image_region_px = (image_region[1] - image_region[0])*(image_region[3] - image_region[2]);
    size_t bias_region_px = (bias_region[1] - bias_region[0])*(bias_region[3] - bias_region[2]);

    // Data cube for processing the flat data
    //        data[0] = frame[0][0,0], data[1] = frame[1][0,0] ... data[num_frames] = frame[0][1,0] etc
    double bias_variance = 0;

    double *data_cube = calloc(num_frames*base->cols*base->rows, sizeof(double));
    double *frame_mean = calloc(num_frames, sizeof(double));
    if (!data_cube || !frame_mean)
        error_jump(processing_error, ret, "Allocation failed");

    for (size_t k = 0; k < num_frames; k++)
    {
        if (verbosity >= 1)
            printf("loading `%s`\n", frame_paths[k]);

        framedata *frame = framedata_load(frame_paths[k]);
        if (!frame)
            error_jump(processing_error, ret, "Error loading frame %s", frame_paths[k]);

        if (frame->rows != base->rows || frame->cols != base->cols)
        {
            framedata_free(frame);
            error_jump(processing_error, ret,
                "Frame %s dimensions mismatch. Expected (%d,%d), was (%d, %d)",
                frame_paths[k], base->rows, base->cols, frame->rows, frame->cols);
        }

        // Subtract dark, normalized to the flat exposure time
        framedata_subtract_bias(frame);
        if (framedata_subtract_normalized(frame, dark))
        {
            framedata_free(frame);
            error_jump(processing_error, ret, "Dark subtraction failed for %s", frame_paths[k]);
        }

        // Store data in cube for processing
        for (size_t j = 0; j < base->rows*base->cols; j++)
        {
            data_cube[num_frames*j + k] = frame->data[j];

            // The frame has been bias subtracted, so the remaining signal
            // in the bias region is caused by readout noise
            if (region_contains(bias_region, j % base->cols, j / base->cols))
                bias_variance += frame->data[j]*frame->data[j];
        }

        // Calculate normalization factor to make the image region 1
        double mean = region_mean(image_region, frame->data, frame->cols);

        // Calculate standard deviation
        double std = 0;
        for (uint16_t j = image_region[2]; j < image_region[3]; j++)
            for (uint16_t i = image_region[0]; i < image_region[1]; i++)
            {
                double temp = frame->data[base->cols*j + i] - mean;
                std += temp*temp;
            }
        std = sqrt(std/image_region_px);

        // Recalculate the mean, excluding outliers at 3 sigma
        frame_mean[k] = 0;
        size_t count = 0;
        for (uint16_t j = image_region[2]; j < image_region[3]; j++)
            for (uint16_t i = image_region[0]; i < image_region[1]; i++)
                if (fabs(frame->data[base->cols*j + i] - mean) < 3*std)
                {
                    frame_mean[k] += frame->data[base->cols*j + i];
                    count++;
                }
        frame_mean[k] /= count;

        framedata_free(frame);
    }

    // Calculate median-mean image for the master-flat
    // Loop over the pixels, sorting the values from each image into increasing order...
    {
        double *cube_slice = calloc(num_frames*base->cols*base->rows, sizeof(double));
        if (!cube_slice)
            error_jump(processing_error, ret, "cube_slice alloc failed");

        for (size_t j = 0; j < base->rows*base->cols; j++)
        {
            for (size_t k = 0; k < num_frames; k++)
                cube_slice[k] = data_cube[num_frames*j + k] / frame_mean[k];

            qsort(cube_slice, num_frames, sizeof(double), compare_double);

            // ...and then average the non-rejected pixels into the output array
            base->data[j] = 0;
            for (int i = minmax; i < num_frames - minmax; i++)
                base->data[j] += cube_slice[i];
            base->data[j] /= (num_frames - 2*minmax);

        }
        free(cube_slice);
    }

    if (bias_region_px)
    {
        double readnoise = sqrt(bias_variance/(num_frames*bias_region_px));

        // Calculate gain from each image
        double *gain = calloc(num_frames, sizeof(double));
        if (!gain)
            error_jump(processing_error, ret, "gain alloc failed");

        // Calculate mean dark level for gain calculation
        double mean_dark = region_mean(image_region, dark->data, dark->cols);
        double mean_gain = 0;

        for (size_t k = 0; k < num_frames; k++)
        {
            // Calculate the variance by taking the difference between
            // the frame and the master-flat
            double var = 0;
            for (uint16_t j = image_region[2]; j < image_region[3]; j++)
                for (uint16_t i = image_region[0]; i < image_region[1]; i++)
                {
                    // Scale master-flat to the same level as the original frame
                    size_t l = base->cols*j + i;
                    double temp = data_cube[num_frames*l + k] - frame_mean[k]*base->data[l];
                    var += temp*temp;
                }

            var /= image_region_px;
            gain[k] = (frame_mean[k] + mean_dark) / (var - readnoise*readnoise);
            mean_gain += gain[k];

            if (verbosity >= 1)
                printf("%zu var: %f mean: %f dark: %f gain: %f\n", k, var, frame_mean[k], mean_dark, gain[k]);
        }
        mean_gain /= num_frames;

        // Find median value
        qsort(gain, num_frames, sizeof(double), compare_double);
        double median_gain = gain[num_frames/2];

        // Set header keys for readout noise and gain
        if (bias_region_px)
        {
            framedata_put_metadata(base, "CCD-READ", FRAME_METADATA_DOUBLE, &readnoise, "Estimated read noise (ADU)");
            framedata_put_metadata(base, "CCD-GAIN", FRAME_METADATA_DOUBLE, &median_gain, "Estimated gain (electrons/ADU)");

            printf("Readnoise: %f\n", readnoise);
            printf("Gain: %f\n", median_gain);
        }

        free(gain);
    }

    // Replace values outside the image region with 1, so overscan survives flatfielding
    if (image_region_px != base->rows*base->cols)
        for (uint16_t j = 0; j < base->rows; j++)
            for (uint16_t i = 0; i < base->cols; i++)
                if (!region_contains(image_region, i, j))
                    base->data[j*base->cols + i] = 1;

    framedata_save(base, outname);

processing_error:
    free(frame_mean);
    free(data_cube);
setup_error:
    if (base)
        framedata_free(base);
    if (dark)
        framedata_free(dark);
insufficient_frames_error:
    free_2d_array(frame_paths, num_frames);
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

    if (num_frames < 2*minmax)
        error_jump(insufficient_frames_error, ret,
            "Insufficient frames. %d found, %d will be discarded", num_frames, 2*minmax);

    framedata *base = framedata_load(frame_paths[0]);
    if (!base)
        error_jump(insufficient_frames_error, ret, "Error loading frame %s", frame_paths[0]);

    // Data cube for processing the flat data
    //        data[0] = frame[0][0,0], data[1] = frame[1][0,0] ... data[num_frames] = frame[0][1,0] etc
    double *data_cube = calloc(num_frames*base->cols*base->rows, sizeof(double));
    if (!data_cube)
        error_jump(setup_error, ret, "data_cube alloc failed");

    for (size_t k = 0; k < num_frames; k++)
    {
        if (verbosity >= 1)
            printf("loading `%s`\n", frame_paths[k]);

        framedata *frame = framedata_load(frame_paths[k]);
        if (!frame)
            error_jump(process_error, ret, "Error loading frame %s", frame_paths[k]);

        if (frame->rows != base->rows || frame->cols != base->cols)
        {
            framedata_free(frame);
            error_jump(process_error, ret,
                "Frame %s dimensions mismatch. Expected (%d,%d), was (%d, %d)",
                frame_paths[k], base->rows, base->cols, frame->rows, frame->cols);
        }

        framedata_subtract_bias(frame);
        for (size_t j = 0; j < base->rows*base->cols; j++)
            data_cube[num_frames*j + k] = frame->data[j];

        framedata_free(frame);
    }

    // Loop over the pixels, sorting the values from each image into increasing order
    for (size_t j = 0; j < base->rows*base->cols; j++)
    {
        qsort(data_cube + num_frames*j, num_frames, sizeof(double), compare_double);

        // then average the non-rejected pixels into the output array
        base->data[j] = 0;
        for (int i = minmax; i < num_frames - minmax; i++)
            base->data[j] += data_cube[num_frames*j + i];
        base->data[j] /= (num_frames - 2*minmax);
    }

    framedata_save(base, outname);
process_error:
    free(data_cube);
setup_error:
    framedata_free(base);
insufficient_frames_error:
    free_2d_array(frame_paths, num_frames);
    return ret;
}

int display_frame(char *data_path, char *frame_name)
{
    int ret = 0;

    datafile *data = datafile_load(data_path);
    if (data == NULL)
        return error("Error opening data file");

    struct observation *obs;
    if (hashmap_get(data->filename_map, frame_name, (void**)(&obs)) == MAP_MISSING)
        error_jump(setup_error, ret, "%s doesn't match any observation in the reduction file", frame_name);

    if (chdir(data->frame_dir))
        error_jump(setup_error, ret, "Invalid frame path: %s", data->frame_dir);

    framedata *frame = framedata_load(obs->filename);
    if (!frame)
        error_jump(setup_error, ret, "Error loading frame %s", obs->filename);

    framedata_subtract_bias(frame);

    if (data->dark_template)
    {
        framedata *dark = framedata_load(data->dark_template);
        if (!dark)
            error_jump(process_error, ret, "Error loading frame %s", data->dark_template);

        if (framedata_subtract(frame, dark))
        {
            framedata_free(dark);
            error_jump(process_error, ret, "Error dark-subtracting frame %s", obs->filename);
        }

        framedata_free(dark);
    }

    if (data->flat_template)
    {
        framedata *flat = framedata_load(data->flat_template);
        if (!flat)
            error_jump(process_error, ret, "Error loading frame %s", data->flat_template);

        if (framedata_divide(frame, flat))
        {
            framedata_free(flat);
            error_jump(process_error, ret, "Error flat-fielding frame %s", obs->filename);
        }

        framedata_free(flat);
    }

    if (init_ds9())
        error_jump(setup_error, ret, "Unable to launch ds9");

    char command[1024];
    snprintf(command, 1024, "xpaset tsreduce array [xdim=%d,ydim=%d,bitpix=-64]", frame->cols, frame->rows);
    if (ts_exec_write(command, frame->data, frame->rows*frame->cols*sizeof(double)))
        error_jump(process_error, ret, "ds9 command failed: %s", command);

    ts_exec_write("xpaset tsreduce scale mode zscale", NULL, 0);
    ts_exec_write("xpaset tsreduce regions delete all", NULL, 0);

    for (size_t i = 0; i < data->target_count; i++)
    {

        double2 xy = obs->pos[i];
        aperture *a = &data->targets[i].aperture;
        snprintf(command, 1024, "xpaset tsreduce regions command '{circle %f %f %f #color=red select=0}'", xy.x + 1, xy.y + 1, a->r);
        ts_exec_write(command, NULL, 0);

        snprintf(command, 1024, "xpaset tsreduce regions command '{annulus %f %f %f %f #select=0}'",
                 xy.x + 1, xy.y + 1, a->s1, a->s2);
        ts_exec_write(command, NULL, 0);
    }

process_error:
    framedata_free(frame);
setup_error:
    datafile_free(data);
    return ret;
}

int print_frame_metadata(char *frame_path)
{
    framedata *frame = framedata_load(frame_path);
    if (!frame)
        return error("Error loading frame %s", frame_path);

    framedata_print_metadata(frame);
    framedata_free(frame);
    return 0;
}

// Load the reduction file at dataPath and reduce any new data
int update_reduction(char *dataPath)
{
    int ret = 0;

    // Read file header
    datafile *data = datafile_load(dataPath);
    if (data == NULL)
        return error("Error opening data file");

    char *datadir = getcwd(NULL, 0);
    if (chdir(data->frame_dir))
        error_jump(data_error, ret, "Invalid frame path: %s", data->frame_dir);

    double readnoise = 0, gain = 1;
    framedata *flat = NULL;

    if (data->flat_template)
    {
        flat = framedata_load(data->flat_template);
        if (!flat)
            error_jump(flat_error, ret, "Error loading frame %s", data->flat_template);

        if (framedata_get_metadata(flat, "CCD-READ", FRAME_METADATA_DOUBLE, &readnoise))
            readnoise = data->ccd_readnoise;

        if (readnoise <= 0)
            error_jump(flat_error, ret, "CCD Read noise unknown. Define CCDReadNoise in %s.", dataPath);

        if (framedata_get_metadata(flat, "CCD-GAIN", FRAME_METADATA_DOUBLE, &gain))
            gain = data->ccd_gain;

        if (gain <= 0)
            error_jump(flat_error, ret, "CCD Gain unknown. Define CCDGain in %s.", dataPath);
    }

    framedata *dark = NULL;
    if (data->dark_template)
    {
        dark = framedata_load(data->dark_template);
        if (!dark)
            error_jump(dark_error, ret, "Error loading frame %s", data->dark_template);
    }

    framedata *reference = framedata_load(data->reference_frame);
    if (!reference)
        error_jump(reference_error, ret, "Error loading reference frame %s", data->reference_frame);

    framedata_subtract_bias(reference);
    if (dark && framedata_subtract_normalized(reference, dark))
        error_jump(reference_error, ret, "Error dark-subtracting reference frame %s", data->reference_frame);

    if (flat && framedata_divide(reference, flat))
        error_jump(reference_error, ret, "Error flat-fielding reference frame %s", data->reference_frame);

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
        if (framedata_get_metadata(frame, "EXPTIME", FRAME_METADATA_DOUBLE, &exptime))
        {
            framedata_free(frame);
            error_jump(process_error, ret, "EXPTIME undefined in %s", frame_paths[i]);
        }

        // Calculate time at the start of the exposure relative to ReferenceTime
        ts_time frame_time;
        if (framedata_start_time(frame, &frame_time))
            error_jump(process_error, ret, "No known time headers found");

        double midtime = ts_difftime(frame_time, data->reference_time) + exptime / 2;

        // Process frame
        framedata_subtract_bias(frame);
        if (dark && framedata_subtract_normalized(frame, dark))
            error_jump(process_error, ret, "Error dark-subtracting frame %s", frame_paths[i]);

        if (flat && framedata_divide(frame, flat))
            error_jump(process_error, ret, "Error flat-fielding frame %s", frame_paths[i]);

        struct observation *obs = datafile_new_observation(data);
        if (!obs)
            error_jump(process_error, ret, "Allocation error");

        // Observation mid time
        obs->time = midtime;

        // Process frame
        double nan = sqrt(-1);

        // Estimate translation from reference frame
        int32_t xt, yt;
        if (framedata_estimate_translation(frame, reference, &xt, &yt))
            error_jump(process_error, ret, "Error calculating frame translation");

        for (size_t i = 0; i < data->target_count; i++)
        {
            // Offset aperture position by frame offset
            aperture a = data->targets[i].aperture;
            a.x += xt;
            a.y += yt;

            obs->star[i] = 0;
            obs->noise[i] = 0;
            obs->sky[i] = 0;
            obs->pos[i] = (double2){0,0};
            obs->fwhm[i] = 0;

            bool failed = false;
            if (!center_aperture(a, frame, &obs->pos[i]))
            {
                double bg = 0;
                if (calculate_background(a, frame, &bg, NULL))
                    bg = 0;

                // Integrate sky over the aperture and normalize per unit time
                obs->sky[i] = bg*M_PI*a.r*a.r / exptime;

                integrate_aperture_and_noise(obs->pos[i], a.r, frame, dark, readnoise, gain, &obs->star[i], &obs->noise[i]);

                obs->star[i] = obs->star[i]/exptime - obs->sky[i];
                obs->noise[i] /= exptime;

                obs->fwhm[i] = estimate_fwhm(frame, obs->pos[i], bg, a.s1);
                if (obs->fwhm[i] < 0)
                    failed = true;
            }
            else
                failed = true;

            if (failed)
            {
                obs->star[i] = nan;
                obs->noise[i] = nan;
                obs->sky[i] = nan;
                obs->pos[i] = (double2){nan,nan};
                obs->fwhm[i] = nan;
            }
        }

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
reference_error:
    framedata_free(reference);
dark_error:
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

    char *filename_fmt = "^%s-[0-9]+.(fits.gz|fits|fit|FIT)";
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
                snprintf(dark_pattern, 1039, filename_fmt, ret);
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
                snprintf(flat_pattern, 1039, filename_fmt, ret);
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

    // Default target prefix to filename
    char *default_prefix = remove_file_suffix(strdup(last_path_component(outname)));
    while (true)
    {
        char namebuf[1039];
        char *ret = prompt_user_input("Enter target prefix", default_prefix);
        snprintf(namebuf, 1039, filename_fmt, ret);
        free(ret);
        data->reference_frame = get_first_matching_file(namebuf);
        if (data->reference_frame)
        {
            data->frame_pattern = strdup(namebuf);
            break;
        }
        printf("No files found matching pattern: %s/%s\n", data->frame_dir, namebuf);
    }
    free(default_prefix);

    framedata *frame = NULL;
    while (true)
    {
        char *ret = prompt_user_input("Enter reference frame", data->reference_frame);
        frame = framedata_load(ret);
        if (frame)
        {
            free(data->reference_frame);
            data->reference_frame = ret;
            break;
        }

        printf("File not found: %s/%s\n", data->frame_dir, ret);
        free(ret);
    }

    framedata_subtract_bias(frame);
    if (framedata_start_time(frame, &data->reference_time))
        error_jump(frameload_error, ret, "No known time headers found");

    if (data->dark_template)
    {
        framedata *dark = framedata_load(data->dark_template);
        if (!dark)
            error_jump(frameload_error, ret, "Error loading frame %s", data->dark_template);

        if (framedata_subtract(frame, dark))
        {
            framedata_free(dark);
            error_jump(frameload_error, ret, "Error dark-subtracting frame %s", data->reference_frame);
        }

        framedata_free(dark);
    }

    if (data->flat_template)
    {
        framedata *flat = framedata_load(data->flat_template);
        if (!flat)
            error_jump(frameload_error, ret, "Error loading frame %s", data->flat_template);

        if (!framedata_has_metadata(flat, "CCD-READ"))
        {
            while (true)
            {
                char *ret = prompt_user_input("Enter CCD Readnoise (ADU):", "3.32");
                data->ccd_readnoise = strtod(ret, NULL);
                free(ret);
                if (data->ccd_readnoise > 0)
                    break;

                printf("Number must be greater than 0\n");
            }
        }

        if (!framedata_has_metadata(flat, "CCD-GAIN"))
        {
            while (true)
            {
                char *ret = prompt_user_input("Enter CCD Gain (ADU):", "2.00");
                data->ccd_gain = strtod(ret, NULL);
                free(ret);

                if (data->ccd_gain > 0)
                    break;

                printf("Number must be greater than 0\n");
            }
        }

        if (framedata_get_metadata(flat, "IM-SCALE", FRAME_METADATA_DOUBLE, &data->ccd_platescale))
        {
            char *ret = prompt_user_input("Enter CCD platescale (arcsec/px):", "0.66");
            data->ccd_platescale = strtod(ret, NULL);
            free(ret);
        }

        if (framedata_divide(frame, flat))
        {
            framedata_free(flat);
            error_jump(frameload_error, ret, "Error flat-fielding frame %s", data->reference_frame);
        }

        framedata_free(flat);
    }

    while (true)
    {
        uint8_t aperture_type = 1;
        double aperture_size = 0;
        while (true)
        {
            printf("Aperture types are  1: 5sigma,  2: 3FWHM,  3: Manual\n");
            char *ret = prompt_user_input("Select aperture type:", "1");
            aperture_type = atoi(ret);
            free(ret);
            if (aperture_type >=1 && aperture_type <= 3)
                break;

            printf("Choice must be between 1 and 3\n");
        }

        if (aperture_type == 3)
        {
            while (true)
            {
                char *ret = prompt_user_input("Enter aperture radius (px):", "8");
                aperture_size = strtod(ret, NULL);
                free(ret);
                if (aperture_size > 0)
                    break;

                printf("Radius must be greater than 0\n");
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

        // Set region mode to annulus
        if (ts_exec_write("xpaset tsreduce regions shape annulus", NULL, 0))
            error_jump(frameload_error, ret, "ds9 command failed");

        printf("Circle the target stars and surrounding sky in ds9\nPress enter in this terminal to continue...");
        getchar();

        // Read zoom from ds9
        char *ds9_zoom;
        if (ts_exec_read("xpaget tsreduce zoom", &ds9_zoom))
            error_jump(frameload_error, ret, "ds9 request zoom failed");
        float zoom = strtod(ds9_zoom, NULL);
        free(ds9_zoom);

        // Sanity-check zoom
        if (zoom < 0.001)
            zoom = 1;

        char *ds9buf;
        if (ts_exec_read("xpaget tsreduce regions", &ds9buf))
            error_jump(frameload_error, ret, "ds9 request regions failed");

        // Parse the region definitions
        data->target_count = 0;
        size_t target_size = 5;
        data->targets = malloc(target_size*sizeof(struct target_data));
        double *sky = malloc(target_size*sizeof(double));
        if (!data->targets || !sky)
            error_jump(target_error, ret, "Target allocation error");

        char *cur = ds9buf;
        double largest_aperture = 0;
        for (; (cur = strstr(cur, "annulus")) != NULL; cur++)
        {
            if (data->target_count >= target_size)
            {
                target_size += 5;
                data->targets = realloc(data->targets, target_size*sizeof(struct target_data));
                sky = realloc(sky, target_size*sizeof(double));
                if (!data->targets || !sky)
                    error_jump(target_error, ret, "Target allocation error");
            }

            // Read aperture coords
            aperture a;
            sscanf(cur, "annulus(%lf,%lf,%lf,%lf)", &a.x, &a.y, &a.s1, &a.s2);

            // ds9 denotes the bottom-left pixel (1,1), tsreduce uses (0,0)
            a.x -= 1;
            a.y -= 1;

            // Estimate initial aperture size as inner sky radius
            a.r = a.s1;

            if (verbosity >= 1)
                printf("Initial aperture xy: (%f,%f) r: %f s:(%f,%f)\n", a.x, a.y, a.r, a.s1, a.s2);

            double2 xy;
            if (center_aperture(a, frame, &xy))
            {
                printf("Centering failed to converge. Removing aperture.\n");
                continue;
            }

            a.x = xy.x;
            a.y = xy.y;

            double sky_intensity, sky_std_dev;
            if (calculate_background(a, frame, &sky_intensity, &sky_std_dev))
            {
                printf("Background calculation failed. Removing aperture.\n");
                continue;
            }
            sky[data->target_count] = sky_intensity;

            switch (aperture_type)
            {
                case 1: // 5*sky sigma
                {
                    // Estimate the radius where the star flux falls to 5 times the std. dev. of the background
                    double lastIntensity = 0;
                    double lastProfile = frame->data[frame->cols*((int)xy.y) + (int)xy.x] - sky_intensity;
                    int maxRadius = (int)a.s2 + 1;
                    for (int radius = 1; radius <= maxRadius; radius++)
                    {
                        double intensity = integrate_aperture(xy, radius, frame) - sky_intensity*M_PI*radius*radius;
                        double profile = (intensity - lastIntensity) / (M_PI*(2*radius-1));
                        if (profile < 5*sky_std_dev)
                        {
                            // Linear interpolate radii to estimate the radius that gives 5*stddev
                            a.r = radius - 1 + (5*sky_std_dev - lastProfile) / (profile - lastProfile);
                            largest_aperture = fmax(largest_aperture, a.r);
                            break;
                        }
                        lastIntensity = intensity;
                        lastProfile = profile;
                    }
                    break;
                }
                case 2: // 3*FWHM
                {
                    double fwhm = estimate_fwhm(frame, xy, sky_intensity, a.s1);
                    if (fwhm < 0)
                    {
                        printf("Invalid fwhm. Removing target\n");
                        continue;
                    }

                    a.r = 3*fwhm/2;
                    largest_aperture = fmax(largest_aperture, a.r);

                    break;
                }
                case 3:
                default:
                    a.r = aperture_size;
                    largest_aperture = aperture_size;
                    break;
            }

            // Set target parameters
            data->targets[data->target_count].aperture = a;
            data->targets[data->target_count].scale = 1.0;

            char label[16];
            if (data->target_count == 0)
                strcpy(label, "Target");
            else
                snprintf(label, 16, "Comparison %zu", data->target_count);
            data->targets[data->target_count].label = strdup(label);

            data->target_count++;
        }
        free(ds9buf);

        // Set aperture radii to the same size, equal to the largest calculated above
        for (size_t i = 0; i < data->target_count; i++)
            data->targets[i].aperture.r = largest_aperture;

        printf("Aperture radius: %.2fpx\n", largest_aperture);

        // Display results in ds9 - errors are non-fatal
        ts_exec_write("xpaset tsreduce regions delete all", NULL, 0);

        for (size_t i = 0; i < data->target_count; i++)
        {
            aperture *t = &data->targets[i].aperture;
            double x = t->x + 1;
            double y = t->y + 1;

            char command[1024];
            snprintf(command, 1024, "xpaset tsreduce regions command '{circle %f %f %f #color=red select=0}'", x, y, t->r);
            ts_exec_write(command, NULL, 0);
            snprintf(command, 1024, "xpaset tsreduce regions command '{annulus %f %f %f %f #background select=0}'", x, y, t->s1, t->s2);
            ts_exec_write(command, NULL, 0);

            double intensity = frame->data[frame->cols*((size_t)t->y) + (size_t)t->x] - sky[i];

            snprintf(command, 1024, "xpaset -p tsreduce regions command '{text %f %f #color=green select=0 text=\"%s\"}'", x, y - t->s2 - 10/zoom, data->targets[i].label);
            ts_exec_write(command, NULL, 0);
            snprintf(command, 1024, "xpaset -p tsreduce regions command '{text %f %f #color=green select=0 text=\"%.0f ADU\"}'", x, y - t->s2 - 25/zoom, intensity);
            ts_exec_write(command, NULL, 0);

            snprintf(command, 1024, "xpaset -p tsreduce regions command '{box %f %f %f %f #color=black select=0 width=%d}'",
                     x, y - t->s2 - 17.5/zoom, 115/zoom, 17.5/zoom, 18);
            ts_exec_write(command, NULL, 0);
        }
        ts_exec_write("xpaset -p tsreduce update now", NULL, 0);

        char *ret = prompt_user_input("Are the displayed apertures correct?:", "y");
        bool done = !strcmp(ret, "y");
        free(ret);
        free(sky);

        if (done)
            break;

        free(data->targets);
    }
    printf("Set %zu targets\n", data->target_count);

    // Save to disk
    if (chdir(datadir))
        error_jump(frameload_error, ret, "Invalid data path: %s", datadir);

    datafile_save(data, outname);
    printf("Saved to %s\n", outname);
    printf("Run `tsreduce update %s` to reduce existing files\n", outname);
    printf("Run `tsreduce plot %s` to preview reduced data\n", outname);

target_error:
frameload_error:
    framedata_free(frame);
create_flat_error:
create_dark_error:
    datafile_free(data);
    return ret;
}

// Update the ds9 preview with a frame and FWHM calculations
int update_preview(char *preview_filename, char *ds9_title, char *autoguide_output)
{
    int ret = 0;
    char ds9_command_buf[1024];

    framedata *frame = framedata_load(preview_filename);
    if (!frame)
        error_jump(frame_error, ret, "Error loading frame %s", preview_filename);
    framedata_subtract_bias(frame);

    double plate_scale = 1;
    if (framedata_get_metadata(frame, "IM-SCALE", FRAME_METADATA_DOUBLE, &plate_scale))
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
    float zoom = strtod(ds9_zoom, NULL);
    free(ds9_zoom);

    // Sanity-check zoom
    if (zoom < 0.001)
        zoom = 1;

    // Set new frame
    snprintf(ds9_command_buf, 1024, "xpaset -p %s file %s", ds9_title, preview_filename);
    ts_exec_write(ds9_command_buf, NULL, 0);

    // Parse the region definitions
    char *cur = ds9_regions;
    bool first = true;
    for (; (cur = strstr(cur, "annulus")) != NULL; cur++)
    {
        // Read aperture coords
        aperture a;
        int select = 1;
        sscanf(cur, "annulus(%lf,%lf,%lf,%lf) # color=red select=%d", &a.x, &a.y, &a.s1, &a.s2, &select);
        if (!select)
            continue;

        // ds9 denotes the bottom-left pixel (1,1), tsreduce uses (0,0)
        a.x -= 1;
        a.y -= 1;

        // Estimate initial aperture size as inner sky radius
        a.r = a.s1;

        double2 xy;
        if (center_aperture(a, frame, &xy))
        {
            error("Aperture converge failed. Removing target");
            continue;
        }
        a.x = xy.x;
        a.y = xy.y;

        // Print centroid of the first target to stdout
        if (autoguide_output && first)
        {
            first = false;

            FILE *guide = fopen(autoguide_output, "a");
            if (!guide)
                error("Failed to open autoguiding output file %s", autoguide_output);
            else
            {
                ts_time start;
                if (framedata_start_time(frame, &start))
                    error("Failed to query frame start time. Skipping centroid output");
                else
                    fprintf(guide, "%lu.%03u %.2f %.2f\n", start.time, start.ms, xy.x, xy.y);

                fclose(guide);
            }
        }

        double sky_intensity, sky_std_dev;
        if (calculate_background(a, frame, &sky_intensity, &sky_std_dev))
        {
            error("Sky calculation failed. Removing target");
            continue;
        }

        double fwhm = estimate_fwhm(frame, xy, sky_intensity, a.s1);
        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{annulus %f %f %f %f}'",
                 ds9_title, a.x + 1, a.y + 1, a.s1, a.s2);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{annulus %f %f %f %f #color=red select=0 background}'",
                 ds9_title, a.x + 1, a.y + 1, a.s1, a.s2);
        ts_exec_write(ds9_command_buf, NULL, 0);

        if (!isnan(fwhm) && fwhm > 0)
        {
            snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{circle %f %f %f #color=red select=0}'",
                     ds9_title, a.x + 1, a.y + 1, fwhm/2);
            ts_exec_write(ds9_command_buf, NULL, 0);
        }

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{text %f %f #color=green select=0 text=\"FWHM: %.2f arcsec\"}'",
                 ds9_title, a.x + 1, a.y + 1 - a.s2 - 10/zoom, fwhm*plate_scale);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{text %f %f #color=green select=0 text=\"Peak: %.0f ADU/px\"}'",
                 ds9_title, a.x + 1, a.y + 1 - a.s2 - 25/zoom, frame->data[frame->cols*((int)xy.y) + (int)xy.x] - sky_intensity);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{text %f %f #color=green select=0 text=\"BG: %.0f ADU/px\"}'",
                 ds9_title, a.x + 1, a.y + 1 - a.s2 - 40/zoom, sky_intensity);
        ts_exec_write(ds9_command_buf, NULL, 0);

        snprintf(ds9_command_buf, 1024, "xpaset -p %s regions command '{box %f %f %f %f #color=black select=0 width=%d}'",
                 ds9_title, a.x + 1, a.y + 1 - a.s2 - 25/zoom, 115/zoom, 25/zoom, 25);
        ts_exec_write(ds9_command_buf, NULL, 0);
    }

    // Display frame time
    double frame_exp = 0;
    char *frame_end = NULL;
    char *frame_date = NULL;
    char *frame_object = NULL;

    // Ignore errors
    if (framedata_get_metadata(frame, "UTC-END", FRAME_METADATA_STRING, &frame_end))
        frame_end = strdup("Unknown");
    if (framedata_get_metadata(frame, "UTC-DATE", FRAME_METADATA_STRING, &frame_date))
        frame_date = strdup("Unknown");
    if (framedata_get_metadata(frame, "OBJECT", FRAME_METADATA_STRING, &frame_object))
        frame_object = strdup("Unknown");
    framedata_get_metadata(frame, "EXPTIME", FRAME_METADATA_DOUBLE, &frame_exp);

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

int frame_translation(const char *frame_path, const char *reference_path, const char *dark_path, const char *flat_path)
{
    int ret = 0;

    framedata *dark = framedata_load(dark_path);
    if (!dark)
        error_jump(load_error, ret, "Error loading frame %s", dark_path);

    framedata *flat = framedata_load(flat_path);
    if (!flat)
        error_jump(load_error, ret, "Error loading frame %s", flat_path);

    framedata *frame = framedata_load(frame_path);
    if (!frame)
        error_jump(load_error, ret, "Error loading frame %s", frame_path);

    framedata *reference = framedata_load(reference_path);
    if (!reference)
        error_jump(load_error, ret, "Error loading frame %s", reference_path);

    framedata_subtract_bias(frame);
    framedata_subtract_bias(reference);

    // Process frames
    framedata_subtract_bias(frame);
    framedata_subtract_bias(reference);
    if (framedata_subtract(frame, dark))
        error_jump(process_error, ret, "Error dark-subtracting frame %s", frame_path);
    if (framedata_subtract(reference, dark))
        error_jump(process_error, ret, "Error dark-subtracting frame %s", reference_path);

    if (framedata_divide(frame, flat))
        error_jump(process_error, ret, "Error flat-fielding frame %s", frame_path);
    if (framedata_divide(reference, flat))
        error_jump(process_error, ret, "Error flat-fielding frame %s", reference_path);

    int32_t xt, yt;
    if (framedata_estimate_translation(frame, reference, &xt, &yt))
        error_jump(process_error, ret, "Error calculating translation between %s and %s", frame_path, reference_path);

    printf("Translation: %d %d\n", xt, yt);

process_error:
load_error:
    framedata_free(frame);
    framedata_free(reference);
    framedata_free(flat);
    framedata_free(dark);
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
    framedata_subtract_bias(frame);

    framedata *dark = framedata_load(data->dark_template);
    if (!dark)
        error_jump(dark_error, ret, "Error loading frame %s", data->dark_template);

    if (framedata_subtract(frame, dark))
        error_jump(flat_error, ret, "Error dark-subtracting frame %s", filename);

    framedata *flat = framedata_load(data->flat_template);
    if (!flat)
        error_jump(flat_error, ret, "Error loading frame %s", data->flat_template);

    if (framedata_divide(frame, flat))
        error_jump(process_error, ret, "Error flat-fielding frame %s", filename);

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
        if (framedata_get_metadata(flat, "CCD-READ", FRAME_METADATA_DOUBLE, &readnoise))
            readnoise = data->ccd_readnoise;

        if (readnoise <= 0)
            error_jump(process_error, ret, "CCD Read noise unknown. Define CCDReadNoise in %s.", dataPath);

        if (framedata_get_metadata(flat, "CCD-GAIN", FRAME_METADATA_DOUBLE, &gain))
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

int process_ccdtime(const char *pattern, const char *date, const char *time)
{
    ts_time start = parse_date_time(date, time);

    char **frame_paths;
    size_t num_frames = get_matching_files(pattern, &frame_paths);
    for (size_t i = 0; i < num_frames; i++)
    {
        framedata *frame = framedata_load(frame_paths[i]);
        
        ts_time frame_time;
        framedata_start_time(frame, &frame_time);
        double gps = ts_difftime(frame_time, start);
        double ccd = 0;
        framedata_get_metadata(frame, "CCD-TIME", FRAME_METADATA_DOUBLE, &ccd);

        printf("%s %f %f %f\n", frame_paths[i], gps, ccd, gps - ccd);
        framedata_free(frame);
    }
    return 0;
}
