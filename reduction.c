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
#include <inttypes.h>

#include "datafile.h"
#include "framedata.h"
#include "helpers.h"
#include "aperture.h"
#include "fit.h"
#include "hashmap.h"

extern int verbosity;

double calculate_readnoise(framedata *master, double *data_cube, size_t num_frames, uint16_t region[4], bool subtract_master)
{
    size_t region_px = (region[1] - region[0])*(region[3] - region[2]);
    double variance = 0;

    for (uint16_t j = region[2]; j < region[3]; j++)
    {
        for (uint16_t i = region[0]; i < region[1]; i++)
        {
            for (size_t k = 0; k < num_frames; k++)
            {
                double diff = data_cube[num_frames*(j*master->cols + i) + k];

                // Bias frames need to have the master bias subtracted to find the difference
                // Flat frames (overscan) have already been bias subtracted by the dark frame
                if (subtract_master)
                    diff -= master->data[j*master->cols + i];

                variance += diff * diff;
            }
        }
    }
    variance /= region_px * num_frames;

    return sqrt(variance);
}

// Create a flat field frame from the frames listed by the command `flatcmd',
// rejecting `minmax' highest and lowest pixel values after subtracting
// the dark frame `masterdark'.
// Save the resulting image to the file `outname'
//
// Also calculates the readnoise and gain if overscan is available
int create_flat(const char *pattern, size_t minmax, const char *masterbias, const char *masterdark, const char *outname)
{
    int ret = 0;

    // Find the filenames that match the specified pattern
    char **frame_paths;
    size_t num_frames = get_matching_files(pattern, &frame_paths);

    // Ensure there are enough frames to discard the requested number of pixels
    if (num_frames <= 2*minmax)
        error_jump(insufficient_frames_error, ret,
            "    Insufficient frames. %d found, %d will be discarded", num_frames, 2*minmax);

    // Load the master bias and dark frames, plus the first flat frame to use as a reference.
    // Frame geometry for all subsequent frames is assumed to match the base frame
    framedata *bias = NULL;
    framedata *dark = NULL;
    framedata *base = NULL;

    if (masterbias)
    {
        bias = framedata_load(masterbias);
        if (!bias)
            error_jump(setup_error, ret, "    Error loading bias frame %s", masterbias);
    }

    dark = framedata_load(masterdark);
    if (!dark)
        error_jump(setup_error, ret, "    Error loading dark frame %s", masterdark);

    base = framedata_load(frame_paths[0]);
    if (!base)
        error_jump(setup_error, ret, "    Error loading base frame %s", frame_paths[0]);

    if (base->rows != dark->rows || base->cols != dark->cols)
        error_jump(setup_error, ret, "    Dark and flat frame sizes don't match");

    uint16_t image_region[4];
    framedata_image_region(base, image_region);
    size_t image_region_px = (image_region[1] - image_region[0])*(image_region[3] - image_region[2]);

    // Data cube for processing the flat data
    //        data[0] = frame[0][0,0], data[1] = frame[1][0,0] ... data[num_frames] = frame[0][1,0] etc
    double *data_cube = calloc(num_frames*base->cols*base->rows, sizeof(double));
    double *frame_mean = calloc(num_frames, sizeof(double));
    if (!data_cube || !frame_mean)
        error_jump(processing_error, ret, "    Allocation failed");

    for (size_t k = 0; k < num_frames; k++)
    {
        if (verbosity >= 1)
            printf("        loading `%s`\n", frame_paths[k]);

        framedata *frame = framedata_load(frame_paths[k]);
        if (!frame)
            error_jump(processing_error, ret, "    Error loading frame %s", frame_paths[k]);

        if (frame->rows != base->rows || frame->cols != base->cols)
        {
            framedata_free(frame);
            error_jump(processing_error, ret,
                "    Frame %s dimensions mismatch. Expected (%d,%d), was (%d, %d)",
                frame_paths[k], base->rows, base->cols, frame->rows, frame->cols);
        }

        // Subtract bias and scaled dark frame
        if (framedata_calibrate(frame, bias, dark, NULL))
        {
            framedata_free(frame);
            error_jump(processing_error, ret, "    Calibration failed for %s", frame_paths[k]);
        }

        // Store data in cube for processing
        for (size_t j = 0; j < base->rows*base->cols; j++)
            data_cube[num_frames*j + k] = frame->data[j];

        // Calculate normalization factor to make the mean image intensity unity
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
            error_jump(processing_error, ret, "    cube_slice alloc failed");

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

    // Calculate gain and (if required) readnoise.
    uint16_t bias_region[4];
    framedata_bias_region(base, bias_region);
    size_t bias_region_px = (bias_region[1] - bias_region[0])*(bias_region[3] - bias_region[2]);

    // Remove existing definitions if they exist (CCD-GAIN is reused by the Puoko-nui acquisition software)
    framedata_remove_metadata(base, "CCD-READ");
    framedata_remove_metadata(base, "CCD-GAIN");

    if (bias || bias_region_px)
    {
        double readnoise;
        if (!bias)
        {
            readnoise = calculate_readnoise(base, data_cube, num_frames, bias_region, false);
            printf("    Calculated CCD-READ: %f\n", readnoise);
        }
        else if (framedata_get_metadata(bias, "CCD-READ", FRAME_METADATA_DOUBLE, &readnoise))
        {
            int unused;
            error_jump(skip_readgain, unused, "    Bias does not specify CCD-READ.");
        }

        // Calculate gain from each image
        double *gain = calloc(num_frames, sizeof(double));
        if (!gain)
            error_jump(processing_error, ret, "    gain alloc failed");

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
                printf("    %zu var: %f mean: %f dark: %f gain: %f\n", k, var, frame_mean[k], mean_dark, gain[k]);
        }
        mean_gain /= num_frames;

        // Find median value
        qsort(gain, num_frames, sizeof(double), compare_double);
        double median_gain = gain[num_frames/2];

        // Set header keys for readout noise and gain
        if (bias || bias_region_px)
        {
            framedata_put_metadata(base, "CCD-READ", FRAME_METADATA_DOUBLE, &readnoise, "Estimated read noise (ADU)");
            framedata_put_metadata(base, "CCD-GAIN", FRAME_METADATA_DOUBLE, &median_gain, "Estimated gain (electrons/ADU)");

            printf("    Calculated CCD-GAIN: %f\n", median_gain);
        }

        free(gain);
    }
skip_readgain:
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
    framedata_free(base);
    framedata_free(dark);
    framedata_free(bias);
insufficient_frames_error:
    free_2d_array(frame_paths, num_frames);

    return ret;
}

int calibration_helper_create_flat(datafile *data, const char *pattern, int discard_minmax)
{
    return create_flat(pattern, discard_minmax, data->bias_template, data->dark_template, data->flat_template);
}

int create_bias_dark_internal(const char *pattern, size_t discard_minmax,
    void (*preprocess)(framedata *, void *), void *preprocess_data,
    void (*postprocess)(framedata *master, double *data_cube, size_t num_frames, void *), void *postprocess_data,
    const char *outname)
{
    int ret = 0;

    char **frame_paths;
    size_t num_frames = get_matching_files(pattern, &frame_paths);

    if (num_frames < 2*discard_minmax)
        error_jump(insufficient_frames_error, ret,
            "Insufficient frames. %d found, %d will be discarded", num_frames, 2*discard_minmax);

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

        preprocess(frame, preprocess_data);

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
        for (int i = discard_minmax; i < num_frames - discard_minmax; i++)
            base->data[j] += data_cube[num_frames*j + i];
        base->data[j] /= (num_frames - 2*discard_minmax);
    }

    postprocess(base, data_cube, num_frames, postprocess_data);
    framedata_save(base, outname);

process_error:
    free(data_cube);
setup_error:
    framedata_free(base);
insufficient_frames_error:
    free_2d_array(frame_paths, num_frames);

    return ret;
}

void preprocess_bias(framedata *frame, void *data)
{
    double *fudge = data;

    // Apply bias offset to allow for proper dark scaling
    for (size_t i = 0; i < frame->cols*frame->rows; i++)
        frame->data[i] += *fudge;
}

void postprocess_bias(framedata *master, double *data_cube, size_t num_frames, void *data)
{
    // Use the full bias frame to calculate the read noise
    uint16_t region[4] = {0, master->cols, 0, master->rows};
    double readnoise = calculate_readnoise(master, data_cube, num_frames, region, true);

    framedata_put_metadata(master, "CCD-READ", FRAME_METADATA_DOUBLE, &readnoise, "Estimated read noise (ADU)");
    printf("    Calculated CCD-READ: %f\n", readnoise);
}

int create_bias(const char *pattern, size_t discard_minmax, double bias_fudge, const char *outname)
{
    return create_bias_dark_internal(pattern, discard_minmax, preprocess_bias, &bias_fudge, postprocess_bias, NULL, outname);
}

int calibration_helper_create_bias(datafile *data, const char *pattern, int discard_minmax)
{
    char *fudge_prompt = prompt_user_input("    Enter bias fudge offset", "0", false);
    double bias_fudge = atof(fudge_prompt);
    free(fudge_prompt);

    return create_bias(pattern, discard_minmax, bias_fudge, data->bias_template);
}

void preprocess_dark(framedata *frame, void *data)
{
    framedata *bias = data;
    if (bias)
        framedata_subtract(frame, bias);
    else
        framedata_subtract_bias(frame);
}

void postprocess_dark(framedata *master, double *data_cube, size_t num_frames, void *data)
{
    // No bias frame
    if (data == NULL)
        return;

    // Set the overscan region to zero -- bias subtraction is covered by the bias frame
    uint16_t bias_region[4];
    framedata_bias_region(master, bias_region);
    size_t bias_region_px = (bias_region[1] - bias_region[0])*(bias_region[3] - bias_region[2]);

    if (bias_region_px)
        for (uint16_t j = bias_region[2]; j < bias_region[3]; j++)
            for (uint16_t i = bias_region[0]; i < bias_region[1]; i++)
                master->data[j*master->cols + i] = 0;
}

int create_dark(const char *pattern, size_t discard_minmax, const char *masterbias, const char *outname)
{
    int ret = 0;
    framedata *bias = NULL;
    if (masterbias)
    {
        bias = framedata_load(masterbias);
        if (!bias)
            error_jump(bias_error, ret, "        Error loading frame %s", masterbias);
    }

    ret = create_bias_dark_internal(pattern, discard_minmax, preprocess_dark, bias, postprocess_dark, (void *)masterbias, outname);

    framedata_free(bias);
bias_error:
    return ret;
}

int calibration_helper_create_dark(datafile *data, const char *pattern, int discard_minmax)
{
    return create_dark(pattern, discard_minmax, data->bias_template, data->dark_template);
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

    if (framedata_calibrate_load(frame, data->bias_template, data->dark_template, data->flat_template))
        error_jump(process_error, ret, "Error processing frame %s", obs->filename);

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

    framedata *bias = NULL;
    if (data->bias_template)
    {
        bias = framedata_load(data->bias_template);
        if (!bias)
            error_jump(bias_error, ret, "Error loading frame %s", data->bias_template);
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

    if (framedata_calibrate(reference, bias, dark, flat))
        error_jump(reference_error, ret, "Error calibrating reference frame %s", data->reference_frame);

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
        if (framedata_calibrate(frame, bias, dark, flat))
            error_jump(process_error, ret, "Error calibrating frame %s", frame_paths[i]);

        struct observation *obs = datafile_new_observation(data);
        if (!obs)
            error_jump(process_error, ret, "Allocation error");

        // Observation mid time
        obs->time = midtime;

        // Process frame
        double nan = sqrt(-1);

        // Estimate translation from reference frame
        int32_t xt, yt;
        bool rotated;
        if (framedata_estimate_translation(frame, reference, &xt, &yt, &rotated))
            error_jump(process_error, ret, "Error calculating frame translation");

        for (size_t i = 0; i < data->target_count; i++)
        {
            // Offset aperture position by frame offset
            aperture a = data->targets[i].aperture;

            a.x += xt;
            a.y += yt;
            if (rotated)
            {
                a.x = frame->cols - a.x - 1;
                a.y = frame->rows - a.y - 1;
            }

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
    framedata_free(dark);
bias_error:
    framedata_free(bias);
flat_error:
    framedata_free(flat);
data_error:
    datafile_free(data);
    return ret;
}

static char *filename_fmt = "^%s(-|.|_)[0-9]+.(fits.gz|fit.gz|fits|fit|FIT)";

int create_master_calibration_file(char *default_name, char *default_prefix, char **out_template, int (*create_calibration_frame)(datafile *data, const char *pattern, int discard_minmax), datafile *data)
{
    char *filename = prompt_user_input("    Enter filename, or ^D to skip", default_name, true);
    if (!filename)
    {
        printf("    Skipping.\n");
        return 0;
    }

    char *path = NULL;
    while (true)
    {
        char *ret = prompt_user_input("    Enter frame directory", ".", false);
        path = canonicalize_path(ret);

        if (!chdir(path))
        {
            free(ret);
            break;
        }

        printf("        Invalid frame directory: %s\n", ret);
        free(path);
        free(ret);
    }

    size_t path_length = strlen(filename) + strlen(path) + 2;
    *out_template = malloc(path_length*sizeof(char));
    snprintf(*out_template, path_length, "%s/%s", path, filename);
    free(filename);

    if (access(*out_template, F_OK) != -1)
    {
        printf("    File exists. Using this file.\n");
        return 0;
    }

    char *pattern;
    int count;
    while (true)
    {
        char *prefix = regex_escape_string(prompt_user_input("    Enter frame prefix", default_prefix, false));
        int len = snprintf(NULL, 0, filename_fmt, prefix);
        pattern = malloc(len + 1);
        sprintf(pattern, filename_fmt, prefix);
        free(prefix);

        char **filenames;
        count = get_matching_files(pattern, &filenames);
        if (count > 0)
        {
            free_2d_array(filenames, count);
            break;
        }

        printf("        No files found matching pattern: %s/%s\n", path, pattern);
        free(pattern);
    }
    free(path);

    int minmax = 0;
    while (true)
    {
        char fallback[32];
        snprintf(fallback, 32, "%d", count / 2);
        char *ret = prompt_user_input("    Enter number of frames around median to average", fallback, false);
        int average = atoi(ret);
        free(ret);

        if (average > 0 && average <= count)
        {
            minmax = (count - average) / 2;
            break;
        }

        printf("        Number must be between 1 and %d\n", count);
    }

    int ret = create_calibration_frame(data, pattern, minmax);
    free(pattern);
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

    printf("Configure master bias frame:\n");
    if (create_master_calibration_file("master-bias.fits.gz", "bias", &data->bias_template, calibration_helper_create_bias, data))
        error_jump(create_bias_error, ret, "ERROR: master bias generation failed");

    chdir(datadir);
    printf("Configure master dark frame:\n");
    if (create_master_calibration_file("master-dark.fits.gz", "dark", &data->dark_template, calibration_helper_create_dark, data))
        error_jump(create_dark_error, ret, "ERROR: master dark generation failed");

    chdir(datadir);
    printf("Configure master flat frame:\n");
    
    if (!data->dark_template)
        printf("    Requires master dark. Skipping.\n");
    else if (create_master_calibration_file("master-flat.fits.gz", "flat", &data->flat_template, calibration_helper_create_flat, data))
        error_jump(create_flat_error, ret, "ERROR: master flat generation failed");

    chdir(datadir);
    printf("Configure reduction:\n");

    while (true)
    {
        char *ret = prompt_user_input("    Enter frame path", ".", false);
        data->frame_dir = canonicalize_path(ret);
        if (!chdir(data->frame_dir))
        {
            free(ret);
            break;
        }
        printf("        Invalid frame path: %s\n", ret);
        free(data->frame_dir);
        free(ret);
    }

    // Default target prefix to filename
    char *default_prefix = remove_file_suffix(strdup(last_path_component(outname)));
    while (true)
    {
        char namebuf[1039];
        char *ret = regex_escape_string(prompt_user_input("    Enter target prefix", default_prefix, false));
        snprintf(namebuf, 1039, filename_fmt, ret);
        free(ret);
        data->reference_frame = get_first_matching_file(namebuf);
        if (data->reference_frame)
        {
            data->frame_pattern = strdup(namebuf);
            break;
        }
        printf("        No files found matching pattern: %s/%s\n", data->frame_dir, namebuf);
    }
    free(default_prefix);

    framedata *frame = NULL;
    while (true)
    {
        char *ret = prompt_user_input("    Enter reference frame", data->reference_frame, false);
        frame = framedata_load(ret);
        if (frame)
        {
            free(data->reference_frame);
            data->reference_frame = ret;
            break;
        }

        printf("        File not found: %s/%s\n", data->frame_dir, ret);
        free(ret);
    }

    if (framedata_calibrate_load(frame, data->bias_template, data->dark_template, data->flat_template))
        error_jump(frameload_error, ret, "    Unable to calibrate reference frame");

    if (framedata_start_time(frame, &data->reference_time))
        error_jump(frameload_error, ret, "    No known time headers found");

    if (data->flat_template)
    {
        framedata *flat = framedata_load(data->flat_template);
        if (!flat)
            error_jump(frameload_error, ret, "    Error loading frame %s", data->flat_template);

        if (!framedata_has_metadata(flat, "CCD-READ"))
        {
            while (true)
            {
                char *ret = prompt_user_input("    Enter CCD Readnoise (ADU):", "3.32", false);
                data->ccd_readnoise = strtod(ret, NULL);
                free(ret);
                if (data->ccd_readnoise > 0)
                    break;

                printf("        Number must be greater than 0\n");
            }
        }

        if (!framedata_has_metadata(flat, "CCD-GAIN"))
        {
            while (true)
            {
                char *ret = prompt_user_input("    Enter CCD Gain (ADU):", "2.00", false);
                data->ccd_gain = strtod(ret, NULL);
                free(ret);

                if (data->ccd_gain > 0)
                    break;

                printf("        Number must be greater than 0\n");
            }
        }

        if (framedata_get_metadata(flat, "IM-SCALE", FRAME_METADATA_DOUBLE, &data->ccd_platescale) && framedata_get_metadata(flat, "SECPIX", FRAME_METADATA_DOUBLE, &data->ccd_platescale))
        {
            char *ret = prompt_user_input("    Enter CCD platescale (arcsec/px):", "0.66", false);
            data->ccd_platescale = strtod(ret, NULL);
            free(ret);
        }

        framedata_free(flat);
    }

    while (true)
    {
        uint8_t aperture_type = 1;
        double aperture_size = 0;
        while (true)
        {
            char *ret = prompt_user_input("    Select aperture type (1: 5sigma,  2: 3FWHM,  3: Manual):", "1", false);
            aperture_type = atoi(ret);
            free(ret);
            if (aperture_type >=1 && aperture_type <= 3)
                break;

            printf("        Choice must be between 1 and 3\n");
        }

        if (aperture_type == 3)
        {
            while (true)
            {
                char *ret = prompt_user_input("Enter aperture radius (px):", "8", false);
                aperture_size = strtod(ret, NULL);
                free(ret);
                if (aperture_size > 0)
                    break;

                printf("        Radius must be greater than 0\n");
            }
        }

        if (init_ds9())
            return error("    Unable to launch ds9");

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

        printf("    Circle the target stars and surrounding sky in ds9\n        Press enter in this terminal to continue...");
        getchar();

        // Read zoom from ds9
        char *ds9_zoom;
        if (ts_exec_read("xpaget tsreduce zoom", &ds9_zoom))
            error_jump(frameload_error, ret, "        ds9 request zoom failed");
        float zoom = strtod(ds9_zoom, NULL);
        free(ds9_zoom);

        // Sanity-check zoom
        if (zoom < 0.001)
            zoom = 1;

        char *ds9buf;
        if (ts_exec_read("xpaget tsreduce regions", &ds9buf))
            error_jump(frameload_error, ret, "        ds9 request regions failed");

        // Parse the region definitions
        data->target_count = 0;
        size_t target_size = 5;
        data->targets = malloc(target_size*sizeof(struct target_data));
        double *sky = malloc(target_size*sizeof(double));
        if (!data->targets || !sky)
            error_jump(target_error, ret, "        Target allocation error");

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
                    error_jump(target_error, ret, "    Target allocation error");
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
                printf("        Initial aperture xy: (%f,%f) r: %f s:(%f,%f)\n", a.x, a.y, a.r, a.s1, a.s2);

            double2 xy;
            if (center_aperture(a, frame, &xy))
            {
                printf("        Centering failed to converge. Removing aperture.\n");
                continue;
            }

            a.x = xy.x;
            a.y = xy.y;

            double sky_intensity, sky_std_dev;
            if (calculate_background(a, frame, &sky_intensity, &sky_std_dev))
            {
                printf("        Background calculation failed. Removing aperture.\n");
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
                        printf("        Invalid fwhm. Removing target\n");
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

        printf("        Aperture radius: %.2fpx\n", largest_aperture);

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

        char *ret = prompt_user_input("    Are the displayed apertures correct?:", "y", false);
        bool done = !strcmp(ret, "y");
        free(ret);
        free(sky);

        if (done)
            break;

        free(data->targets);
    }
    printf("        Set %zu targets\n", data->target_count);

    // Save to disk
    if (chdir(datadir))
        error_jump(frameload_error, ret, "    Invalid data path: %s", datadir);

    datafile_save(data, outname);
    printf("Saved to %s\n", outname);
    printf("Run `tsreduce update %s` to reduce existing files\n", outname);
    printf("Run `tsreduce plot %s` to preview reduced data\n", outname);

target_error:
frameload_error:
    framedata_free(frame);
create_flat_error:
create_dark_error:
create_bias_error:
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
                    fprintf(guide, "%" PRIuMAX ".%03" PRIu16 " %.2f %.2f\n", (intmax_t)start.time, start.ms, xy.x, xy.y);

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
int calculate_bjd(char *date, char *time, char *ra_string, char *dec_string)
{
    // Convert ra from HH:MM:SS to radians
    double a,b,c;
    sscanf(ra_string, "%lf:%lf:%lf", &a, &b, &c);
    double ra = (a + b/60 + c/3600)*M_PI/12;

    // Convert dec from DD:'':"" to radians
    sscanf(dec_string, "%lf:%lf:%lf", &a, &b, &c);
    double dec = copysign((fabs(a) + b/60 + c/3600)*M_PI/180, a);

    printf("%.8Lf\n", ts_time_to_bjd(parse_date_time(date, time), ra, dec));
    return 0;
}

/*
 * Create a timeseries file from a list of datafiles.
 * Times are converted to BJD after the specified reference
 *
 * If comparison_magnitude != 0 then output absolute magnitudes (relative to comparison)
 * instead of mmi, and the polynomial subtraction is not done.
 */
int create_ts(char *reference_date, char *reference_time, char **filenames, size_t num_datafiles, char *ts_filename, bool use_ratio, float comparison_magnitude)
{
    bool output_mag = comparison_magnitude != 0.0;
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

    if (!datafiles[0]->coord_ra || !datafiles[0]->coord_dec)
        error_jump(coord_error, ret, "Datafile %s doesn't specify star coordinates", filenames[0]);

    // Convert ra from HH:MM:SS to radians
    double a,b,c;
    sscanf(datafiles[0]->coord_ra, "%lf:%lf:%lf", &a, &b, &c);
    double ra = (a + b/60 + c/3600)*M_PI/12;

    // Convert dec from DD:'':"" to radians
    sscanf(datafiles[0]->coord_dec, "%lf:%lf:%lf", &a, &b, &c);
    double dec = copysign((fabs(a) + b/60 + c/3600)*M_PI/180, a);

    double reference_bjd = ts_time_to_bjd(parse_date_time(reference_date, reference_time), ra, dec);
    printf("Reference BJD: %f\n", reference_bjd);

    FILE *out = fopen(ts_filename, "w+");
    if (!out)
        error_jump(output_error, ret, "Error opening file %s", ts_filename);

    // Print file header
    if (use_ratio)
    {
        if (output_mag)
        {
            fprintf(out, "# tsreduce create-mag output file\n");
            fprintf(out, "# Comparison Magnitude: %.3f\n", comparison_magnitude);        
        }
        else
            fprintf(out, "# tsreduce create-ratio output file\n");
    }
    else
        fprintf(out, "# tsreduce create-ts output file\n");

    fprintf(out, "# Reference time: %s %s UTC; %.8f BJD\n", reference_date, reference_time, reference_bjd);
    fprintf(out, "# Files:\n");
    for (size_t i = 0; i < num_datafiles; i++)
        fprintf(out, "#   %s\n", filenames[i]);

    // Convert data to BJD relative to the reference time
    size_t num_saved = 0;
    for (size_t i = 0; i < num_datafiles; i++)
    {
        long double start_bjd = ts_time_to_bjd(datafiles[i]->reference_time, ra, dec);

        // Calculate precessed RA and DEC at the start of each night
        // This is already far more accurate than we need
        printf("%s start BJD: %Lf\n", filenames[i], start_bjd);

        struct photometry_data *pd = datafile_generate_photometry(datafiles[i]);
        if (!pd)
            error_jump(processing_error, ret, "Error generating photometry data for data %s", filenames[i]);

        for (size_t j = 0; j < pd->filtered_count; j++)
        {
            ts_time obstime = datafiles[i]->reference_time;
            obstime.time += (time_t)(pd->time[j]);
            obstime.ms += round(1000*fmod(pd->time[j], 1));

            double intensity = use_ratio ? pd->ratio[j] : pd->mmi[j];
            double error = use_ratio ? pd->ratio_noise[j] : pd->mmi_noise[j];            
            if (output_mag && use_ratio)
            {
                intensity = comparison_magnitude - 2.5 * log10(intensity);
                error = 2.5 / log(10) * error / intensity;
            }
            
            fprintf(out, "%.8Lf %f %f\n", ts_time_to_bjd(obstime, ra, dec) - reference_bjd, intensity, error);
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

int frame_translation(const char *frame_path, const char *reference_path, const char *bias_path, const char *dark_path, const char *flat_path)
{
    int ret = 0;

    framedata *frame = framedata_load(frame_path);
    if (!frame)
        error_jump(frame_error, ret, "Error loading frame %s", frame_path);

    framedata *reference = framedata_load(reference_path);
    if (!reference)
        error_jump(reference_error, ret, "Error loading frame %s", reference_path);

    if (framedata_calibrate_load(frame, bias_path, dark_path, flat_path))
        error_jump(process_error, ret, "Error dark-subtracting frame %s", frame_path);

    if (framedata_calibrate_load(reference, bias_path, dark_path, flat_path))
        error_jump(process_error, ret, "Error dark-subtracting frame %s", reference_path);

    int32_t xt, yt;
    bool rotated;
    if (framedata_estimate_translation(frame, reference, &xt, &yt, &rotated))
        error_jump(process_error, ret, "Error calculating translation between %s and %s", frame_path, reference_path);

    printf("Translation: %d %d%s\n", xt, yt, rotated ? " (after 180 deg rotation)" : "");

process_error:
    framedata_free(reference);
reference_error:
    framedata_free(frame);
frame_error:
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
