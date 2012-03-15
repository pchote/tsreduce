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

#include "tsreduce.h"
#include "framedata.h"
#include "helpers.h"
#include "aperture.h"


// Calculate the standard deviation of pixels in a data cube
// The standard deviation is calculated as the sum of differences from the per-pixel
// mean, over all frames (the deviation is not calculated per-frame, and then averaged)
//
//    Data is a pointer to a data cube
//        data[0] = frame[0][0,0], data[1] = frame[1][0,0] ... data[num_frames] = frame[0][1,0] etc
//    num_frames is the number of frames worth of data
//    num_pixels is the number of pixels in each frame
double calculate_cube_stddev(double *data, int num_frames, int num_pixels)
{
    // Calculate a mean master-bias
    double *mean = (double *)malloc(num_pixels*sizeof(double));
    for (int j = 0; j < num_pixels; j++)
    {
        mean[j] = 0;
        for (int i = 0; i < num_frames; i++)
            mean[j] += data[num_frames*j + i];
        mean[j] /= num_frames;
    }

    // Calculate and return standard deviation from pixel mean over the cube
    double std = 0;
    for (int j = 0; j < num_pixels; j++)
        for (int i = 0; i < num_frames; i++)
            std += (data[num_frames*j + i] - mean[j])*(data[num_frames*j + i] - mean[j]);
    return sqrt(std/(num_pixels*num_frames));
}

// Calculate the CCD read noise from the overscan bias in a collection of frames
// Read noise is given by the standard deviation away from the mean bias level
int ccd_readnoise(const char *framePattern)
{
    char **frames;
    int num_frames = get_matching_files(framePattern, &frames);

    // TODO: Check fits headers for an appropriate header
    int x_min = 525;
    int x_max = 535;
    int y_min = 5;
    int y_max = 508;

    int num_pixels = (y_max - y_min)*(x_max - x_min);
    double *data = (double *)malloc(num_pixels*num_frames*sizeof(double));
    int p = 0;

    // Load overscan region into a data cube
    for (int k = 0; k < num_frames; k++)
    {
        framedata frame = framedata_new(frames[k], FRAMEDATA_DBL);

        // Copy overscan pixels
        for (int j = y_min; j < y_max; j++)
            for (int i = x_min; i < x_max; i++)
                data[p++] = frame.dbl_data[j*frame.cols + i];

        framedata_free(frame);
    }

    double readnoise = calculate_cube_stddev(data, num_frames, num_pixels);
    printf("Read noise from %d frames: %f ADU\n", num_frames, readnoise);

    free(data);
    free_2d_array(frames, num_frames);

    return 0;
}

// Convenience function for calculating the mean signal in a sub-region of a frame
// Assumes that the frame type is double, and that the region is inside the frame
double mean_in_region(framedata *frame, int x_min, int x_max, int y_min, int y_max)
{
    int num_px = (x_max - x_min)*(y_max - y_min);
    double mean = 0;
    for (int j = y_min; j < y_max; j++)
        for (int i = x_min; i < x_max; i++)
            mean += frame->dbl_data[j*frame->cols + i]/num_px;

    return mean;
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

    // Calculate and subtract bias if the frame has overscan
    // TODO: have acquisition software save regions into a header key
    // and use this for determining overscan
    bool use_overscan = (flat->rows != flat->cols);
    int bias_rgn[4] = {525, 535, 5, 508};
    int image_rgn[4] = {0, 512, 0, 512};
    int num_image_px = (image_rgn[1] - image_rgn[0])*(image_rgn[3] - image_rgn[2]);

    double mean_bias = 0;
    if (use_overscan)
        mean_bias = mean_in_region(flat, bias_rgn[0], bias_rgn[1], bias_rgn[2], bias_rgn[3]);

    // Dark subtract frame
    for (int i = 0; i < flat->rows*flat->cols; i++)
        flat->dbl_data[i] -= flatexp*1.0/darkexp*dark->dbl_data[i] + mean_bias;

    // Calculate mean in image area only
    double mean = mean_in_region(flat, image_rgn[0], image_rgn[1], image_rgn[2], image_rgn[3]);

    // Calculate standard deviation
    double std = 0;
    for (int j = image_rgn[2]; j < image_rgn[3]; j++)
        for (int i = image_rgn[0]; i < image_rgn[1]; i++)
        {
            double temp = flat->dbl_data[flat->cols*j + i] - mean;
            std += temp*temp;
        }
    std = sqrt(std/num_image_px);

    // Recalculate the mean, excluding outliers at 3 sigma
    double mean_new = 0;
    int count = 0;
    for (int j = image_rgn[2]; j < image_rgn[3]; j++)
        for (int i = image_rgn[0]; i < image_rgn[1]; i++)
            if (fabs(flat->dbl_data[flat->cols*j + i] - mean) < 3*std)
            {
                mean_new += flat->dbl_data[flat->cols*j + i];
                count++;
            }
    mean_new /= count;

    // Normalize flat so the image region mean is 1
    for (int i = 0; i < flat->rows*flat->cols; i++)
        flat->dbl_data[i] /= mean_new;

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
    int return_status = 0;

    // Find the filenames that match the specified pattern
    char **frame_paths;
    int num_frames = get_matching_files(pattern, &frame_paths);

    // Ensure there are enough frames to discard the requested number of pixels
    if (num_frames <= 2*minmax)
    {
        return_status = error("Insufficient frames. %d found, %d will be discarded", num_frames, 2*minmax);
        goto insufficient_frames;
    }

    // Load the master dark frame
    // Frame geometry for all subsequent frames is assumed to match the master-dark
    framedata dark = framedata_new(masterdark, FRAMEDATA_DBL);

    // Extract image and bias regions
    // TODO: Acquisition code should save this into the fits header
    int bias_rgn[4] = {525, 535, 5, 508};
    int image_rgn[4] = {0, 512, 0, 512};
    int num_image_px = (image_rgn[1] - image_rgn[0])*(image_rgn[3] - image_rgn[2]);

    int num_bias_px = 0;
    bool use_overscan = (dark.rows != dark.cols);
    if (use_overscan)
    {
        num_bias_px = (bias_rgn[1] - bias_rgn[0])*(bias_rgn[3] - bias_rgn[2]);

        // Subtract bias from master-dark
        double mean_bias = mean_in_region(&dark, bias_rgn[0], bias_rgn[1], bias_rgn[2], bias_rgn[3]);
        for (int j = 0; j < dark.rows*dark.cols; j++)
            dark.dbl_data[j] -= mean_bias;
    }

    // Calculate mean dark level
    double mean_dark = mean_in_region(&dark, image_rgn[0], image_rgn[1], image_rgn[2], image_rgn[3]);

    // Data cube for processing the flat data
    //        data[0] = frame[0][0,0], data[1] = frame[1][0,0] ... data[num_frames] = frame[0][1,0] etc
    double *data_cube = (double *)malloc(num_frames*dark.cols*dark.rows*sizeof(double));
    if (data_cube == NULL)
    {
        return_status = error("malloc failed");
        goto datacube_failed;
    }

    // Mean intensity of each flat frame
    double *mean_flat = (double *)malloc(num_frames*sizeof(double));
    if (mean_flat == NULL)
    {
        return_status = error("malloc failed");
        goto meanflat_failed;
    }

    // We need to iterate over the individual frames multiple times to calculate gain
    framedata *frames = (framedata *)malloc(num_frames*sizeof(framedata));
    if (frames == NULL)
    {
        return_status = error("malloc failed");
        goto frames_failed;
    }

    // Load flat field frames into a data cube for the image data, and a cube for the bias data
    // data[0] = flat[0][0,0], data[1] = flat[1][0,0] ... data[num_frames] = flat[0][0,1] etc
    for(int i = 0; i < num_frames; i++)
    {
        printf("loading `%s`\n", frame_paths[i]);
        frames[i] = framedata_new(frame_paths[i], FRAMEDATA_DBL);

        if (frames[i].rows != dark.rows || frames[i].cols != dark.cols)
        {
            // Cleanup the frames that we have allocated
            for (int j = 0; j <= i; j++)
                framedata_free(frames[j]);

            return_status = error("Frame %s dimensions mismatch. Expected (%d,%d), was (%d, %d)", frames[i], dark.rows, dark.cols, frames[i].rows, frames[i].cols);
            goto loadflat_failed;
        }

        // Dark-subtract and normalize the frame, recording the mean level for calculating the gain
        mean_flat[i] = prepare_flat(&frames[i], &dark);

        // Store normalized data in data cube
        for (int j = 0; j < dark.rows*dark.cols; j++)
            data_cube[num_frames*j+i] = frames[i].dbl_data[j];
    }

    // Calculate read noise from overscan region
    // Subtracting master-dark removed the mean level, so just add
    double readnoise = 0;
    if (use_overscan)
    {
        double var = 0;
        for (int j = bias_rgn[2]; j < bias_rgn[3]; j++)
            for (int i = bias_rgn[0]; i < bias_rgn[1]; i++)
                for (int k = 0; k < num_frames; k++)
                {
                    // Need to multiply by the mean level for the frame to give the variance in ADU
                    double temp = mean_flat[k]*data_cube[num_frames*(dark.cols*j + i) + k];
                    var += temp*temp;
                }
        readnoise = sqrt(var/(num_bias_px*num_frames));
    }

    // Calculate median image for the master-flat
    // Loop over the pixels, sorting the values from each image into increasing order
    double *median_flat = (double *)malloc(dark.rows*dark.cols*sizeof(double));
    if (median_flat == NULL)
    {
        return_status = error("malloc failed");
        goto median_failed;
    }

    for (int j = 0; j < dark.rows*dark.cols; j++)
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
    {
        return_status = error("malloc failed");
        goto gain_failed;
    }

    if (use_overscan)
    {
        for(int k = 0; k < num_frames; k++)
        {
            // Calculate the variance by taking the difference between
            // the (normalized) frame and the (normalized) master-flat
            double var = 0;
            for (int j = image_rgn[2]; j < image_rgn[3]; j++)
                for (int i = image_rgn[0]; i < image_rgn[1]; i++)
                {
                    double temp = mean_flat[k]*(frames[k].dbl_data[dark.cols*j + i] - median_flat[dark.cols*j + i]);
                    var += temp*temp;
                }
            var /= num_image_px;
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

    // Replace values in the overscan region with 1, so overscan survives flatfielding
    if (use_overscan)
        for (int j = bias_rgn[2]; j < bias_rgn[3]; j++)
            for (int i = bias_rgn[0]; i < bias_rgn[1]; i++)
                median_flat[j*dark.cols + i] = 1;

    // Create a new fits file
    fitsfile *out;
    int status = 0;
    char *outbuf;
    asprintf(&outbuf, "!%s", outname);
    fits_create_file(&out, outbuf, &status);
    free(outbuf);

    // Create the primary array image (16-bit short integer pixels
    fits_create_img(out, DOUBLE_IMG, 2, (long []){dark.cols, dark.rows}, &status);

    // Set header keys for readout noise and gain
    if (use_overscan)
    {
        fits_update_key(out, TDOUBLE, "CCD-READ", &readnoise, "Estimated read noise (ADU)", &status);
        fits_update_key(out, TDOUBLE, "CCD-GAIN", &median_gain, "Estimated gain (electrons/ADU)", &status);

        printf("Readnoise: %f\n", readnoise);
        printf("Gain mean: %f median: %f\n", mean_gain, median_gain);
    }

    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, dark.rows*dark.cols, median_flat, &status))
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

    return return_status;
}

// Create a darkframe from the frames listed by the command `darkcmd',
// rejecting `minmax' highest and lowest pixel values.
// Save the resulting image to the file `outname'
int create_dark(const char *pattern, int minmax, const char *outname)
{
    char **frames;
    int numMatched = get_matching_files(pattern, &frames);
    if (numMatched < 2*minmax)
    {
        free_2d_array(frames, numMatched);
        return error("Insufficient frames. %d found, %d will be discarded", numMatched, 2*minmax);
    }
    
    framedata base = framedata_new(frames[0], FRAMEDATA_DBL);
    int exptime = framedata_get_header_int(&base, "EXPTIME");
    double *dark = (double *)malloc(base.rows*base.cols*sizeof(double));
    if (dark == NULL)
        return error("malloc failed");
    
    // Load the flat frames, discarding the 5 outermost pixels for each
    if (load_reject_minmax( (const char **)frames, numMatched, base.rows, base.cols, minmax, minmax, dark, NULL, NULL))
    {
        free(dark);
        framedata_free(base);
        free_2d_array(frames, numMatched);
        return error("frame loading failed");
    }
    
    // Create a new fits file
    fitsfile *out;
    int status = 0;
    char outbuf[2048];
    sprintf(outbuf, "!%s", outname);
    
    fits_create_file(&out, outbuf, &status);
    
    // Create the primary array image (16-bit short integer pixels
    fits_create_img(out, DOUBLE_IMG, 2, (long []){base.cols, base.rows}, &status);
    fits_update_key(out, TINT, "EXPTIME", &exptime, "Actual integration time (sec)", &status);
    
    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, base.rows*base.cols, dark, &status))
    {
        fits_close_file(out, &status);
        free(dark);
        framedata_free(base);
        free_2d_array(frames, numMatched);
        return error("fits_write_img failed with status %d", status);
    }
    
    fits_close_file(out, &status);
    free(dark);
    framedata_free(base);
    free_2d_array(frames, numMatched);
    
    return 0;
}

// Dark-subtract and flatfield the frame at `framePath' using the dark
// at `darkPath' and the flatfield at `flatPath'.
// Save the resulting image as a (double) floating-point fits file `outname'
int reduce_single_frame(char *framePath, char *darkPath, char *flatPath, char *outPath)
{
    framedata base = framedata_new(framePath, FRAMEDATA_DBL);
    framedata dark = framedata_new(darkPath, FRAMEDATA_DBL);
    framedata flat = framedata_new(flatPath, FRAMEDATA_DBL);
    
    // Subtract dark counts
    if (dark.dbl_data != NULL)
        framedata_subtract(&base, &dark);
    
    // Flat field image
    if (flat.dbl_data != NULL)
        framedata_divide(&base, &flat);
    
    // Create a new fits file
    fitsfile *out;
    int status = 0;
    char outbuf[2048];
    sprintf(outbuf, "!%s", outPath);
    fits_create_file(&out, outbuf, &status);
    
    // Create the primary array image
    fits_create_img(out, DOUBLE_IMG, 2, (long []){base.cols, base.rows}, &status);
    
    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, base.rows*base.cols, base.dbl_data, &status))
    {
        fits_close_file(out, &status);
        framedata_free(flat);
        framedata_free(dark);
        framedata_free(base);
        return error("fits_write_img failed with status %d", status);
    }
    
    fits_close_file(out, &status);
    framedata_free(flat);
    framedata_free(dark);
    framedata_free(base);
    return 0;
}


// Load the reduction file at dataPath and reduce any new data
int update_reduction(char *dataPath)
{
    // Processed dark and flat templates
    framedata dark;
    dark.data = NULL;
    dark.dbl_data = NULL;
    framedata flat;
    flat.dbl_data = NULL;
    
    // Read file header
    datafile data = read_data_header(dataPath);
    int start_obs = data.num_obs;

    if (data.file == NULL)
        return error("Error opening data file");
    
    if (data.version < 3)
        return error("Invalid data file version `%d'. Requires version >= 3", data.version);
    
    chdir(data.frame_dir);
    
    if (data.flat_template != NULL)
        flat = framedata_new(data.flat_template, FRAMEDATA_DBL);
    
    if (data.dark_template != NULL)
        dark = framedata_new(data.dark_template, FRAMEDATA_DBL);
    
    // Compile the filepattern into a regex
    regex_t regex;
    int regerr = 0;
    if ((regerr = regcomp(&regex, data.frame_pattern, REG_EXTENDED | REG_NOSUB)))
    {
        char errbuf[1024];
        regerror(regerr, &regex, errbuf, 1024);
        regfree(&regex);
        fclose(data.file);
        return error("Error compiling `%s` into a regular expression: %s", data.frame_pattern, errbuf);
    }
    
    // Iterate through the files in the directory
    struct dirent **matched;
    char filename[NAME_MAX];
    int numMatched = scandir(".", &matched, 0, alphasort);
    for (int i = 0; i < numMatched; i++)
    {
        strncpy(filename, matched[i]->d_name, NAME_MAX);
        filename[NAME_MAX-1] = '\0';
        free(matched[i]);
        
        // Ignore files that don't match the regex
        if (regexec(&regex, filename, 0, NULL, 0))
            continue;
        
        // Check whether the frame has been processed
        int processed = FALSE;
        for (int i = 0; i < data.num_obs; i++)
            if (strcmp(filename, data.obs[i].filename) == 0)
            {
                processed = TRUE;
                break;
            }

        if (processed)
            continue;
        
        printf("Reducing %s\n", filename);

        framedata frame = framedata_new(filename, FRAMEDATA_DBL);
        
        int exptime = framedata_get_header_int(&frame, "EXPTIME");
        
        // Calculate time at the *start* of the exposure relative to ReferenceTime
        double starttime = 0;
        
        char datebuf[128], timebuf[128], datetimebuf[257];
        struct tm t;
        if (framedata_has_header_string(&frame, "UTC-BEG"))
        {
            framedata_get_header_string(&frame, "UTC-DATE", datebuf);
            framedata_get_header_string(&frame, "UTC-BEG", timebuf);
            
            sprintf(datetimebuf, "%s %s", datebuf, timebuf);
        }
        else if (framedata_has_header_string(&frame, "GPSTIME"))
            framedata_get_header_string(&frame, "GPSTIME", datetimebuf);
        
        strptime(datetimebuf, "%Y-%m-%d %H:%M:%S", &t);
        time_t frame_time = timegm(&t);
        starttime = difftime(frame_time, data.reference_time);
        
        // Subtract dark counts
        if (dark.dbl_data != NULL)
            framedata_subtract(&frame, &dark);
        
        // Flat field image
        if (flat.dbl_data != NULL)
            framedata_divide(&frame, &flat);
        
        //
        // Process frame
        //
        
        // Observation start time
        fprintf(data.file, "%.1f ", starttime);
        data.obs[data.num_obs].time = starttime;
        
        // Target stars
        double comparisonIntensity = 0;
        double targetIntensity = 0;
        for (int i = 0; i < data.num_targets; i++)
        {
            // Use the aperture position from the previous frame
            // as a starting point if it is valid
            target t = data.targets[i];
            if (data.num_obs > 0)
            {
                double2 last = data.obs[data.num_obs-1].pos[i];
                if (last.x > 0 && last.x < frame.cols)
                    t.x = last.x;
                if (last.y > 0 && last.y < frame.rows)
                    t.y = last.y;
            }
            
            double2 xy = center_aperture(t, &frame);
            double sky = 0;
            double intensity = 0;
            
            if (xy.x < 0) // center_aperture returns negative on error
            {
                xy.x = 0;
                xy.y = 0;
            }
            else
            {
                double r = t.r;
                double2 bg = calculate_background(t, &frame);

                sky = bg.x*M_PI*r*r / exptime;
                intensity = integrate_aperture(xy, r, &frame) / exptime - sky;
            }

            fprintf(data.file, "%.2f ", intensity); // intensity (ADU/s)
            fprintf(data.file, "%.2f ", sky); // sky intensity (ADU/s)
            fprintf(data.file, "%.2f %.2f ", xy.x, xy.y); // Aperture center
            
            data.obs[data.num_obs].star[i] = intensity;
            data.obs[data.num_obs].sky[i] = sky;
            data.obs[data.num_obs].pos[i].x = xy.x;
            data.obs[data.num_obs].pos[i].y = xy.y;
            
            if (i == 0)
                targetIntensity = intensity;
            else
                comparisonIntensity += intensity;
        }
        
        // Ratio
        double ratio = comparisonIntensity > 0 ? targetIntensity / comparisonIntensity : 0;
        fprintf(data.file, "%.3e ", ratio);
        data.obs[data.num_obs].ratio = ratio;
        
        // Filename
        fprintf(data.file, "%s\n", filename);
        strncpy(data.obs[data.num_obs].filename, filename, sizeof(filename));
        data.num_obs++;
        
        framedata_free(frame);
    }
    if (dark.data != NULL)
        framedata_free(dark);
    
    free(matched);
    regfree(&regex);
    fclose(data.file);
    printf("Reduced %d observations\n", data.num_obs - start_obs);
    return 0;
}


// Create a reduction file at filePath, with frames from the dir framePath
// matching the regex framePattern, with dark and flat frames given by
// darkTemplate and flatTemplate
int create_reduction_file(char *framePath, char *framePattern, char *darkTemplate, char *flatTemplate, char *filePath)
{
    FILE *data = fopen(filePath, "wx");
    if (data == NULL)
        return error("Unable to create data file: %s. Does it already exist?", filePath);
    
    char pathBuf[PATH_MAX];
    realpath(framePath, pathBuf);
    
    if (chdir(pathBuf))
    {
        fclose(data);
        return error("Invalid frame path: %s", pathBuf);
    }
    
    if (!init_ds9("tsreduce"))
    {
        fclose(data);
        return error("Unable to launch ds9");
    }
    
    char filenamebuf[NAME_MAX];
    if (!get_first_matching_file(framePattern, filenamebuf, NAME_MAX))
    {
        fclose(data);
        return error("No matching files found");
    }
    
    // Open the file to find the reference time
    framedata frame = framedata_new(filenamebuf, FRAMEDATA_DBL);
    char datetimebuf[257];
    if (framedata_has_header_string(&frame, "UTC-BEG"))
    {
        char datebuf[128], timebuf[128];
        framedata_get_header_string(&frame, "UTC-DATE", datebuf);
        framedata_get_header_string(&frame, "UTC-BEG", timebuf);
        
        sprintf(datetimebuf, "%s %s", datebuf, timebuf);
    }
    else if (framedata_has_header_string(&frame, "GPSTIME"))
        framedata_get_header_string(&frame, "GPSTIME", datetimebuf);
    
    if (darkTemplate != NULL)
    {
        framedata dark = framedata_new(darkTemplate, FRAMEDATA_DBL);
        framedata_subtract(&frame, &dark);
        framedata_free(dark);
    }
    
    if (flatTemplate != NULL)
    {
        framedata flat = framedata_new(flatTemplate, FRAMEDATA_DBL);
        framedata_divide(&frame, &flat);
        framedata_free(flat);
    }
    
    char command[128];
    snprintf(command, 128, "array [xdim=%d,ydim=%d,bitpix=-64]", frame.cols, frame.rows);
    if (tell_ds9("tsreduce", command, frame.dbl_data, frame.rows*frame.cols*sizeof(double)))
    {
        framedata_free(frame);
        fclose(data);
        return error("ds9 command failed: %s", command);
    }
    
    // Set scaling mode
    if (tell_ds9("tsreduce", "scale mode 99.5", NULL, 0))
    {
        framedata_free(frame);
        fclose(data);
        return error("ds9 command failed: scale mode 99.5");
    }
    
    // Flip X axis
    if (tell_ds9("tsreduce", "orient x", NULL, 0))
    {
        framedata_free(frame);
        fclose(data);
        return error("ds9 command failed: orient x");
    }
    
    printf("Circle the target stars and surrounding sky in ds9 then press enter to continue...\n");
    getchar();
    
    char *ds9buf;
    if (ask_ds9("tsreduce", "regions", &ds9buf) || ds9buf == NULL)
    {
        framedata_free(frame);
        fclose(data);
        return error("ds9 request regions failed");
    }
    
    // Parse the region definitions
    target targets[MAX_TARGETS];
    int num_targets = 0;
    char *cur = ds9buf;
    while ((cur = strstr(cur, "circle")) != NULL)
    {
        if (num_targets == MAX_TARGETS)
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
        
        int n = 0;
        double move = 0;
        target last;
        // Converge on the best aperture size and position
        do
        {
            if (n++ == 20)
            {
                printf("WARNING: Aperture centering did not converge");
                break;
            }
            
            last = t;
            // Calculate a rough center and background: Estimates will improve as we converge
            double2 xy = center_aperture(t, &frame);
            if (xy.x < 0)
            {
                free(ds9buf);
                framedata_free(frame);
                fclose(data);
                return error("Aperture did not converge");
            }
            t.x = xy.x; t.y = xy.y;
            double2 bg = calculate_background(t, &frame);
            
            double lastIntensity = 0;
            double lastProfile = frame.dbl_data[frame.cols*((int)xy.y) + (int)xy.x];
            
            // Estimate the radius where the star flux falls to 5 times the std. dev. of the background
            int maxRadius = (int)t.s2 + 1;
            for (int radius = 1; radius <= maxRadius; radius++)
            {
                double intensity = integrate_aperture(xy, radius, &frame) - bg.x*M_PI*radius*radius;
                double profile = (intensity - lastIntensity) / (M_PI*(2*radius-1));
                if (profile < 5*bg.y)
                {
                    t.r = radius - 1 + (5*bg.y - lastProfile) / (profile - lastProfile);
                    break;
                }
                lastIntensity = intensity;
                lastProfile = profile;
            }
            
            printf("Iteration %d: Center: (%f, %f) Radius: %f Background: %f +/- %f\n", n, t.x, t.y, t.r, bg.x, bg.y);
            move = (xy.x-last.x)*(xy.x-last.x) + (xy.y-last.y)*(xy.y-last.y);
        } while (move >= 0.00390625);
        
        // Set target parameters
        targets[num_targets++] = t;
        
        // Increment by one char so strstr will find the next instance
        cur++;
    }
    free(ds9buf);
    
    // Set all target radii equal to the largest one
    // Temporary workaround (hopefully) until we can determine why different aperture sizes ratio badly
    double largest = 0;
    for (int i = 0; i < num_targets; i++)
        largest = fmax(largest, targets[i].r);
    for (int i = 0; i < num_targets; i++)
        targets[i].r = largest;
    
    printf("Founds %d targets\n", num_targets);
    
    // Write the file
    fprintf(data, "# Puoko-nui Online reduction output\n");
    fprintf(data, "# Version: 4\n");
    fprintf(data, "# FrameDir: %s\n", pathBuf);
    fprintf(data, "# FramePattern: %s\n", framePattern);
    fprintf(data, "# DarkTemplate: %s\n", darkTemplate);
    fprintf(data, "# FlatTemplate: %s\n", flatTemplate);
    fprintf(data, "# ReferenceTime: %s\n", datetimebuf);
    for (int i = 0; i < num_targets; i++)
        fprintf(data, "# Target: (%f, %f, %f, %f, %f) [1.0]\n", targets[i].x, targets[i].y, targets[i].r, targets[i].s1, targets[i].s2);
    
    // Display results in ds9
    snprintf(command, 128, "regions delete all");
    if (tell_ds9("tsreduce", command, NULL, 0))
    {
        framedata_free(frame);
        fclose(data);
        return error("ds9 command failed: %s", command);
    }
    
    for (int i = 0; i < num_targets; i++)
    {
        double x = targets[i].x + 1;
        double y = targets[i].y + 1;
        
        snprintf(command, 128, "regions command {circle %f %f %f #color=red}", x, y, targets[i].r);
        if (tell_ds9("tsreduce", command, NULL, 0))
            fprintf(stderr, "ds9 command failed: %s\n", command);
        snprintf(command, 128, "regions command {circle %f %f %f #background}", x, y, targets[i].s1);
        if (tell_ds9("tsreduce", command, NULL, 0))
            fprintf(stderr, "ds9 command failed: %s\n", command);
        snprintf(command, 128, "regions command {circle %f %f %f #background}", x, y, targets[i].s2);
        if (tell_ds9("tsreduce", command, NULL, 0))
            fprintf(stderr, "ds9 command failed: %s\n", command);
    }
    
    framedata_free(frame);
    fclose(data);
    return 0;
}
