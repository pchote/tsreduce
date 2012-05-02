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
#include "fit.h"

static time_t get_frame_time(framedata *frame)
{
    char datebuf[128], timebuf[128], datetimebuf[257];
    struct tm t;
    if (framedata_has_header_string(frame, "UTC-BEG"))
    {
        framedata_get_header_string(frame, "UTC-DATE", datebuf);
        framedata_get_header_string(frame, "UTC-BEG", timebuf);
        sprintf(datetimebuf, "%s %s", datebuf, timebuf);
    }
    else if (framedata_has_header_string(frame, "GPSTIME"))
        framedata_get_header_string(frame, "GPSTIME", datetimebuf);

    strptime(datetimebuf, "%Y-%m-%d %H:%M:%S", &t);
    return timegm(&t);
}

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

    if (!num_frames)
        return error("No matching frames");

    // Load first frame to get frame info
    framedata frame = framedata_new(frames[0], FRAMEDATA_DBL);
    int *br = frame.regions.bias_region;
    double *data = (double *)malloc(frame.regions.bias_px*num_frames*sizeof(double));
    framedata_free(frame);
    int p = 0;

    // Load overscan region into a data cube
    for (int k = 0; k < num_frames; k++)
    {
        frame = framedata_new(frames[k], FRAMEDATA_DBL);

        // Copy overscan pixels
        for (int j = br[2]; j < br[3]; j++)
            for (int i = br[0]; i < br[1]; i++)
                data[p++] = frame.dbl_data[j*frame.cols + i];

        framedata_free(frame);
    }

    double readnoise = calculate_cube_stddev(data, num_frames, frame.regions.bias_px);
    printf("Read noise from %d frames: %f ADU\n", num_frames, readnoise);

    free(data);
    free_2d_array(frames, num_frames);

    return 0;
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
        flat->dbl_data[i] -= exp_ratio*dark->dbl_data[i];

    // Calculate mean in image area only
    int *ir = flat->regions.image_region;
    double mean = mean_in_region(flat, ir);

    // Calculate standard deviation
    double std = 0;
    for (int j = ir[2]; j < ir[3]; j++)
        for (int i = ir[0]; i < ir[1]; i++)
        {
            double temp = flat->dbl_data[flat->cols*j + i] - mean;
            std += temp*temp;
        }
    std = sqrt(std/flat->regions.image_px);

    // Recalculate the mean, excluding outliers at 3 sigma
    double mean_new = 0;
    int count = 0;
    for (int j = ir[2]; j < ir[3]; j++)
        for (int i = ir[0]; i < ir[1]; i++)
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
    if (dark.regions.has_overscan)
    {
        int *br = dark.regions.bias_region;
        double var = 0;
        for (int j = br[2]; j < br[3]; j++)
            for (int i = br[0]; i < br[1]; i++)
                for (int k = 0; k < num_frames; k++)
                {
                    // Need to multiply by the mean level for the frame to give the variance in ADU
                    double temp = mean_flat[k]*data_cube[num_frames*(dark.cols*j + i) + k];
                    var += temp*temp;
                }
        readnoise = sqrt(var/(dark.regions.bias_px*num_frames));
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

    if (dark.regions.has_overscan)
    {
        // Calculate mean dark level for gain calculation
        int *ir = dark.regions.image_region;
        double mean_dark = mean_in_region(&dark, ir);

        for(int k = 0; k < num_frames; k++)
        {
            // Calculate the variance by taking the difference between
            // the (normalized) frame and the (normalized) master-flat
            double var = 0;
            for (int j = ir[2]; j < ir[3]; j++)
                for (int i = ir[0]; i < ir[1]; i++)
                {
                    double temp = mean_flat[k]*(frames[k].dbl_data[dark.cols*j + i] - median_flat[dark.cols*j + i]);
                    var += temp*temp;
                }
            var /= dark.regions.image_px;
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
    if (dark.regions.has_overscan)
    {
        int *br = dark.regions.image_region;
        for (int k = 0; k < dark.rows*dark.cols; k++)
        {
            int x = k % dark.cols;
            int y = k / dark.cols;
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
    fits_create_img(out, DOUBLE_IMG, 2, (long []){dark.cols, dark.rows}, &status);

    // Set header keys for readout noise and gain
    if (dark.regions.has_overscan)
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
    int return_status = 0;

    char **frame_paths;
    int num_frames = get_matching_files(pattern, &frame_paths);
    if (num_frames < 2*minmax)
    {
        return_status = error("Insufficient frames. %d found, %d will be discarded", num_frames, 2*minmax);
        goto insufficient_frames;
    }
    
    framedata base = framedata_new(frame_paths[0], FRAMEDATA_DBL);
    int exptime = framedata_get_header_int(&base, "EXPTIME");
    double *median_dark = (double *)malloc(base.rows*base.cols*sizeof(double));
    if (median_dark == NULL)
    {
        return_status = error("malloc failed");
        goto dark_failed;
    }
    
    // Data cube for processing the flat data
    //        data[0] = frame[0][0,0], data[1] = frame[1][0,0] ... data[num_frames] = frame[0][1,0] etc
    double *data_cube = (double *)malloc(num_frames*base.cols*base.rows*sizeof(double));
    if (data_cube == NULL)
    {
        return_status = error("malloc failed");
        goto datacube_failed;
    }

    for( int i = 0; i < num_frames; i++)
    {
        printf("loading `%s`\n", frame_paths[i]);
        framedata f = framedata_new(frame_paths[i], FRAMEDATA_DBL);

        if (f.rows != base.rows || f.cols != base.cols)
        {
            framedata_free(f);
            return_status = error("Frame %s dimensions mismatch. Expected (%d,%d), was (%d, %d)", frame_paths[i], base.rows, base.cols, f.rows, f.cols);
            goto loaddark_failed;
        }

        subtract_bias(&f);
        for (int j = 0; j < base.rows*base.cols; j++)
            data_cube[num_frames*j+i] = f.dbl_data[j];

        framedata_free(f);
    }
    
    // Loop over the pixels, sorting the values from each image into increasing order
    for (int j = 0; j < base.rows*base.cols; j++)
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
    fits_create_img(out, DOUBLE_IMG, 2, (long []){base.cols, base.rows}, &status);
    fits_update_key(out, TINT, "EXPTIME", &exptime, "Actual integration time (sec)", &status);
    
    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, base.rows*base.cols, median_dark, &status))
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

    return return_status;
}

// Dark-subtract and flatfield the frame at `framePath' using the dark
// at `darkPath' and the flatfield at `flatPath'.
// Save the resulting image as a (double) floating-point fits file `outname'
int reduce_single_frame(char *framePath, char *darkPath, char *flatPath, char *outPath)
{
    framedata base = framedata_new(framePath, FRAMEDATA_DBL);
    framedata dark = framedata_new(darkPath, FRAMEDATA_DBL);
    framedata flat = framedata_new(flatPath, FRAMEDATA_DBL);
    
    // Subtract bias
    subtract_bias(&base);

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
        
        // Calculate time at the start of the exposure relative to ReferenceTime
        time_t frame_time = get_frame_time(&frame);
        double starttime = difftime(frame_time, data.reference_time);
        
        // Subtract bias
        subtract_bias(&frame);

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
            
            double sky = 0;
            double intensity = 0;
            double2 xy = {0,0};

            if (!center_aperture(t, &frame, &xy))
            {
                if (calculate_background(t, &frame, &sky, NULL))
                    sky = 0;

                // Integrate sky over the aperture and normalize per unit time
                sky *= M_PI*t.r*t.r / exptime;
                intensity = integrate_aperture(xy, t.r, &frame) / exptime - sky;
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
    subtract_bias(&frame);
    time_t frame_time = get_frame_time(&frame);

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
            double2 xy;
            if (center_aperture(t, &frame, &xy))
            {
                free(ds9buf);
                framedata_free(frame);
                fclose(data);
                return error("Aperture did not converge");
            }

            t.x = xy.x; t.y = xy.y;
            double sky_intensity, sky_std_dev;
            if (calculate_background(t, &frame, &sky_intensity, &sky_std_dev))
            {
                free(ds9buf);
                framedata_free(frame);
                fclose(data);
                return error("Background calculation failed");
            }

            double lastIntensity = 0;
            double lastProfile = frame.dbl_data[frame.cols*((int)xy.y) + (int)xy.x];
            
            // Estimate the radius where the star flux falls to 5 times the std. dev. of the background
            int maxRadius = (int)t.s2 + 1;
            for (int radius = 1; radius <= maxRadius; radius++)
            {
                double intensity = integrate_aperture(xy, radius, &frame) - sky_intensity*M_PI*radius*radius;
                double profile = (intensity - lastIntensity) / (M_PI*(2*radius-1));
                if (profile < 5*sky_std_dev)
                {
                    t.r = radius - 1 + (5*sky_std_dev - lastProfile) / (profile - lastProfile);
                    break;
                }
                lastIntensity = intensity;
                lastProfile = profile;
            }
            
            printf("Iteration %d: Center: (%f, %f) Radius: %f Background: %f +/- %f\n", n, t.x, t.y, t.r, sky_intensity, sky_std_dev);
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

    char datetimebuf[20];
    strftime(datetimebuf, 20, "%F %H %M %S\n", gmtime(&frame_time));
    fprintf(data, "# ReferenceTime: %s", datetimebuf);
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

int create_mmi(char *dataPath)
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
    double *time = malloc(data.num_obs*sizeof(double));
    if (time == NULL)
        return error("malloc failed");

    double *ratio = malloc(data.num_obs*sizeof(double));
    if (ratio == NULL)
        return error("malloc failed");

    double *polyfit = malloc(data.num_obs*sizeof(double));
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
        double target = data.obs[i].star[0];
        double comparison = 0;
        for (int j = 1; j < data.num_targets; j++)
            comparison += data.obs[i].star[j];

        ratio[i] = target/comparison;
        ratio_mean += ratio[i];
    }

    ratio_mean /= data.num_obs;

    // Calculate standard deviation
    double ratio_std = 0;
    for (int i = 0; i < data.num_obs; i++)
        ratio_std += (ratio[i] - ratio_mean)*(ratio[i] - ratio_mean);
    ratio_std = sqrt(ratio_std/data.num_obs);

    if (fit_polynomial_d(time, ratio, data.num_obs, coeffs, data.plot_fit_degree))
    {
        free(coeffs);
        free(ratio);
        free(time);
        fclose(data.file);
        return error("Fit failed");
    }

    double *mmi = malloc(data.num_obs*sizeof(double));
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

    printf("# points = %d; dt = %f\n", data.num_obs, (time[data.num_obs-1] - time[0])/3600.0);
    printf("# tgt = %s\n", data.frame_pattern);
    char buf[25];
    strftime(buf, 25, "# UT start = %H %M %S\n", gmtime(&data.reference_time));
    printf("%s", buf);
    printf("#\n");

    for (int i = 0; i < data.num_obs; i++)
        printf("%f %f\n", time[i]/3600.0, mmi[i]);

    free(coeffs);
    free(mmi);
    free(ratio);
    free(time);
    fclose(data.file);

    return 0;
}

//
// Iterate a range of aperture sizes over a data set to find the best signal-to-noise ratio
//
int evaluate_aperture_snr(char *dataPath, double minAperture, double maxAperture, int numApertures)
{
    // Read file header
    datafile data = read_data_header(dataPath);

    if (data.file == NULL)
        return error("Error opening data file");

    if (data.version < 3)
        return error("Invalid data file version `%d'. Requires version >= 3", data.version);

    chdir(data.frame_dir);

    framedata flat = framedata_new(data.flat_template, FRAMEDATA_DBL);
    framedata dark = framedata_new(data.dark_template, FRAMEDATA_DBL);
    double readnoise = framedata_get_header_dbl(&flat, "CCD-READ");
    double gain = framedata_get_header_dbl(&flat, "CCD-GAIN");

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

    struct dirent **matched;
    char filename[NAME_MAX];
    int numMatched = scandir(".", &matched, 0, alphasort);

    double *radii = malloc(numApertures*sizeof(double));
    double *snr = malloc(numApertures*sizeof(double));
    for (int k = 0; k < numApertures; k++)
    {
        radii[k] = minAperture + k*(maxAperture - minAperture)/(numApertures - 1);
        for (int i = 0; i < data.num_targets; i++)
            data.targets[i].r = radii[k];

        double totalSignal = 0;
        double totalNoise = 0;
        for (int j = 0; j < numMatched; j++)
        {
            strncpy(filename, matched[j]->d_name, NAME_MAX);
            filename[NAME_MAX-1] = '\0';

            // Ignore files that don't match the regex
            if (regexec(&regex, filename, 0, NULL, 0))
                continue;

            framedata frame = framedata_new(filename, FRAMEDATA_DBL);
            int exptime = framedata_get_header_int(&frame, "EXPTIME");

            // Calculate time at the start of the exposure relative to ReferenceTime
            double starttime = difftime(get_frame_time(&frame), data.reference_time);

            bool skip = false;
            for (int i = 0; i < data.num_blocked_ranges; i++)
                if (starttime >= data.blocked_ranges[i].x && starttime <= data.blocked_ranges[i].y)
                {
                    printf("Ignoring time %f\n", starttime);
                    skip = true;
                    break;
                }

            if (!skip)
            {
                // Preprocess frame
                subtract_bias(&frame);
                framedata_subtract(&frame, &dark);
                framedata_divide(&frame, &flat);

                // Process frame
                double comparisonIntensity = 0;
                double targetIntensity = 0;
                double comparisonNoise = 0;
                double targetNoise = 0;
                for (int i = 0; i < data.num_targets; i++)
                {
                    // Use the aperture position from the previous frame
                    // as a starting point if it is valid
                    target t = data.targets[i];
                    if (j > 0)
                    {
                        double2 last = data.obs[j-1].pos[i];
                        if (last.x > 0 && last.x < frame.cols)
                            t.x = last.x;
                        if (last.y > 0 && last.y < frame.rows)
                            t.y = last.y;
                    }

                    double sky = 0;
                    double intensity = 0;
                    double noise = 0;
                    double2 xy = {0,0};

                    if (!center_aperture(t, &frame, &xy))
                    {
                        if (calculate_background(t, &frame, &sky, NULL))
                            sky = 0;

                        // Integrate sky over the aperture and normalize per unit time
                        sky *= M_PI*t.r*t.r / exptime;

                        integrate_aperture_and_noise(xy, t.r, &frame, &dark, readnoise, gain, &intensity, &noise);
                        intensity = intensity/exptime - sky;
                        noise /= exptime;
                    }

                    data.obs[j].pos[i] = xy;

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
                double ratioNoise = (targetNoise/targetIntensity + comparisonNoise/comparisonIntensity)*ratio;
                data.num_obs++;
                totalSignal += ratio;
                totalNoise += ratioNoise;
            }
            framedata_free(frame);
        }

        snr[k] = totalSignal/totalNoise;
        printf("Radius: %.2f SNR: %f\n", radii[k], snr[k]);

    }

    framedata_free(dark);
    for (int i = 0; i < numMatched; i++)
        free(matched[i]);
    free(matched);
    regfree(&regex);
    fclose(data.file);

    printf("# Radius SNR\n");
    for (int i = 0; i < numApertures; i++)
        printf("%.2f %f\n", radii[i], snr[i]);
    free(radii);
    free(snr);
    return 0;
}
