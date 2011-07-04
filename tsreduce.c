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
#include <fitsio.h>

#include "tsreduce.h"
#include "framedata.h"
#include "helpers.h"
#include "aperture.h"


// Subtracts the dark count, then normalizes the frame to average to unity
// (not that we actually care about the total photometric counts)
void normalize_flat(framedata *flat, void *data)
{
    framedata *dark = (framedata *)data;
    
    if (dark->dtype != FRAMEDATA_DBL || flat->dtype != FRAMEDATA_DBL)
        error("normalize_flat frames must be type DBL");
    
    if (dark->rows != flat->rows || dark->cols != flat->rows)
        error("normalize_flat frames must have same size");
    
    int flatexp = framedata_get_header_int(flat, "EXPTIME");
    int darkexp = framedata_get_header_int(dark, "EXPTIME");
    int n = flat->rows*flat->cols;
    
    // Calculate mean
    double mean = 0;
    for (int i = 0; i < n; i++)
    {
        flat->dbl_data[i] -= flatexp*1.0/darkexp*dark->dbl_data[i];
        mean += flat->dbl_data[i];
    }
    mean /= n;
    
    // Calculate standard deviation
    double std = 0;
    for (int i = 0; i < n; i++)
        std += (flat->dbl_data[i] - mean)*(flat->dbl_data[i] - mean);
    std = sqrt(std/n);
    
    // Recalculate the mean, excluding outliers at 3 sigma
    double mean_new = 0;
    for (int i = 0; i < n; i++)
        if (fabs(flat->dbl_data[i] - mean) < 3*std)
            mean_new += flat->dbl_data[i];
    mean_new /= n;
    
    // Normalize flat to unity counts
    for (int i = 0; i < n; i++)
        flat->dbl_data[i] /= mean_new;
}

// Create a flat field frame from the frames listed by the command `flatcmd',
// rejecting `minmax' highest and lowest pixel values after subtracting
// the dark frame `masterdark'.
// Save the resulting image to the file `outname'
int create_flat(const char *flatcmd, int minmax, const char *masterdark, const char *outname)
{
    char **frames = (char **)malloc(MAX_FRAMES*sizeof(char*));
    for (int i = 0; i < MAX_FRAMES; i++)
        frames[i] = (char *)malloc(PATH_MAX*sizeof(char));
    
    int numflats = get_matching_files(flatcmd, frames, MAX_FRAMES);
    if (numflats < 2*minmax)
        return error("Insufficient frames. %d found, %d will be discarded", numflats, 2*minmax);
    
    double *flat = (double *)malloc(512*512*sizeof(double));
    
    framedata dark = framedata_new(masterdark, FRAMEDATA_DBL);
    
    // Load the flat frames, discarding the 5 outermost pixels for each
    load_reject_minmax( (const char **)frames, numflats, dark.rows, dark.cols, minmax, minmax, flat, normalize_flat, (void *)&dark);
    framedata_free(dark);
    
    // Create a new fits file
    fitsfile *out;
    int status = 0;
    char outbuf[2048];
    sprintf(outbuf, "!%s", outname);
    fits_create_file(&out, outbuf, &status);
    
    /* Create the primary array image (16-bit short integer pixels */
	long size[2] = { 512, 512 };
	fits_create_img(out, DOUBLE_IMG, 2, size, &status);
    
    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, 512*512, flat, &status))
        return error("fits_write_img failed with status %d", status);
    
    fits_close_file(out, &status);
    free(flat);
    
    for (int i = 0; i < MAX_FRAMES; i++)
        free(frames[i]);
    free(frames);
    
    return 0;
}

// Create a darkframe from the frames listed by the command `darkcmd',
// rejecting `minmax' highest and lowest pixel values.
// Save the resulting image to the file `outname'
int create_dark(const char *darkcmd, int minmax, const char *outname)
{
    char **frames = (char **)malloc(MAX_FRAMES*sizeof(char*));
    for (int i = 0; i < MAX_FRAMES; i++)
        frames[i] = (char *)malloc(PATH_MAX*sizeof(char));
    
    int numdarks = get_matching_files(darkcmd, frames, MAX_FRAMES);
    if (numdarks < 2*minmax)
        return error("Insufficient frames. %d found, at least %d are required", numdarks, 2*minmax);
    
    framedata base = framedata_new(frames[0], FRAMEDATA_INT);
    int exptime = framedata_get_header_int(&base, "EXPTIME");
    double *dark = (double *)malloc(base.rows*base.cols*sizeof(double));
    
    
    // Load the flat frames, discarding the 5 outermost pixels for each
    load_reject_minmax( (const char **)frames, numdarks, base.rows, base.cols, minmax, minmax, dark, NULL, NULL);
    
    // Create a new fits file
    fitsfile *out;
    int status = 0;
    char outbuf[2048];
    sprintf(outbuf, "!%s", outname);
    
    fits_create_file(&out, outbuf, &status);
    
    /* Create the primary array image (16-bit short integer pixels */
	long size[2] = { 512, 512 };
	fits_create_img(out, DOUBLE_IMG, 2, size, &status);
    
    fits_update_key(out, TINT, "EXPTIME", &exptime, "Actual integration time (sec)", &status);
    
    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, base.rows*base.cols, dark, &status))
        return error("fits_write_img failed with status %d", status);
    
    fits_close_file(out, &status);
    free(dark);
    
    framedata_free(base);
    
    for (int i = 0; i < MAX_FRAMES; i++)
        free(frames[i]);
    free(frames);
    
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
    
    /* Create the primary array image */
    long size[2] = { 512, 512 };
    fits_create_img(out, DOUBLE_IMG, 2, size, &status);
    
    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, base.rows*base.cols, base.dbl_data, &status))
        return error("fits_write_img failed with status %d", status);
    
    fits_close_file(out, &status);
    
    return 0;
}

// Read the header section from the data file
// defined as all the lines at the top of the file
// that start with #
datafile read_data_header(char *dataFile)
{
    datafile h;
    h.file = NULL;
    h.frame_dir[0] = '\0';
    h.frame_pattern[0] = '\0';
    h.dark_template[0] = '\0';
    h.flat_template[0] = '\0';
    
    // Open the data file (created with `tsreduce init`)
    h.file = fopen(dataFile, "r+");
    if (h.file == NULL)
        return h;

    char linebuf[1024];
    while (fgets(linebuf, sizeof(linebuf)-1, h.file) != NULL)
    {
        //
        // Header
        //
        if (!strncmp(linebuf,"# FramePattern:", 15))
            sscanf(linebuf, "# FramePattern: %128s\n", h.frame_pattern);

        else if (!strncmp(linebuf,"# FrameDir:", 11))
            sscanf(linebuf, "# FrameDir: %1024s\n", h.frame_dir);

        else if (!strncmp(linebuf,"# DarkTemplate:", 15))
            sscanf(linebuf, "# DarkTemplate: %128s\n", h.dark_template);

        else if (!strncmp(linebuf,"# FlatTemplate:", 15))
            sscanf(linebuf, "# FlatTemplate: %128s\n", h.flat_template);

        else if (!strncmp(linebuf,"# Target:", 9))
        {
            sscanf(linebuf, "# Target: (%lf, %lf, %lf, %lf, %lf)\n",
                   &h.targets[h.num_targets].x,
                   &h.targets[h.num_targets].y,
                   &h.targets[h.num_targets].r,
                   &h.targets[h.num_targets].s1,
                   &h.targets[h.num_targets].s2);
            h.num_targets++;
        }
        else if (!strncmp(linebuf,"# ReferenceTime:", 16))
        {
            struct tm t;
            strptime(linebuf, "# ReferenceTime: %Y-%m-%d %H:%M:%S\n", &t);
            h.reference_time = timegm(&t);
        }
        else if (!strncmp(linebuf,"# Version:", 10))
            sscanf(linebuf, "# Version: %d\n", &h.version);
        
        if (linebuf[0] == '#')
            continue;
        
        //
        // Observations
        //
        
        // Time
        char *ctx;
        h.obs[h.num_obs].time = atof(strtok_r(linebuf, " ", &ctx));
        
        // Target intensity / sky / aperture x / aperture y
        for (int i = 0; i < 3; i++)
        {
            h.obs[h.num_obs].star[i] = atof(strtok_r(NULL, " ", &ctx));
            h.obs[h.num_obs].sky[i] = atof(strtok_r(NULL, " ", &ctx));
            h.obs[h.num_obs].pos[i].x = atof(strtok_r(NULL, " ", &ctx));
            h.obs[h.num_obs].pos[i].y = atof(strtok_r(NULL, " ", &ctx));
        }
        
        // Ratio
        h.obs[h.num_obs].ratio = atof(strtok_r(NULL, " ", &ctx));
        
        // Filename
        strncpy(h.obs[h.num_obs].filename, strtok_r(NULL, " ", &ctx), sizeof(h.obs[h.num_obs].filename));
        
        // Strip newline
        h.obs[h.num_obs].filename[strlen(h.obs[h.num_obs].filename)-1] = '\0';
        
        h.num_obs++;

    }
    return h;
}

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
    
    if (data.file == NULL)
        return error("Error opening data file");

    if (data.version != 2)
        return error("Invalid data file version `%d'. Requires version `%d'", data.version, 2);

    chdir(data.frame_dir);
    
    if (data.flat_template != NULL)
        flat = framedata_new(data.flat_template, FRAMEDATA_DBL);
    
    if (data.dark_template != NULL)
        dark = framedata_new(data.dark_template, FRAMEDATA_DBL);

    // Iterate through the list of files matching the filepattern
    char filenamebuf[PATH_MAX+8];
    sprintf(filenamebuf, "/bin/ls %s", data.frame_pattern);
    FILE *ls = popen(filenamebuf, "r");
    if (ls == NULL)
        return error("failed to list directory");
    
    while (fgets(filenamebuf, sizeof(filenamebuf)-1, ls) != NULL)
    {
        // Strip the newline character from the end of the filename
        filenamebuf[strlen(filenamebuf)-1] = '\0';
        char *filename = filenamebuf;
        
        // Check whether the frame has been processed
        int processed = FALSE;
        for (int i = 0; i < data.num_obs; i++)
            if (strcmp(filename, data.obs[i].filename) == 0)
            {
                processed = TRUE;
                break;
            }
        printf("%s: processed %d\n", filename, processed);
        
        if (processed)
            continue;
        
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

        // Target stars
        double comparison = 0;
        double target = 0;
        for (int i = 0; i < data.num_targets; i++)
        {
            double2 xy = converge_aperture(data.targets[i], &frame);
            double sky = 0;
            double intensity = 0;

            if (xy.x > 0) // converge_aperture returns negative on error
            {
                double r = data.targets[i].r;
                double2 bg = calculate_background(data.targets[i], &frame);

                sky = bg.x*M_PI*r*r / exptime;
                intensity = integrate_aperture(xy, r, &frame) / exptime - sky;
            }

            fprintf(data.file, "%f ", intensity); // intensity (ADU/s)
            fprintf(data.file, "%f ", sky); // sky intensity (ADU/s)
            fprintf(data.file, "%f %f ", xy.x, xy.y); // Aperture center

            if (i == 0)
                target = intensity;
            else
                comparison += intensity;
        }

        // Ratio, Filename
        fprintf(data.file, "%f %s\n", target / comparison, filename);

        framedata_free(frame);
    }
    if (dark.data != NULL)
        framedata_free(dark);
    
    pclose(ls);
    fclose(data.file);
    
    return 0;
}

int main( int argc, char *argv[] )
{
    // `tsreduce create-flat "/bin/ls dome-*.fits.gz" 5 master-dark.fits.gz master-dome.fits.gz`
    if (argc == 6 && strncmp(argv[1], "create-flat", 11) == 0)
        return create_flat(argv[2], atoi(argv[3]), argv[4], argv[5]);
    
    // `tsreduce create-dark "/bin/ls dark-*.fits.gz" 5 master-dark.fits.gz`
    else if (argc == 5 && strncmp(argv[1], "create-dark", 11) == 0)
        return create_dark(argv[2], atoi(argv[3]), argv[4]);
    
    // `tsreduce reduce rawframe.fits.gz master-dark.fits.gz master-flat.fits.gz reduced.fits.gz`
    else if (argc == 6 && strncmp(argv[1], "reduce", 6) == 0)
        return reduce_single_frame(argv[2], argv[3], argv[4], argv[5]);
    
    // `tsreduce update ~/data/20110704/gwlib.dat`
    else if (argc == 3 && strncmp(argv[1], "update", 6) == 0)
        return update_reduction(argv[2]);
    
    else
        error("Invalid args");
    return 0;
}
