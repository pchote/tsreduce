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
#include <fitsio.h>
#include <xpa.h>
#include <cpgplot.h>

#include "tsreduce.h"
#include "framedata.h"
#include "helpers.h"
#include "aperture.h"

// Subtracts the dark count, then normalizes the frame to average to unity
// (not that we actually care about the total photometric counts)
int normalize_flat(framedata *flat, void *data)
{
    framedata *dark = (framedata *)data;

    if (dark->dtype != FRAMEDATA_DBL || flat->dtype != FRAMEDATA_DBL)
        return error("normalize_flat frames must be type DBL");

    if (dark->rows != flat->rows || dark->cols != flat->rows)
        return error("normalize_flat frames must have same size");

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

    return 0;
}

// Create a flat field frame from the frames listed by the command `flatcmd',
// rejecting `minmax' highest and lowest pixel values after subtracting
// the dark frame `masterdark'.
// Save the resulting image to the file `outname'
int create_flat(const char *pattern, int minmax, const char *masterdark, const char *outname)
{
    char **frames;
    int numMatched = get_matching_files(pattern, &frames);

    if (numMatched <= 2*minmax)
    {
        free_2d_array(frames, numMatched);
        return error("Insufficient frames. %d found, %d will be discarded", numMatched, 2*minmax);
    }

    framedata dark = framedata_new(masterdark, FRAMEDATA_DBL);
    double *flat = (double *)malloc(dark.rows*dark.cols*sizeof(double));

    // Load the flat frames, discarding the 5 outermost pixels for each
    if (load_reject_minmax( (const char **)frames, numMatched, dark.rows, dark.cols, minmax, minmax, flat, normalize_flat, (void *)&dark))
    {
        framedata_free(dark);
        free(flat);
        free_2d_array(frames, numMatched);
        return error("frame loading failed");
    }
    framedata_free(dark);

    // Create a new fits file
    fitsfile *out;
    int status = 0;
    char outbuf[2048];
    sprintf(outbuf, "!%s", outname);
    fits_create_file(&out, outbuf, &status);

    /* Create the primary array image (16-bit short integer pixels */
	long size[2] = { dark.rows, dark.cols };
	fits_create_img(out, DOUBLE_IMG, 2, size, &status);

    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, dark.rows*dark.cols, flat, &status))
    {
        fits_close_file(out, &status);
        free(flat);
        free_2d_array(frames, numMatched);
        return error("fits_write_img failed with status %d", status);
    }

    fits_close_file(out, &status);
    free(flat);
    free_2d_array(frames, numMatched);
    return 0;
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

    framedata base = framedata_new(frames[0], FRAMEDATA_INT);
    int exptime = framedata_get_header_int(&base, "EXPTIME");
    double *dark = (double *)malloc(base.rows*base.cols*sizeof(double));

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
	long size[2] = { base.rows, base.cols };
	fits_create_img(out, DOUBLE_IMG, 2, size, &status);

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
    long size[2] = { base.rows, base.cols };
    fits_create_img(out, DOUBLE_IMG, 2, size, &status);

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

        // Skip header / comment lines
        if (linebuf[0] == '#')
            continue;

        // Skip empty lines
        if (linebuf[0] == '\n')
            continue;

        //
        // Observations
        //

        // Time
        char *ctx;
        h.obs[h.num_obs].time = atof(strtok_r(linebuf, " ", &ctx));

        // Target intensity / sky / aperture x / aperture y
        for (int i = 0; i < h.num_targets; i++)
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

    if (data.file == NULL)
        return error("Error opening data file");

    if (data.version != 3)
        return error("Invalid data file version `%d'. Requires version `%d'", data.version, 3);

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

            printf("using aperture (%f,%f) orig (%f,%f)\n", t.x, t.y, data.targets[i].x, data.targets[i].y);
            double2 xy = converge_aperture(t, &frame);
            double sky = 0;
            double intensity = 0;

            if (xy.x > 0) // converge_aperture returns negative or nan
            {
                double r = t.r;
                double2 bg = calculate_background(t, &frame);

                sky = bg.x*M_PI*r*r / exptime;
                intensity = integrate_aperture(xy, r, &frame) / exptime - sky;
            }
            else
            {
                xy.x = 0;
                xy.y = 0;
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

    return 0;
}

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
    snprintf(command, 128, "array [xdim=%d,ydim=%d,bitpix=-64]", frame.rows, frame.cols);
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
            double2 xy = converge_aperture(t, &frame);
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
    fprintf(data, "# Version: 3\n");
    fprintf(data, "# FrameDir: %s\n", pathBuf);
    fprintf(data, "# FramePattern: %s\n", framePattern);
    fprintf(data, "# DarkTemplate: %s\n", darkTemplate);
    fprintf(data, "# FlatTemplate: %s\n", flatTemplate);
    fprintf(data, "# ReferenceTime: %s\n", datetimebuf);
    for (int i = 0; i < num_targets; i++)
        fprintf(data, "# Target: (%f, %f, %f, %f, %f)\n", targets[i].x, targets[i].y, targets[i].r, targets[i].s1, targets[i].s2);

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


// Takes a ratio value, normalizes it by the pre-calculated
//  polynomial fit, and returns as a percentage change
float normalize_ratio(int i, float d)
{
    float fit = 9.21E-11*i*i*i - 1.90E-07*i*i + 1.47E-04*i + 3.41E-01;
    return 100*(d - fit)/d;
}

int plot_fits(char *dataPath)
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

    // Load raw lightcurve and ratio data
    float time[MAX_OBS];
    float sky[MAX_OBS];
    float ratio[MAX_OBS];
    float raw1[MAX_OBS], raw2[MAX_OBS], raw3[MAX_OBS];

    float dftinput[MAX_OBS];
    float dft[500];
    float freq[500];

    double mean = 0;
    for (int i = 0; i < data.num_obs; i++)
    {
        time[i] = data.obs[i].time;
        ratio[i] = normalize_ratio(i, data.obs[i].ratio);
        sky[i] = data.obs[i].sky[0];

        raw1[i] = data.obs[i].star[0];
        raw2[i] = data.obs[i].star[1];
        raw3[i] = data.obs[i].star[2];

        dftinput[i] = ratio[i];
        mean += ratio[i];
    }
    mean /= data.num_obs;

    // Calculate standard deviation
    double std = 0;
    for (int i = 0; i < data.num_obs; i++)
        std += (dftinput[i] - mean)*(dftinput[i] - mean);
    std = sqrt(std/data.num_obs);

    double newmean = 0;
    int meancount = 0;

    // Discard outliers and recalculate mean
    for (int i = 0; i < data.num_obs; i++)
    {
        if (fabs(dftinput[i]-mean) > 3*std)
            dftinput[i] = 0;
        else
        {
            newmean += dftinput[i];
            meancount++;
        }
    }
    newmean /= meancount;

    // Calculate mmi intensities
    for (int i = 0; i < data.num_obs; i++)
        dftinput[i] = 1e-3*(dftinput[i] - newmean)/newmean;

    // Ratio limits
    float ratioMin = -12;
    float ratioMax = 13;

    float rawMin = 0.1;
    float rawMax = 3999;


    if (cpgopen("7/xs") <= 0)
        return error("Unable to open PGPLOT window");

    // 800 x 480
    cpgpap(9.41, 0.6);
    cpgask(0);
    cpgslw(3);
    cpgsfs(2);
    cpgscf(2);

    // Ratio
    cpgsvp(0.1, 0.9, 0.65, 0.92);
    cpgmtxt("t", 2, 0.5, 0.5, "UTC Hour");
    cpgmtxt("l", 2.5, 0.5, 0.5, "% Change");
    cpgswin(0, time[data.num_obs-1], ratioMin, ratioMax);
    cpgline(data.num_obs, time, ratio);

    // Rescale time axis to hours
    cpgswin(13.275, 13.275+time[data.num_obs-1]/3600, ratioMin, ratioMax);
    cpgbox("bcstm", 1, 4, "bcstn", 10, 2);

    // Raw Data
    cpgsvp(0.1, 0.9, 0.3, 0.65);
    cpgmtxt("l", 2.5, 0.5, 0.5, "Counts Per Second");
    cpgswin(0,time[data.num_obs-1], rawMin, rawMax);
    cpgsci(4);
    cpgline(data.num_obs, time, raw1);
    cpgsci(2);
    cpgline(data.num_obs, time, raw2);
    cpgsci(8);
    cpgline(data.num_obs, time, raw3);
    cpgsci(1);

    // Rescale time axis to hours
    cpgswin(13.275, 13.275+time[data.num_obs-1]/3600, rawMin, rawMax);
    cpgbox("bcst", 1, 4, "bcstn", 1000, 2);

    // DFT
    cpgsvp(0.1, 0.9, 0.075, 0.3);

    float fmax = 0.006;
    double pmax = 39.9e-3;
    calculate_amplitude_spectrum(0, fmax, time, dftinput, data.num_obs, freq, dft, 500);
    cpgswin(0, fmax, 0, pmax);
    cpgsci(12);
    cpgline(500, freq, dft);
    cpgsci(1);

    cpgswin(0, fmax/1e-6, 0, pmax*1000);
    cpgbox("bcstn", 0, 0, "bcnst", 10, 2);
    cpgmtxt("b", 2.5, 0.5, 0.5, "Frequency (\\gmHz)");
    cpgmtxt("l", 2, 0.5, 0.5, "Amplitude (mma)");

    cpgend();

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

    // `tsreduce display-targets ~/data/20110704/gwlib.dat`
    else if ((argc == 3 || argc == 4) && strncmp(argv[1], "display-targets", 15) == 0)
    {
        int obs = argc == 4 ? atoi(argv[3]) : -1;
        return display_targets(argv[2], obs);
    }

    // `tsreduce create ~/data/20110307 ec04207-[0-9]+.fits.gz master-dark.fits master-skyflat-20110305.fits 20110307.dat`
    else if (argc == 7 && strncmp(argv[1], "create", 6) == 0)
        return create_reduction_file(argv[2], argv[3], argv[4], argv[5], argv[6]);

    else if (argc == 5 && strncmp(argv[1], "profile", 7) == 0)
        return calculate_profile(argv[2], atoi(argv[3]), atoi(argv[4]));

    else if (argc == 3 && strncmp(argv[1], "repeats", 7) == 0)
        return detect_repeats(argv[2]);

    else if (argc == 3 && strncmp(argv[1], "plot", 4) == 0)
        plot_fits(argv[2]);
    else
        error("Invalid args");
    return 0;
}
