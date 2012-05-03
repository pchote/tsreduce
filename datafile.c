/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */


#include <sys/time.h>
#include <time.h>
#include <string.h>

#include "datafile.h"

#define PLOT_FIT_DEGREE_DEFAULT 4
#define PLOT_MAX_RAW_DEFAULT 5000
#define PLOT_NUM_UHZ_DEFAULT 1000
#define PLOT_MIN_UHZ_DEFAULT 0
#define PLOT_MAX_UHZ_DEFAULT 10000

/*
 * Allocate a datafile on the heap and set default values
 */
datafile *datafile_alloc()
{
    datafile *dp = malloc(sizeof(datafile));
    dp->file = NULL;
    dp->version = 0;
    dp->frame_dir = NULL;
    dp->frame_pattern = NULL;
    dp->dark_template = NULL;
    dp->flat_template = NULL;
    dp->plot_fit_degree = PLOT_FIT_DEGREE_DEFAULT;
    dp->plot_max_raw = PLOT_MAX_RAW_DEFAULT;
    dp->plot_num_uhz = PLOT_NUM_UHZ_DEFAULT;
    dp->plot_min_uhz = PLOT_MIN_UHZ_DEFAULT;
    dp->plot_max_uhz = PLOT_MAX_UHZ_DEFAULT;
    dp->num_obs = 0;
    dp->num_targets = 0;
    dp->num_blocked_ranges = 0;
    return dp;
}

/*
 * Allocate a new datafile and load with data from the specified file
 */
datafile* datafile_load(char *filename)
{
    datafile *dp = datafile_alloc();

    // Open the data file (created with `tsreduce init`)
    dp->file = fopen(filename, "r+");
    if (dp->file == NULL)
        return NULL;

    char linebuf[1024];
    char stringbuf[1024];
    while (fgets(linebuf, sizeof(linebuf)-1, dp->file) != NULL)
    {
        //
        // Header
        //
        if (!strncmp(linebuf,"# FramePattern:", 15))
        {
            sscanf(linebuf, "# FramePattern: %1024s\n", stringbuf);
            dp->frame_pattern = strndup(stringbuf, 1024);
        }
        else if (!strncmp(linebuf,"# FrameDir:", 11))
        {
            sscanf(linebuf, "# FrameDir: %1024s\n", stringbuf);
            dp->frame_dir = strndup(stringbuf, 1024);
        }
        else if (!strncmp(linebuf,"# DarkTemplate:", 15))
        {
            sscanf(linebuf, "# DarkTemplate: %1024s\n", stringbuf);
            dp->dark_template = strndup(stringbuf, 1024);
        }
        else if (!strncmp(linebuf,"# FlatTemplate:", 15))
        {
            sscanf(linebuf, "# FlatTemplate: %1024s\n", stringbuf);
            dp->flat_template = strndup(stringbuf, 1024);
        }
        else if (!strncmp(linebuf,"# Target:", 9))
        {
            if (dp->version == 3)
            {
                sscanf(linebuf, "# Target: (%lf, %lf, %lf, %lf, %lf)\n",
                       &dp->targets[dp->num_targets].x,
                       &dp->targets[dp->num_targets].y,
                       &dp->targets[dp->num_targets].r,
                       &dp->targets[dp->num_targets].s1,
                       &dp->targets[dp->num_targets].s2);
                dp->targets[dp->num_targets].plot_scale = 1;
                dp->num_targets++;
            }
            else
            {
                sscanf(linebuf, "# Target: (%lf, %lf, %lf, %lf, %lf) [%lf]\n",
                       &dp->targets[dp->num_targets].x,
                       &dp->targets[dp->num_targets].y,
                       &dp->targets[dp->num_targets].r,
                       &dp->targets[dp->num_targets].s1,
                       &dp->targets[dp->num_targets].s2,
                       &dp->targets[dp->num_targets].plot_scale);
                dp->num_targets++;
            }
        }
        else if (!strncmp(linebuf,"# ReferenceTime:", 16))
        {
            struct tm t;
            strptime(linebuf, "# ReferenceTime: %Y-%m-%d %H:%M:%S\n", &t);
            dp->reference_time = timegm(&t);
        }
        else if (!strncmp(linebuf,"# Version:", 10))
            sscanf(linebuf, "# Version: %d\n", &dp->version);
        else if (!strncmp(linebuf,"# PlotFitDegree:", 16))
            sscanf(linebuf, "# PlotFitDegree: %d\n", &dp->plot_fit_degree);
        else if (!strncmp(linebuf,"# PlotMaxRaw:", 13))
            sscanf(linebuf, "# PlotMaxRaw: %lf\n", &dp->plot_max_raw);
        else if (!strncmp(linebuf,"# PlotMinUhz:", 13))
            sscanf(linebuf, "# PlotMinUhz: %lf\n", &dp->plot_min_uhz);
        else if (!strncmp(linebuf,"# PlotMaxUhz:", 13))
            sscanf(linebuf, "# PlotMaxUhz: %lf\n", &dp->plot_max_uhz);
        else if (!strncmp(linebuf,"# PlotNumUhz:", 13))
            sscanf(linebuf, "# PlotNumUhz: %d\n", &dp->plot_num_uhz);
        else if (!strncmp(linebuf,"# BlockRange:", 13))
        {
            sscanf(linebuf, "# BlockRange: (%lf, %lf)\n",
                   &dp->blocked_ranges[dp->num_blocked_ranges].x,
                   &dp->blocked_ranges[dp->num_blocked_ranges].y);
            dp->num_blocked_ranges++;
        }

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
        dp->obs[dp->num_obs].time = atof(strtok_r(linebuf, " ", &ctx));

        // Target intensity / sky / aperture x / aperture y
        for (int i = 0; i < dp->num_targets; i++)
        {
            dp->obs[dp->num_obs].star[i] = atof(strtok_r(NULL, " ", &ctx));
            dp->obs[dp->num_obs].sky[i] = atof(strtok_r(NULL, " ", &ctx));
            dp->obs[dp->num_obs].pos[i].x = atof(strtok_r(NULL, " ", &ctx));
            dp->obs[dp->num_obs].pos[i].y = atof(strtok_r(NULL, " ", &ctx));
        }

        // Ratio
        dp->obs[dp->num_obs].ratio = atof(strtok_r(NULL, " ", &ctx));

        // Filename
        strncpy(dp->obs[dp->num_obs].filename, strtok_r(NULL, " ", &ctx),
                sizeof(dp->obs[dp->num_obs].filename));

        // Strip newline
        dp->obs[dp->num_obs].filename[strlen(dp->obs[dp->num_obs].filename)-1] = '\0';

        dp->num_obs++;

    }
    return dp;
}

void datafile_free(datafile *data)
{
    if (data->file)
        fclose(data->file);
    if (data->frame_dir)
        free(data->frame_dir);
    if (data->frame_pattern)
        free(data->frame_pattern);
    if (data->dark_template)
        free(data->dark_template);
    if (data->flat_template)
        free(data->flat_template);

    free(data);
}
}
