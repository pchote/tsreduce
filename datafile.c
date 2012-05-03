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

// Read a reduced data file into a struct datafile
datafile load_reduced_data(char *dataFile)
{
    datafile h;
    h.file = NULL;
    h.version = 0;
    h.frame_dir = NULL;
    h.frame_pattern = NULL;
    h.dark_template = NULL;
    h.flat_template = NULL;
    h.num_obs = 0;
    h.num_targets = 0;
    h.num_blocked_ranges = 0;
    h.plot_fit_degree = PLOT_FIT_DEGREE_DEFAULT;
    h.plot_max_raw = PLOT_MAX_RAW_DEFAULT;
    h.plot_num_uhz = PLOT_NUM_UHZ_DEFAULT;
    h.plot_min_uhz = PLOT_MIN_UHZ_DEFAULT;
    h.plot_max_uhz = PLOT_MAX_UHZ_DEFAULT;

    // Open the data file (created with `tsreduce init`)
    h.file = fopen(dataFile, "r+");
    if (h.file == NULL)
        return h;

    char linebuf[1024];
    char stringbuf[1024];
    while (fgets(linebuf, sizeof(linebuf)-1, h.file) != NULL)
    {
        //
        // Header
        //
        if (!strncmp(linebuf,"# FramePattern:", 15))
        {
            sscanf(linebuf, "# FramePattern: %1024s\n", stringbuf);
            h.frame_pattern = strndup(stringbuf, 1024);
        }
        else if (!strncmp(linebuf,"# FrameDir:", 11))
        {
            sscanf(linebuf, "# FrameDir: %1024s\n", stringbuf);
            h.frame_dir = strndup(stringbuf, 1024);
        }
        else if (!strncmp(linebuf,"# DarkTemplate:", 15))
        {
            sscanf(linebuf, "# DarkTemplate: %1024s\n", stringbuf);
            h.dark_template = strndup(stringbuf, 1024);
        }
        else if (!strncmp(linebuf,"# FlatTemplate:", 15))
        {
            sscanf(linebuf, "# FlatTemplate: %1024s\n", stringbuf);
            h.flat_template = strndup(stringbuf, 1024);
        }
        else if (!strncmp(linebuf,"# Target:", 9))
        {
            if (h.version == 3)
            {
                sscanf(linebuf, "# Target: (%lf, %lf, %lf, %lf, %lf)\n",
                       &h.targets[h.num_targets].x,
                       &h.targets[h.num_targets].y,
                       &h.targets[h.num_targets].r,
                       &h.targets[h.num_targets].s1,
                       &h.targets[h.num_targets].s2);
                h.targets[h.num_targets].plot_scale = 1;
                h.num_targets++;
            }
            else
            {
                sscanf(linebuf, "# Target: (%lf, %lf, %lf, %lf, %lf) [%lf]\n",
                       &h.targets[h.num_targets].x,
                       &h.targets[h.num_targets].y,
                       &h.targets[h.num_targets].r,
                       &h.targets[h.num_targets].s1,
                       &h.targets[h.num_targets].s2,
                       &h.targets[h.num_targets].plot_scale);
                h.num_targets++;
            }
        }
        else if (!strncmp(linebuf,"# ReferenceTime:", 16))
        {
            struct tm t;
            strptime(linebuf, "# ReferenceTime: %Y-%m-%d %H:%M:%S\n", &t);
            h.reference_time = timegm(&t);
        }
        else if (!strncmp(linebuf,"# Version:", 10))
            sscanf(linebuf, "# Version: %d\n", &h.version);
        else if (!strncmp(linebuf,"# PlotFitDegree:", 16))
            sscanf(linebuf, "# PlotFitDegree: %d\n", &h.plot_fit_degree);
        else if (!strncmp(linebuf,"# PlotMaxRaw:", 13))
            sscanf(linebuf, "# PlotMaxRaw: %lf\n", &h.plot_max_raw);
        else if (!strncmp(linebuf,"# PlotMinUhz:", 13))
            sscanf(linebuf, "# PlotMinUhz: %lf\n", &h.plot_min_uhz);
        else if (!strncmp(linebuf,"# PlotMaxUhz:", 13))
            sscanf(linebuf, "# PlotMaxUhz: %lf\n", &h.plot_max_uhz);
        else if (!strncmp(linebuf,"# PlotNumUhz:", 13))
            sscanf(linebuf, "# PlotNumUhz: %d\n", &h.plot_num_uhz);
        else if (!strncmp(linebuf,"# BlockRange:", 13))
        {
            sscanf(linebuf, "# BlockRange: (%lf, %lf)\n",
                   &h.blocked_ranges[h.num_blocked_ranges].x,
                   &h.blocked_ranges[h.num_blocked_ranges].y);
            h.num_blocked_ranges++;
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
}
