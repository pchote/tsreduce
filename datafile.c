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

// Read a reduced data file into a struct datafile
datafile load_reduced_data(char *dataFile)
{
    datafile h;
    h.file = NULL;
    h.version = 0;
    h.frame_dir[0] = '\0';
    h.frame_pattern[0] = '\0';
    h.dark_template[0] = '\0';
    h.flat_template[0] = '\0';
    h.num_obs = 0;
    h.num_targets = 0;
    h.num_blocked_ranges = 0;
    h.plot_fit_degree = 4;
    h.plot_max_raw = 5000;
    h.plot_num_uhz = 1000;
    h.plot_min_uhz = 0;
    h.plot_max_uhz = 10000;

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
