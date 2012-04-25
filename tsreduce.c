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
#include "fit.h"

#include "dft_analysis.h"
#include "reduction.h"
#include "ts_analysis.h"

// Read the header section from the data file
// defined as all the lines at the top of the file
// that start with #
datafile read_data_header(char *dataFile)
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

    else if ((argc >= 3 && argc <= 5) && strncmp(argv[1], "plot", 4) == 0)
        return plot_fits(argv[2], argc > 3 ? argv[3] : NULL, argc > 4 ? argv[4] : NULL);

    // `tsreduce model july2011_run2.ts dftfreq.dat 0.0678829 6.3125 0.0001 fit.dat [residuals.dat]`
    else if ((argc == 8 || argc == 9) && strncmp(argv[1], "model", 5) == 0)
        return model_fit(argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), argv[7], (argc == 9) ? argv[8] : NULL);

    // `tsreduce dft july2011_run2.ts 100 10000 0.01 dft.dat [dftfreq.dat]`
    else if ((argc == 7 || argc == 8) && strncmp(argv[1], "dft", 3) == 0)
        return dft_bjd(argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]), argv[6], (argc == 8) ? argv[7] : NULL);

    // `tsreduce findfreq july2011_run2.ts dftfreq.dat 100 10000 1`
    else if (argc == 7 && strncmp(argv[1], "findfreq", 8) == 0)
        return find_max_freq(argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]));

    // `tsreduce window july2011_run2.ts 1000 800 1200 0.01 window.dat`
    else if (argc == 8 && strncmp(argv[1], "window", 6) == 0)
        return dft_window(argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]), argv[7]);

    // `tsreduce optimisefreq july2011_run2.ts dftfreq.dat`
    else if (argc == 4 && strncmp(argv[1], "optimizefreqs", 13) == 0)
        return nonlinear_fit(argv[2], argv[3]);

    else if (argc == 3 && strncmp(argv[1], "fittime", 7) == 0)
        return fit_time(argv[2]);

    else if (argc == 4 && strncmp(argv[1], "offsettime", 10) == 0)
        return offset_time(argv[2], atof(argv[3]));
    else if (argc == 3 && strncmp(argv[1], "readnoise", 9) == 0)
        return ccd_readnoise(argv[2]);
    else if (argc == 3 && strncmp(argv[1], "mmi", 3) == 0)
        return create_mmi(argv[2]);

    else
        error("Invalid args");
    return 0;
}
