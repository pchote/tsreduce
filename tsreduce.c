/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "datafile.h"
#include "helpers.h"
#include "dft_analysis.h"
#include "reduction.h"
#include "ts_analysis.h"

int verbosity = 0;

int main( int argc, char *argv[] )
{
    // `tsreduce create-flat "flat-[0-9]+.fits.gz" 5 master-dark.fits.gz master-flat.fits.gz`
    if (argc == 6 && strcmp(argv[1], "create-flat") == 0)
        return create_flat(argv[2], atoi(argv[3]), argv[4], argv[5]);

    // `tsreduce create-dark "dark-[0-9]+.fits.gz" 5 master-dark.fits.gz`
    else if (argc == 5 && strcmp(argv[1], "create-dark") == 0)
        return create_dark(argv[2], atoi(argv[3]), argv[4]);

    // `tsreduce update ~/data/20110704/gwlib.dat`
    else if (argc == 3 && strcmp(argv[1], "update") == 0)
        return update_reduction(argv[2]);

    // `tsreduce tracer test.dat`
    else if (argc == 3 && strcmp(argv[1], "tracer") == 0)
        return display_tracer(argv[2]);

    // `tsreduce display test.dat l19-2-0660.fits.gz`
    else if (argc == 4 && strcmp(argv[1], "display") == 0)
        return display_frame(argv[2], argv[3]);

    // `tsreduce display test.dat l19-2-0660.fits.gz`
    else if (argc == 3 && strcmp(argv[1], "metadata") == 0)
        return print_frame_metadata(argv[2]);

    // `tsreduce model july2011_run2.ts dftfreq.dat 0.0678829 6.3125 0.0001 fit.dat [residuals.dat]`
    else if ((argc == 8 || argc == 9) && strcmp(argv[1], "model") == 0)
        return model_fit(argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), argv[7], (argc == 9) ? argv[8] : NULL);

    // `tsreduce dft july2011_run2.ts 100 10000 0.01 dft.dat [dftfreq.dat]`
    else if ((argc == 7 || argc == 8) && strcmp(argv[1], "dft") == 0)
        return dft_bjd(argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]), argv[6], (argc == 8) ? argv[7] : NULL);

    // `tsreduce findfreq july2011_run2.ts dftfreq.dat 100 10000 1`
    else if (argc == 7 && strcmp(argv[1], "findfreq") == 0)
        return find_max_freq(argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]));

    // `tsreduce window july2011_run2.ts 1000 800 1200 0.01 window.dat`
    else if (argc == 8 && strcmp(argv[1], "window") == 0)
        return dft_window(argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]), argv[7]);

    // `tsreduce optimisefreq july2011_run2.ts dftfreq.dat`
    else if (argc == 4 && strcmp(argv[1], "optimizefreqs") == 0)
        return nonlinear_fit(argv[2], argv[3]);

    // `tsreduce reduce-range ec05221.dat 2 15 0.5 ec05221-range`
    else if (argc == 7 && strcmp(argv[1], "reduce-range") == 0)
        return reduce_aperture_range(argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]), argv[6]);

    // `tsreduce plot-range ec05221-range-[0-9]+.dat`
    else if (argc == 3 && strcmp(argv[1], "plot-range") == 0)
        return plot_range(argv[2]);

    // `tsreduce plot ec04207.dat [ts.ps/cps 10 dft.ps/cps 10]
    else if ((argc == 3 || argc == 7) && strcmp(argv[1], "plot") == 0)
    {
        if (argc == 3)
        {
#if (defined _WIN32)
            double size = 10;
            char *ts_device = "ts.gif/gif";
            char *dft_device = "dft.gif/gif";
#else
            double size = 9.41;
            char *ts_device = "5/xs";
            char *dft_device = "6/xs";
#endif
            return online_plot(argv[2], ts_device, size, dft_device, size);
        }
        else
            return online_plot(argv[2], argv[3], atof(argv[4]), argv[5], atof(argv[6]));
    }

    // `tsreduce playback ec04207.dat 100 1 [ts.ps/cps 10 dft.ps/cps 10]
    else if ((argc == 5 || argc == 9) && strcmp(argv[1], "playback") == 0)
    {
        if (argc == 5)
        {
#if (defined _WIN32)
            double size = 10;
            char *ts_device = "ts.gif/gif";
            char *dft_device = "dft.gif/gif";
#else
            double size = 9.41;
            char *ts_device = "5/xs";
            char *dft_device = "6/xs";
#endif
            return playback_reduction(argv[2], atoi(argv[3]), atoi(argv[4]), ts_device, size, dft_device, size);
        }
        else
            return playback_reduction(argv[2], atoi(argv[3]), atoi(argv[4]), argv[5], atof(argv[6]), argv[7], atof(argv[8]));
    }

    // `tsreduce create-ts 2011-07-27 13:00:00 20110727.dat [...] july2011_run2_rereduce.ts`
    else if (argc >= 3 && strcmp(argv[1], "create-ts") == 0)
        return create_ts(argv[2], argv[3], &argv[4], argc - 5, argv[argc - 1]);

    // `tsreduce shuffle-dft july2011_run2_rereduce.ts july2011_run2_rereduce_paperfreqs.dat 100 10000 1 rereduce.ran 1000`
    else if (argc == 9 && strncmp(argv[1], "shuffle-dft", 11) == 0)
        return shuffle_dft(argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), argv[7], atoi(argv[8]));

    // `tsreduce create ~/data/20110307 ec04207-[0-9]+.fits.gz master-dark.fits master-skyflat-20110305.fits 20110307.dat`
    else if (argc == 3 && strcmp(argv[1], "create") == 0)
        return create_reduction_file(argv[2]);

    else if (argc == 4 && strcmp(argv[1], "preview") == 0)
        return update_preview(argv[2], argv[3]);

    // Misc one-off utility functions
    else if (argc == 5 && strcmp(argv[1], "profile") == 0)
        return calculate_profile(argv[2], atoi(argv[3]), atoi(argv[4]));
    else if (argc == 3 && strcmp(argv[1], "repeats") == 0)
        return detect_repeats(argv[2]);
    else if (argc == 7 && strcmp(argv[1], "bjd") == 0)
        return calculate_bjd(argv[2], argv[3], argv[4], argv[5], atoi(argv[6]));
    else
        error("Invalid args");
    return 0;
}
