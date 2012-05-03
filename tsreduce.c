/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#include <string.h>

#include "datafile.h"
#include "helpers.h"
#include "dft_analysis.h"
#include "reduction.h"
#include "ts_analysis.h"

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
    else if (argc == 6 && strncmp(argv[1], "snr", 3) == 0)
        return evaluate_aperture_snr(argv[2], atof(argv[3]), atof(argv[4]), atoi(argv[5]));

    else
        error("Invalid args");
    return 0;
}
