/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cpgplot.h>

#include "helpers.h"
#include "fit.h"


/*
 * Allocates memory and loads timeseries data, returning pointers via the passed args
 * Returns the number of data points loaded, or -1 on error
 */
static int load_tsfile(char *tsFile, double **time, double **mmi)
{
    if (*time != NULL || *mmi != NULL)
    {
        error("time or mmi array is non-null");
        return -1;
    }

    char linebuf[1024];
    FILE *file = fopen(tsFile, "r+");
    if (file == NULL)
    {
        error("Unable to open file: %s", tsFile);
        return -1;
    }

    // Count the number of entries to allocate
    int total_obs = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total_obs++;
    rewind(file);

    int num_obs = 0;
    *time = (double *)malloc(total_obs*sizeof(double));
    *mmi = (double *)malloc(total_obs*sizeof(double));
    if (*time == NULL || *mmi == NULL)
    {
        error("Unable to allocate memory");
        return -1;
    }

    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num_obs < total_obs)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        sscanf(linebuf, "%lf %lf\n", &(*time)[num_obs], &(*mmi)[num_obs]);

        // Convert to seconds
        (*time)[num_obs] *= 86400;
        num_obs++;
    }

    fclose(file);
    return num_obs;
}

/*
 * Allocates memory and loads frequency data, returning a pointers via the passed arg
 * Returns the number of frequencies loaded, or -1 on error
 */
static int load_freqfile(char *freqFile, double **freqs)
{
    char linebuf[1024];
    FILE *file = fopen(freqFile, "r+");
    if (file == NULL)
    {
        error("Unable to open file: %s", freqFile);
        return -1;
    }

    int total_freqs = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total_freqs++;
    rewind(file);

    int num_freqs = 0;
    *freqs = (double *)malloc(total_freqs*sizeof(double));
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num_freqs < total_freqs)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        int use;
        sscanf(linebuf, "%*d %lf %d %*s\n", &(*freqs)[num_freqs], &use);

        // Convert to Hz
        (*freqs)[num_freqs] *= 1e-6;
        if (use)
            num_freqs++;
    }
    fclose(file);
    return num_freqs;
}

static void print_freq_table(double *freqs, double *amplitudes, int numFreqs)
{
    printf("Freq (uHz) Period (s) Amp (mma) Phase (deg)\n");
    for (int i = 0; i < numFreqs; i++)
    {
        double a = amplitudes[2*i+1];
        double b = amplitudes[2*i];
        double amp = sqrt(a*a + b*b);
        double phase = atan2(b, a)*180/M_PI;
        while (phase > 360) phase -= 360;
        while (phase < 0) phase += 360;
        printf("%10.2f %10.2f %9.2f %11.2f\n", 1e6*freqs[i], 1/freqs[i], amp, phase);
    }
}

static double calculate_chi2(double *freqs, double *amplitudes, int numFreqs, double *mmi, double *time, int numObs)
{
    double chi2 = 0;
    for (int i = 0; i < numObs; i++)
    {
        double model = 0;
        for (int j = 0; j < numFreqs; j++)
        {
            double phase = 2*M_PI*freqs[j]*time[i];
            model += amplitudes[2*j]*cos(phase);
            model += amplitudes[2*j+1]*sin(phase);
        }
        chi2 += (mmi[i] - model)*(mmi[i]-model);
    }
    return chi2;
}

/*
 * Fit the frequencies defined in freqFile to the data in tsFile
 * and output a model lightcurve between startTime and endTime with increments of dt to modelFile
 * Residuals are saved to residualsFile if it is non-NULL
 */
int model_fit(char *tsFile, char *freqFile, double startTime, double endTime, double dt, char *modelFile, char *residualsFile)
{
    double *time = NULL, *mmi = NULL, *fit_freqs = NULL;
    int num_obs = load_tsfile(tsFile, &time, &mmi);
    if (num_obs <= 0)
        return error("ts load failed");
    else if (num_obs == 0)
    {
        free(time);
        free(mmi);
        return error("No observations in ts file");
    }
    printf("Read %d observations\n", num_obs);

    // Load freqs
    int num_freqs = load_freqfile(freqFile, &fit_freqs);

    if (num_freqs < 0)
    {
        free(time);
        free(mmi);
        return error("freq load failed");
    }
    printf("Read %d freqs\n", num_freqs);

    if (num_freqs == 0)
    {
        free(fit_freqs);
        free(time);
        free(mmi);
        return error("No frequencies found");
    }

    // Fit amplitudes for each freq
    double *fit_amplitudes = (double *)malloc(2*num_freqs*sizeof(double));
    if (fit_sinusoids(time, mmi, num_obs, fit_freqs, num_freqs, fit_amplitudes))
    {
        free(fit_amplitudes);
        free(fit_freqs);
        free(time);
        free(mmi);
        return error("Fit failed");
    }

    print_freq_table(fit_freqs, fit_amplitudes, num_freqs);
    double chi2 = calculate_chi2(fit_freqs, fit_amplitudes, num_freqs, mmi, time, num_obs);
    printf("Chi^2: %f\n", chi2);

    // Output model curve
    FILE *file = fopen(modelFile, "w");
    if (file == NULL)
    {
        free(fit_amplitudes);
        free(fit_freqs);
        free(time);
        free(mmi);
        return error("Unable to open file: %s", modelFile);
    }

    for (double t = startTime; t <= endTime; t += dt)
    {
        double fit = 0;
        for (int j = 0; j < num_freqs; j++)
        {
            // convert time from BJD to seconds
            double phase = 2*M_PI*fit_freqs[j]*t*86400;
            fit += fit_amplitudes[2*j]*cos(phase);
            fit += fit_amplitudes[2*j+1]*sin(phase);
        }
        fprintf(file, "%f %f\n", t, fit);
    }
    fclose(file);

    // Output residuals
    if (residualsFile != NULL)
    {
        file = fopen(residualsFile, "w");
        if (file == NULL)
        {
            free(fit_amplitudes);
            free(fit_freqs);
            free(time);
            free(mmi);
            return error("Unable to open file: %s", residualsFile);
        }

        for (int i = 0; i < num_obs; i++)
        {
            double model = 0;
            for (int j = 0; j < num_freqs; j++)
            {
                double phase = 2*M_PI*fit_freqs[j]*time[i];
                model += fit_amplitudes[2*j]*cos(phase);
                model += fit_amplitudes[2*j+1]*sin(phase);
            }
            fprintf(file, "%f %f\n", time[i]/86400, mmi[i] - model);
        }
        fclose(file);
    }

    free(fit_amplitudes);
    free(fit_freqs);
    free(time);
    free(mmi);

    return 0;
}

/*
 * Calculate the DFT of the BJD/mma data in tsFile between minUHz and maxUHz in increments dUHz
 * Output is saved to outFile with frequencies in uHz.
 * If freqFile is non-NULL, the data is first prewhitened with the contained frequencies
 */
int dft_bjd(char *tsFile, double minUHz, double maxUHz, double dUHz, char *outFile, char *freqFile)
{
    double *time = NULL, *mmi = NULL;
    int num_obs = load_tsfile(tsFile, &time, &mmi);
    if (num_obs <= 0)
        return error("ts load failed");
    else if (num_obs == 0)
    {
        free(time);
        free(mmi);
        return error("No observations in ts file");
    }
    printf("Read %d observations\n", num_obs);

    // Load freqs
    if (freqFile)
    {
        // Load freqs
        double *fit_freqs;
        int num_freqs = load_freqfile(freqFile, &fit_freqs);

        if (num_freqs < 0)
        {
            free(time);
            free(mmi);
            return error("freq load failed");
        }
        printf("Read %d freqs\n", num_freqs);

        if (num_freqs == 0)
        {
            free(fit_freqs);
            free(time);
            free(mmi);
            return error("No frequencies found");
        }

        // Fit amplitudes for each freq
        double *fit_amplitudes = (double *)malloc(2*num_freqs*sizeof(double));
        if (fit_sinusoids(time, mmi, num_obs, fit_freqs, num_freqs, fit_amplitudes))
        {
            free(fit_amplitudes);
            free(fit_freqs);
            free(time);
            free(mmi);
            return error("Fit failed");
        }

        print_freq_table(fit_freqs, fit_amplitudes, num_freqs);
        double chi2 = calculate_chi2(fit_freqs, fit_amplitudes, num_freqs, mmi, time, num_obs);
        printf("Chi^2: %f\n", chi2);


        // Prewhiten
        for (int i = 0; i < num_obs; i++)
        {
            double model = 0;
            for (int j = 0; j < num_freqs; j++)
            {
                double phase = 2*M_PI*fit_freqs[j]*time[i];
                model += fit_amplitudes[2*j]*cos(phase);
                model += fit_amplitudes[2*j+1]*sin(phase);
            }
            mmi[i] -= model;
        }
        free(fit_amplitudes);
        free(fit_freqs);
    }

    // Calculate DFT
    int num_uhz = (int)((maxUHz - minUHz)/dUHz);
    double *freq = (double *)malloc(num_uhz*sizeof(double));
    double *ampl = (double *)malloc(num_uhz*sizeof(double));

    calculate_amplitude_spectrum(minUHz*1e-6, maxUHz*1e-6, time, mmi, num_obs, freq, ampl, num_uhz);

    // Save output
    FILE *file = fopen(outFile, "w");
    if (file == NULL)
    {
        free(ampl);
        free(freq);
        free(time);
        free(mmi);
        return error("Unable to open file: %s", outFile);
    }

    for (int i = 0; i < num_uhz; i++)
        fprintf(file, "%f %f\n", 1e6*freq[i], ampl[i]);

    fclose(file);

    free(ampl);
    free(freq);
    free(time);
    free(mmi);

    return 0;
}

/*
 * Calculate the DFT window at freq (in uHz) for the BJD/mma data in tsFile between minUHz and maxUHz in increments dUHz
 * Output is saved to outFile with frequencies in uHz.
 */
int dft_window(char *tsFile, double windowFreq, double minUHz, double maxUHz, double dUHz, char *outFile)
{
    double *time = NULL, *mmi = NULL;
    int num_obs = load_tsfile(tsFile, &time, &mmi);
    if (num_obs <= 0)
        return error("ts load failed");
    else if (num_obs == 0)
    {
        free(time);
        free(mmi);
        return error("No observations in ts file");
    }
    printf("Read %d observations\n", num_obs);


    // Calculate DFT
    int num_uhz = (int)((maxUHz - minUHz)/dUHz);
    double *freq = (double *)malloc(num_uhz*sizeof(double));
    double *ampl = (double *)malloc(num_uhz*sizeof(double));
    calculate_amplitude_spectrum(minUHz*1e-6, maxUHz*1e-6, time, mmi, num_obs, freq, ampl, num_uhz);

    // Save output
    FILE *file = fopen(outFile, "w");
    if (file != NULL)
    {
        for (int i = 0; i < num_uhz; i++)
            fprintf(file, "%f %f\n", 1e6*freq[i], ampl[i]);
        fclose(file);
    }

    free(ampl);
    free(freq);
    free(time);
    free(mmi);
    return file != NULL ? 0 : error("Unable to open file: %s", outFile);
}

/*
 * Calculate the DFT of the BJD/mma data in tsFile between minUHz and maxUHz in increments dUHz
 * and find the frequency with the largest amplitude
 * Results are plotted in a pgplot window
 * If freqFile is non-NULL, the data is first prewhitened with the contained frequencies
 */
int find_max_freq(char *tsFile, char *freqFile, double minUHz, double maxUHz, double dUHz)
{
    double *time = NULL, *mmi = NULL, *fit_freqs = NULL;
    int num_obs = load_tsfile(tsFile, &time, &mmi);
    if (num_obs <= 0)
        return error("ts load failed");
    else if (num_obs == 0)
    {
        free(time);
        free(mmi);
        return error("No observations in ts file");
    }
    printf("Read %d observations\n", num_obs);

    // Load freqs
    int num_freqs = load_freqfile(freqFile, &fit_freqs);

    if (num_freqs < 0)
    {
        free(time);
        free(mmi);
        return error("freq load failed");
    }
    printf("Read %d freqs\n", num_freqs);

    // Fit amplitudes for each freq
    double *fit_amplitudes = (double *)malloc(2*num_freqs*sizeof(double));
    if (fit_sinusoids(time, mmi, num_obs, fit_freqs, num_freqs, fit_amplitudes))
    {
        free(fit_amplitudes);
        free(fit_freqs);
        free(time);
        free(mmi);
        return error("Fit failed");
    }

    print_freq_table(fit_freqs, fit_amplitudes, num_freqs);
    double chi2 = calculate_chi2(fit_freqs, fit_amplitudes, num_freqs, mmi, time, num_obs);
    printf("Chi^2: %f\n", chi2);

    // Prewhiten
    for (int i = 0; i < num_obs; i++)
    {
        double model = 0;
        for (int j = 0; j < num_freqs; j++)
        {
            double phase = 2*M_PI*fit_freqs[j]*time[i];
            model += fit_amplitudes[2*j]*cos(phase);
            model += fit_amplitudes[2*j+1]*sin(phase);
        }
        mmi[i] -= model;
    }
    free(fit_amplitudes);
    free(fit_freqs);

    // Calculate DFT
    int num_uhz = (int)((maxUHz - minUHz)/dUHz);
    double *freq = (double *)malloc(num_uhz*sizeof(double));
    double *ampl = (double *)malloc(num_uhz*sizeof(double));

    calculate_amplitude_spectrum(minUHz*1e-6, maxUHz*1e-6, time, mmi, num_obs, freq, ampl, num_uhz);

    // Determine max dft ampl and convert freq,ampl arrays to floats for PGPLOT
    double max_ampl = 0;
    double max_ampl_freq = 0;
    float *pgfreq = (float *)malloc(num_uhz*sizeof(float));
    float *pgampl = (float *)malloc(num_uhz*sizeof(float));
    for (int i = 0; i < num_uhz; i++)
    {
        if (ampl[i] > max_ampl)
        {
            max_ampl = ampl[i];
            max_ampl_freq = freq[i];
        }

        pgfreq[i] = freq[i];
        pgampl[i] = ampl[i];
    }
    printf("Max amplitude of %f at %f uHz\n", max_ampl, max_ampl_freq*1e6);

    // Display
    if (cpgopen("6/xs") <= 0)
    {
        free(pgampl);
        free(pgfreq);
        free(ampl);
        free(freq);
        free(time);
        free(mmi);
        return error("Unable to open PGPLOT window");
    }

    // 800 x 480
    cpgpap(9.41, 0.6);
    cpgask(0);
    cpgslw(3);
    cpgsfs(2);
    cpgscf(2);

    // DFT
    cpgsvp(0.1, 0.9, 0.075, 0.9);
    cpgswin(minUHz, maxUHz, 0, max_ampl*1.1);
    cpgbox("bcstn", 0, 0, "bcnst", 0, 0);

    cpgswin(minUHz*1e-6, maxUHz*1e-6, 0, max_ampl*1.1);

    cpgsci(12);
    cpgline(num_uhz, pgfreq, pgampl);

    cpgsci(2);
    cpgmove(max_ampl_freq, max_ampl);
    cpgdraw(max_ampl_freq, 0);
    cpgsci(1);

    cpgmtxt("b", 2.5, 0.5, 0.5, "Frequency (\\gmHz)");
    cpgmtxt("l", 2, 0.5, 0.5, "Amplitude (mma)");

    cpgsch(1.25);
    cpgmtxt("t", 1.5, 0.5, 0.5, "Amplitude Spectrum");


    free(pgampl);
    free(pgfreq);
    free(ampl);
    free(freq);
    free(time);
    free(mmi);

    return 0;
}

static int step_freq(double *fit_freqs, double *fit_amplitudes, int num_freqs,
                     double *mmi, double *time, int num_obs,
                     int vary_freq, double centerFreq, double stepSize, int numSteps,
                     double *best_freq, double *best_chi2)
{
    for (int j = -numSteps; j <= numSteps; j++)
    {
        fit_freqs[vary_freq] = centerFreq + j*stepSize;

        if (fit_sinusoids(time, mmi, num_obs, fit_freqs, num_freqs, fit_amplitudes))
            continue;

        double chi2 = calculate_chi2(fit_freqs, fit_amplitudes, num_freqs, mmi, time, num_obs);
        if (chi2 < *best_chi2)
        {
            *best_chi2 = chi2;
            *best_freq = fit_freqs[vary_freq];
        }
    }
    return 0;
}

int nonlinear_fit(char *tsFile, char *freqFile)
{
    double *time = NULL, *mmi = NULL, *init_freqs = NULL;
    int num_obs = load_tsfile(tsFile, &time, &mmi);

    if (num_obs <= 0)
        return error("ts load failed");

    else if (num_obs == 0)
    {
        free(time);
        free(mmi);
        return error("No observations in ts file");
    }
    printf("Read %d observations\n", num_obs);

    // Load freqs
    int num_freqs = load_freqfile(freqFile, &init_freqs);
    if (num_freqs < 0)
    {
        free(time);
        free(mmi);
        return error("freq load failed");
    }
    printf("Read %d freqs\n", num_freqs);

    // Fit amplitudes for each freq
    double *fit_amplitudes = (double *)malloc(2*num_freqs*sizeof(double));
    if (fit_sinusoids(time, mmi, num_obs, init_freqs, num_freqs, fit_amplitudes))
    {
        free(fit_amplitudes);
        free(init_freqs);
        free(time);
        free(mmi);
        return error("Fit failed");
    }

    // Calculate amplitude and phase for each freq
    printf("Initial settings:\n");
    print_freq_table(init_freqs, fit_amplitudes, num_freqs);
    double chi2 = calculate_chi2(init_freqs, fit_amplitudes, num_freqs, mmi, time, num_obs);
    printf("Chi^2: %f\n", chi2);

    double *fit_freqs = (double *)malloc(num_freqs*sizeof(double));
    for (int j = 0; j < num_freqs; j++)
        fit_freqs[j] = init_freqs[j];

    // Vary frequency between x-20 .. x + 20 in 0.01 steps
    double coarse_step = 1e-6;
    int coarse_step_count = 2;
    double fine_step = 0.01e-6;
    int fine_step_count = 20;

    double initialchi2 = chi2;
    double lastouterchi2 = chi2;
    double lastchi2 = chi2;
    
    do
    {
        lastouterchi2 = chi2;
        for (int i = 0; i < num_freqs; i++)
        {
            double best_chi2 = chi2;
            double best_freq = fit_freqs[i];

            // Start with a rough step
            printf("Freq %10.2f\n", 1e6*init_freqs[i]);
            int iterate_mode = 2;
            double step = coarse_step;
            int num_steps = coarse_step_count;
            printf("\tUsing coarse steps\n");
            do
            {
                double last_best_freq = best_freq;
                if (step_freq(fit_freqs, fit_amplitudes, num_freqs,
                              mmi, time, num_obs,
                              i, best_freq, step, num_steps,
                              &best_freq, &best_chi2))
                {
                    free(fit_amplitudes);
                    free(init_freqs);
                    free(fit_freqs);
                    free(time);
                    free(mmi);
                    return error("Variation failed");
                }

                if (best_chi2 < lastchi2)
                {
                    printf("\t%10.2f -> %10.2f dchi2: %f\n", 1e6*last_best_freq, 1e6*best_freq, lastchi2 - best_chi2);
                    lastchi2 = best_chi2;
                }
                else
                {
                    // Switch to fine iteration
                    if (--iterate_mode)
                    {
                        printf("\tUsing fine steps\n");
                        step = fine_step;
                        num_steps = fine_step_count;
                    }
                }
            }
            while (iterate_mode);

            fit_freqs[i] = best_freq;
            chi2 = lastchi2;
        }
    } while (chi2 < lastouterchi2);

    printf("%f -> %f (%f)\n", initialchi2, chi2, initialchi2 - chi2);
    for (int i = 0; i < num_freqs; i++)
        printf("%d %10.2f 1\n", i, 1e6*fit_freqs[i]);

    free(fit_freqs);
    free(fit_amplitudes);
    free(init_freqs);
    free(time);
    free(mmi);

    return 0;
}

/*
 * Shift 2011 first july run observations forward and backward
 * in time to find the best fit to the 2nd july run model
 */
int fit_time(char *tsFile)
{
    // Second july run fit parameters
    int num_freqs = 14;
    double freqs[] =
    {
        0.0022362400,
        0.0044724500,
        0.0023613300,
        0.0029726500,
        0.0067087500,
        0.0052088500,
        0.0046091000,
        0.0016810000,
        0.0022336800,
        0.0074450700,
        0.0029182400,
        0.0089449700,
        0.0037938300,
        0.0023581800
    };
    double amplitudes[] =
    {
        -32.735226, -20.399597,
        11.016796, -1.155572,
         8.461928, -1.198719,
        -4.585997,  3.090083,
        -3.266993,  2.921126,
         0.534985, -3.280352,
         1.121596,  2.381564,
         1.657419, -1.233287,
        -0.241483, -2.256900,
         0.901674,  1.725242,
        -1.622693, -0.596785,
         0.864272, -1.399492,
        -0.009832, -1.506967,
         1.560148, -1.033600
    };

    double *time = NULL, *mmi = NULL;
    int num_obs = load_tsfile(tsFile, &time, &mmi);

    if (num_obs <= 0)
        return error("ts load failed");

    else if (num_obs == 0)
    {
        free(time);
        free(mmi);
        return error("No observations in ts file");
    }
    printf("Read %d observations\n", num_obs);

    double best_chi2 = -1;
    double best_offset = 0;
    int num_steps = 40;
    double step_size = 20;
    for(int i = -num_steps; i <= num_steps; i++)
    {
        double offset = i*step_size;

        double chi2 = 0;
        for (int i = 0; i < num_obs; i++)
        {
            double model = 0;
            for (int j = 0; j < num_freqs; j++)
            {
                double phase = 2*M_PI*freqs[j]*(time[i]-offset);
                model += amplitudes[2*j]*cos(phase);
                model += amplitudes[2*j+1]*sin(phase);
            }
            chi2 += (mmi[i] - model)*(mmi[i]-model);
        }
        printf("%f: %f\n", offset, chi2);
        if (chi2 < best_chi2 || best_chi2 < 0)
        {
            best_chi2 = chi2;
            best_offset = offset;
        }
    }
    printf("Best offset is %f (%d frames) with chi2 %f\n", best_offset, (int)best_offset/20, best_chi2);
    free(time);
    free(mmi);
    return 0;
}
