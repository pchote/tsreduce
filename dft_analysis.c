/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cpgplot.h>

#include "helpers.h"
#include "fit.h"
#include "random.h"

/*
 * Allocates memory and loads timeseries data, returning pointers via the passed args
 * Returns the number of data points loaded, or -1 on error
 */
static int load_tsfile(char *tsFile, double **time, double **mma, double **err, size_t *loaded)
{
    int ret = 0;
    if (!time || !mma || !err)
        return error("Invalid input array");

    char linebuf[1024];
    FILE *file = fopen(tsFile, "r+");
    if (!file)
        return error("Unable to open file: %s", tsFile);

    // Count the number of entries to allocate
    int total_obs = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total_obs++;
    rewind(file);

    *time = (double *)malloc(total_obs*sizeof(double));
    if (!*time)
        error_jump(time_alloc_error, ret, "Error allocating time array");

    *mma = (double *)malloc(total_obs*sizeof(double));
    if (!*mma)
        error_jump(mma_alloc_error, ret, "Error allocating mma array");

    *err = (double *)malloc(total_obs*sizeof(double));
    if (!*err)
        error_jump(err_alloc_error, ret, "Error allocating err array");

    int num_obs = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num_obs < total_obs)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        int read = sscanf(linebuf, "%lf %lf %lf\n", &(*time)[num_obs], &(*mma)[num_obs], &(*err)[num_obs]);

        // No error defined; set to unity
        if (read == 2)
            (*err)[num_obs] = 1;

        // Convert to seconds
        (*time)[num_obs] *= 86400;
        num_obs++;
    }
    *loaded = num_obs;
    fclose(file);
    return 0;

err_alloc_error:
    free(*mma);
mma_alloc_error:
    free(*time);
time_alloc_error:
    fclose(file);
    return ret;
}

/*
 * Allocates memory and loads frequency data, returning a pointers via the passed args
 */
static int load_freqfile(char *freqFile, char ***labels, double **freqs, int **mode, size_t *loaded)
{
    int ret = 0;
    if (!labels || !freqs || !mode)
        return error("Invalid input array");

    FILE *file = fopen(freqFile, "r+");
    if (file == NULL)
        return error("Unable to open file: %s", freqFile);

    int total_freqs = 0;
    char linebuf[1024];
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total_freqs++;
    rewind(file);

    int num_freqs = 0;
    *labels = (char **)malloc(total_freqs*sizeof(char **));
    if (!labels)
        error_jump(labels_alloc_error, ret, "Error allocating labels array");

    *freqs = (double *)malloc(total_freqs*sizeof(double));
    if (!freqs)
        error_jump(freqs_alloc_error, ret, "Error allocating freqs array");

    *mode = (int *)malloc(total_freqs*sizeof(int));
    if (!mode)
        error_jump(mode_alloc_error, ret, "Error allocating mode array");

    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num_freqs < total_freqs)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        int m;
        char label[1024];
        sscanf(linebuf, "%1024s %lf %d %*s\n", label, &(*freqs)[num_freqs], &m);

        (*freqs)[num_freqs] *= 1e-6;
        (*mode)[num_freqs] = m;
        if (m)
        {
            (*labels)[num_freqs] = strdup(label);
            num_freqs++;
        }
    }

    fclose(file);
    *loaded = num_freqs;
    return 0;

mode_alloc_error:
    free(*freqs);
freqs_alloc_error:
    free(*labels);
labels_alloc_error:
    fclose(file);
    return ret;
}

static void print_freq_table(char **labels, double *freqs, double *amplitudes, size_t numFreqs)
{
    printf("ID     Freq    Period    Amp    Phase \n");
    printf("       (uHz)     (s)    (mma)   (deg)\n");
    for (size_t i = 0; i < numFreqs; i++)
    {
        double a = amplitudes[2*i+1];
        double b = amplitudes[2*i];
        double amp = sqrt(a*a + b*b);
        double phase = atan2(b, a)*180/M_PI;
        while (phase > 360) phase -= 360;
        while (phase < 0) phase += 360;
        printf("%-5s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n", labels[i], 1e6*freqs[i], 1/freqs[i], amp, phase, 1/freqs[i]/60, 1/freqs[i]/3600);
    }
}

static double calculate_chi2(double *freqs, double *amplitudes, size_t numFreqs,
                             double *time, double *mma, double *err, size_t numObs)
{
    double chi2 = 0;
    for (size_t i = 0; i < numObs; i++)
    {
        double model = 0;
        for (size_t j = 0; j < numFreqs; j++)
        {
            double phase = 2*M_PI*freqs[j]*time[i];
            model += amplitudes[2*j]*cos(phase);
            model += amplitudes[2*j+1]*sin(phase);
        }
        double t = (mma[i] - model)/err[i];
        chi2 += t*t;
    }
    return chi2;
}

/*
 * Load a tsfile and prewhiten by the frequencies in freqFile.
 * If freqFile is null, the file will be loaded without prewhitening
 * and the freq arrays will be set to NULL
 */
static int load_and_fit_freqs(char *tsFile, char *freqFile,
                              double **time, double **mma, double **err, size_t *num_obs,
                              char ***freq_labels, double **freqs, double **freq_amplitudes, int **freq_modes, size_t *num_freqs)
{
    int ret = 0;

    if (load_tsfile(tsFile, time, mma, err, num_obs))
        return error("ts load failed");

    if (*num_obs == 0)
        error_jump(no_obs_error, ret, "No observations in ts file");

    printf("Read %zu observations\n", *num_obs);

    if (freqFile)
    {
        if (load_freqfile(freqFile, freq_labels, freqs, freq_modes, num_freqs))
            error_jump(freq_load_error, ret, "Freq load failed");

        if (num_freqs == 0)
            error_jump(no_freqs_error, ret, "No frequencies specified");

        printf("Read %zu freqs\n", *num_freqs);

        // Fit amplitudes for each freq
        *freq_amplitudes = (double *)malloc(2*(*num_freqs)*sizeof(double));
        if (!(*freq_amplitudes))
            error_jump(freq_amplitude_alloc_failed, ret, "Error allocating freq_amplitude array");

        if (fit_sinusoids(*time, *mma, *err, *num_obs, *freqs, *num_freqs, *freq_amplitudes))
            error_jump(fit_failed_error, ret, "Sinusoid fit failed");

        print_freq_table(*freq_labels, *freqs, *freq_amplitudes, *num_freqs);
        double chi2 = calculate_chi2(*freqs, *freq_amplitudes, *num_freqs, *time, *mma, *err, *num_obs);
        int dof = *num_obs - 3*(*num_freqs);
        printf("Chi^2: %f / DOF: %d = %f\n", chi2, dof, chi2/dof);
    }
    else
    {
        *freq_labels = NULL;
        *freqs = NULL;
        *freq_amplitudes = NULL;
        *freq_modes = 0;
        *num_freqs = 0;
    }

    return ret;

fit_failed_error:
    free(*freq_amplitudes);
freq_amplitude_alloc_failed:
no_freqs_error:
    free(*freq_labels);
    free(*freqs);
    free(*freq_modes);
freq_load_error:
no_obs_error:
    free(*time);
    free(*mma);
    free(*err);
    return ret;
}

/*
 * Fit the frequencies defined in freqFile to the data in tsFile
 * and output a model lightcurve between startTime and endTime with increments of dt to modelFile
 * Residuals are saved to residualsFile if it is non-NULL
 */
int model_fit(char *tsFile, char *freqFile, double startTime, double endTime, double dt, char *modelFile, char *residualsFile)
{
    int ret = 0;

    double *time, *mma, *err, *freqs, *freq_amplitudes;
    char **freq_labels;
    int *freq_modes;
    size_t num_obs, num_freqs;

    if (load_and_fit_freqs(tsFile, freqFile,
        &time, &mma, &err, &num_obs,
        &freq_labels, &freqs, &freq_amplitudes, &freq_modes, &num_freqs))
    {
        error_jump(load_failed_error, ret, "Error processing data");
    }

    if (!freqs)
        error_jump(no_freqs_error, ret, "freqfile is required");

    // Output model curve
    FILE *file = fopen(modelFile, "w");
    if (!file)
        error_jump(outfile_open_error, ret, "Error opening output file %s", modelFile);

    for (double t = startTime; t <= endTime; t += dt)
    {
        double fit = 0;
        for (size_t j = 0; j < num_freqs; j++)
        {
            // convert time from BJD to seconds
            double phase = 2*M_PI*freqs[j]*t*86400;
            fit += freq_amplitudes[2*j]*cos(phase);
            fit += freq_amplitudes[2*j+1]*sin(phase);
        }
        fprintf(file, "%f %f\n", t, fit);
    }
    fclose(file);

    // Output residuals
    if (residualsFile != NULL)
    {
        file = fopen(residualsFile, "w");
        if (!file)
            error_jump(outfile_open_error, ret, "Error opening output file %s", modelFile);

        for (size_t i = 0; i < num_obs; i++)
        {
            double model = 0;
            for (size_t j = 0; j < num_freqs; j++)
            {
                double phase = 2*M_PI*freqs[j]*time[i];
                model += freq_amplitudes[2*j]*cos(phase);
                model += freq_amplitudes[2*j+1]*sin(phase);
            }
            fprintf(file, "%f %f\n", time[i]/86400, mma[i] - model);
        }
        fclose(file);
    }

outfile_open_error:
    free(freq_labels);
    free(freqs);
    free(freq_amplitudes);
    free(freq_modes);
no_freqs_error:
    free(time);
    free(mma);
    free(err);
load_failed_error:
    return ret;
}

/*
 * Calculate the DFT of the BJD/mma data in tsFile between minUHz and maxUHz in increments dUHz
 * Output is saved to outFile with frequencies in uHz.
 * If freqFile is non-NULL, the data is first prewhitened with the contained frequencies
 */
int dft_bjd(char *tsFile, double minUHz, double maxUHz, double dUHz, char *outFile, char *freqFile)
{
    int ret = 0;

    double *time, *mma, *err;
    size_t num_obs;

    // Prewhiten if a frequency file is given
    if (freqFile)
    {
        double *freqs, *freq_amplitudes;
        char **freq_labels;
        int *freq_modes;
        size_t num_freqs;
        if (load_and_fit_freqs(tsFile, freqFile,
                           &time, &mma, &err, &num_obs,
                           &freq_labels, &freqs, &freq_amplitudes, &freq_modes, &num_freqs))
        {
            error_jump(load_failed_error, ret, "Error processing data");
        }

        for (size_t i = 0; i < num_obs; i++)
            for (size_t j = 0; j < num_freqs; j++)
            {
                double phase = 2*M_PI*freqs[j]*time[i];
                mma[i] -= freq_amplitudes[2*j]*cos(phase);
                mma[i] -= freq_amplitudes[2*j+1]*sin(phase);
            }

        free(freq_labels);
        free(freqs);
        free(freq_amplitudes);
        free(freq_modes);
    }
    else if (load_tsfile(tsFile, &time, &mma, &err, &num_obs))
        error_jump(load_failed_error, ret, "Error processing data");

    // Calculate DFT
    int num_uhz = (int)((maxUHz - minUHz)/dUHz);
    double *dftfreq = (double *)malloc(num_uhz*sizeof(double));
    if (!dftfreq)
        error_jump(dftfreq_alloc_error, ret, "Error allocating dftfreq");

    double *dftampl = (double *)malloc(num_uhz*sizeof(double));
    if (!dftampl)
        error_jump(dftampl_alloc_error, ret, "Error allocating dftfreq");

    calculate_amplitude_spectrum(time, mma, num_obs, minUHz*1e-6, maxUHz*1e-6, dftfreq, dftampl, num_uhz);

    // Save output
    FILE *file = fopen(outFile, "w");
    if (!file)
        error_jump(outfile_open_error, ret, "Error opening output file %s", outFile);

    for (size_t i = 0; i < num_uhz; i++)
        fprintf(file, "%f %f\n", 1e6*dftfreq[i], dftampl[i]);

    fclose(file);
outfile_open_error:
    free(dftampl);
dftampl_alloc_error:
    free(dftfreq);
dftfreq_alloc_error:
    free(time);
    free(mma);
    free(err);
load_failed_error:
    return ret;
}

/*
 * Calculate the DFT window at freq (in uHz) for the BJD/mma data in tsFile between minUHz and maxUHz in increments dUHz
 * Output is saved to outFile with frequencies in uHz.
 */
int dft_window(char *tsFile, double windowFreq, double minUHz, double maxUHz, double dUHz, char *outFile)
{
    int ret = 0;

    double *time, *mma, *err;
    size_t num_obs;

    if (load_tsfile(tsFile, &time, &mma, &err, &num_obs))
        error_jump(load_failed_error, ret, "Error loading data");

    // Generate sinusoid
    for (size_t i = 0; i < num_obs; i++)
        mma[i] = sin(2*M_PI*windowFreq*1e-6*time[i]);

    // Calculate DFT
    int num_uhz = (int)((maxUHz - minUHz)/dUHz);
    double *dftfreq = (double *)malloc(num_uhz*sizeof(double));
    if (!dftfreq)
        error_jump(dftfreq_alloc_error, ret, "Error allocating dftfreq");

    double *dftampl = (double *)malloc(num_uhz*sizeof(double));
    if (!dftampl)
        error_jump(dftampl_alloc_error, ret, "Error allocating dftfreq");

    calculate_amplitude_spectrum(time, mma, num_obs, minUHz*1e-6, maxUHz*1e-6, dftfreq, dftampl, num_uhz);

    FILE *file = fopen(outFile, "w");
    if (!file)
        error_jump(outfile_open_error, ret, "Error opening output file %s", outFile);

    for (size_t i = 0; i < num_uhz; i++)
        fprintf(file, "%f %f\n", 1e6*dftfreq[i], dftampl[i]);

    fclose(file);

outfile_open_error:
    free(dftampl);
dftampl_alloc_error:
    free(dftfreq);
dftfreq_alloc_error:
    free(time);
    free(mma);
    free(err);
load_failed_error:
    return ret;
}

/*
 * Calculate the DFT of the BJD/mma data in tsFile between minUHz and maxUHz in increments dUHz
 * and find the frequency with the largest amplitude
 * Results are plotted in a pgplot window
 * If freqFile is non-NULL, the data is first prewhitened with the contained frequencies
 */
int find_max_freq(char *tsFile, char *freqFile, double minUHz, double maxUHz, double dUHz)
{
    int ret = 0;

    double *time, *mma, *err, *freqs, *freq_amplitudes;
    char **freq_labels;
    int *freq_modes;
    size_t num_obs, num_freqs;

    if (load_and_fit_freqs(tsFile, freqFile,
                           &time, &mma, &err, &num_obs,
                           &freq_labels, &freqs, &freq_amplitudes, &freq_modes, &num_freqs))
    {
        error_jump(load_failed_error, ret, "Error processing data");
    }

    // Prewhiten
    for (size_t i = 0; i < num_obs; i++)
        for (size_t j = 0; j < num_freqs; j++)
        {
            double phase = 2*M_PI*freqs[j]*time[i];
            mma[i] -= freq_amplitudes[2*j]*cos(phase);
            mma[i] -= freq_amplitudes[2*j+1]*sin(phase);
        }

    // Calculate DFT
    int num_uhz = (int)((maxUHz - minUHz)/dUHz);
    double *dftfreq = (double *)malloc(num_uhz*sizeof(double));
    if (!dftfreq)
        error_jump(dftfreq_alloc_error, ret, "Error allocating dftfreq");

    double *dftampl = (double *)malloc(num_uhz*sizeof(double));
    if (!dftampl)
        error_jump(dftampl_alloc_error, ret, "Error allocating dftfreq");

    calculate_amplitude_spectrum(time, mma, num_obs, minUHz*1e-6, maxUHz*1e-6, dftfreq, dftampl, num_uhz);

    double max_ampl = 0;
    double max_ampl_freq = 0;
    for (size_t i = 0; i < num_uhz; i++)
        if (dftampl[i] > max_ampl)
        {
            max_ampl = dftampl[i];
            max_ampl_freq = dftfreq[i];
        }

    printf("Max amplitude of %f at %f uHz\n", max_ampl, max_ampl_freq*1e6);

    float *pgfreq = cast_double_array_to_float(dftfreq, num_uhz);
    float *pgampl = cast_double_array_to_float(dftampl, num_uhz);

    // Display
    if (cpgopen("6/xs") <= 0)
        error_jump(pgplot_open_error, ret, "Unable to open PGPLOT window");

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

    cpgend();

pgplot_open_error:
    free(dftampl);
dftampl_alloc_error:
    free(dftfreq);
dftfreq_alloc_error:
    if (freqs)
    {
        free(freq_labels);
        free(freqs);
        free(freq_amplitudes);
        free(freq_modes);
    }
    free(time);
    free(mma);
    free(err);
load_failed_error:
    return ret;
}

static int step_freq(double *fit_freqs, double *fit_amplitudes, int num_freqs,
                     double *time, double *mma, double *err, int num_obs,
                     int vary_freq, double centerFreq, double stepSize, int numSteps,
                     double *best_freq, double *best_chi2)
{
    for (int j = -numSteps; j <= numSteps; j++)
    {
        fit_freqs[vary_freq] = centerFreq + j*stepSize;

        if (fit_sinusoids(time, mma, err, num_obs, fit_freqs, num_freqs, fit_amplitudes))
            continue;

        double chi2 = calculate_chi2(fit_freqs, fit_amplitudes, num_freqs, time, mma, err, num_obs);
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
    int ret = 0;

    double *time, *mma, *err, *init_freqs, *fit_amplitudes;
    char **freq_labels;
    int *freq_modes;
    size_t num_obs, num_freqs;

    if (load_and_fit_freqs(tsFile, freqFile,
                           &time, &mma, &err, &num_obs,
                           &freq_labels, &init_freqs, &fit_amplitudes, &freq_modes, &num_freqs))
    {
        error_jump(load_failed_error, ret, "Error processing data");
    }

    double *fit_freqs = (double *)malloc(num_freqs*sizeof(double));
    if (!fit_freqs)
        error_jump(fit_freqs_alloc_error, ret, "Error allocating fit_freqs");

    for (size_t j = 0; j < num_freqs; j++)
        fit_freqs[j] = init_freqs[j];

    // Vary frequency between x-20 .. x + 20 in 0.01 steps
    double coarse_step = 1.0/(24*3600);
    int coarse_step_count = 5;
    double fine_step = 0.01e-6;
    int fine_step_count = 10;

    double chi2 = calculate_chi2(init_freqs, fit_amplitudes, num_freqs, time, mma, err, num_obs);
    double initialchi2 = chi2;
    double lastouterchi2 = chi2;
    double lastchi2 = chi2;
    do
    {
        lastouterchi2 = chi2;
        for (size_t i = 0; i < num_freqs; i++)
        {
            // Only fit freqs with 3rd column >= 2
            if (freq_modes[i] < 2)
                continue;

            double best_chi2 = chi2;
            double best_freq = fit_freqs[i];

            // Start with a rough step
            printf("Freq %10.2f\n", 1e6*init_freqs[i]);
            int iterate_mode = 2;
            double step = coarse_step;
            int num_steps = coarse_step_count;
            printf("\tUsing %d coarse steps\n", num_steps);
            do
            {
                double last_best_freq = best_freq;
                if (step_freq(fit_freqs, fit_amplitudes, num_freqs,
                              time, mma, err, num_obs,
                              i, best_freq, step, num_steps,
                              &best_freq, &best_chi2))
                {
                    error_jump(variation_failed_error , ret, "Variation failed");
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
    for (size_t i = 0; i < num_freqs; i++)
        printf("%-4s %7.2f %d\n", freq_labels[i], 1e6*fit_freqs[i], freq_modes[i]);

variation_failed_error:
    free(fit_freqs);
fit_freqs_alloc_error:
    free(freq_labels);
    free(init_freqs);
    free(fit_amplitudes);
    free(freq_modes);
    free(time);
    free(mma);
    free(err);
load_failed_error:
    return ret;
}

int shuffle_dft(char *tsFile, char *freqFile, double minUHz, double maxUHz, double dUHz, char *outFile, size_t repeats)
{
    int ret = 0;

    double *time, *mma, *err, *freqs, *freq_amplitudes;
    size_t num_obs, num_freqs;
    char **freq_labels;
    int *freq_modes;

    if (load_and_fit_freqs(tsFile, freqFile,
                           &time, &mma, &err, &num_obs,
                           &freq_labels, &freqs, &freq_amplitudes, &freq_modes, &num_freqs))
    {
        error_jump(load_failed_error, ret, "Error processing data");
    }

    // Prewhiten
    for (size_t i = 0; i < num_obs; i++)
        for (size_t j = 0; j < num_freqs; j++)
        {
            double phase = 2*M_PI*freqs[j]*time[i];
            mma[i] -= freq_amplitudes[2*j]*cos(phase);
            mma[i] -= freq_amplitudes[2*j+1]*sin(phase);
        }

    int num_uhz = (int)((maxUHz - minUHz)/dUHz);
    double *dftfreq = (double *)malloc(num_uhz*sizeof(double));
    if (!dftfreq)
        error_jump(dftfreq_alloc_error, ret, "Error allocating dftfreq");

    double *dftampl = (double *)malloc(num_uhz*sizeof(double));
    if (!dftampl)
        error_jump(dftampl_alloc_error, ret, "Error allocating dftfreq");

    // TODO: Seed correctly
    uint32_t seed = 19937;
    random_generator *rand = random_create(seed);
    if (!rand)
        error_jump(rand_alloc_error, ret, "Error allocating random generator");

    FILE *file = fopen(outFile, "w");
    if (!file)
        error_jump(outfile_open_error, ret, "Error opening output file %s", outFile);

    // File header
    fprintf(file, "# Prewhitened and shuffled max DFT amplitudes\n");
    fprintf(file, "# Data: %s\n", tsFile);
    fprintf(file, "# Freqs: %s\n", freqFile);
    fprintf(file, "# Seed: %u\n", seed);
    fprintf(file, "# Freq, Max Ampl, Mean Ampl\n");

    // Calculate output
    for (size_t i = 0; i < repeats; i++)
    {
        printf("%zu of %zu\n", i + 1, repeats);
        shuffle_double_array(time, num_obs, rand);
        calculate_amplitude_spectrum(time, mma, num_obs, minUHz*1e-6, maxUHz*1e-6, dftfreq, dftampl, num_uhz);

        // Find mean and max intensity
        double maxfreq = 0;
        double meanampl  = 0;
        double maxampl = 0;
        for (size_t j = 0; j < num_uhz; j++)
        {
            meanampl += dftampl[j];
            if (dftampl[j] > maxampl)
            {
                maxampl = dftampl[j];
                maxfreq = dftfreq[j];
            }
        }
        meanampl /= num_uhz;
        fprintf(file, "%f %f %f\n", maxfreq, maxampl, meanampl);
    }
    fclose(file);

outfile_open_error:
    random_free(rand);
rand_alloc_error:
    free(dftampl);
dftampl_alloc_error:
    free(dftfreq);
dftfreq_alloc_error:
    free(time);
    free(mma);
    free(err);
    free(freq_labels);
    free(freqs);
    free(freq_amplitudes);
    free(freq_modes);
load_failed_error:
    return ret;
}

int monitor_phase_amplitude(char *ts_file, double base_uhz, size_t harmonic_count, double window_width)
{
    int ret = 0;

    double *time, *mma, *err;
    size_t num_obs;

    if (load_tsfile(ts_file, &time, &mma, &err, &num_obs))
        error_jump(load_failed_error, ret, "Error processing data");

    size_t num_freqs = harmonic_count + 1;
    double *freqs = calloc(num_freqs, sizeof(double));
    if (!freqs)
        error_jump(load_failed_error, ret, "Error allocating freqs array");

    double *freq_amplitudes = calloc(2*num_freqs, sizeof(double));
    if (!freq_amplitudes)
        error_jump(load_failed_error, ret, "Error allocating freq_amplitude array");

    double *phase_start = calloc(num_freqs, sizeof(double));
    if (!phase_start)
        error_jump(load_failed_error, ret, "Error allocating phase_start array");

    for (size_t i = 0; i <= harmonic_count; i++)
        freqs[i] = (i+1)*base_uhz*1e-6;

    for (size_t i = 0; i < num_obs-1; i++)
    {
        if (time[i] + window_width > time[num_obs-1])
            break;

        // Find the last observation within the time window
        size_t j = 1;
        while (j + 1 < num_obs - i && time[i + j] <= time[i] + window_width)
            j++;

        if (j < 100)
            continue;

        if (fit_sinusoids(&time[i], &mma[i], &err[i], j, freqs, num_freqs, freq_amplitudes))
            error_jump(fit_failed_error, ret, "Sinusoid fit failed");

        double mean_time =  (time[i + j - 1] + time[i])/86400/2;

        double fit = 0;
        for (size_t l = 0; l < num_freqs; l++)
        {
            // convert time from BJD to seconds
            double phase = 2*M_PI*freqs[l]*mean_time*86400;
            fit += freq_amplitudes[2*l]*cos(phase);
            fit += freq_amplitudes[2*l+1]*sin(phase);
        }

        printf("%f %f ", mean_time, fit);

        for (size_t k = 0; k < num_freqs; k++)
        {
            double a = freq_amplitudes[2*k+1];
            double b = freq_amplitudes[2*k];
            double amp = sqrt(a*a + b*b);
            double phase = atan2(b, a)/(2*M_PI) - phase_start[k];
            while (phase > 1) phase -= 1;
            while (phase < 0) phase += 1;

            if (phase > 0.5)
                phase -= 1;

            if (i == 0)
            {
                phase_start[k] = phase;
                phase = 0;
            }

            double time_offset = phase/(freqs[k]*60);
            printf("%f %f %f ", amp, phase, time_offset);
        }

        printf("\n");
    }
fit_failed_error:
load_failed_error:
    return ret;
}

typedef struct
{
    double min;
    double max;
    double step;
    uint8_t color;
} animated_window_search_region;

int animated_window(char *tsfile)
{
    int ret = 0;

    double display_freq_min = 500;
    double display_freq_max = 4000;
    double display_freq_step = 5;
    size_t window_length = 3600*3/30; // data points -> 3h
    size_t window_increment = 2; // Number of data points to advance the window each increment

    size_t display_freq_count = (size_t)((display_freq_max - display_freq_min)/display_freq_step);

    animated_window_search_region peak_freqs[] =
    {
        {850, 950, 0.5, 2},
        {850*2, 950*2, 1, 3},
        {850*3, 950*3, 1, 4}
    };
    size_t peak_freq_count = 3;

    // Determine the maximum number of frequency steps so we can preallocate
    // a single block of memory for all calculations
    size_t peak_freq_working_size = 0;
    for (size_t i = 0; i < peak_freq_count; i++)
        peak_freq_working_size = fmax(peak_freq_working_size,
                (size_t)((peak_freqs[i].max - peak_freqs[i].min)/peak_freqs[i].step));

    double *time, *mmi, *err;
    size_t num_obs;

    if (load_tsfile(tsfile, &time, &mmi, &err, &num_obs))
        error_jump(load_failed_error, ret, "Error processing data");

    // Number of time steps to slide the window through
    size_t num_steps = (num_obs - window_length - 1)/window_increment;

    double *display_dftfreq = (double *)malloc(display_freq_count*sizeof(double));
    if (!display_dftfreq)
        error_jump(display_dftfreq_alloc_error, ret, "Error allocating display_dftfreq");

    double *display_dftampl = (double *)malloc(display_freq_count*sizeof(double));
    if (!display_dftampl)
        error_jump(display_dftampl_alloc_error, ret, "Error allocating display_dftampl");

    double *peak_dftfreq = (double *)malloc(peak_freq_working_size*sizeof(double));
    if (!peak_dftfreq)
        error_jump(peak_dftfreq_alloc_error, ret, "Error allocating peak_dftfreq");

    double *peak_dftampl = (double *)malloc(peak_freq_working_size*sizeof(double));
    if (!peak_dftampl)
        error_jump(peak_dftampl_alloc_error, ret, "Error allocating peak_dftampl");

    float *peak_time = (float *)malloc(num_steps*sizeof(float));
    if (!peak_time)
        error_jump(peak_time_alloc_error, ret, "Error allocating peak_time");

    float *peak_freq = (float *)malloc(peak_freq_count*num_steps*sizeof(float));
    if (!peak_freq)
        error_jump(peak_freq_alloc_error, ret, "Error allocating peak_freq");

    double *mmi_windowed = (double *)malloc(window_length*sizeof(double));
    if (!mmi_windowed)
        error_jump(mmi_windowed_alloc_error, ret, "Error allocating mmi_windowed");

    if (cpgopen("6/xs") <= 0)
        error_jump(pgplot_open_error, ret, "Unable to open PGPLOT window");

    // 800 x 480
    cpgpap(9.41, 0.6);
    cpgask(0);
    cpgslw(3);
    cpgsfs(2);
    cpgscf(2);

    // Draw Labels
    cpgsvp(0.1, 0.9, 0.55, 0.9);
    cpgmtxt("t", 2, 0.5, 0.5, "Frequency (\\gmHz)");
    cpgmtxt("l", 2, 0.5, 0.5, "Amplitude (mma)");

    cpgsvp(0.1, 0.9, 0.075, 0.5);
    cpgmtxt("l", 2, 0.5, 0.5, "Frequency (\\gmHz)");
    cpgmtxt("b", 2, 0.5, 0.5, "Time (Hours)");

    for (size_t i = 0; i < num_steps; i++)
    {
        // Copy amplitudes into mmi_window
        for (size_t j = 0; j < window_length; j++)
            mmi_windowed[j] = mmi[window_increment*i+j];

        // Taper start and end
        size_t taper_count = window_length/2;
        for (size_t j = 0; j < taper_count; j++)
        {
            mmi_windowed[j] *= (j*1.0/taper_count);
            mmi_windowed[window_length - j - 1] *= (j*1.0/taper_count);
        }

        cpgbbuf();

        // Erase top and bottom displays
		cpgsfs(1);
		cpgsci(0);
        cpgsvp(0.1, 0.9, 0.5, 0.9);
        cpgswin(0, 1, 0, 1);
		cpgrect(0, 1, 0, 1);

        cpgsvp(0.1, 0.9, 0.075, 0.5);
        cpgswin(0, 1, 0, 1);
		cpgrect(0, 1, 0, 1);
        cpgsfs(2);

        // Draw borders
        cpgsci(1);
        cpgsvp(0.1, 0.9, 0.5, 0.9);
        cpgswin(display_freq_min, display_freq_max, 0, 60);
        cpgbox("bcstm", 0, 0, "bcnst", 0, 0);

        cpgsvp(0.1, 0.9, 0.075, 0.5);
        cpgswin(0, 11, peak_freqs[0].min, peak_freqs[0].max);
        cpgbox("bcstn", 0, 0, "bcnst", 0, 0);

        // Calculate and draw display region
        calculate_amplitude_spectrum(&time[window_increment*i], mmi_windowed, window_length,
                                     display_freq_min*1e-6, display_freq_max*1e-6,
                                     display_dftfreq, display_dftampl, display_freq_count);
        cpgsci(1);

        float *pgfreq = cast_double_array_to_float(display_dftfreq, display_freq_count);
        float *pgampl = cast_double_array_to_float(display_dftampl, display_freq_count);

        cpgsvp(0.1, 0.9, 0.5, 0.9);
        cpgswin(display_freq_min*1e-6, display_freq_max*1e-6, 0, 60);
        cpgline(display_freq_count, pgfreq, pgampl);

        // Calculate peak frequencies
        peak_time[i] = time[window_increment*i + window_length/2]/3600;
        for (size_t j = 0; j < peak_freq_count; j++)
        {
            cpgsci(peak_freqs[j].color);

            size_t count = (size_t)((peak_freqs[j].max - peak_freqs[j].min)/peak_freqs[j].step);
            // Find max intensity
            calculate_amplitude_spectrum(&time[window_increment*i], mmi_windowed, window_length,
                                         peak_freqs[j].min*1e-6, peak_freqs[j].max*1e-6,
                                         peak_dftfreq, peak_dftampl, count);

            double maxfreq = 0;
            double maxampl = 0;
            for (size_t k = 0; k < count; k++)
            {
                if (peak_dftampl[k] > maxampl)
                {
                    maxampl = peak_dftampl[k];
                    maxfreq = peak_dftfreq[k];
                }
            }

            peak_freq[j*num_steps + i] = maxfreq*1e6;

            cpgsvp(0.1, 0.9, 0.5, 0.9);
            cpgswin(display_freq_min*1e-6, display_freq_max*1e-6, 0, 60);
            cpgmove(maxfreq, 60);
            cpgdraw(maxfreq, 0);

            cpgsvp(0.1, 0.9, 0.075, 0.5);
            cpgswin(0, 11, peak_freqs[j].min, peak_freqs[j].max);
            cpgline(i, peak_time, &peak_freq[j*num_steps]);
        }
        cpgebuf();
    }
    cpgend();

pgplot_open_error:
    free(mmi_windowed);
mmi_windowed_alloc_error:
    free(peak_freq);
peak_freq_alloc_error:
    free(peak_time);
peak_time_alloc_error:
    free(peak_dftampl);
peak_dftampl_alloc_error:
    free(peak_dftfreq);
peak_dftfreq_alloc_error:
    free(display_dftampl);
display_dftampl_alloc_error:
    free(display_dftfreq);
display_dftfreq_alloc_error:
    free(time);
    free(mmi);
    free(err);
load_failed_error:
    return ret;
}
