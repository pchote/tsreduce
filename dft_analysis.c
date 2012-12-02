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
