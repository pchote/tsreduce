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
 * Fit the frequencies defined in freqFile to the data in tsFile
 * and output a model lightcurve between startTime and endTime with increments of dt to modelFile
 * Residuals are saved to residualsFile if it is non-NULL
 */
int model_fit(char *tsFile, char *freqFile, double startTime, double endTime, double dt, char *modelFile, char *residualsFile)
{
    char linebuf[1024];
    FILE *file = fopen(tsFile, "r+");
    if (file == NULL)
        return error("Unable to open file: %s", tsFile);

    // Count the number of entries to allocate
    int total_obs = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total_obs++;
    rewind(file);

    int num_obs = 0;
    double *time = (double *)malloc(total_obs*sizeof(double));
    double *mmi = (double *)malloc(total_obs*sizeof(double));
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num_obs < total_obs)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        sscanf(linebuf, "%lf %lf\n", &time[num_obs], &mmi[num_obs]);

        // Convert to seconds
        time[num_obs] *= 86400;
        num_obs++;
    }
    fclose(file);
    printf("Read %d observations\n", num_obs);

    if (num_obs == 0)
    {
        free(time);
        free(mmi);
        return error("No observations found");
    }

    // Load freqs
    file = fopen(freqFile, "r+");
    if (file == NULL)
    {
        free(time);
        free(mmi);
        return error("Unable to open file: %s", freqFile);
    }

    int total_freqs = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total_freqs++;
    rewind(file);

    int num_freqs = 0;
    double *fit_freqs = (double *)malloc(total_freqs*sizeof(double));
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num_freqs < total_freqs)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        int use;
        sscanf(linebuf, "%*d %lf %d\n", &fit_freqs[num_freqs], &use);

        // Convert to Hz
        fit_freqs[num_freqs] *= 1e-6;
        if (use)
            num_freqs++;
    }
    fclose(file);

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

    // Calculate amplitude and phase for each freq
    printf("Freq (uHz) Period (s) Amp (mma) Phase (deg)\n");
    for (int i = 0; i < num_freqs; i++)
    {
        double a = fit_amplitudes[2*i+1];
        double b = fit_amplitudes[2*i];
        double amp = sqrt(a*a + b*b);
        double phase = atan2(b, a)*180/M_PI;
        while (phase > 360) phase -= 360;
        while (phase < 0) phase += 360;

        printf("%10.2f %10.2f %9.2f %11.2f\n", 1e6*fit_freqs[i], 1/fit_freqs[i], amp, phase);
    }

    // Output model curve
    file = fopen(modelFile, "w");
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
    char linebuf[1024];
    FILE *file = fopen(tsFile, "r+");
    if (file == NULL)
        return error("Unable to open file: %s", tsFile);

    // Count the number of entries to allocate
    int total_obs = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total_obs++;
    rewind(file);

    int num_obs = 0;
    double *time = (double *)malloc(total_obs*sizeof(double));
    double *mmi = (double *)malloc(total_obs*sizeof(double));
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num_obs < total_obs)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        sscanf(linebuf, "%lf %lf\n", &time[num_obs], &mmi[num_obs]);

        // Convert to seconds
        time[num_obs] *= 86400;
        num_obs++;
    }
    fclose(file);
    printf("Read %d observations\n", num_obs);

    if (num_obs == 0)
    {
        free(time);
        free(mmi);
        return error("No observations found");
    }

    // Load freqs
    if (freqFile)
    {
        file = fopen(freqFile, "r+");
        if (file == NULL)
        {
            free(time);
            free(mmi);
            return error("Unable to open file: %s", freqFile);
        }

        int total_freqs = 0;
        while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
            if (linebuf[0] != '#' && linebuf[0] != '\n')
                total_freqs++;
        rewind(file);

        int num_freqs = 0;
        double *fit_freqs = (double *)malloc(total_freqs*sizeof(double));
        while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num_freqs < total_freqs)
        {
            // Skip comment / empty lines
            if (linebuf[0] == '#' || linebuf[0] == '\n')
                continue;

            int use;
            sscanf(linebuf, "%*d %lf %d\n", &fit_freqs[num_freqs], &use);

            // Convert to Hz
            fit_freqs[num_freqs] *= 1e-6;
            if (use)
                num_freqs++;
        }
        fclose(file);

        printf("Using %d of %d freqs\n", num_freqs, total_freqs);
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

        // Calculate amplitude and phase for each freq
        printf("Freq (uHz) Period (s) Amp (mma) Phase (deg)\n");
        for (int i = 0; i < num_freqs; i++)
        {
            double a = fit_amplitudes[2*i+1];
            double b = fit_amplitudes[2*i];
            double amp = sqrt(a*a + b*b);
            double phase = atan2(b, a)*180/M_PI;
            while (phase > 360) phase -= 360;
            while (phase < 0) phase += 360;

            printf("%10.2f %10.2f %9.2f %11.2f\n", 1e6*fit_freqs[i], 1/fit_freqs[i], amp, phase);
        }

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
        fclose(file);
        free(fit_amplitudes);
        free(fit_freqs);
    }

    // Calculate DFT
    int num_uhz = (int)((maxUHz - minUHz)/dUHz);
    double *freq = (double *)malloc(num_uhz*sizeof(double));
    double *ampl = (double *)malloc(num_uhz*sizeof(double));

    calculate_amplitude_spectrum(minUHz*1e-6, maxUHz*1e-6, time, mmi, num_obs, freq, ampl, num_uhz);

    // Save output
    file = fopen(outFile, "w");
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
 * Calculate the DFT of the BJD/mma data in tsFile between minUHz and maxUHz in increments dUHz
 * and find the frequency with the largest amplitude
 * Results are plotted in a pgplot window
 * If freqFile is non-NULL, the data is first prewhitened with the contained frequencies
 */
int find_max_freq(char *tsFile, char *freqFile, double minUHz, double maxUHz, double dUHz)
{
    char linebuf[1024];
    FILE *file = fopen(tsFile, "r+");
    if (file == NULL)
        return error("Unable to open file: %s", tsFile);

    // Count the number of entries to allocate
    int total_obs = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total_obs++;
    rewind(file);

    int num_obs = 0;
    double *time = (double *)malloc(total_obs*sizeof(double));
    double *mmi = (double *)malloc(total_obs*sizeof(double));
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num_obs < total_obs)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        sscanf(linebuf, "%lf %lf\n", &time[num_obs], &mmi[num_obs]);

        // Convert to seconds
        time[num_obs] *= 86400;
        num_obs++;
    }
    fclose(file);
    printf("Read %d observations\n", num_obs);

    if (num_obs == 0)
    {
        free(time);
        free(mmi);
        return error("No observations found");
    }

    // Load freqs

    file = fopen(freqFile, "r+");
    if (file == NULL)
    {
        free(time);
        free(mmi);
        return error("Unable to open file: %s", freqFile);
    }

    int total_freqs = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total_freqs++;
    rewind(file);

    int num_freqs = 0;
    double *fit_freqs = (double *)malloc(total_freqs*sizeof(double));
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num_freqs < total_freqs)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        int use;
        sscanf(linebuf, "%*d %lf %d %*s\n", &fit_freqs[num_freqs], &use);

        // Convert to Hz
        fit_freqs[num_freqs] *= 1e-6;
        if (use)
            num_freqs++;
    }
    fclose(file);

    printf("Using %d of %d freqs\n", num_freqs, total_freqs);

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

    // Calculate amplitude and phase for each freq
    printf("Freq (uHz) Period (s) Amp (mma) Phase (deg)\n");
    for (int i = 0; i < num_freqs; i++)
    {
        double a = fit_amplitudes[2*i+1];
        double b = fit_amplitudes[2*i];
        double amp = sqrt(a*a + b*b);
        double phase = atan2(b, a)*180/M_PI;
        while (phase > 360) phase -= 360;
        while (phase < 0) phase += 360;

        printf("%10.2f %10.2f %9.2f %11.2f\n", 1e6*fit_freqs[i], 1/fit_freqs[i], amp, phase);
    }

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
    fclose(file);
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
