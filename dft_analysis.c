/*
 * Copyright 2010, 2011, 2012 Paul Chote
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

struct ts_data
{
    // Timeseries data
    double *time;
    double *mma;
    double *err;
    size_t obs_count;

    // Frequency data
    double *freq;
    char **freq_label;
    int *freq_mode;
    size_t freq_count;

    // Fit amplitudes (2*freq_count elements)
    double *freq_amplitude;
};

static void ts_data_free(struct ts_data *data)
{
    free(data->time);
    free(data->mma);
    free(data->err);

    free(data->freq);
    free_2d_array(data->freq_label, data->freq_count);
    free(data->freq_mode);
    free(data->freq_amplitude);
    free(data);
}

static int ts_data_fit_sinusoids(struct ts_data *data)
{
    return fit_sinusoids(data->time, data->mma, data->err, data->obs_count, data->freq, data->freq_count, data->freq_amplitude);
}

static int load_tsfile(const char *ts_path, struct ts_data *data)
{
    int ret = 0;
    char linebuf[1024];

    FILE *file = fopen(ts_path, "r+");
    if (!file)
        error_jump(file_error, ret, "Unable to load %s", ts_path);

    size_t total = 0;
    while (fgets(linebuf, sizeof(linebuf) - 1, file))
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total++;
    rewind(file);

    data->time = calloc(total, sizeof(double));
    data->mma = calloc(total, sizeof(double));
    data->err = calloc(total, sizeof(double));
    if (!data->time || !data->mma || !data->err)
        error_jump(allocation_error, ret, "Allocation error");

    size_t count = 0;
    while (fgets(linebuf, sizeof(linebuf) - 1, file) && count < total)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        int read = sscanf(linebuf, "%lf %lf %lf\n", &data->time[count], &data->mma[count], &data->err[count]);

        // No error defined; set to unity
        if (read == 2)
            data->err[count] = 1;

        // Convert to seconds
        data->time[count] *= 86400;
        count++;
    }
    data->obs_count = count;
    fclose(file);

    return ret;

allocation_error:
    free(data->time); data->time = NULL;
    free(data->mma); data->mma = NULL;
    free(data->err); data->err = NULL;

file_error:
    data->obs_count = 0;
    return ret;
}

static int load_freqfile(const char *freq_path, struct ts_data *data)
{
    int ret = 0;
    char linebuf[1024];

    // Load freqs
    FILE *file = fopen(freq_path, "r+");
    if (!file)
        error_jump(file_error, ret, "Unable to load %s", freq_path);

    size_t total = 0;
    while (fgets(linebuf, sizeof(linebuf) - 1, file))
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total++;
    rewind(file);

    size_t count = 0;
    data->freq = calloc(total, sizeof(double));
    data->freq_amplitude = calloc(2*total, sizeof(double));
    data->freq_label = calloc(total, sizeof(char **));
    data->freq_mode = calloc(total, sizeof(int));
    if (!data->freq || !data->freq_amplitude || !data->freq_label || !data->freq_mode)
        error_jump(allocation_error, ret, "Allocation error");

    while (fgets(linebuf, sizeof(linebuf) - 1, file) && count < total)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        char label[1024];
        sscanf(linebuf, "%1024s %lf %d %*s\n", label, &data->freq[count], &data->freq_mode[count]);

        data->freq[count] *= 1e-6;
        if (data->freq_mode[count] != 0)
        {
            data->freq_label[count] = strdup(label);
            count++;
        }
    }

    data->freq_count = count;

    if (ts_data_fit_sinusoids(data))
        error_jump(fit_error, ret, "Amplitude fit failed");

    return ret;
fit_error:
allocation_error:
    free(data->freq); data->freq = NULL;
    free(data->freq_amplitude); data->freq_amplitude = NULL;
    free_2d_array(data->freq_label, total); data->freq_label = NULL;
    free(data->freq_mode); data->freq_mode = NULL;
file_error:
    data->freq_count = 1;
    return ret;
}

static double ts_data_chi2(struct ts_data *data)
{
    double chi2 = 0;
    for (size_t i = 0; i < data->obs_count; i++)
    {
        double model = 0;
        for (size_t j = 0; j < data->freq_count; j++)
        {
            double phase = 2*M_PI*data->freq[j]*data->time[i];
            model += data->freq_amplitude[2*j]*cos(phase);
            model += data->freq_amplitude[2*j+1]*sin(phase);
        }
        double t = (data->mma[i] - model)/data->err[i];
        chi2 += t*t;
    }
    return chi2;
}

static void ts_data_prewhiten(struct ts_data *data)
{
    for (size_t i = 0; i < data->obs_count; i++)
        for (size_t j = 0; j < data->freq_count; j++)
        {
            double phase = 2*M_PI*data->freq[j]*data->time[i];
            data->mma[i] -= data->freq_amplitude[2*j]*cos(phase);
            data->mma[i] -= data->freq_amplitude[2*j+1]*sin(phase);
        }
}

static void ts_data_print_table(struct ts_data *data)
{
    printf("# ID     Freq    Period    Amp    Phase \n");
    printf("#        (uHz)     (s)    (mma)   (deg)\n");
    for (size_t i = 0; i < data->freq_count; i++)
    {
        double a = data->freq_amplitude[2*i+1];
        double b = data->freq_amplitude[2*i];
        double amp = sqrt(a*a + b*b);
        double phase = atan2(b, a)*180/M_PI;
        while (phase > 360) phase -= 360;
        while (phase < 0) phase += 360;
        printf("# %-5s %7.2f %7.2f %7.2f %7.2f\n",
               data->freq_label[i], 1e6*data->freq[i], 1/data->freq[i], amp, phase);
    }

    double chi2 = ts_data_chi2(data);
    int dof = data->obs_count - 3*data->freq_count;
    printf("# Chi^2: %f / DOF: %d = %f\n", chi2, dof, chi2/dof);
}

struct ts_data *ts_data_load(const char *ts_path, const char *freq_path)
{
    int ret = 0;
    struct ts_data *data = calloc(1, sizeof(struct ts_data));
    if (!data)
        return NULL;

    if (load_tsfile(ts_path, data))
        error_jump(load_error, ret, "Error loading timeseries data", freq_path);

    if (freq_path && load_freqfile(freq_path, data))
        error_jump(load_error, ret, "Error loading frequency data", freq_path);

    ts_data_print_table(data);

    return data;
load_error:
    ts_data_free(data);
    return NULL;
}

/*
 * Fit the frequencies defined in freqFile to the data in tsFile
 * and output a model lightcurve between startTime and endTime with increments of dt to modelFile
 * Residuals are saved to residualsFile if it is non-NULL
 */
int model_fit(char *ts_path, char *freq_path, double start_time, double end_time, double dt, char *model_path, char *residuals_path)
{
    int ret = 0;

    struct ts_data *data = ts_data_load(ts_path, freq_path);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    // Output model curve
    FILE *file = fopen(model_path, "w");
    if (!file)
        error_jump(outfile_open_error, ret, "Error opening output file %s", model_path);

    for (double t = start_time; t <= end_time; t += dt)
    {
        double fit = 0;
        for (size_t j = 0; j < data->freq_count; j++)
        {
            // convert time from BJD to seconds
            double phase = 2*M_PI*data->freq[j]*t*86400;
            fit += data->freq_amplitude[2*j]*cos(phase);
            fit += data->freq_amplitude[2*j+1]*sin(phase);
        }
        fprintf(file, "%f %f\n", t, fit);
    }
    fclose(file);

    // Output residuals
    if (residuals_path)
    {
        file = fopen(residuals_path, "w");
        if (!file)
            error_jump(outfile_open_error, ret, "Error opening output file %s", residuals_path);

        for (size_t i = 0; i < data->obs_count; i++)
        {
            double model = 0;
            for (size_t j = 0; j < data->freq_count; j++)
            {
                double phase = 2*M_PI*data->freq[j]*data->time[i];
                model += data->freq_amplitude[2*j]*cos(phase);
                model += data->freq_amplitude[2*j+1]*sin(phase);
            }
            fprintf(file, "%f %f\n", data->time[i]/86400, data->mma[i] - model);
        }
        fclose(file);
    }

outfile_open_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}

/*
 * Calculate the DFT of the BJD/mma data in tsFile between minUHz and maxUHz in increments dUHz
 * Output is saved to outFile with frequencies in uHz.
 * If freqFile is non-NULL, the data is first prewhitened with the contained frequencies
 */
int dft_bjd(char *ts_path, double min_uhz, double max_uhz, double d_uhz, char *out_path, char *freq_path)
{
    int ret = 0;
    struct ts_data *data = ts_data_load(ts_path, freq_path);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    // Prewhiten if a frequency file is given
    if (data->freq_count)
        ts_data_prewhiten(data);

    // Calculate DFT
    size_t dft_count = (size_t)((max_uhz - min_uhz)/d_uhz);
    double *dft_freq = calloc(dft_count, sizeof(double));
    if (!dft_freq)
        error_jump(dftfreq_alloc_error, ret, "Error allocating dft_freq");

    double *dft_ampl = calloc(dft_count, sizeof(double));
    if (!dft_ampl)
        error_jump(dftampl_alloc_error, ret, "Error allocating dft_ampl");

    calculate_amplitude_spectrum(data->time, data->mma, data->obs_count, min_uhz*1e-6, max_uhz*1e-6, dft_freq, dft_ampl, dft_count);

    // Save output
    FILE *file = fopen(out_path, "w");
    if (!file)
        error_jump(outfile_open_error, ret, "Error opening output file %s", out_path);

    for (size_t i = 0; i < dft_count; i++)
        fprintf(file, "%f %f\n", 1e6*dft_freq[i], dft_ampl[i]);

    fclose(file);
outfile_open_error:
    free(dft_ampl);
dftampl_alloc_error:
    free(dft_freq);
dftfreq_alloc_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}

/*
 * Calculate the DFT window at freq (in uHz) for the BJD/mma data in tsFile between minUHz and maxUHz in increments dUHz
 * Output is saved to outFile with frequencies in uHz.
 */
int dft_window(char *ts_path, double window_freq, double min_uhz, double max_uhz, double d_uhz, char *out_path)
{
    int ret = 0;
    struct ts_data *data = ts_data_load(ts_path, NULL);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    // Generate sinusoid
    for (size_t i = 0; i < data->obs_count; i++)
        data->mma[i] = sin(2*M_PI*window_freq*1e-6*data->time[i]);

    // Calculate DFT
    size_t dft_count = (size_t)((max_uhz - min_uhz)/d_uhz);
    double *dft_freq = calloc(dft_count, sizeof(double));
    if (!dft_freq)
        error_jump(dftfreq_alloc_error, ret, "Error allocating dft_freq");

    double *dft_ampl = calloc(dft_count, sizeof(double));
    if (!dft_ampl)
        error_jump(dftampl_alloc_error, ret, "Error allocating dft_ampl");

    calculate_amplitude_spectrum(data->time, data->mma, data->obs_count, min_uhz*1e-6, max_uhz*1e-6, dft_freq, dft_ampl, dft_count);

    FILE *file = fopen(out_path, "w");
    if (!file)
        error_jump(outfile_open_error, ret, "Error opening output file %s", out_path);

    for (size_t i = 0; i < dft_count; i++)
        fprintf(file, "%f %f\n", 1e6*dft_freq[i], dft_ampl[i]);

    fclose(file);

outfile_open_error:
    free(dft_ampl);
dftampl_alloc_error:
    free(dft_freq);
dftfreq_alloc_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}

/*
 * Calculate the DFT of the BJD/mma data in tsFile between minUHz and maxUHz in increments dUHz
 * and find the frequency with the largest amplitude
 * Results are plotted in a pgplot window
 * If freqFile is non-NULL, the data is first prewhitened with the contained frequencies
 */
int find_max_freq(char *ts_path, char *freq_path, double min_uhz, double max_uhz, double d_uhz)
{
    int ret = 0;
    struct ts_data *data = ts_data_load(ts_path, freq_path);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    // Prewhiten
    ts_data_prewhiten(data);

    // Calculate DFT
    size_t dft_count = (size_t)((max_uhz - min_uhz)/d_uhz);
    double *dft_freq = calloc(dft_count, sizeof(double));
    if (!dft_freq)
        error_jump(dftfreq_alloc_error, ret, "Error allocating dft_freq");

    double *dft_ampl = calloc(dft_count, sizeof(double));
    if (!dft_ampl)
        error_jump(dftampl_alloc_error, ret, "Error allocating dft_ampl");

    calculate_amplitude_spectrum(data->time, data->mma, data->obs_count, min_uhz*1e-6, max_uhz*1e-6, dft_freq, dft_ampl, dft_count);

    double ampl_max = 0;
    double freq = 0;
    for (size_t i = 0; i < dft_count; i++)
        if (dft_ampl[i] > ampl_max)
        {
            ampl_max = dft_ampl[i];
            freq = dft_freq[i];
        }

    printf("Max amplitude of %f at %f uHz\n", ampl_max, freq*1e6);

    float *pgfreq = cast_double_array_to_float(dft_freq, dft_count);
    float *pgampl = cast_double_array_to_float(dft_ampl, dft_count);

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
    cpgswin(min_uhz, max_uhz, 0, ampl_max*1.1);
    cpgbox("bcstn", 0, 0, "bcnst", 0, 0);

    cpgswin(min_uhz*1e-6, max_uhz*1e-6, 0, ampl_max*1.1);

    cpgsci(12);
    cpgline(dft_count, pgfreq, pgampl);

    cpgsci(2);
    cpgmove(freq, ampl_max);
    cpgdraw(freq, 0);
    cpgsci(1);

    cpgmtxt("b", 2.5, 0.5, 0.5, "Frequency (\\gmHz)");
    cpgmtxt("l", 2, 0.5, 0.5, "Amplitude (mma)");

    cpgsch(1.25);
    cpgmtxt("t", 1.5, 0.5, 0.5, "Amplitude Spectrum");

    cpgend();

pgplot_open_error:
    free(dft_ampl);
dftampl_alloc_error:
    free(dft_freq);
dftfreq_alloc_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}

int nonlinear_fit(char *ts_path, char *freq_path)
{
    int ret = 0;
    struct ts_data *data = ts_data_load(ts_path, freq_path);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    double *initial_freq = calloc(data->freq_count, sizeof(double));
    if (!initial_freq)
        error_jump(allocation_error, ret, "Error allocating fit_freqs");

    for (size_t j = 0; j < data->freq_count; j++)
        initial_freq[j] = data->freq[j];

    // Vary frequency between x-20 .. x + 20 in 0.01 steps
    double coarse_step = 1.0/(24*3600);
    int32_t coarse_step_count = 5;
    double fine_step = 0.01e-6;
    int32_t fine_step_count = 10;

    double chi2 = ts_data_chi2(data);
    double initialchi2 = chi2;
    double lastouterchi2 = chi2;
    double lastchi2 = chi2;
    do
    {
        lastouterchi2 = chi2;
        for (size_t i = 0; i < data->freq_count; i++)
        {
            printf("Varying freq %zu %g\n", i, data->freq[i]);
            // Only fit freqs with 3rd column >= 2
            if (data->freq_mode[i] < 2)
                continue;

            double best_chi2 = chi2;
            double best_freq = data->freq[i];

            // Start with a rough step
            printf("Freq %10.2f\n", 1e6*data->freq[i]);
            uint8_t iterate_mode = 2;
            double step = coarse_step;
            int32_t step_count = coarse_step_count;
            printf("\tUsing %d coarse steps\n", step_count);
            do
            {
                double last_best_freq = best_freq;
                for (int32_t j = -step_count; j <= step_count; j++)
                {
                    data->freq[i] = best_freq + j*step;
                    if (ts_data_fit_sinusoids(data))
                        continue;

                    double chi2 = ts_data_chi2(data);
                    if (chi2 < best_chi2)
                    {
                        best_chi2 = chi2;
                        best_freq = data->freq[i];
                    }
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
                        step_count = fine_step_count;
                    }
                }
            }
            while (iterate_mode);

            data->freq[i] = best_freq;
            chi2 = lastchi2;
        }
    } while (chi2 < lastouterchi2);

    printf("%f -> %f (%f)\n", initialchi2, chi2, initialchi2 - chi2);
    for (size_t i = 0; i < data->freq_count; i++)
        printf("%-4s %7.2f %d\n", data->freq_label[i], 1e6*data->freq[i], data->freq_mode[i]);

    free(initial_freq);
allocation_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}

int shuffle_dft(char *ts_path, char *freq_path, double min_uhz, double max_uhz, double d_uhz, char *out_path, size_t repeats)
{
    int ret = 0;
    struct ts_data *data = ts_data_load(ts_path, freq_path);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    // Prewhiten
    ts_data_prewhiten(data);

    size_t dft_count = (size_t)((max_uhz - min_uhz)/d_uhz);
    double *dft_freq = calloc(dft_count, sizeof(double));
    if (!dft_freq)
        error_jump(dftfreq_alloc_error, ret, "Error allocating dft_freq");

    double *dft_ampl = calloc(dft_count, sizeof(double));
    if (!dft_ampl)
        error_jump(dftampl_alloc_error, ret, "Error allocating dft_freq");

    // TODO: Seed correctly
    uint32_t seed = 19937;
    random_generator *rand = random_create(seed);
    if (!rand)
        error_jump(rand_alloc_error, ret, "Error allocating random generator");

    FILE *file = fopen(out_path, "w");
    if (!file)
        error_jump(outfile_open_error, ret, "Error opening output file %s", out_path);

    // File header
    fprintf(file, "# Prewhitened and shuffled max DFT amplitudes\n");
    fprintf(file, "# Data: %s\n", ts_path);
    fprintf(file, "# Freqs: %s\n", freq_path);
    fprintf(file, "# Seed: %u\n", seed);
    fprintf(file, "# Freq, Max Ampl, Mean Ampl\n");

    // Calculate output
    for (size_t i = 0; i < repeats; i++)
    {
        printf("%zu of %zu\n", i + 1, repeats);
        random_shuffle_double_array(rand, data->time, data->obs_count);
        calculate_amplitude_spectrum(data->time, data->mma, data->obs_count, min_uhz*1e-6, max_uhz*1e-6, dft_freq, dft_ampl, dft_count);

        // Find mean and max intensity
        double freq = 0;
        double ampl_mean = 0;
        double ampl_max = 0;
        for (size_t j = 0; j < dft_count; j++)
        {
            ampl_mean += dft_ampl[j];
            if (dft_ampl[j] > ampl_max)
            {
                ampl_max = dft_ampl[j];
                freq = dft_freq[j];
            }
        }
        ampl_mean /= dft_count;
        fprintf(file, "%f %f %f\n", freq, ampl_max, ampl_mean);
    }
    fclose(file);

outfile_open_error:
    random_free(rand);
rand_alloc_error:
    free(dft_ampl);
dftampl_alloc_error:
    free(dft_freq);
dftfreq_alloc_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}

int monitor_phase_amplitude(char *ts_path, double base_uhz, size_t harmonic_count, double window_width)
{
    int ret = 0;
    struct ts_data *data = calloc(1, sizeof(struct ts_data));
    if (!data)
        return 1;

    if (load_tsfile(ts_path, data))
        error_jump(load_failed_error, ret, "Error loading timeseries data");

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

    for (size_t i = 0; i < data->obs_count-1; i++)
    {
        if (data->time[i] + window_width > data->time[data->obs_count-1])
            break;

        // Find the last observation within the time window
        size_t j = 1;
        while (j + 1 < data->obs_count - i && data->time[i + j] <= data->time[i] + window_width)
            j++;

        if (j < 100)
            continue;

        if (fit_sinusoids(&data->time[i], &data->mma[i], &data->err[i], j, freqs, num_freqs, freq_amplitudes))
            error_jump(fit_failed_error, ret, "Sinusoid fit failed");

        double mean_time =  (data->time[i + j - 1] + data->time[i])/86400/2;

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

int animated_window(char *ts_path)
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

    struct ts_data *data = calloc(1, sizeof(struct ts_data));
    if (!data)
        return 1;

    // Number of time steps to slide the window through
    size_t num_steps = (data->obs_count - window_length - 1)/window_increment;

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
            mmi_windowed[j] = data->mma[window_increment*i+j];

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
        calculate_amplitude_spectrum(&data->time[window_increment*i], mmi_windowed, window_length,
                                     display_freq_min*1e-6, display_freq_max*1e-6,
                                     display_dftfreq, display_dftampl, display_freq_count);
        cpgsci(1);

        float *pgfreq = cast_double_array_to_float(display_dftfreq, display_freq_count);
        float *pgampl = cast_double_array_to_float(display_dftampl, display_freq_count);

        cpgsvp(0.1, 0.9, 0.5, 0.9);
        cpgswin(display_freq_min*1e-6, display_freq_max*1e-6, 0, 60);
        cpgline(display_freq_count, pgfreq, pgampl);

        // Calculate peak frequencies
        peak_time[i] = data->time[window_increment*i + window_length/2]/3600;
        for (size_t j = 0; j < peak_freq_count; j++)
        {
            cpgsci(peak_freqs[j].color);

            size_t count = (size_t)((peak_freqs[j].max - peak_freqs[j].min)/peak_freqs[j].step);
            // Find max intensity
            calculate_amplitude_spectrum(&data->time[window_increment*i], mmi_windowed, window_length,
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
    ts_data_free(data);
    return ret;
}
