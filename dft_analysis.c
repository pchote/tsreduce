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
#include <float.h>

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
    data->freq_count = 0;
    return ret;
}

static int load_freq_harmonics(double base_freq, size_t freq_count, struct ts_data *data)
{
    int ret = 0;

    data->freq = calloc(freq_count, sizeof(double));
    data->freq_amplitude = calloc(2*freq_count, sizeof(double));
    data->freq_label = calloc(freq_count, sizeof(char **));
    data->freq_mode = calloc(freq_count, sizeof(int));
    if (!data->freq || !data->freq_amplitude || !data->freq_label || !data->freq_mode)
        error_jump(allocation_error, ret, "Allocation error");

    for (size_t i = 0; i < freq_count; i++)
    {
        char buf[32];
        snprintf(buf, 32, "%zu", i);
        data->freq[i] = (i+1)*base_freq;
        data->freq_label[i] = strdup(buf);
        data->freq_mode[i] = 1;
    }
    data->freq_count = freq_count;

    if (ts_data_fit_sinusoids(data))
        error_jump(fit_error, ret, "Amplitude fit failed");

    return ret;
fit_error:
allocation_error:
    free(data->freq); data->freq = NULL;
    free(data->freq_amplitude); data->freq_amplitude = NULL;
    free_2d_array(data->freq_label, freq_count); data->freq_label = NULL;
    free(data->freq_mode); data->freq_mode = NULL;
    data->freq_count = 0;
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
        error_jump(load_error, ret, "Error loading timeseries data", ts_path);

    if (freq_path && load_freqfile(freq_path, data))
        error_jump(load_error, ret, "Error loading frequency data", freq_path);

    ts_data_print_table(data);

    return data;
load_error:
    ts_data_free(data);
    return NULL;
}

struct ts_data *ts_data_load_harmonics(const char *ts_path, double base_freq, size_t freq_count)
{
    int ret = 0;
    struct ts_data *data = calloc(1, sizeof(struct ts_data));
    if (!data)
        return NULL;

    if (load_tsfile(ts_path, data))
        error_jump(load_error, ret, "Error loading timeseries data", ts_path);

    if (load_freq_harmonics(base_freq, freq_count, data))
        error_jump(load_error, ret, "Error loading frequency data (base: %g, count: %zu)", base_freq, freq_count);

    ts_data_print_table(data);

    return data;
load_error:
    ts_data_free(data);
    return NULL;
}

struct histogram_data
{
    double *center;
    size_t *count;
    size_t bin_count;

    double minimum;
    double maximum;
    double bin_width;
    size_t excluded_count;
};

void histogram_free(struct histogram_data *hd)
{
    free(hd->center);
    free(hd->count);
    free(hd);
}

struct histogram_data *histogram_create(double *data, size_t count, double minimum, double maximum, size_t bin_count)
{
    struct histogram_data *hd = calloc(1, sizeof(struct histogram_data));
    if (!hd)
        return NULL;

    hd->minimum = minimum;
    hd->maximum = maximum;
    hd->bin_count = bin_count;
    hd->bin_width = (maximum - minimum)/bin_count;

    hd->center = calloc(bin_count, sizeof(double));
    hd->count = calloc(bin_count, sizeof(size_t));
    if (!hd->center || !hd->count)
    {
        histogram_free(hd);
        return NULL;
    }

    for (size_t i = 0; i < bin_count; i++)
        hd->center[i] = minimum + (i + 0.5)*hd->bin_width;

    for (size_t i = 0; i < count; i++)
    {
        double bin = (data[i] - minimum)/hd->bin_width;
        if (bin < 0 || bin >= bin_count)
            hd->excluded_count++;
        else
            hd->count[(size_t)bin]++;
    }

    return hd;
}

int histogram_fit_gaussian(struct histogram_data *hd, double fit_min, double fit_max, double params[3])
{
    int ret = 0;

    double *center = calloc(hd->bin_count, sizeof(double));
    double *count = calloc(hd->bin_count, sizeof(double));
    if (!center || !count)
        error_jump(error, ret, "Histogram fit allocation failed");

    printf("Source histogram:\n");
    size_t n = 0;
    for (size_t i = 0; i < hd->bin_count; i++)
    {
        printf("%g %zu\n", hd->center[i], hd->count[i]);
        // Discard zero values (which break the fit) and bins outside the requested fit range
        if (hd->count[i] == 0 || hd->center[i] - hd->bin_width/2 < fit_min || hd->center[i] + hd->bin_width/2 > fit_max)
            continue;

        // Transform center of range to 0
        center[n] = hd->center[i] - (fit_min + fit_max) / 2;
        count[n] = hd->count[i];
        n++;
    }

    printf("Fit histogram:\n");
    for (size_t i = 0; i < n; i++)
        printf("%g %g\n", center[i], count[i]);

    if (fit_gaussian(center, count, n, params))
        error_jump(error, ret, "Fit failed");

    // Correct offset
    params[1] += (fit_min + fit_max) / 2;
error:
    free(count);
    free(center);
    return ret;
}

int histogram_export(struct histogram_data *hd, double gaussian_params[3], const char *path)
{
    int ret = 0;

    FILE *output = fopen(path, "w");
    if (!output)
        error_jump(file_error, ret, "Unable to open histogram for writing: %s", path);

    fprintf(output, "# Fit parameters: ampl: %g mu: %g sigma: %g\n", gaussian_params[0], gaussian_params[1], gaussian_params[2]);
    for (size_t i = 0; i < hd->bin_count; i++)
    {
        double arg = (hd->center[i] - gaussian_params[1])/gaussian_params[2];
        double fit = gaussian_params[0]*exp(-0.5*arg*arg);
        fprintf(output, "%g %g %zu %g %g\n", hd->center[i] - hd->bin_width/2, hd->center[i] + hd->bin_width/2, hd->count[i], hd->center[i], fit);
    }

    fclose(output);
file_error:
    return ret;
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

int monitor_phase_amplitude(char *ts_path, double base_uhz, size_t freq_count, double window_width)
{
    int ret = 0;
    struct ts_data *data = ts_data_load_harmonics(ts_path, base_uhz*1e-6, freq_count);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    double *phase_start = calloc(data->freq_count, sizeof(double));
    if (!phase_start)
        error_jump(load_failed_error, ret, "Error allocating phase_start array");

    for (size_t i = 0; i < data->freq_count; i++)
        data->freq[i] = (i + 1)*base_uhz*1e-6;

    double *orig_time = data->time;
    double *orig_mma = data->mma;
    double *orig_err = data->err;
    size_t orig_obs_count = data->obs_count;

    for (size_t i = 0; i < orig_obs_count - 1; i++)
    {
        if (orig_time[i] + window_width > orig_time[orig_obs_count - 1])
            break;

        // Find the last observation within the time window
        size_t j = 1;
        while (j + 1 < orig_obs_count - i && orig_time[i + j] <= orig_time[i] + window_width)
            j++;

        if (j < 100)
            continue;

        data->time = &orig_time[i];
        data->mma = &orig_mma[i];
        data->err = &orig_err[i];
        data->obs_count = j;

        if (ts_data_fit_sinusoids(data))
            error_jump(fit_failed_error, ret, "Sinusoid fit failed");

        double mean_time =  (orig_time[i + j - 1] + orig_time[i])/86400/2;

        double fit = 0;
        for (size_t l = 0; l < data->freq_count; l++)
        {
            // convert time from BJD to seconds
            double phase = 2*M_PI*data->freq[l]*mean_time*86400;
            fit += data->freq_amplitude[2*l]*cos(phase);
            fit += data->freq_amplitude[2*l+1]*sin(phase);
        }

        printf("%f %f ", mean_time, fit);

        for (size_t k = 0; k < data->freq_count; k++)
        {
            double a = data->freq_amplitude[2*k+1];
            double b = data->freq_amplitude[2*k];
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

            double time_offset = phase/(data->freq[k]*60);
            printf("%f %f %f ", amp, phase, time_offset);
        }

        printf("\n");
    }

    data->time = orig_time;
    data->mma = orig_mma;
    data->err = orig_err;
    data->obs_count = orig_obs_count;
fit_failed_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}

int fit_baseline_polynomial(char *ts_path, size_t poly_degree)
{
    int ret = 0;
    struct ts_data *data = ts_data_load(ts_path, NULL);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    double *coeffs = calloc(poly_degree + 1, sizeof(double));
    if (!coeffs)
        error_jump(coeffs_alloc_error, ret, "Error allocating coeffs array");

    if (fit_polynomial(data->time, data->mma, data->err, data->obs_count, coeffs, poly_degree))
        error_jump(fit_failed_error, ret, "Polynomial fit failed");

    for (size_t i = 0; i < data->obs_count; i++)
    {
        double model = evaluate_polynomial(coeffs, poly_degree, data->time[i]);
        printf("%f %f %f\n", data->time[i]/86400, data->mma[i] - model, data->err[i]);
        fprintf(stderr, "%f %f\n", data->time[i]/86400, model);
    }
fit_failed_error:
    free(coeffs);
coeffs_alloc_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}

int fit_gwlib_freqshift(char *ts_path, double first_uhz, double second_uhz, size_t harmonic_count)
{
    int ret = 0;
    struct ts_data *data = ts_data_load_harmonics(ts_path, first_uhz*1e-6, harmonic_count);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    if (data->obs_count == 0)
        error_jump(no_obs_error, ret, "No observations in ts file");

    size_t orig_obs_count = data->obs_count;
    double *orig_time = data->time;
    double *orig_mma = data->mma;
    double *orig_err = data->err;
    for (size_t i = 2*data->freq_count; i < orig_obs_count - 2*data->freq_count; i++)
    {
        // Fit first half
        for (size_t j = 0; j <= harmonic_count; j++)
            data->freq[j] = (j+1)*first_uhz*1e-6;

        data->time = orig_time;
        data->mma = orig_mma;
        data->err = orig_err;
        data->obs_count = i;

        if (ts_data_fit_sinusoids(data))
        {
            fprintf(stderr, "%zu: first fit failed\n", i);
            continue;
        }

        data->obs_count = i;
        double first = ts_data_chi2(data);

        // Fit second half
        for (size_t j = 0; j <= harmonic_count; j++)
            data->freq[j] = (j+1)*second_uhz;

        data->time = orig_time + i;
        data->mma = orig_mma + i;
        data->err = orig_err + i;
        data->obs_count = orig_obs_count - i;

        if (ts_data_fit_sinusoids(data))
        {
            fprintf(stderr, "%zu: second fit failed\n", i);
            continue;
        }

        double second = ts_data_chi2(data);
        printf("%f %f %f %f\n", orig_time[i]/86400, first, second, first+second);
    }

    data->time = orig_time;
    data->mma = orig_mma;
    data->err = orig_err;
    data->obs_count = orig_obs_count;

no_obs_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}

static void step_freq_fit(struct ts_data *data, double freq_min, double freq_max, double freq_step,
            double *best_freq, double *best_chi2)
{

    // Step through fit frequencies
    double fit_freq = freq_min;

    *best_freq = 0;
    *best_chi2 = DBL_MAX;
    while (fit_freq <= freq_max)
    {
        for (size_t j = 0; j < data->freq_count; j++)
            data->freq[j] = 1e-6*(j+1)*fit_freq;

        if (ts_data_fit_sinusoids(data))
        {
            fprintf(stderr, "fit failed\n");
            continue;
        }
        double chi2 = ts_data_chi2(data);
        if (chi2 < *best_chi2)
        {
            *best_freq = fit_freq;
            *best_chi2 = chi2;
        }
        fit_freq += freq_step;
    }
}

int noise_histogram(const char *ts_path, double base_uhz, size_t freq_count,
                    double min_mma, double max_mma, size_t bin_count,
                    double fit_min_mma, double fit_max_mma, size_t randomize_count,
                    const char *output_prefix)
{
    int ret = 0;
    size_t output_path_len = strlen(output_prefix) + 10;
    char *output_path_buffer = malloc(output_path_len);
    if (!output_path_buffer)
        error_jump(load_failed_error, ret, "Allocation error");


    double freq_search_min = base_uhz - 10;
    double freq_search_max = base_uhz + 10;
    double freq_search_step = 0.05;

    struct ts_data *data = ts_data_load_harmonics(ts_path, 1e-6*base_uhz, freq_count);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    // Optimize fit
    {
        double old_chi2 = ts_data_chi2(data);
        double best_freq, best_chi2;
        step_freq_fit(data, freq_search_min, freq_search_max, freq_search_step, &best_freq, &best_chi2);
        printf("Optimized base freq from %g to %g: delta-chi2: %g\n", base_uhz, best_freq, old_chi2 - best_chi2);

        for (size_t j = 0; j < freq_count; j++)
            data->freq[j] = 1e-6*(j+1)*best_freq;

        if (ts_data_fit_sinusoids(data))
            error_jump(alloc_failed, ret, "Amplitude fit failed");
    }

    ts_data_prewhiten(data);


    // DFT Residuals
    {
        double max_uhz = 5000;
        double min_uhz = 0;
        double d_uhz = 1;

        size_t dft_count = (size_t)((max_uhz - min_uhz)/d_uhz);
        double *dft_freq = calloc(dft_count, sizeof(double));
        //if (!dft_freq)
        //    error_jump(dftfreq_alloc_error, ret, "Error allocating dft_freq");

        double *dft_ampl = calloc(dft_count, sizeof(double));
        //if (!dft_ampl)
        //    error_jump(dftampl_alloc_error, ret, "Error allocating dft_ampl");

        calculate_amplitude_spectrum(data->time, data->mma, data->obs_count, min_uhz*1e-6, max_uhz*1e-6, dft_freq, dft_ampl, dft_count);

        // Save output
        snprintf(output_path_buffer, output_path_len, "%s.dft", output_prefix);
        FILE *file = fopen(output_path_buffer, "w");
        //if (!file)
        //    error_jump(outfile_open_error, ret, "Error opening output file %s", dft_path);

        for (size_t i = 0; i < dft_count; i++)
            fprintf(file, "%f %f\n", 1e6*dft_freq[i], dft_ampl[i]);

        fclose(file);
    }

    // Calculate gaussian noise spectrum
    struct histogram_data *noise_histogram = histogram_create(data->mma, data->obs_count, min_mma, max_mma, bin_count);
    if (!noise_histogram)
        error_jump(noise_histogram_failed, ret, "Noise histogram failed");

    double noise_gaussian[3];
    if (histogram_fit_gaussian(noise_histogram, fit_min_mma, fit_max_mma, noise_gaussian))
        error_jump(noise_histogram_processing_failed, ret, "Noise histogram fit failed");

    snprintf(output_path_buffer, output_path_len, "%s.hist", output_prefix);
    if (histogram_export(noise_histogram, noise_gaussian, output_path_buffer))
        error_jump(noise_histogram_processing_failed, ret, "Noise histogram fit failed");

    snprintf(output_path_buffer, output_path_len, "%s.rand", output_prefix);
    FILE *randomized = fopen(output_path_buffer, "w");
    if (!randomized)
        error_jump(randomized_failed, ret, "Unable to open randomized data for writing: %s", output_path_buffer);

    uint32_t seed = time(NULL);
    random_generator *rand = random_create(seed);
    if (!rand)
        error_jump(rand_error, ret, "Error creating random generator");

    // Set timeseries noise column to its mean value
    double err_mean = 0;
    for (size_t i = 0; i < data->obs_count; i++)
        err_mean += data->err[i];
    err_mean /= data->obs_count;
    for (size_t i = 0; i < data->obs_count; i++)
        data->err[i] = err_mean;

    // Copy original freqs for generating synthesized data
    double *orig_freq = calloc(data->freq_count, sizeof(double));
    double *orig_ampl = calloc(2*data->freq_count, sizeof(double));
    if (!orig_freq || !orig_ampl)
        error_jump(freq_alloc_error, ret, "Error copying freqs");

    for (size_t j = 0; j < data->freq_count; j++)
    {
        orig_freq[j] = data->freq[j];
        orig_ampl[2*j] = data->freq_amplitude[2*j];
        orig_ampl[2*j+1] = data->freq_amplitude[2*j+1];
    }

    // Synthesize artificial data with random noise and find best-fit frequency
    double *fitted_freq = calloc(randomize_count, sizeof(double));
    if (!fitted_freq)
        error_jump(fitted_alloc_failed, ret, "Allocation failed");

    fprintf(randomized, "# Seed: %u\n", seed);
    for (size_t k = 0; k < randomize_count; k++)
    {
        if (k > 0 && !(k % 10))
            printf("%zu...\n", k);

        for (size_t i = 0; i < data->obs_count; i++)
        {
            data->mma[i] = random_normal(rand, noise_gaussian[1], noise_gaussian[2]);
            for (size_t j = 0; j < data->freq_count; j++)
            {
                double phase = 2*M_PI*orig_freq[j]*data->time[i];
                data->mma[i] += orig_ampl[2*j]*cos(phase);
                data->mma[i] += orig_ampl[2*j+1]*sin(phase);
            }
        }

        // Step through fit frequencies
        double best_freq, best_chi2;
        step_freq_fit(data, freq_search_min, freq_search_max, freq_search_step, &best_freq, &best_chi2);
        fprintf(randomized, "%zu %.2f %g\n", k, best_freq, best_chi2);
        fflush(randomized);

        fitted_freq[k] = best_freq;
    }

    struct histogram_data *freq_histogram = histogram_create(fitted_freq, randomize_count, freq_search_min, freq_search_max, (freq_search_max - freq_search_min) / freq_search_step);
    if (!freq_histogram)
        error_jump(freq_histogram_failed, ret, "Freq histogram failed");

    double freq_gaussian[3];
    if (histogram_fit_gaussian(freq_histogram, freq_search_min, freq_search_max, freq_gaussian))
        error_jump(freq_histogram_processing_failed, ret, "Freq histogram fit failed");

    snprintf(output_path_buffer, output_path_len, "%s.freqhist", output_prefix);
    if (histogram_export(freq_histogram, freq_gaussian, output_path_buffer))
        error_jump(freq_histogram_processing_failed, ret, "Freq histogram fit failed");

freq_histogram_processing_failed:
    histogram_free(freq_histogram);
freq_histogram_failed:
fitted_alloc_failed:
    free(orig_freq);
    free(orig_ampl);
freq_alloc_error:
    random_free(rand);
rand_error:
    fclose(randomized);
randomized_failed:
noise_histogram_processing_failed:
    histogram_free(noise_histogram);
noise_histogram_failed:
alloc_failed:
    ts_data_free(data);
load_failed_error:
    return ret;

}
