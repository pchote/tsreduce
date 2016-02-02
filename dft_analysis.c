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
	double bjda, bjdb;
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
	{
		// Load time header
		if (linebuf[0] == '#' && strncmp(linebuf, "# Reference time:", 17) == 0)
		{
			int bjda, bjdb;
	        int r = sscanf(linebuf, "# Reference time: %*d-%*d-%*d %*d:%*d:%*f UTC; %d.%d BJD\n", &bjda, &bjdb);
			if (r != 2)
		        error_jump(allocation_error, ret, "Invalid time header: `%s`", linebuf);

			data->bjda = bjda;
			data->bjdb = bjdb > 0 ? copysign(bjdb / pow(10, (int)log10(bjdb) + 1), bjda) : 0;
		}

        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total++;
	}
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

    // Normalize data to zero-mean
    double mean = 0;
    for (size_t i = 0; i < count; i++)
        mean += data->mma[i];

    mean /= count;
    for (size_t i = 0; i < count; i++)
        data->mma[i] -= mean;

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
    printf("# ID     Freq    Period      Amp        Phase \n");
    printf("#        (uHz)     (s)      (mma)       (rad)\n");
    for (size_t i = 0; i < data->freq_count; i++)
    {
        double a = data->freq_amplitude[2*i+1];
        double b = data->freq_amplitude[2*i];
        double amp = sqrt(a*a + b*b);
        double phase = atan2(b, a);
        while (phase > 2*M_PI) phase -= 2*M_PI;
        while (phase < 0) phase += 2*M_PI;
        printf("# %-5s %7.2f %7.2f %11.6f %11.6f\n",
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

int histogram_fit_gaussian(struct histogram_data *hd, bool fit_center, double step_sigma, double *center, double *amplitude, double *sigma)
{
    int ret = 0;

    double *x = calloc(hd->bin_count, sizeof(double));
    double *y = calloc(hd->bin_count, sizeof(double));
    if (!x || !y)
        error_jump(error, ret, "Histogram fit allocation failed");

    size_t n = 0;
    for (size_t i = 0; i < hd->bin_count; i++)
    {
        // Discard zero values (which break the fit) and bins outside the requested fit range
        if (hd->count[i] == 0)
            continue;

        x[n] = hd->center[i];
        y[n] = hd->count[i];
        n++;
    }

    // Fit with an analytic gaussian as a rough initial guess
    double params[3];
    if (fit_analytical_gaussian(x, y, n, params))
        error_jump(error, ret, "Analytic Fit failed");

    // Center histogram
    if (fit_center)
    {
        *center = params[1];
        for (size_t i = 0; i < n; i++)
            x[i] -= params[1];
    }
    else
        *center = 0;

    double min_sigma = fmax((int)(hd->bin_width / step_sigma), 1) * step_sigma;        
    double max_sigma = (int)(params[2] * 4 / step_sigma) * step_sigma;

    if (fit_gaussian(x, y, n, min_sigma, max_sigma, step_sigma, sigma, amplitude))
        error_jump(error, ret, "Fit failed");

error:
    free(x);
    free(y);
    return ret;
}

int histogram_export(struct histogram_data *hd, double center, double amplitude, double sigma, const char *path)
{
    int ret = 0;

    FILE *output = fopen(path, "w");
    if (!output)
        error_jump(file_error, ret, "Unable to open histogram for writing: %s", path);

    fprintf(output, "# Fit parameters: ampl: %g mu: %g sigma: %g\n", amplitude, center, sigma);

    for (size_t i = 0; i < hd->bin_count; i++)
    {
        double arg = (hd->center[i] - center)/sigma;
        double fit = amplitude*exp(-0.5*arg*arg);
        fprintf(output, "%g %g %zu %g %g\n", hd->center[i] - hd->bin_width/2, hd->center[i] + hd->bin_width/2, hd->count[i], hd->center[i], fit);
    }

    fclose(output);
file_error:
    return ret;
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

int optimize_frequency_fit(struct ts_data *data, double *initialchi2, double *finalchi2)
{
    // Vary frequency between x-20 .. x + 20 in 0.01 steps
    double coarse_step = 0.001e-6;
    int32_t coarse_step_count = 5;
    double fine_step = 0.0001e-6;
    int32_t fine_step_count = 10;

    double chi2 = ts_data_chi2(data);
    *initialchi2 = chi2;
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

	*finalchi2 = chi2;
	return 0;
}

int nonlinear_fit(char *ts_path, char *freq_path)
{
    int ret = 0;
    struct ts_data *data = ts_data_load(ts_path, freq_path);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

	double initialchi2, finalchi2;
	optimize_frequency_fit(data, &initialchi2, &finalchi2);

    printf("%f -> %f (%f)\n", initialchi2, finalchi2, initialchi2 - finalchi2);
    for (size_t i = 0; i < data->freq_count; i++)
        printf("%-4s %7.5f %d\n", data->freq_label[i], 1e6*data->freq[i], data->freq_mode[i]);

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

int noise_histogram(const char *ts_path, const char *freq_path,
                    double min_mma, double max_mma, size_t bin_count,
                    double fit_min_mma, double fit_max_mma, size_t randomize_count,
                    const char *output_prefix)
{
    int ret = 0;
    size_t output_path_len = strlen(output_prefix) + 20;
    char *output_path_buffer = malloc(output_path_len);
    if (!output_path_buffer)
        error_jump(load_failed_error, ret, "Allocation error");

    struct ts_data *data = ts_data_load(ts_path, freq_path);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    double max_uhz = 10000;
    double min_uhz = 0;
    double d_uhz = 1;

	printf("Calculating initial DFT...\n");

    // Base DFT
    {
        size_t dft_count = (size_t)((max_uhz - min_uhz)/d_uhz);
        double *dft_freq = calloc(dft_count, sizeof(double));
        double *dft_ampl = calloc(dft_count, sizeof(double));

        calculate_amplitude_spectrum(data->time, data->mma, data->obs_count, min_uhz*1e-6, max_uhz*1e-6, dft_freq, dft_ampl, dft_count);

        // Save output
        snprintf(output_path_buffer, output_path_len, "%s.dft", output_prefix);
        FILE *file = fopen(output_path_buffer, "w");

        for (size_t i = 0; i < dft_count; i++)
            fprintf(file, "%f %f\n", 1e6*dft_freq[i], dft_ampl[i]);

        fclose(file);
    }

    ts_data_prewhiten(data);
	printf("Calculating prewhitened DFT...\n");

    // Prewhitened DFT
    {
        size_t dft_count = (size_t)((max_uhz - min_uhz)/d_uhz);
        double *dft_freq = calloc(dft_count, sizeof(double));
        double *dft_ampl = calloc(dft_count, sizeof(double));

        calculate_amplitude_spectrum(data->time, data->mma, data->obs_count, min_uhz*1e-6, max_uhz*1e-6, dft_freq, dft_ampl, dft_count);

        // Save output
        snprintf(output_path_buffer, output_path_len, "%s.prewhitened", output_prefix);
        FILE *file = fopen(output_path_buffer, "w");

        for (size_t i = 0; i < dft_count; i++)
            fprintf(file, "%f %f\n", 1e6*dft_freq[i], dft_ampl[i]);

        fclose(file);
    }

    double *noise_samples = calloc(data->obs_count, sizeof(double));
    for (size_t i = 0; i < data->obs_count; i++)
        noise_samples[i] = data->mma[i];

    {
    	printf("Calculating noise...\n");
        // Calculate gaussian noise spectrum
        struct histogram_data *noise_histogram = histogram_create(noise_samples, data->obs_count, min_mma, max_mma, bin_count);
        if (!noise_histogram)
            error_jump(noise_histogram_failed, ret, "Noise histogram failed");

        double noise_center, noise_amplitude, noise_sigma;
        if (histogram_fit_gaussian(noise_histogram, false, 0.01, &noise_center, &noise_amplitude, &noise_sigma))
            error_jump(noise_histogram_processing_failed, ret, "Noise histogram fit failed");

        snprintf(output_path_buffer, output_path_len, "%s.hist", output_prefix);
        if (histogram_export(noise_histogram, noise_center, noise_amplitude, noise_sigma, output_path_buffer))
            error_jump(noise_histogram_processing_failed, ret, "Noise histogram fit failed");

noise_histogram_processing_failed:
            histogram_free(noise_histogram);
    }

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

    // Copy original amplitudes for generating synthesized data
    double *orig_ampl = calloc(2*data->freq_count, sizeof(double));
    if (!orig_ampl)
        error_jump(ampl_alloc_error, ret, "Error copying amplitudes");

	// Generate a data table for all the results
	double *calculated_amplitude = calloc(data->freq_count*randomize_count, sizeof(double));
	double *calculated_phase = calloc(data->freq_count*randomize_count, sizeof(double));

    // Synthesize artificial data with random noise and find best-fit phase
    fprintf(randomized, "# Files: %s %s\n", ts_path, freq_path);
    fprintf(randomized, "# Seed: %u\n", seed);
	fprintf(randomized, "#");
    for (size_t j = 0; j < data->freq_count; j++)
    {
        orig_ampl[2*j] = data->freq_amplitude[2*j];
        orig_ampl[2*j+1] = data->freq_amplitude[2*j+1];

		if (data->freq_mode[j] == 2)
        	fprintf(randomized, "      %8.3f    ", 1e6*data->freq[j]);
    }

	fprintf(randomized, "\n");

    for (size_t k = 0; k < randomize_count; k++)
    {
        if (k > 0 && !(k % 10))
            printf("%zu...\n", k);

        for (size_t i = 0; i < data->obs_count; i++)
        {
            data->mma[i] = noise_samples[random_uint32_max(rand, data->obs_count)];
            for (size_t j = 0; j < data->freq_count; j++)
            {
                double phase = 2*M_PI*data->freq[j]*data->time[i];
                data->mma[i] += orig_ampl[2*j]*cos(phase);
                data->mma[i] += orig_ampl[2*j+1]*sin(phase);
            }
        }

	    if (ts_data_fit_sinusoids(data))
	        error_jump(fit_error, ret, "Amplitude fit failed");

		// Print results as we progress
    	fprintf(randomized, " ");
	    for (size_t j = 0; j < data->freq_count; j++)
		{
			if (data->freq_mode[j] == 2)
			{
		        double a = data->freq_amplitude[2*j+1];
		        double b = data->freq_amplitude[2*j];
		        double phase = atan2(b, a);
		        double ampl = sqrt(a*a + b*b);
				if (phase < 0)
					phase += 2*M_PI;

				if (k > 0)
				{
					double last = calculated_phase[j*randomize_count + k - 1];
					double delta = phase - last;
					if (delta > M_PI)
						phase -= 2*M_PI;
					if (delta < -M_PI)
						phase += 2*M_PI;
				}

				calculated_amplitude[j*randomize_count + k] = ampl;
				calculated_phase[j*randomize_count + k] = phase;
	        	fprintf(randomized, " %8.5f %8.5f", ampl, phase);
			}
		}

		fprintf(randomized, "\n");
        fflush(randomized);
    }

	fprintf(randomized, "#");
    for (size_t j = 0; j < data->freq_count; j++)
	{
		// Use the guess to calculate a proper fit
	    struct histogram_data *phase_histogram = histogram_create(&calculated_phase[j*randomize_count], randomize_count, 0, 2*M_PI, 1000);
	    if (!phase_histogram)
			printf("Noise histogram failed\n");

        double phase_center, phase_amplitude, phase_sigma;
        if (histogram_fit_gaussian(phase_histogram, true, 0.0001, &phase_center, &phase_amplitude, &phase_sigma))
            error_jump(noise_histogram_processing_failed, ret, "Noise histogram fit failed");

        char foo[1024];        
        snprintf(foo, 1024, "%s.hist.%s.lo", output_prefix, data->freq_label[j]);
        histogram_export(phase_histogram, phase_center, phase_amplitude, phase_sigma, foo);
        histogram_free(phase_histogram);
    	fprintf(randomized, " %8.5f %8.5f", phase_sigma, 0.0f);
		
		// Output values to collated data files
		double value = atan2(orig_ampl[2*j], orig_ampl[2*j + 1]);
        while (value > 2*M_PI) value -= 2*M_PI;
        while (value < 0) value += 2*M_PI;

		char outname[1024];
		snprintf(outname, 1024, "%s.%s", data->freq_label[j], "phase");
	    FILE *out = fopen(outname, "a");
		fprintf(out, "%.8Lf %10.5f %10.5f %10.5f %s %s\n", (long double)data->bjda + (long double)data->bjdb, value, phase_center, phase_sigma, ts_path, freq_path);
		fclose(out);
	}
	fprintf(randomized, "\n");

fit_error:
    free(orig_ampl);
ampl_alloc_error:
    random_free(rand);
rand_error:
    fclose(randomized);
randomized_failed:
noise_histogram_failed:
    ts_data_free(data);
load_failed_error:
    return ret;
}

int o_minus_c(const char *data, const char *str_ref_bjd, long double period, long double phase_offset)
{
    int ret = 0;

    FILE *file = fopen(data, "r");
    if (!file)
        return error("Unable to open file: %s", data);

    // Count the number of entries to allocate
    char linebuf[1024];
    size_t total = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL)
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total++;

    rewind(file);

    long double *bjd = calloc(total, sizeof(long double));
    if (!bjd)
        error_jump(bjd_alloc_error, ret, "Error allocating bjdb array");

    double *phase = calloc(total, sizeof(double));
    if (!phase)
        error_jump(phase_alloc_error, ret, "Error allocating phase array");

    size_t num = 0;
    while (fgets(linebuf, sizeof(linebuf)-1, file) != NULL && num < total)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        int read = sscanf(linebuf, "%Lf %lf %*s\n", &bjd[num], &phase[num]);
        if (read != 2)
            error_jump(invalid_data_error, ret, "Invalid data line (%d entries)", read);

        num++;
    }
    fclose(file);

    long double ref_bjd;
    sscanf(str_ref_bjd, "%Lf", &ref_bjd);

    printf("# Reference BJD: %Lf\n", ref_bjd);
    printf("# Reference period: %Lf\n", period);
    printf("#        dBJD          cycles           dphase           dsec\n");
    for (size_t i = 0; i < total; i++)
    {
        long double dbjd = bjd[i] - ref_bjd;
        long double cycles = dbjd*86400.0L / period;

        //long double cyclesa = (dbjda * (long double)86400.0 / period);
        //long double cyclesb = (dbjdb * (long double)8.64e-5 / period);
        long double offset = cycles*2*M_PI - phase[i] + phase_offset;
        while (offset > M_PI)
            offset -= 2*M_PI;
        while (offset < -M_PI)
            offset += 2*M_PI;

        printf("%15.6Lf %15.6Lf %15.6Lf %15.6Lf\n", dbjd, cycles, offset, offset * period / (2 * M_PI));
    }

invalid_data_error:
    free(phase);
phase_alloc_error:
    free(bjd);
bjd_alloc_error:
    return ret;
}

int monitor_phase_amplitude(char *ts_path, double base_uhz, size_t freq_count, double window_width)
{
    int ret = 0;
    struct ts_data *data = ts_data_load_harmonics(ts_path, base_uhz*1e-6, freq_count);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    for (size_t i = 0; i < data->freq_count; i++)
        data->freq[i] = (i + 1)*base_uhz*1e-6;

    double *orig_time = data->time;
    double *orig_mma = data->mma;
    double *orig_err = data->err;
    size_t orig_obs_count = data->obs_count;

    struct mode
    {
        double amplitude;
        double phase_cycles;
        double phase_time;
        double dphase;
    };

    struct calculation
    {
        double time;
        double amplitude;
        struct mode modes[10];
    };

    struct calculation *calculation = calloc(orig_obs_count - 1, sizeof(struct calculation));
    size_t i = 0;
    for (; i < orig_obs_count - 1; i++)
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

        double mean_time = (orig_time[i + j - 1] + orig_time[i])/86400/2;

        double fit = 0;
        for (size_t l = 0; l < data->freq_count; l++)
        {
            // convert time from BJD to seconds
            double phase = 2*M_PI*data->freq[l]*mean_time*86400;
            fit += data->freq_amplitude[2*l]*cos(phase);
            fit += data->freq_amplitude[2*l+1]*sin(phase);
        }

        calculation[i].time = mean_time;
        calculation[i].amplitude = fit;

        for (size_t k = 0; k < data->freq_count; k++)
        {
            double a = data->freq_amplitude[2*k+1];
            double b = data->freq_amplitude[2*k];
            double amp = sqrt(a*a + b*b);
            double phase = atan2(b, a)/(2*M_PI);

            if (i > 0)
            {
                while (phase - calculation[i-1].modes[k].phase_cycles > 0.5) phase -= 1;
                while (phase - calculation[i-1].modes[k].phase_cycles < 0.5) phase += 1;
            }

            if (phase > 0.5)
                phase -= 1;

            if (phase > 0 && mean_time > 0.376733 && k == 3)
                phase -= 1;

            calculation[i].modes[k].amplitude = amp;
            calculation[i].modes[k].phase_cycles = phase;
            calculation[i].modes[k].phase_time = phase/(data->freq[k]*60);
        }
    }

    for (size_t j = 0; j < i; j++)
    {
        printf("%f %f ", calculation[j].time, calculation[j].amplitude);
        for (size_t k = 0; k < data->freq_count; k++)
            printf("%f %f %f ", calculation[j].modes[k].amplitude,
                   calculation[j].modes[k].phase_cycles - calculation[0].modes[k].phase_cycles,
                   calculation[j].modes[k].phase_time - calculation[0].modes[k].phase_time);
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

int print_run_data(const char *ts_path, double exptime)
{
    int ret = 0;
    struct ts_data *data = ts_data_load(ts_path, NULL);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");

    long double bjd = (long double)data->bjda + (long double)data->bjdb;
    double run_length = data->time[data->obs_count - 1] - data->time[0];
    printf("%s & %.8Lf & %.2f & %zu & %.0f\\\\\n", ts_path, bjd, run_length / 86400, data->obs_count, data->obs_count * exptime / run_length * 100);
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
            double *best_freq, double *best_ampl, double *best_chi2)
{

    // Step through fit frequencies
    double fit_freq = freq_min;

    *best_freq = 0;
    *best_chi2 = DBL_MAX;
    while (fit_freq <= freq_max)
    {
        for (size_t j = 0; j < data->freq_count; j++)
            if (data->freq_mode[j] == 2)
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
            *best_ampl = sqrt(data->freq_amplitude[0]*data->freq_amplitude[0] + data->freq_amplitude[1]*data->freq_amplitude[1]);
            *best_chi2 = chi2;
        }
        fit_freq += freq_step;
    }
}

int gwlib_noise_histogram(const char *ts_path, const char *freq_path,
                    double min_mma, double max_mma, size_t bin_count, size_t randomize_count,
                    const char *output_prefix)
{
    int ret = 0;
    size_t output_path_len = strlen(output_prefix) + 20;
    char *output_path_buffer = malloc(output_path_len);
    if (!output_path_buffer)
        error_jump(load_failed_error, ret, "Allocation error");


    struct ts_data *data = ts_data_load(ts_path, freq_path);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");
    double base_uhz = data->freq[0] * 1e6;
    double freq_search_min = base_uhz - 5;
    double freq_search_max = base_uhz + 5;
    double freq_search_step = 0.1;

    // Optimize fit
    /*
    {
        double old_chi2 = ts_data_chi2(data);
        double best_freq, best_chi2, best_ampl;
        step_freq_fit(data, freq_search_min, freq_search_max, freq_search_step, &best_freq, &best_ampl, &best_chi2);
        printf("Optimized base freq from %g to %g: delta-chi2: %g\n", base_uhz, best_freq, old_chi2 - best_chi2);

        for (size_t j = 0; j < freq_count; j++)
            data->freq[j] = 1e-6*(j+1)*best_freq;

        if (ts_data_fit_sinusoids(data))
            error_jump(alloc_failed, ret, "Amplitude fit failed");
    }
    */

    double initial_freq = data->freq[0] * 1e6;
    double initial_ampl = sqrt(data->freq_amplitude[0]*data->freq_amplitude[0] + data->freq_amplitude[1]*data->freq_amplitude[1]);
    printf("Initial freq: %f\n", initial_freq);

    double ampl_search_min = initial_ampl - 5;
    double ampl_search_max = initial_ampl + 5;
    double ampl_search_step = 0.25;
    ts_data_print_table(data);

    // DFT
    {
        double max_uhz = 5000;
        double min_uhz = 0;
        double d_uhz = 1;

        size_t dft_count = (size_t)((max_uhz - min_uhz)/d_uhz);
        double *dft_freq = calloc(dft_count, sizeof(double));
        double *dft_ampl = calloc(dft_count, sizeof(double));

        calculate_amplitude_spectrum(data->time, data->mma, data->obs_count, min_uhz*1e-6, max_uhz*1e-6, dft_freq, dft_ampl, dft_count);

        // Save output
        snprintf(output_path_buffer, output_path_len, "%s.dft", output_prefix);
        FILE *file = fopen(output_path_buffer, "w");

        for (size_t i = 0; i < dft_count; i++)
            fprintf(file, "%f %f\n", 1e6*dft_freq[i], dft_ampl[i]);

        fclose(file);
    }

    ts_data_prewhiten(data);


    // DFT Residuals
    {
        double max_uhz = 5000;
        double min_uhz = 0;
        double d_uhz = 1;

        size_t dft_count = (size_t)((max_uhz - min_uhz)/d_uhz);
        double *dft_freq = calloc(dft_count, sizeof(double));
        double *dft_ampl = calloc(dft_count, sizeof(double));

        calculate_amplitude_spectrum(data->time, data->mma, data->obs_count, min_uhz*1e-6, max_uhz*1e-6, dft_freq, dft_ampl, dft_count);

        // Save output
        snprintf(output_path_buffer, output_path_len, "%s.prewhitened", output_prefix);
        FILE *file = fopen(output_path_buffer, "w");

        for (size_t i = 0; i < dft_count; i++)
            fprintf(file, "%f %f\n", 1e6*dft_freq[i], dft_ampl[i]);

        fclose(file);
    }

    // Calculate gaussian noise spectrum
    struct histogram_data *noise_histogram = histogram_create(data->mma, data->obs_count, min_mma, max_mma, bin_count);
    if (!noise_histogram)
        error_jump(noise_histogram_failed, ret, "Noise histogram failed");

    double noise_center, noise_amplitude, noise_sigma;
    if (histogram_fit_gaussian(noise_histogram, false, 0.01, &noise_center, &noise_amplitude, &noise_sigma))
        error_jump(noise_histogram_processing_failed, ret, "Noise histogram fit failed");

    snprintf(output_path_buffer, output_path_len, "%s.hist", output_prefix);
    if (histogram_export(noise_histogram, noise_center, noise_amplitude, noise_sigma, output_path_buffer))
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

    // Synthesize artificial data with random noise and find best-fit frequency
    double *fitted_ampl = calloc(randomize_count, sizeof(double));
    if (!fitted_ampl)
        error_jump(fitted_alloc_failed, ret, "Allocation failed");

    fprintf(randomized, "# Seed: %u\n", seed);
    for (size_t k = 0; k < randomize_count; k++)
    {
        if (k > 0 && !(k % 10))
            printf("%zu...\n", k);

        for (size_t i = 0; i < data->obs_count; i++)
        {
            data->mma[i] = random_normal(rand, noise_center, noise_sigma);
            for (size_t j = 0; j < data->freq_count; j++)
            {
                double phase = 2*M_PI*orig_freq[j]*data->time[i];
                data->mma[i] += orig_ampl[2*j]*cos(phase);
                data->mma[i] += orig_ampl[2*j+1]*sin(phase);
            }
        }

        // Step through fit frequencies
        double best_freq, best_chi2, best_ampl;
        step_freq_fit(data, freq_search_min, freq_search_max, freq_search_step, &best_freq, &best_ampl, &best_chi2);
        fprintf(randomized, "%zu %.2f %.2f %g\n", k, best_freq, best_ampl, best_chi2);
        fflush(randomized);

        fitted_freq[k] = best_freq - initial_freq;
        fitted_ampl[k] = best_ampl - initial_ampl;
    }

    struct histogram_data *freq_histogram = histogram_create(fitted_freq, randomize_count, freq_search_min - initial_freq, freq_search_max - initial_freq, (freq_search_max - freq_search_min) / freq_search_step);
    if (!freq_histogram)
        error_jump(freq_histogram_failed, ret, "Freq histogram failed");

    double freq_center, freq_amplitude, freq_sigma;
    if (histogram_fit_gaussian(freq_histogram, true, 0.01, &freq_center, &freq_amplitude, &freq_sigma))
        error_jump(freq_histogram_processing_failed, ret, "Freq histogram fit failed");

    snprintf(output_path_buffer, output_path_len, "%s.freqhist", output_prefix);
    if (histogram_export(freq_histogram, freq_center, freq_amplitude, freq_sigma, output_path_buffer))
        error_jump(freq_histogram_processing_failed, ret, "Freq histogram fit failed");

    FILE *output = fopen(output_path_buffer, "a+");
    double base_period = 1e6/base_uhz;
    double up_period = 1e6/(base_uhz - freq_sigma);
    double down_period = 1e6/(base_uhz + freq_sigma);
        
    fprintf(output, "# Freq: %.2f \\pm %.2f\n", base_period, (up_period - down_period) / 2);
    fclose(output);

    struct histogram_data *ampl_histogram = histogram_create(fitted_ampl, randomize_count, ampl_search_min - initial_ampl, ampl_search_max - initial_ampl, (ampl_search_max - ampl_search_min) / ampl_search_step);
    if (!ampl_histogram)
        error_jump(freq_histogram_failed, ret, "ampl histogram failed");

    double ampl_center, ampl_amplitude, ampl_sigma;
    if (histogram_fit_gaussian(ampl_histogram, true, 0.01, &ampl_center, &ampl_amplitude, &ampl_sigma))
        error_jump(freq_histogram_processing_failed, ret, "ampl histogram fit failed");

    snprintf(output_path_buffer, output_path_len, "%s.amplhist", output_prefix);
    if (histogram_export(ampl_histogram, ampl_center, ampl_amplitude, ampl_sigma, output_path_buffer))
        error_jump(freq_histogram_processing_failed, ret, "ampl histogram fit failed");

    {
        // Save output
        snprintf(output_path_buffer, output_path_len, "%s.data", output_prefix);
        FILE *file = fopen(output_path_buffer, "w");

        fprintf(file, "%s %.2f $\\pm$ %.2f & %f $\\pm$ %2f\n", output_prefix, initial_freq, freq_sigma, initial_ampl, ampl_sigma);

        fclose(file);
    }
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

int shuffle_dft_harmonics(char *ts_path, double base_uhz, size_t freq_count, double min_uhz, double max_uhz, double d_uhz, char *out_path, size_t repeats)
{
    int ret = 0;

    double freq_search_min = base_uhz - 10;
    double freq_search_max = base_uhz + 10;
    double freq_search_step = 0.05;

    struct ts_data *data = ts_data_load_harmonics(ts_path, 1e-6*base_uhz, freq_count);
    if (!data)
        error_jump(load_failed_error, ret, "Error processing data");
    
    // Optimize fit
    double best_freq;
    {
        double old_chi2 = ts_data_chi2(data);
        double best_chi2, best_ampl;
        step_freq_fit(data, freq_search_min, freq_search_max, freq_search_step, &best_freq, &best_ampl, &best_chi2);
        printf("Optimized base freq from %g to %g: delta-chi2: %g\n", base_uhz, best_freq, old_chi2 - best_chi2);

        for (size_t j = 0; j < freq_count; j++)
            data->freq[j] = 1e-6*(j+1)*best_freq;

        if (ts_data_fit_sinusoids(data))
            error_jump(dftfreq_alloc_error, ret, "Amplitude fit failed");
    }

    ts_data_print_table(data);
    ts_data_prewhiten(data);

    size_t dft_count = (size_t)((max_uhz - min_uhz)/d_uhz);
    double *dft_freq = calloc(dft_count, sizeof(double));
    if (!dft_freq)
        error_jump(dftfreq_alloc_error, ret, "Error allocating dft_freq");

    double *dft_ampl = calloc(dft_count, sizeof(double));
    if (!dft_ampl)
        error_jump(dftampl_alloc_error, ret, "Error allocating dft_freq");

    uint32_t seed = time(NULL);
    random_generator *rand = random_create(seed);
    if (!rand)
        error_jump(rand_alloc_error, ret, "Error allocating random generator");

    FILE *file = fopen(out_path, "w");
    if (!file)
        error_jump(outfile_open_error, ret, "Error opening output file %s", out_path);

    // File header
    fprintf(file, "# Prewhitened and shuffled max DFT amplitudes\n");
    fprintf(file, "# Data: %s\n", ts_path);
    fprintf(file, "# Freqs: %g", best_freq);
    for (size_t j = 1; j < freq_count; j++)
        fprintf(file, ", %g", (j+1)*best_freq);
    fprintf(file, "\n");
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
        fflush(file);
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


int test()
{
    char *paths[6] =
    {
        "20120517.fixed",
        "20120518.fixed",
        "20120519.fixed",
        "20120520.fixed",
        "20120521.fixed",
        "20120522.fixed",
    };

    for (uint8_t i = 0; i < 6; i++)
    {
        struct ts_data *data = ts_data_load(paths[i], NULL);

        for (size_t j = 0; j < data->obs_count; j++)
            printf("%g %g %g\n", i + data->time[j]/86400, data->mma[j], data->err[j]);

        ts_data_free(data);
    }
    return 0;
}
