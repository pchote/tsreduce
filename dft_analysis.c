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
#include <stdbool.h>

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

static void set_color_table()
{
    float l[9] = {0.0, 0.33, 0.66, 1.0};
    float r[9] = {1.0, 0.0, 0.0, 1.0};
    float g[9] = {1.0, 0.0, 1.0, 0.0};
    float b[9] = {1.0, 1.0, 0.0, 0.0};
    cpgctab(l, r, g, b, 4, 1.0, 0.5);
}

static void plot_subrun(double start, double end, const char *label, const char *label2)
{
    double y = -80;
    double texty = -280;
    texty = -250;
    double mid = (start + end) / 2;
    
    cpgsch(1.2);
    cpgptxt(mid * 86400, texty, 0, 0.5f, label);
    cpgptxt(mid * 86400, texty - 250, 0, 0.5f, label2);
    cpgsch(1.6);
}

static void plot_tick(double y, const char *label, double min, double max)
{
    cpgmove(0.9925, y);
    cpgdraw(1, y);
    cpgmtxt("RV", 0.5, (y - min) / (max - min), 0, label);
}

int colorplot(const char *ts_path)
{
    int ret = 0;

    struct ts_data *data = calloc(1, sizeof(struct ts_data));
    if (!data)
        return 1;

    if (load_tsfile(ts_path, data))
        error_jump(load_failed_error, ret, "Error loading timeseries data");

    double freq_min = 0;
    double freq_max = 4250;
    size_t freq_steps = 500;
    size_t time_steps = 500;
    double window_extent = 3600 * 1.5;
    
    double time_min = data->time[0];
    double time_max = data->time[data->obs_count - 1];

    double *mmi_windowed = (double *)malloc(data->obs_count*sizeof(double));
    if (!mmi_windowed)
        error_jump(mmi_windowed_alloc_error, ret, "Error allocating mmi_windowed");

    if (cpgopen("plot.ps/vcps") <= 0)
//    if (cpgopen("2/xs") <= 0)
        error_jump(pgplot_open_error, ret, "Unable to open PGPLOT window");

    // 800 x 480
    cpgpap(6, 0.3);
    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(2);
    cpgsch(1.5);
    

    printf("%zu\n", data->obs_count);

	float x_scale = (time_max-time_min)/(time_steps-1);
	float y_scale = (freq_max-freq_min)/(freq_steps-1);
	float tr[] = {time_min-x_scale, x_scale, 0, freq_min-y_scale, 0, y_scale};

	float min_amplitude = 0;
	float max_amplitude = 50;
	float *amplitude = calloc(time_steps*freq_steps, sizeof(float));
    if (!amplitude)
        error_jump(amplitude_alloc_error, ret, "Error allocating amplitude");

    double *temp_ampl = calloc(freq_steps, sizeof(double));
    if (!temp_ampl)
        error_jump(temp_ampl_alloc_error, ret, "Error allocating temp_ampl");
    
    double *temp_freq = calloc(freq_steps, sizeof(double));
    if (!temp_freq)
        error_jump(temp_freq_alloc_error, ret, "Error allocating temp_freq");

    double *temp_time = calloc(data->obs_count, sizeof(double));
    if (!temp_time)
        error_jump(temp_time_alloc_error, ret, "Error allocating temp_time");
    
    double *temp_mmi = calloc(data->obs_count, sizeof(double));
    if (!temp_mmi)
        error_jump(temp_mmi_alloc_error, ret, "Error allocating temp_mmi");

    double dt = (time_max - time_min)/time_steps;

	set_color_table();

    for (size_t k = 2; k >= 1; k--)
    {
        double freq_diff = freq_max - freq_min;
        double mid_freq = k * (freq_max + freq_min) / 2;
        double local_freq_min = mid_freq - freq_diff/2;
        double local_freq_max = mid_freq + freq_diff/2;
        
        if (k == 2)
        {
            mid_freq = 3000;
            local_freq_min = mid_freq-500;
            local_freq_max = mid_freq+500;
        }

    	for (size_t i = 0; i < time_steps; i++)
    	{
            double mid_time = time_min + i*dt;
            size_t n = 0;
        	for (size_t j = 0; j < data->obs_count; j++)
            {
                if ((data->time[j] >= mid_time - window_extent) && (data->time[j] <= mid_time + window_extent))
                {
                    temp_time[n] = data->time[j];
                    temp_mmi[n] = data->mma[j];
                    n++;
                }
            }

            // Require the window to be at least half full
            if (n < window_extent / 30)
                continue;

            // Window function
            if (k == 2)
            {
                // Generate sinusoid
            	for (size_t j = 0; j < n; j++)
                    temp_mmi[j] = max_amplitude * sin(2*M_PI*mid_freq*1e-6*temp_time[j]);
            }

            calculate_amplitude_spectrum(temp_time, temp_mmi, n, local_freq_min*1e-6, local_freq_max*1e-6, temp_freq, temp_ampl, freq_steps);

        	for (size_t j = 0; j < freq_steps; j++)
                amplitude[j*time_steps + i] = temp_ampl[j];
    	}

        if (k == 2)
        {
            cpgsvp(0.1, 0.85, 0.8, 0.95);
            cpgmtxt("l", 4, 0.5, 0.5, "Window");
//            cpgmtxt("t", 0.5, 0.5, 0.5, ts_path);
            cpgswin(time_min, time_max, freq_min, freq_max);
            
        }
        else
        {
            cpgsvp(0.1, 0.85, 0.11, 0.79);
            cpgmtxt("l", 4, 0.5, 0.5, "Frequency (\\gmHz)");
            cpgswin(time_min, time_max, freq_min-700, freq_max);
      
            plot_subrun(-0.031250, 0.206961,"2 March", "2011");
            plot_subrun(0.269461, 0.405001, "4 March", "2011");
            plot_subrun(0.467501, 0.741213, "1 July ", "2011");
            plot_subrun(0.803713, 1.058433, "2 July ", "2011");
            plot_subrun(1.120933, 1.493017, "4 July ", "2011");
            plot_subrun(1.555517, 1.738032, "6 July ", "2011");
            plot_subrun(1.800532, 2.111979, "27 July ", "2011");
            plot_subrun(2.174479, 2.474109, "1 Aug", "2011");
            plot_subrun(2.536609, 2.744235, "2 Aug", "2011");
            plot_subrun(2.806735, 2.972541, "24 March", "2012");
            plot_subrun(3.035041, 3.317176, "25 March", "2012");
            plot_subrun(3.379676, 3.735184, "23 April", "2012");  
            plot_subrun(3.82904969, 4.02669433, "13 March", "2013");  
        }

    	//cpgsitf(2); // Use a sqrt mapping between value and colour
    	cpgimag(amplitude, time_steps, freq_steps, 1, time_steps, 1, freq_steps, min_amplitude, max_amplitude, tr);

        float label_stride = 1000;
        if (k == 2)
        {
            local_freq_min = -500;
            local_freq_max = 500;
            label_stride = 500;
        }

        double offset = (2456064.75658279 - 2456047.50638806) * 86400;
        printf("Times: %f %f\n", time_min + offset, time_max + offset);
        if (k == 2)
            cpgswin(time_min + offset, time_max + offset, local_freq_min, local_freq_max);
        else
            cpgswin(time_min + offset, time_max + offset, local_freq_min-700, local_freq_max);
        
        
        if (k == 2)
            cpgtbox("bcst", 12*3600, 4, "bstnv", label_stride, 4);
        else
            cpgtbox("bcstnZHXY", 12*3600, 4, "bstnv", label_stride, 4);

        cpgbox("0", 0, 0, "c", 0, 0);

        if (k != 2)
        {
            cpgswin(0, 1, local_freq_min-700, local_freq_max);
            plot_tick(4000, "250", local_freq_min - 700, local_freq_max);
            plot_tick(3333, "300", local_freq_min - 700, local_freq_max);
            plot_tick(2500, "400", local_freq_min - 700, local_freq_max);
            plot_tick(2000, "500", local_freq_min - 700, local_freq_max);
            plot_tick(1666, "600", local_freq_min - 700, local_freq_max);
            plot_tick(1000, "1000", local_freq_min - 700, local_freq_max);
            plot_tick(500,  "2000", local_freq_min - 700, local_freq_max);
            plot_tick(100, "10000", local_freq_min - 700, local_freq_max);
            
            cpgmtxt("r", 4, 0.5, 0.5, "Period (s)");
        }
    }
    

    cpgsvp(0.1, 0.85, 0.11, 0.95);
    cpgmtxt("B", 2.5, 0.5, 0.5, "Time (Hours)");
    /*
    for (size_t v = 0; v < 2; v++)
    {
        size_t k = 2;
        size_t l = k;//5 - k;
        double m = (k == 2 ? -0.02 : 0);
        cpgsvp(0.1, 0.89, n + (l-1)*o - m, n + l*o - m);

        char buf[1024]; strcpy(buf, ts_path);
        buf[strlen(buf)-3] = '\0';
        cpgmtxt("t", 0.5, 0.5, 0.5, buf);
        
        cpgsvp(0.1, 0.89, n, n + 3*o);
        cpgswin(time_min / 86400, time_max / 86400, freq_min, freq_max);
        cpgbox("bstn", 0, 0, "0", 0, 0);
    }
    */
    cpgsvp(0.91, 0.92, 0.11, 0.95);

    size_t amplitude_steps = 50;
    float *ampl_scale_data = calloc(amplitude_steps, sizeof(float));
    float ampl_scale = (max_amplitude-min_amplitude)/(amplitude_steps-1);
    for (size_t i = 0; i < amplitude_steps; i++)
        ampl_scale_data[i] = min_amplitude + i * ampl_scale;

    float atr[] = {-0.5, 1, 0, min_amplitude-ampl_scale, 0, ampl_scale, 0, 1};
    cpgswin(0, 1, min_amplitude, max_amplitude);
    cpgimag(ampl_scale_data, 1, amplitude_steps, 1, 1, 1, amplitude_steps, min_amplitude, max_amplitude, atr);

    cpgbox("bc", 0, 0, "bcsmv", 0, 0);
    cpgmtxt("r", 3.2, 0.5, 0.5, "Amplitude (mma)");

    free(temp_mmi);
temp_mmi_alloc_error:
    free(temp_time);
temp_time_alloc_error:
    free(temp_freq);
temp_freq_alloc_error:
    free(temp_ampl);
temp_ampl_alloc_error:
    free(amplitude);
amplitude_alloc_error:
    cpgend();
pgplot_open_error:
    free(mmi_windowed);
mmi_windowed_alloc_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}
