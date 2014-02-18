/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef DFT_ANALYSIS_H
#define DFT_ANALYSIS_H

int model_fit(char *tsFile, char *freqFile, double startTime, double endTime, double dt, char *modelFile, char *residualsFile);
int dft_bjd(char *tsFile, double minUHz, double maxUHz, double dUHz, char *outFile, char *freqFile);
int dft_window(char *tsFile, double freq, double minUHz, double maxUHz, double dUHz, char *outFile);
int find_max_freq(char *tsFile, char *freqFile, double minUHz, double maxUHz, double dUHz);
int nonlinear_fit(char *tsFile, char *freqFile);
int shuffle_dft(char *tsFile, char *freqFile, double minUHz, double maxUHz, double dUHz, char *outFile, size_t repeats);
int noise_histogram(const char *ts_path, const char *freq_path,
                    double min_mma, double max_mma, size_t bin_count,
                    double fit_min_mma, double fit_max_mma, size_t randomize_count,
                    const char *output_prefix);
int o_minus_c(const char *data, const char *ref_bjd, long double period, long double phase_offset);
int monitor_phase_amplitude(char *ts_file, double base_uhz, size_t harmonic_count, double window_width);
int fit_baseline_polynomial(char *tsfile, size_t poly_degree);

int fit_gwlib_freqshift(char *ts_file, double first_uhz, double second_uhz, size_t harmonic_count);

int gwlib_noise_histogram(const char *ts_path, double base_uhz, size_t freq_count,
                    double min_mma, double max_mma, size_t bin_count, size_t randomize_count,
                    const char *output_prefix);
int shuffle_dft_harmonics(char *ts_path, double base_uhz, size_t freq_count, double min_uhz, double max_uhz, double d_uhz, char *out_path, size_t repeats);

#endif
