/*
 * Copyright 2010, 2011 Paul Chote
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

int monitor_phase_amplitude(char *ts_file, double base_uhz, size_t harmonic_count, double window_width);
int animated_window(char *tsfile);
int prewhiten_variable_freqs(char *tsfile, char *freqfile);

#endif
