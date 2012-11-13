/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef TS_ANALYSIS_H
#define TS_ANALYSIS_H

int display_tracer(char *dataPath);
int calculate_profile(char *dataPath, int obsIndex, int targetIndex);
int detect_repeats(char *dataPath);
int online_plot(char *dataPath, char *tsDevice, double tsSize, char *dftDevice, double dftSize);
int playback_reduction(char *dataPath, int delay, int step, char *tsDevice, double tsSize, char *dftDevice, double dftSize);

int reduce_aperture_range(char *base_name, double min, double max, double step, char *prefix);
int plot_range(char *datafile_pattern);
int report_time(char *dataPath);

#endif
