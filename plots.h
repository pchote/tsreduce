/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef PLOTS_H
#define PLOTS_H
int online_focus_plot(char *data_path, const char *device, double size);
int online_plot(char *dataPath, char *tsDevice, double tsSize, char *dftDevice, double dftSize);
int playback_reduction(char *dataPath, int delay, int step, char *tsDevice, double tsSize, char *dftDevice, double dftSize);
int plot_range(char *datafile_pattern);

#endif
