/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef PLOTS_H
#define PLOTS_H

int online_focus_plot(char *data_path, const char *device, double size);
int online_plot(char *data_path, char *ts_device, char *dft_device, double size);
int playback_reduction(char *data_path, int delay, int step, char *ts_device, char *dft_device, double size);
int plot_range(char *datafile_pattern);

#endif
