/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef FIT_H
#define FIT_H

int fit_sinusoids(double *x, double *y, double *e, size_t c, double *freqs, size_t numFreqs, double *amplitudes);
int fit_polynomial(double *x, double *y, double *e, size_t c, double *coeffs, size_t degree);
int fit_gaussian(double *x, double *y, size_t n, double min_sigma, double max_sigma, size_t count, double *sigma, double *ampl);

#endif
