/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef FIT_H
#define FIT_H

int fit_sinusoids(double *x, double *y, int c, double *freqs, int numFreqs, double *amplitudes);
int fit_polynomial(float *x, float *y, int c, double *coeffs, int degree);
int fit_polynomial_d(double *x, double *y, int c, double *coeffs, int degree);
int fit_polynomial_with_errors_d(double *x, double *y, double *e, int c, double *coeffs, int degree);
#endif