/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef ASTRO_CONVERT_H
#define ASTRO_CONVERT_H

void getbary(double rar, double decr, double djd, double *bjd, double *baryc, double *helioc);
void precess(double ra0, double dec0, double *ra, double *dec, double t0, double t1);
double jd(double day, double month, double year, double ut);

#endif