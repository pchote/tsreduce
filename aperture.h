/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#include "framedata.h"

#ifndef APERTURE_H
#define APERTURE_H

// Compatability tweaks
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

/* Represents an aquired frame */
typedef struct
{
    double x;
    double y;
} double2;

typedef struct
{
    double x;
    double y;
    double r;
    double s1;
    double s2;
    double plot_scale;
} target;

int center_aperture(target r, framedata *frame, double2 *center);
int calculate_background(target r, framedata *frame, double *sky_mode, double *sky_std_dev);
double integrate_aperture(double2 xy, double r, framedata *frame);

#endif
