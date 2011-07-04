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
} target;

typedef struct
{
    double star[3];
    double sky[3];
    double2 pos[2];
    double time;
    double ratio;
    char filename[64];
} record;

double2 center_aperture(target reg, double2 bg2, framedata *frame);
double2 converge_aperture(target r, framedata *frame);
double2 calculate_background(target r, framedata *frame);
double integrate_aperture(double2 xy, double r, framedata *frame);

#endif
