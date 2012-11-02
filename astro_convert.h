/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <time.h>
#include "aperture.h"

#ifndef ASTRO_CONVERT_H
#define ASTRO_CONVERT_H

double jdtobjd(double jd, double2 coords);
double2 precess(double2 coords, double t0, double t1);
double tmtojd(struct tm *t);
double tmtoyear(struct tm *t);
double utcttoffset(time_t ut);

#endif