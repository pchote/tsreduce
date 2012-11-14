/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#ifndef FRAMEDATA_H
#define FRAMEDATA_H

#include <stdbool.h>
#include "helpers.h"

typedef struct
{
    bool has_overscan;
    int image_region[4];
    int image_px;
    int bias_region[4];
    int bias_px;
} frameregions;

typedef struct
{
    void *fptr;
    int rows;
    int cols;
    double *data;
    frameregions regions;
} framedata;

framedata *framedata_load(const char *filename);
void framedata_free(framedata *frame);

int framedata_get_header_long(framedata *this, const char *key, long *value);
int framedata_get_header_dbl(framedata *this, const char *key, double *value);
char *framedata_get_header_string(framedata *this, const char *key);
void framedata_subtract(framedata *this, framedata *other);
void framedata_add(framedata *this, framedata *other);
void framedata_subtract(framedata *this, framedata *other);
void framedata_multiply(framedata *this, int div);
void framedata_divide_const(framedata *this, int div);
void framedata_divide(framedata *this, framedata *div);
ts_time framedata_start_time(framedata *frame);

double mean_in_region(framedata *frame, int rgn[4]);
void subtract_bias(framedata *frame);

#endif