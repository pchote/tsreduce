/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#include <fitsio.h>

#ifndef FRAMEDATA_H
#define FRAMEDATA_H

typedef enum
{
    FRAMEDATA_INT,
    FRAMEDATA_DBL
} framedata_type;

typedef struct
{
    fitsfile *_fptr;
    int rows;
    int cols;
    framedata_type dtype;
    int *data;
    double *dbl_data;
} framedata;


framedata framedata_new(const char *filename, framedata_type dtype);
int framedata_get_header_int(framedata *this, const char *key);
int framedata_has_header_string(framedata *this, const char *key);
int framedata_has_header_string(framedata *this, const char *key);
void framedata_get_header_string(framedata *this, const char *key, char *ret);
void framedata_subtract(framedata *this, framedata *other);
void framedata_add(framedata *this, framedata *other);
void framedata_subtract(framedata *this, framedata *other);
void framedata_multiply(framedata *this, int div);
void framedata_divide_const(framedata *this, int div);
void framedata_divide(framedata *this, framedata *div);
void framedata_free(framedata this);

int compare_double(const void *a, const void *b);
int load_reject_minmax( const char **frames, int numFrames, int rows, int cols, int rejectHigh, int rejectLow, double *outFrame, int (*preprocess_func)(framedata*, void*), void *preprocess_data);

#endif