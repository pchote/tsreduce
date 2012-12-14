/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#ifndef FRAMEDATA_H
#define FRAMEDATA_H

#include <stdbool.h>
#include <stdint.h>
#include "helpers.h"
#include "hashmap.h"

#define FRAME_METADATA_STRING 0
#define FRAME_METADATA_INT 1
#define FRAME_METADATA_DOUBLE 2
#define FRAME_METADATA_BOOL 3

#define FRAME_METADATA_OK 0
#define FRAME_METADATA_MISSING 1
#define FRAME_METADATA_INVALID_TYPE 2
struct frame_metadata
{
    char *key;
    char *comment;

    uint8_t type;
    union
    {
        char *s;
        bool b;
        double d;
        int64_t i;
    } value;

    struct frame_metadata *prev;
    struct frame_metadata *next;
};

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
    int rows;
    int cols;
    double *data;
    frameregions regions;

    struct frame_metadata *metadata_start;
    struct frame_metadata *metadata_end;
    map_t metadata_map;
} framedata;

framedata *framedata_load(const char *filename);
int framedata_save(framedata *fd, const char *path);
void framedata_free(framedata *frame);
void framedata_print_metadata(framedata *fd);

bool framedata_has_metadata(framedata *fd, const char *key);
int framedata_get_metadata(framedata *fd, const char *key, int type, void *data);
int framedata_put_metadata(framedata *fd, const char *key, int type, void *data, const char *comment);
int framedata_remove_metadata(framedata *fd, const char *key);

int framedata_subtract(framedata *fd, framedata *other);
int framedata_divide(framedata *fd, framedata *other);
int framedata_start_time(framedata *frame, ts_time *time);

double mean_in_region(framedata *frame, int rgn[4]);
void subtract_bias(framedata *frame);

#endif
