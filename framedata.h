/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef FRAMEDATA_H
#define FRAMEDATA_H

#include <stdbool.h>
#include <stdint.h>
#include <wcslib/wcslib.h>
#include "helpers.h"
#include "hashmap.h"

#define FRAME_METADATA_STRING 0
#define FRAME_METADATA_INT 1
#define FRAME_METADATA_DOUBLE 2
#define FRAME_METADATA_BOOL 3

#define FRAME_METADATA_OK 0
#define FRAME_METADATA_MISSING 1
#define FRAME_METADATA_INVALID_TYPE 2

typedef struct
{
    int rows;
    int cols;
    double *data;

    struct frame_metadata *metadata_start;
    struct frame_metadata *metadata_end;
    map_t metadata_map;
    
    int nwcs;
    struct wcsprm *wcs;
} framedata;

typedef struct
{
    // WCS conversion
    struct wcsprm *reference;
    struct wcsprm *frame;

    // Fallback translation
    int32_t dx;
    int32_t dy;
} framedata_translation;

framedata *framedata_load(const char *filename);
int framedata_save(framedata *fd, const char *path);
void framedata_free(framedata *frame);
void framedata_print_metadata(framedata *fd);

bool framedata_has_metadata(framedata *fd, const char *key);
int framedata_get_metadata(framedata *fd, const char *key, int type, void *data);
int framedata_put_metadata(framedata *fd, const char *key, int type, void *data, const char *comment);
int framedata_remove_metadata(framedata *fd, const char *key);

int framedata_subtract(framedata *fd, framedata *other);
int framedata_subtract_normalized(framedata *fd, framedata *other);
void framedata_subtract_bias(framedata *fd);
int framedata_calibrate(framedata *frame, framedata *bias, framedata *dark, framedata *flat);
int framedata_calibrate_load(framedata *frame, const char *bias_path, const char *dark_path, const char *flat_path);
int framedata_divide(framedata *fd, framedata *other);
int framedata_start_time(framedata *frame, ts_time *time);

framedata_translation *framedata_init_translation(framedata *frame, framedata *reference);
void framedata_get_translation(framedata_translation *ft, int32_t x, int32_t y, int32_t *dx, int32_t *dy);
void framedata_free_translation(framedata_translation *ft);

int framedata_image_region(framedata *frame, uint16_t region[4]);
int framedata_bias_region(framedata *frame, uint16_t region[4]);

#endif
