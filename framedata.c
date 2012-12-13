/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#include <fitsio2.h>
#include <string.h>

#include "framedata.h"
#include "helpers.h"

static int list_metadata_entry(any_t unused, any_t _meta)
{
    struct frame_metadata *metadata = _meta;
    printf("%s = ", metadata->key);
    switch(metadata->type)
    {
        case FRAME_METADATA_STRING:
            printf("%s (string) / ", metadata->value.s);
            break;
        case FRAME_METADATA_INT:
            printf("%lld (int) / ", metadata->value.i);
            break;
        case FRAME_METADATA_DOUBLE:
            printf("%f (double) / ", metadata->value.d);
            break;
        case FRAME_METADATA_BOOL:
            printf("%d (string) / ", metadata->value.b);
            break;
    }
    printf("%s\n", metadata->comment);
    return MAP_OK;
}

static int free_metadata_entry(any_t unused, any_t _meta)
{
    struct frame_metadata *metadata = _meta;
    free(metadata->key);
    free(metadata->comment);

    if (metadata->type == FRAME_METADATA_STRING)
        free(metadata->value.s);

    free(metadata);
    return MAP_OK;
}

framedata *framedata_load(const char *filename)
{
    int ret = 0;
	int status = 0;
    fitsfile *input;

    framedata *fd = calloc(1, sizeof(framedata));
    if (!fd)
        return NULL;

    fd->metadata_map = hashmap_new();

    if (fits_open_image(&input, filename, READONLY, &status))
    {
        char fitserr[128];
        while (fits_read_errmsg(fitserr))
            error("%s\n", fitserr);

        error_jump(error, ret, "fits_open_image failed with error %d; %s", status, filename);
    }

    // Query the image size
    fits_read_key(input, TINT, "NAXIS1", &fd->cols, NULL, &status);
    fits_read_key(input, TINT, "NAXIS2", &fd->rows, NULL, &status);
    if (status)
        error_jump(error, ret, "querying NAXIS failed");

    fd->data = (double *)malloc(fd->cols*fd->rows*sizeof(double));
    if (!fd->data)
        error_jump(error, ret, "Error allocating frame data");

    if (fits_read_pix(input, TDOUBLE, (long []){1, 1}, fd->cols*fd->rows, 0, fd->data, NULL, &status))
        error_jump(error, ret, "fits_read_pix failed");

    // Load header keys
    int keyword_count = 0;
    if (fits_get_hdrspace(input, &keyword_count, NULL, &status))
        error_jump(error, ret, "fits_get_hdrspace failed");

    for (size_t i = 0; i < keyword_count; i++)
    {
        char card[FLEN_CARD];
        char key[FLEN_KEYWORD];
        char value[FLEN_VALUE];
        char comment[FLEN_COMMENT];

        struct frame_metadata *metadata = calloc(1, sizeof(struct frame_metadata));
        if (!metadata)
            error_jump(key_read_error, ret, "Error allocating key %zu", i);

        // We're only interested in user keys
        if (fits_read_record(input, i + 1, card, &status))
            error_jump(key_read_error, ret, "Error reading card %zu", i);

        if (fits_get_keyclass(card) != TYP_USER_KEY)
            continue;

        if (fits_read_keyn(input, i + 1, key, value, comment, &status))
            error_jump(key_read_error, ret, "Error reading key %zu", i);

        // Parse value
        char type;
        if (fits_get_keytype(value, &type, &status))
            error_jump(key_read_error, ret, "Error determining type for '%s'", value);

        switch (type)
        {
            case 'I':
            {
                LONGLONG val;
                if (ffc2jj(value, &val, &status))
                    error_jump(key_read_error, ret, "Error parsing '%s' as integer", value);

                metadata->type = FRAME_METADATA_INT;
                metadata->value.i = val;
                break;
            }
            case 'F':
            {
                if (ffc2dd(value, &metadata->value.d, &status))
                    error_jump(key_read_error, ret, "Error parsing '%s' as double", value);

                metadata->type = FRAME_METADATA_DOUBLE;
                break;
            }
            case 'L':
            {
                int val;
                if (ffc2ll(value, &val, &status))
                    error_jump(key_read_error, ret, "Error parsing '%s' as boolean", value);

                metadata->type = FRAME_METADATA_BOOL;
                metadata->value.b = val;
                break;
            }
            default:
            {
                char val[FLEN_VALUE];
                if (ffc2s(value, val, &status))
                    error_jump(key_read_error, ret, "Error parsing '%s' as type '%c'", value, type);

                metadata->type = FRAME_METADATA_STRING;
                metadata->value.s = strdup(val);
                if (!metadata->value.s)
                    error_jump(key_read_error, ret, "Error parsing '%s' as type '%c'", value, type);
                break;
            }
        }

        metadata->key = strdup(key);
        metadata->comment = strdup(comment);
        if (!metadata->key || !metadata->comment)
            error_jump(key_read_error, ret, "Error reading key %zu", i);

        hashmap_put(fd->metadata_map, metadata->key, metadata);
    }

    // Load image regions
    int *ir = fd->regions.image_region;
    int *br = fd->regions.bias_region;
    ir[0] = 0; ir[1] = fd->cols; ir[2] = 0; ir[3] = fd->rows;
    br[0] = br[1] = br[2] = br[3] = 0;

    struct frame_metadata *m;
    if (hashmap_get(fd->metadata_map, "BIAS-RGN", (void**)(&m)) == MAP_OK)
    {
        fd->regions.has_overscan = true;
        sscanf(m->value.s, "[%d, %d, %d, %d]", &br[0], &br[1], &br[2], &br[3]);
    }

    if (hashmap_get(fd->metadata_map, "IMAG-RGN", (void**)(&m)) == MAP_OK)
        sscanf(m->value.s, "[%d, %d, %d, %d]", &ir[0], &ir[1], &ir[2], &ir[3]);

    fd->regions.image_px = (ir[1] - ir[0])*(ir[3] - ir[2]);
    fd->regions.bias_px = (br[1] - br[0])*(br[3] - br[2]);

    fits_close_file(input, &status);

    return fd;

key_read_error:
    hashmap_iterate(fd->metadata_map, free_metadata_entry, NULL);
error:
    framedata_free(fd);

    fits_close_file(input, &status);
    return NULL;
}

bool framedata_has_metadata(framedata *fd, const char *key)
{
    struct frame_metadata *metadata;
    return hashmap_get(fd->metadata_map, key, (void**)(&metadata)) == MAP_OK;
}

int framedata_get_metadata(framedata *fd, const char *key, int type, void *data)
{
    struct frame_metadata *metadata;
    if (hashmap_get(fd->metadata_map, key, (void**)(&metadata)) == MAP_MISSING)
        return FRAME_METADATA_MISSING;

    switch (type)
    {
        case FRAME_METADATA_STRING:
        {
            // Allow any type to be returned as string
            char buf[100];
            char *val;
            switch (metadata->type)
            {
                case FRAME_METADATA_STRING: val = metadata->value.s; break;
                case FRAME_METADATA_BOOL: val = metadata->value.b ? "true" : "false"; break;
                case FRAME_METADATA_DOUBLE:
                {
                    snprintf(buf, 100, "%f", metadata->value.d);
                    val = buf;
                    break;
                }
                case FRAME_METADATA_INT:
                {
                    snprintf(buf, 100, "%lld", metadata->value.i);
                    val = buf;
                    break;
                }
            }

            *(char **)data = strdup(val);
            return FRAME_METADATA_OK;
        }
        case FRAME_METADATA_INT:
        {
            if (metadata->type != FRAME_METADATA_INT)
                return FRAME_METADATA_INVALID_TYPE;

            *(int64_t *)data = metadata->value.i;
            return FRAME_METADATA_OK;
        }
        case FRAME_METADATA_DOUBLE:
        {
            // Allow int to be returned as double
            if (metadata->type != FRAME_METADATA_DOUBLE &&
                metadata->type != FRAME_METADATA_INT)
                return FRAME_METADATA_INVALID_TYPE;

            *(double *)data = metadata->type == FRAME_METADATA_DOUBLE ?
                metadata->value.d : (double)metadata->value.i;
            return FRAME_METADATA_OK;
        }
        case FRAME_METADATA_BOOL:
        {
            // Allow int to be returned as bool
            if (metadata->type != FRAME_METADATA_BOOL &&
                metadata->type != FRAME_METADATA_INT)
                return FRAME_METADATA_INVALID_TYPE;

            *(bool *)data = metadata->type == FRAME_METADATA_BOOL ?
                metadata->value.b : (bool)metadata->value.i;
            return FRAME_METADATA_OK;
        }
    }

    return FRAME_METADATA_OK;
}

int framedata_subtract(framedata *fd, framedata *other)
{
    if (fd->cols != other->cols || fd->rows != other->rows)
        return error("Frame size mismatch");

    for (size_t i = 0; i < fd->cols*fd->rows; i++)
        fd->data[i] -= other->data[i];

    return 0;
}

int framedata_divide(framedata *fd, framedata *other)
{
    if (fd->cols != other->cols || fd->rows != other->rows)
        return error("Frame size mismatch");

    for (size_t i = 0; i < fd->cols*fd->rows; i++)
        fd->data[i] /= other->data[i];

    return 0;
}

void framedata_free(framedata *frame)
{
    if (!frame)
        return;

    if (frame->data)
        free(frame->data);

    hashmap_iterate(frame->metadata_map, free_metadata_entry, NULL);
    hashmap_free(frame->metadata_map);

    free(frame);
}

int framedata_start_time(framedata *fd, ts_time *out_time)
{
    struct frame_metadata *date, *time;
    hashmap_get(fd->metadata_map, "UTC-DATE", (void **)(&date));
    hashmap_get(fd->metadata_map, "UTC-BEG", (void **)(&time));

    if (date && date->type == FRAME_METADATA_STRING &&
        time && time->type == FRAME_METADATA_STRING)
    {
        *out_time = parse_date_time(date->value.s, time->value.s);
        return 0;
    }

    // Legacy keywords
    hashmap_get(fd->metadata_map, "GPSTIME", (void **)(&date));
    if (date && date->type == FRAME_METADATA_STRING)
    {
        *out_time = parse_time(date->value.s);
        return 0;
    }

    return 1;
}

// Convenience function for calculating the mean signal in a sub-region of a frame
// Assumes that the frame type is double, and that the region is inside the frame
double mean_in_region(framedata *frame, int rgn[4])
{
    int num_px = (rgn[1] - rgn[0])*(rgn[3] - rgn[2]);
    double mean = 0;
    for (int j = rgn[2]; j < rgn[3]; j++)
        for (int i = rgn[0]; i < rgn[1]; i++)
            mean += frame->data[j*frame->cols + i]/num_px;
    return mean;
}

// Calculate and subtract the mean bias level from a frame
void subtract_bias(framedata *frame)
{
    // Calculate and subtract bias if the frame has overscan
    if (!frame->regions.has_overscan)
        return;

    double mean_bias = mean_in_region(frame, frame->regions.bias_region);
    for (int i = 0; i < frame->rows*frame->cols; i++)
        frame->data[i] -= mean_bias;
}
