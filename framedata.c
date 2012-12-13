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

static framedata *framedata_alloc()
{
    framedata *fp = malloc(sizeof(framedata));
    if (!fp)
        return NULL;

    fp->rows = 0;
    fp->cols = 0;
    fp->data = NULL;
    fp->regions.has_overscan = false;
    fp->metadata_map = hashmap_new();

    return fp;
}

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

    framedata *fp = framedata_alloc();
    if (!fp)
    {
        error("framedata_alloc failed");
        return NULL;
    }

	int status = 0;
    fitsfile *input;

    if (fits_open_image(&input, filename, READONLY, &status))
    {
        char fitserr[128];
        while (fits_read_errmsg(fitserr))
            error("%s\n", fitserr);

        error_jump(error, ret, "fits_open_image failed with error %d; %s", status, filename);
    }

    // Query the image size
    fits_read_key(input, TINT, "NAXIS1", &fp->cols, NULL, &status);
    fits_read_key(input, TINT, "NAXIS2", &fp->rows, NULL, &status);
    if (status)
        error_jump(error, ret, "querying NAXIS failed");

    fp->data = (double *)malloc(fp->cols*fp->rows*sizeof(double));
    if (!fp->data)
        error_jump(error, ret, "Error allocating frame data");

    if (fits_read_pix(input, TDOUBLE, (long []){1, 1}, fp->cols*fp->rows, 0, fp->data, NULL, &status))
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

        hashmap_put(fp->metadata_map, metadata->key, metadata);
    }

    // Load image regions
    int *ir = fp->regions.image_region;
    int *br = fp->regions.bias_region;
    ir[0] = 0; ir[1] = fp->cols; ir[2] = 0; ir[3] = fp->rows;
    br[0] = br[1] = br[2] = br[3] = 0;

    char *bias_region_str = framedata_get_header_string(fp, "BIAS-RGN");
    if (bias_region_str)
    {
        fp->regions.has_overscan = true;
        sscanf(bias_region_str, "[%d, %d, %d, %d]", &br[0], &br[1], &br[2], &br[3]);
    }

    char *image_region_str = framedata_get_header_string(fp, "IMAG-RGN");
    if (image_region_str)
        sscanf(image_region_str, "[%d, %d, %d, %d]", &ir[0], &ir[1], &ir[2], &ir[3]);

    fp->regions.image_px = (ir[1] - ir[0])*(ir[3] - ir[2]);
    fp->regions.bias_px = (br[1] - br[0])*(br[3] - br[2]);

    fits_close_file(input, &status);

    return fp;

key_read_error:
    hashmap_iterate(fp->metadata_map, free_metadata_entry, NULL);
error:
    framedata_free(fp);

    fits_close_file(input, &status);
    return NULL;
}

struct frame_metadata *framedata_metadata(framedata *fd, char *key)
{
    struct frame_metadata *metadata;
    if (hashmap_get(fd->metadata_map, key, (void**)(&metadata)) == MAP_MISSING)
        return NULL;

    return metadata;
}

int framedata_get_header_long(framedata *fd, const char *key, long *value)
{
    struct frame_metadata *metadata;
    if (hashmap_get(fd->metadata_map, key, (void**)(&metadata)) == MAP_MISSING)
        return FRAME_METADATA_MISSING;

    if (metadata->type != FRAME_METADATA_INT)
        return FRAME_METADATA_INVALID_TYPE;

    *value = (long)metadata->value.i;
    return FRAME_METADATA_OK;
}

int framedata_get_header_dbl(framedata *fd, const char *key, double *value)
{
    struct frame_metadata *metadata;
    if (hashmap_get(fd->metadata_map, key, (void**)(&metadata)) == MAP_MISSING)
        return FRAME_METADATA_MISSING;

    if (metadata->type == FRAME_METADATA_DOUBLE)
        *value = metadata->value.d;
    else if (metadata->type == FRAME_METADATA_INT)
        *value = (double)metadata->value.i;
    else
        return FRAME_METADATA_INVALID_TYPE;

    return FRAME_METADATA_OK;
}

char *framedata_get_header_string(framedata *fd, const char *key)
{
    struct frame_metadata *metadata;
    if (hashmap_get(fd->metadata_map, key, (void**)(&metadata)) == MAP_MISSING)
        return NULL;

    if (metadata->type != FRAME_METADATA_STRING)
        return NULL;

    return strdup(metadata->value.s);
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

int framedata_start_time(framedata *frame, ts_time *out_time)
{
    char *date = framedata_get_header_string(frame, "UTC-DATE");
    char *time = framedata_get_header_string(frame, "UTC-BEG");
    if (date && time)
    {
        *out_time = parse_date_time(date, time);
        free(date);
        free(time);
        return 0;
    }

    // Legacy keywords
    char *datetime = framedata_get_header_string(frame, "GPSTIME");
    if (datetime)
    {
        *out_time = parse_time(datetime);
        free(datetime);
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
