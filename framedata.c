/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <fitsio2.h>
#include <string.h>
#include <math.h>

#include "framedata.h"
#include "helpers.h"

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

        // We're only interested in user keys
        if (fits_read_record(input, i + 1, card, &status))
            error_jump(error, ret, "Error reading card %zu", i);

        // Ignore keywords we aren't interested in
        int class = fits_get_keyclass(card);
        if (class < TYP_WCS_KEY || class == TYP_COMM_KEY)
            continue;

        if (fits_read_keyn(input, i + 1, key, value, comment, &status))
            error_jump(error, ret, "Error reading key %zu", i);

        // Parse value
        char type;
        if (fits_get_keytype(value, &type, &status))
        {
            error("Unable to determine type for key '%s' - skipping.", key);
            continue;
        }

        switch (type)
        {
            case 'I':
            {
                LONGLONG val;
                if (ffc2jj(value, &val, &status))
                    error_jump(error, ret, "Error parsing '%s' as integer", value);

                framedata_put_metadata(fd, key, FRAME_METADATA_INT, &(int64_t){val}, comment);
                break;
            }
            case 'F':
            {
                double val;
                if (ffc2dd(value, &val, &status))
                    error_jump(error, ret, "Error parsing '%s' as double", value);

                framedata_put_metadata(fd, key, FRAME_METADATA_DOUBLE, &val, comment);
                break;
            }
            case 'L':
            {
                int val;
                if (ffc2ll(value, &val, &status))
                    error_jump(error, ret, "Error parsing '%s' as boolean", value);

                framedata_put_metadata(fd, key, FRAME_METADATA_BOOL, &(bool){val}, comment);
                break;
            }
            default:
            {
                char val[FLEN_VALUE];
                if (ffc2s(value, val, &status))
                    error_jump(error, ret, "Error parsing '%s' as type '%c'", value, type);

                framedata_put_metadata(fd, key, FRAME_METADATA_STRING, val, comment);
                break;
            }
        }

    }

    fits_close_file(input, &status);

    return fd;

error:
    framedata_free(fd);

    fits_close_file(input, &status);
    return NULL;
}

int framedata_save(framedata *fd, const char *path)
{
    fitsfile *out;
    int status = 0;

    // Prepending path with a '!' tells cfitsio to
    // overwrite any existing file.
    size_t filepath_len = strlen(path) + 2;
    char *filepath = malloc(filepath_len*sizeof(char));
    snprintf(filepath, filepath_len, "!%s", path);

    fits_create_file(&out, filepath, &status);
    free(filepath);

    // Create the primary array image (16-bit short integer pixels
    fits_create_img(out, DOUBLE_IMG, 2, (long []){fd->cols, fd->rows}, &status);

    // Set header keys
    for (struct frame_metadata *m = fd->metadata_start; m; m = m->next)
    {
        switch (m->type)
        {
            case FRAME_METADATA_STRING:
                fits_update_key(out, TSTRING, m->key, m->value.s, m->comment, &status);
                break;
            case FRAME_METADATA_INT:
                fits_update_key(out, TLONG, m->key, &(long){m->value.i}, m->comment, &status);
                break;
            case FRAME_METADATA_DOUBLE:
                fits_update_key(out, TDOUBLE, m->key, &m->value.d, m->comment, &status);
                break;
            case FRAME_METADATA_BOOL:
                fits_update_key(out, TLOGICAL, m->key, &(int){m->value.b}, m->comment, &status);
                break;
        }

    }

    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, fd->rows*fd->cols, fd->data, &status))
        error("fits_write_img failed with status %d", status);

    fits_close_file(out, &status);
    return 0;
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

int framedata_put_metadata(framedata *fd, const char *key, int type, void *data, const char *comment)
{
    struct frame_metadata *m = calloc(1, sizeof(struct frame_metadata));
    if (!m)
        return error("Metadata allocation failed");

    // Check for and remove existing key
    struct frame_metadata *find;
    if (hashmap_get(fd->metadata_map, key, (void**)(&find)) == MAP_OK)
        framedata_remove_metadata(fd, key);

    switch (type)
    {
        case FRAME_METADATA_INT: m->value.i = *(int64_t *)data; break;
        case FRAME_METADATA_DOUBLE: m->value.d = *(double *)data; break;
        case FRAME_METADATA_BOOL: m->value.b = *(bool *)data; break;
        case FRAME_METADATA_STRING: m->value.s = strdup((char *)data); break;
        default:
            free(m);
            return error("Unknown metadata type: %d", type);
    }

    m->type = type;
    m->key = strdup(key);
    m->comment = strdup(comment);

    if (!m->key || !m->comment ||
        (type == FRAME_METADATA_STRING && !m->value.s))
    {
        if (type == FRAME_METADATA_STRING)
            free(m->value.s);
        free(m);
        return error("Null string when inserting metadata");
    }

    // Add to list and hashmap
    m->prev = fd->metadata_end;
    m->next = NULL;
    if (fd->metadata_end)
        fd->metadata_end->next = m;
    else
        fd->metadata_start = m;
    fd->metadata_end = m;

    hashmap_put(fd->metadata_map, m->key, m);

    return 0;
}

int framedata_remove_metadata(framedata *fd, const char *key)
{
    struct frame_metadata *m;
    if (hashmap_get(fd->metadata_map, key, (void**)(&m)) == MAP_MISSING)
        return FRAME_METADATA_MISSING;

    if (m->prev)
        m->prev->next = m->next;
    if (m->next)
        m->next->prev = m->prev;

    hashmap_remove(fd->metadata_map, key);
    free_metadata_entry(NULL, m);

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

void framedata_subtract_bias(framedata *fd)
{
    struct frame_metadata *m;
    if (hashmap_get(fd->metadata_map, "BIAS-RGN", (void**)(&m)) == MAP_MISSING)
        return;

    uint16_t br[4] = {0, 0, 0, 0};
    sscanf(m->value.s, "[%hu, %hu, %hu, %hu]", &br[0], &br[1], &br[2], &br[3]);

    double mean_bias = region_mean(br, fd->data, fd->cols);
    for (size_t i = 0; i < fd->rows*fd->cols; i++)
        fd->data[i] -= mean_bias;
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

void framedata_print_metadata(framedata *fd)
{
    for (struct frame_metadata *m = fd->metadata_start; m; m = m->next)
    {
        printf("%s = ", m->key);
        switch (m->type)
        {
            case FRAME_METADATA_STRING:
                printf("%s (string) / ", m->value.s);
                break;
            case FRAME_METADATA_INT:
                printf("%lld (int) / ", m->value.i);
                break;
            case FRAME_METADATA_DOUBLE:
                printf("%f (double) / ", m->value.d);
                break;
            case FRAME_METADATA_BOOL:
                printf("%s (bool) / ", m->value.b ? "T" : "F");
                break;
        }
        printf("%s\n", m->comment);
    }
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

    // SBIG CCD-OPS timestamp
    hashmap_get(fd->metadata_map, "DATE-OBS", (void **)(&date));
    if (date && date->type == FRAME_METADATA_STRING)
    {
        *out_time = parse_time_ccdops(date->value.s);
        return 0;
    }

    return 1;
}

int framedata_image_region(framedata *frame, uint16_t region[4])
{
    uint16_t r[4] = {0, frame->cols, 0, frame->rows};
    char *str;
    if (framedata_get_metadata(frame, "IMAG-RGN", FRAME_METADATA_STRING, &str) == FRAME_METADATA_OK)
        sscanf(str, "[%hu, %hu, %hu, %hu]", &r[0], &r[1],
                                            &r[2], &r[3]);
    memcpy(region, r, 4*sizeof(uint16_t));
    return 0;
}

static int sum_into_axes(framedata *frame, uint16_t region[4], double **x, double **y)
{
    int ret = 0;
    uint16_t xs = region[1] - region[0];
    uint16_t ys = region[3] - region[2];

    double *xa = calloc(xs, sizeof(double));
    double *ya = calloc(ys, sizeof(double));
    if (!xa || !ya)
        error_jump(error, ret, "allocation_failed");

    // Sum into x and y axes
    for (uint16_t j = 0; j < ys; j++)
        for (uint16_t i = 0; i < xs; i++)
        {
            double px = frame->data[frame->cols*(region[2] + j) + (i + region[0])];
            xa[i] += px;
            ya[j] += px;
        }

    // Subtracting mean level puts the baseline below zero, so clamp at zero
    double xm = mean_exclude_sigma(xa, xs, 1);
    for (uint16_t i = 0; i < xs; i++)
        xa[i] = fmax(xa[i] - xm, 0);

    double ym = mean_exclude_sigma(ya, ys, 1);
    for (uint16_t j = 0; j < ys; j++)
        ya[j] = fmax(ya[j] - ym, 0);

    *x = xa;
    *y = ya;

    return 0;
error:
    free(xa);
    free(ya);
    x = NULL;
    y = NULL;
    return ret;
}

int32_t find_max_correlatation(double *a, double *b, uint16_t n)
{
    int32_t best_idx = 0;
    double best = -1;
    for (int32_t j = -n; j < n; j++)
    {
        double ret = 0;
        for (uint16_t i = 0; i < n; i++)
            ret += (i+j > 0) ? a[i]*b[i+j] : 0;

        if (ret > best)
        {
            best = ret;
            best_idx = j;
        }
    }
    return best_idx;
}

int framedata_estimate_translation(framedata *frame, framedata *reference, int32_t *xt, int32_t *yt)
{
    int ret = 0;

    uint16_t fr[4], rr[4];
    if (framedata_image_region(frame, fr))
        error_jump(error, ret, "Image region failed");
    if (framedata_image_region(frame, rr))
        error_jump(error, ret, "Image region failed");

    if (fr[0] != rr[0] || fr[1] != rr[1] || fr[2] != rr[2] || fr[3] != rr[3])
        error_jump(error, ret, "frame sizes don't match");

    double *fx, *fy, *rx, *ry;
    if (sum_into_axes(frame, fr, &fx, &fy))
        error_jump(error, ret, "Summing into axes failed");
    if (sum_into_axes(reference, rr, &rx, &ry))
        error_jump(error, ret, "Summing into axes failed");

    *xt = find_max_correlatation(rx, fx, rr[1] - rr[0]);
    *yt = find_max_correlatation(ry, fy, rr[3] - rr[2]);
    return ret;
error:
    *xt = 0;
    *yt = 0;
    return ret;
}

