/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <fitsio2.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

#include "framedata.h"
#include "helpers.h"
#include "fit.h"

extern int verbosity;

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

static int load_and_fix_padding(const char *filename, uint8_t **framedata, size_t *length)
{
    *framedata = NULL;
    *length = 0;

    FILE *frame = fopen(filename, "r");
    if (!frame)
        return error("Unable to open frame %s", filename);

    // Query file size
    fseek(frame, 0L, SEEK_END);
    size_t file_size = (size_t)ftell(frame);
    fseek(frame, 0L, SEEK_SET);

    size_t mem_size = 2880 * (file_size / 2880);
    if (mem_size != file_size)
    {
        error("Warning: Truncated file detected: %s", filename);

        // Add an extra block to hold the remainder plus padding
        mem_size += 2880;
    }

    uint8_t *data = calloc(mem_size, sizeof(uint8_t));
    if (!data)
    {
        fclose(frame);
        return error("Unable to allocate %zu bytes for %s", mem_size, filename);
    }

    if (fread(data, sizeof(uint8_t), file_size, frame) != file_size)
    {
        free(data);
        fclose(frame);
        return error("Unable to read %zu bytes from %s", file_size, filename);
    }

    fclose(frame);

    *framedata = data;
    *length = mem_size;
    return 0;
}

static void print_fits_error()
{
    char fitserr[128];
    while (fits_read_errmsg(fitserr))
        error("        %s\n", fitserr);
}

framedata *framedata_load(const char *filename)
{
    int ret = 0;
    int status = 0;
    fitsfile *input;
    uint8_t *padded_data = NULL;

    // Suppress unused variable warning
    (void)ret;

    framedata *fd = calloc(1, sizeof(framedata));
    if (!fd)
        return NULL;

    fd->metadata_map = hashmap_new();

    if (fits_open_image(&input, filename, READONLY, &status))
    {
        if (verbosity >= 1)
            print_fits_error();

        error_jump(error, ret, "        fits_open_image failed with error %d; %s", status, filename);
    }

    // The FITS specification requires files be an integer number of 2880 blocks.
    // If the length doesn't match then the file may be corrupted or saved by
    // software that violates the specification (e.g. Quilt).
    //
    // It is assumed that any uncompressed files that aren't an integer block size have
    // been saved using Quilt, and are padded to the required length and re-opened.
    // This may produce unexpected results if the frame was truncated for other reasons
    // (e.g. interrupted transfer).

    char type[20]; // MAX_PREFIX_LEN from cfitsio's cfileio.c
    fits_url_type(input, type, &status);

    if (strcmp("file://", type) == 0 && (input->Fptr->filesize % 2880) != 0)
    {
        size_t length = 0;
        if (load_and_fix_padding(filename, &padded_data, &length))
            error_jump(error, ret, "load_and_fix_padding failed for frame: %s", filename);

        fits_close_file(input, &status);
        if (fits_open_memfile(&input, filename, READONLY, (void **)(&padded_data), &length, 0, NULL, &status))
        {
            print_fits_error();
            error_jump(error, ret, "fits_open_memfile failed with error %d; %s", status, filename);
        }
    }

    // Query the image size
    fits_read_key(input, TINT, "NAXIS1", &fd->cols, NULL, &status);
    fits_read_key(input, TINT, "NAXIS2", &fd->rows, NULL, &status);
    if (status)
    {
        print_fits_error();
        error_jump(error, ret, "querying NAXIS failed");
    }

    fd->data = (double *)malloc(fd->cols*fd->rows*sizeof(double));
    if (!fd->data)
        error_jump(error, ret, "Error allocating frame data");

    if (fits_read_pix(input, TDOUBLE, (long []){1, 1}, fd->cols*fd->rows, 0, fd->data, NULL, &status))
        error_jump(error, ret, "fits_read_pix failed");

    // Load header keys
    int keyword_count = 0;
    if (fits_get_hdrspace(input, &keyword_count, NULL, &status))
    {
        print_fits_error();
        error_jump(error, ret, "fits_get_hdrspace failed");
    }

    for (size_t i = 0; i < keyword_count; i++)
    {
        char card[FLEN_CARD];
        char key[FLEN_KEYWORD];
        char value[FLEN_VALUE];
        char comment[FLEN_COMMENT];

        // We're only interested in user keys
        if (fits_read_record(input, i + 1, card, &status))
        {
            print_fits_error();
            error_jump(error, ret, "Error reading card %zu", i);
        }

        // Ignore keywords we aren't interested in
        int class = fits_get_keyclass(card);
        if (class < TYP_WCS_KEY || class == TYP_COMM_KEY)
            continue;

        if (fits_read_keyn(input, i + 1, key, value, comment, &status))
        {
            print_fits_error();
            error_jump(error, ret, "Error reading key %zu", i);
        }

        // Parse value
        char type;
        if (fits_get_keytype(value, &type, &status))
        {
            print_fits_error();
            error("Unable to determine type for key '%s' - skipping.", key);
            continue;
        }

        switch (type)
        {
            case 'I':
            {
                LONGLONG val;
                if (ffc2jj(value, &val, &status))
                {
                    print_fits_error();
                    error_jump(error, ret, "Error parsing '%s' as integer", value);
                }

                framedata_put_metadata(fd, key, FRAME_METADATA_INT, &(int64_t){val}, comment);
                break;
            }
            case 'F':
            {
                double val;
                if (ffc2dd(value, &val, &status))
                {
                    print_fits_error();
                    error_jump(error, ret, "Error parsing '%s' as double", value);
                }

                framedata_put_metadata(fd, key, FRAME_METADATA_DOUBLE, &val, comment);
                break;
            }
            case 'L':
            {
                int val;
                if (ffc2ll(value, &val, &status))
                {
                    print_fits_error();
                    error_jump(error, ret, "Error parsing '%s' as boolean", value);
                }

                framedata_put_metadata(fd, key, FRAME_METADATA_BOOL, &(bool){val}, comment);
                break;
            }
            default:
            {
                char val[FLEN_VALUE];
                if (ffc2s(value, val, &status))
                {
                    print_fits_error();
                    error_jump(error, ret, "Error parsing '%s' as type '%c'", value, type);
                }

                framedata_put_metadata(fd, key, FRAME_METADATA_STRING, val, comment);
                break;
            }
        }
    }

    fits_close_file(input, &status);
    if (padded_data)
        free(padded_data);

    return fd;

error:
    fits_close_file(input, &status);
    framedata_free(fd);
    if (padded_data)
        free(padded_data);

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
    {
        print_fits_error();
        error("fits_write_img failed with status %d", status);
    }

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
                    snprintf(buf, 100, "%" PRIi64, metadata->value.i);
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
            // Allow int and string to be returned as double
            switch (metadata->type)
            {
                case FRAME_METADATA_DOUBLE: *(double *)data = metadata->value.d; break;
                case FRAME_METADATA_INT: *(double *)data = (double)metadata->value.i; break;
                case FRAME_METADATA_STRING:
                {
                    // Empty strings should be treated as a missing key
                    if (strlen(metadata->value.s) == 0)
                        return FRAME_METADATA_MISSING;

                    *(double *)data = atof(metadata->value.s);
                    break;
                }
            }

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

int framedata_subtract_normalized(framedata *fd, framedata *other)
{
    if (fd->cols != other->cols || fd->rows != other->rows)
        return error("Frame size mismatch");

    double exp;
    if (framedata_get_metadata(fd, "EXPTIME", FRAME_METADATA_DOUBLE, &exp))
        return error("Unable to fetch EXPTIME");

    double other_exp = 0;
    if (framedata_get_metadata(other, "EXPTIME", FRAME_METADATA_DOUBLE, &other_exp))
        return error("Unable to fetch other EXPTIME");

    for (size_t i = 0; i < fd->cols*fd->rows; i++)
        fd->data[i] -= exp/other_exp*other->data[i];

    return 0;
}

void framedata_subtract_bias(framedata *fd)
{
    uint16_t br[4];
    framedata_bias_region(fd, br);
    double mean_bias = region_mean(br, fd->data, fd->cols);
    for (size_t i = 0; i < fd->rows*fd->cols; i++)
        fd->data[i] -= mean_bias;
}

int framedata_calibrate(framedata *frame, framedata *bias, framedata *dark, framedata *flat)
{
    if (bias)
        framedata_subtract(frame, bias);
    else
        framedata_subtract_bias(frame);

    if (dark)
        if (framedata_subtract_normalized(frame, dark))
            return error("Error dark-subtracting frame");

    if (flat)
        if (framedata_divide(frame, flat))
            return error("Error flat-fielding frame");

    return 0;
}

int framedata_calibrate_load(framedata *frame, const char *bias_path, const char *dark_path, const char *flat_path)
{
    int ret = 0;
    framedata *bias = NULL;
    framedata *dark = NULL;
    framedata *flat = NULL;

    if (bias_path)
    {
        bias = framedata_load(bias_path);
        if (!bias)
            error_jump(process_error, ret, "Error loading frame %s", bias_path);
    }

    if (dark_path)
    {
        dark = framedata_load(dark_path);
        if (!dark)
            error_jump(process_error, ret, "Error loading frame %s", dark_path);
    }

    if (flat_path)
    {
        flat = framedata_load(flat_path);
        if (!flat)
            error_jump(process_error, ret, "Error loading frame %s", flat_path);
    }

    ret = framedata_calibrate(frame, bias, dark, flat);

process_error:
    framedata_free(flat);
    framedata_free(dark);

    return ret;
}

int framedata_divide(framedata *fd, framedata *other)
{
    if (fd->cols != other->cols || fd->rows != other->rows)
        return error("Frame size mismatch");

    for (size_t i = 0; i < fd->cols*fd->rows; i++)
        fd->data[i] /= other->data[i];

    return 0;
}

static double cubic_interpolate(double* p, double x)
{
    return p[1] + 0.5f*x*(p[2] - p[0] + x*(2*p[0] - 5*p[1] + 4*p[2] - p[3] + x*(3*(p[1] - p[2]) + p[3] - p[0])));
}

static double bicubic_interpolate(double p[16], double x, double y)
{
    // interpolate along the x axis
    double py[4] = {
        cubic_interpolate(&p[0], x),
        cubic_interpolate(&p[4], x),
        cubic_interpolate(&p[8], x),
        cubic_interpolate(&p[12], x)
    };

    // interpolate the interpolated control points along the y axis
    return cubic_interpolate(py, y);
}

#define clamp(x, min, max) ((x) < (min) ? (min) : (x) > (max) ? (max) : (x))

void set_external_background_map(size_t i, double bg, void *data)
{
    ((double *)data)[i] = bg;
}

void subtract_frame_background(size_t i, double bg, void *data)
{
    ((framedata *)data)->data[i] -= bg;
}

int calculate_background_map_internal(framedata *frame, uint16_t min_tile_size, void (*px_worker)(size_t, double, void *), void *worker_data)
{
    uint16_t rr[4];
    if (framedata_image_region(frame, rr))
        return error("Image region failed");

    uint16_t width = rr[1] - rr[0];
    uint16_t height = rr[3] - rr[2];

    // Extend tile width/height by a few px if necessary
    // to reduce the margin of unused pixels on the right/top edge
    uint16_t x_tile_count = width / min_tile_size;
    uint16_t y_tile_count = height / min_tile_size;
    uint16_t tile_width = width / x_tile_count;
    uint16_t tile_height = height / y_tile_count;

    // Extra margin accounts for the edge tiles that are extended to fill the remaining gap
    double *tile_buf = calloc(tile_width * tile_height, sizeof(double));
    double *tile_median = calloc(x_tile_count * y_tile_count, sizeof(double));

    if (!tile_buf || !tile_median)
    {
        free(tile_buf);
        free(tile_median);
        return error("Failed to allocate background map buffers");
    }

    // Calculate tile medians
    for (uint16_t ty = 0; ty < y_tile_count; ty++)
    {
        for (uint16_t tx = 0; tx < x_tile_count; tx++)
        {
            uint16_t x_start = rr[0] + tx * tile_width;
            uint16_t y_start = rr[2] + ty * tile_height;

            size_t i = 0;
            for (uint16_t dy = 0; dy < tile_height; dy++)
                for (uint16_t dx = 0; dx < tile_width; dx++)
                    tile_buf[i++] = frame->data[(y_start + dy) * frame->cols + x_start + dx];

            qsort(tile_buf, tile_width * tile_height, sizeof(double), compare_double);

            tile_median[ty * x_tile_count + tx] = tile_buf[tile_width * tile_height / 2];
        }
    }

    // Bicubic interpolate between the tiles to generate the background map
    // The coarse points are in the middle of the tiles, so we must enumerate
    // over an extra set of tiles on each border in order to evaluate the outer
    // half tile edges of the image.  The pixels outside the image are filtered
    // when interpolating
    for (int32_t ty = -1; ty < y_tile_count; ty++)
    {
        for (int32_t tx = -1; tx < x_tile_count; tx++)
        {
            // Calculate the 16 control points for interpolating over this tile
            double p[16];
            for (uint8_t j = 0; j < 4; j++)
                for (uint8_t i = 0; i < 4; i++)
                    p[j * 4 + i] = tile_median[clamp(ty - 1 + j, 0, y_tile_count - 1) * x_tile_count + clamp(tx - 1 + i, 0, x_tile_count - 1)];

            // Interpolate over all the pixels inside the image area inside this tile
            for (uint16_t dy = 0; dy < tile_height; dy++)
            {
                int32_t y = ty * tile_height + tile_height / 2 + dy;
                if (y < 0 || y >= height)
                    continue;

                for (uint16_t dx = 0; dx < tile_width; dx++)
                {
                    int32_t x = tx * tile_width + tile_width / 2 + dx;
                    if (x < 0 || x >= width)
                        continue;

                    double bg = bicubic_interpolate(p, dx * 1.0 / tile_width, dy * 1.0 / tile_height);
                    px_worker((rr[2] + y) * frame->cols + rr[0] + x, bg, worker_data);
                }
            }
        }
    }

    free(tile_buf);
    free(tile_median);
    return 0;
}

int framedata_calculate_background_map(framedata *frame, uint16_t min_tile_size, double *background_map)
{
    return calculate_background_map_internal(frame, min_tile_size, set_external_background_map, background_map);
}

int framedata_subtract_background_map(framedata *frame, uint16_t min_tile_size)
{
    return calculate_background_map_internal(frame, min_tile_size, subtract_frame_background, frame);
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
                printf("%" PRIi64 " (int) / ", m->value.i);
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

    // Puoko-nui
    hashmap_get(fd->metadata_map, "UTC-DATE", (void **)(&date));
    hashmap_get(fd->metadata_map, "UTC-BEG", (void **)(&time));
    if (date && date->type == FRAME_METADATA_STRING &&
        time && time->type == FRAME_METADATA_STRING)
    {
        *out_time = parse_date_time(date->value.s, time->value.s);
        return 0;
    }

    // Puoko-nui (legacy)
    hashmap_get(fd->metadata_map, "GPSTIME", (void **)(&date));
    if (date && date->type == FRAME_METADATA_STRING)
    {
        *out_time = parse_time(date->value.s);
        return 0;
    }

    // PROMPT
    hashmap_get(fd->metadata_map, "DATE", (void **)(&date));
    hashmap_get(fd->metadata_map, "TIME-OBS", (void **)(&time));
    if (date && date->type == FRAME_METADATA_STRING &&
        time && time->type == FRAME_METADATA_STRING)
    {
        *out_time = parse_date_time(date->value.s, time->value.s);
        return 0;
    }

    // Quilt
    hashmap_get(fd->metadata_map, "DATE-OBS", (void **)(&date));
    hashmap_get(fd->metadata_map, "UTC", (void **)(&time));
    if (date && date->type == FRAME_METADATA_STRING &&
        time && time->type == FRAME_METADATA_STRING)
    {
        *out_time = parse_date_time(date->value.s, time->value.s);
        return 0;
    }

    // SBIG CCD-OPS
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

int framedata_bias_region(framedata *frame, uint16_t region[4])
{
    uint16_t r[4] = {0, 0, 0, 0};
    char *str;
    if (framedata_get_metadata(frame, "BIAS-RGN", FRAME_METADATA_STRING, &str) == FRAME_METADATA_OK)
        sscanf(str, "[%hu, %hu, %hu, %hu]", &r[0], &r[1],
                                            &r[2], &r[3]);
    memcpy(region, r, 4*sizeof(uint16_t));
    return 0;
}

static int sum_into_axes(framedata *frame, uint16_t region[4], double **x, double **y, bool rotate)
{
    int ret = 0;
    uint16_t xs = region[1] - region[0];
    uint16_t ys = region[3] - region[2];

    double *xi = calloc(xs, sizeof(double));
    double *xa = calloc(xs, sizeof(double));
    double *yi = calloc(ys, sizeof(double));
    double *ya = calloc(ys, sizeof(double));
    if (!xi || !xa || !yi || !ya)
        error_jump(error, ret, "allocation_failed");

    // Sum into x and y axes
    for (uint16_t j = 0; j < ys; j++)
    {
        yi[j] = j;
        uint16_t jj = rotate ? ys - j - 1 : j;
        for (uint16_t i = 0; i < xs; i++)
        {
            uint16_t ii = rotate ? xs - i - 1 : i;
            double px = frame->data[frame->cols*(region[2] + jj) + (ii + region[0])];
            xa[i] += px;
            ya[j] += px;
        }
    }

    for (uint16_t i = 0; i < xs; i++)
        xi[i] = i;

    // Subtract the mean background level then clamp to zero to
    // discard the majority of the sky background
    double xm = mean_exclude_sigma(xa, xs, 1);
    for (uint16_t i = 0; i < xs; i++)
        xa[i] = fmax(xa[i] - xm, 0);

    double ym = mean_exclude_sigma(ya, ys, 1);
    for (uint16_t j = 0; j < ys; j++)
        ya[j] = fmax(ya[j] - ym, 0);

    *x = xa;
    *y = ya;

    free(xi);
    free(yi);

    return 0;
error:
    free(xi);
    free(xa);
    free(yi);
    free(ya);
    x = NULL;
    y = NULL;
    return ret;
}

double find_max_correlatation(double *a, double *b, uint16_t n)
{
    int32_t best_idx = 0;
    double best = -1;
    for (int32_t j = -n; j < n; j++)
    {
        double ret = 0;
        for (uint16_t i = 0; i < n; i++)
        {
            // Assume data outside the array bounds are zero
            uint16_t k = i + j;
            ret += (k > 0 && k < n) ? a[i]*b[k] : 0;
        }

        if (ret > best)
        {
            best = ret;
            best_idx = j;
        }
    }

    // Best fit is at the edge of the correlation.
    // Can't improve our estimate!
    if (best_idx == -n || best_idx == n - 1)
        return best_idx;
    
    // Improve estimate by doing a quadratic fit of
    // the three points around the peak
    double offset[3] = { -1, 0, 1 };
    double corr[3] = { 0, 1, 0 };
    double p[3];
    for (uint16_t i = 0; i < n; i++)
    {
        corr[0] += a[i]*b[i + best_idx - 1];
        corr[2] += a[i]*b[i + best_idx + 1];
    }

    corr[0] /= best;
    corr[2] /= best;

    if (fit_polynomial(offset, corr, NULL, 3, p, 2))
        return best_idx;

    return best_idx - p[1] / (2 * p[2]);
}

int framedata_estimate_translation(framedata *frame, framedata *reference, double *xt, double *yt, bool *rotated)
{
    int ret = 0;

    uint16_t fr[4], rr[4];
    if (framedata_image_region(frame, fr))
        error_jump(error, ret, "Image region failed");
    if (framedata_image_region(frame, rr))
        error_jump(error, ret, "Image region failed");

    if (fr[0] != rr[0] || fr[1] != rr[1] || fr[2] != rr[2] || fr[3] != rr[3])
        error_jump(error, ret, "frame sizes don't match");

    // The PROMPT telescopes can rotate their frames during a run
    double cd11, cd22, refcd11, refcd22;
    bool checkrotation = !framedata_get_metadata(reference, "CD1_1", FRAME_METADATA_DOUBLE, &refcd11) &&
        !framedata_get_metadata(reference, "CD2_2", FRAME_METADATA_DOUBLE, &refcd22) &&
        !framedata_get_metadata(frame, "CD1_1", FRAME_METADATA_DOUBLE, &cd11) &&
        !framedata_get_metadata(frame, "CD2_2", FRAME_METADATA_DOUBLE, &cd22);

    if (checkrotation && signbit(cd11) != signbit(refcd11) && signbit(cd22) != signbit(refcd22))
        *rotated = true;

    double *fx, *fy, *rx, *ry;
    if (sum_into_axes(frame, fr, &fx, &fy, *rotated))
        error_jump(error, ret, "Summing into axes failed");
    if (sum_into_axes(reference, rr, &rx, &ry, false))
        error_jump(error, ret, "Summing into axes failed");

    *xt = find_max_correlatation(rx, fx, rr[1] - rr[0]);
    *yt = find_max_correlatation(ry, fy, rr[3] - rr[2]);

    free(fx);
    free(fy);
    free(rx);
    free(ry);

    return ret;
error:
    *xt = 0;
    *yt = 0;
    return ret;
}

