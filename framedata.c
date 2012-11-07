/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#include <fitsio.h>

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
    return fp;
}

framedata *framedata_load(const char *filename)
{
    framedata *fp = framedata_alloc();
    if (!fp)
    {
        error("framedata_alloc failed");
        return NULL;
    }

	int status = 0;
    if (fits_open_image((fitsfile **)(&fp->fptr), filename, READONLY, &status))
    {
        char fitserr[128];
        while (fits_read_errmsg(fitserr))
            error("%s\n", fitserr);

        error("fits_open_image failed with error %d; %s", status, filename);
        framedata_free(fp);
        return NULL;
    }

    // Query the image size
    fits_read_key(fp->fptr, TINT, "NAXIS1", &fp->cols, NULL, &status);
    fits_read_key(fp->fptr, TINT, "NAXIS2", &fp->rows, NULL, &status);
    if (status)
    {
        error("querying NAXIS failed");
        framedata_free(fp);
        return NULL;
    }

    fp->data = (double *)malloc(fp->cols*fp->rows*sizeof(double));
    if (fp->data == NULL)
    {
        error("malloc failed");
        framedata_free(fp);
        return NULL;
    }

    if (fits_read_pix(fp->fptr, TDOUBLE, (long []){1, 1}, fp->cols*fp->rows, 0, fp->data, NULL, &status))
    {
        error("fits_read_pix failed");
        framedata_free(fp);
        return NULL;
    }

    // Load image regions
    // TODO: have acquisition software save regions into a header key
    fp->regions.has_overscan = (fp->cols != fp->rows);
    int *ir = fp->regions.image_region;
    int *br = fp->regions.bias_region;
    if (fp->regions.has_overscan)
    {
        ir[0] = 0; ir[1] = 512; ir[2] = 0; ir[3] = 512;
        br[0] = 525; br[1] = 535; br[2] = 5; br[3] = 508;
    }
    else
    {
        ir[0] = 0; ir[1] = fp->cols; ir[2] = 0; ir[3] = fp->rows;
        br[0] = br[1] = br[2] = br[3] = 0;
    }

    fp->regions.image_px = (ir[1] - ir[0])*(ir[3] - ir[2]);
    fp->regions.bias_px = (br[1] - br[0])*(br[3] - br[2]);
    return fp;
}

int framedata_get_header_long(framedata *this, const char *key, long *value)
{
    int ret, status = 0;
    if (fits_read_key(this->fptr, TLONG, key, &ret, NULL, &status))
        return 1;

    *value = ret;
    return 0;
}

int framedata_get_header_dbl(framedata *this, const char *key, double *value)
{
    double ret;
    int status = 0;
    if (fits_read_key(this->fptr, TDOUBLE, key, &ret, NULL, &status))
        return 1;

    *value = ret;
    return 0;
}

int framedata_has_header_string(framedata *this, const char *key)
{
    int status = 0;
    char buf[128];
    fits_read_key(this->fptr, TSTRING, key, &buf, NULL, &status);
    return status != KEY_NO_EXIST;
}

void framedata_get_header_string(framedata *this, const char *key, char *ret)
{
    int status = 0;
    if (fits_read_key(this->fptr, TSTRING, key, ret, NULL, &status))
        die("framedata_get_header_string failed for key %s", key);
}

void framedata_subtract(framedata *this, framedata *other)
{
    if (this->cols != other->cols || this->rows != other->rows)
        die("Attempting to subtract frame with different size");

    for (int i = 0; i < this->cols*this->rows; i++)
        this->data[i] -= other->data[i];
}

void framedata_add(framedata *this, framedata *other)
{
    if (this->cols != other->cols || this->rows != other->rows)
        die("Attempting to add frame with different size");

    for (int i = 0; i < this->cols*this->rows; i++)
        this->data[i] += other->data[i];
}

void framedata_multiply(framedata *this, int mul)
{
    for (int i = 0; i < this->cols*this->rows; i++)
        this->data[i] *= mul;
}

void framedata_divide_const(framedata *this, int div)
{
    for (int i = 0; i < this->cols*this->rows; i++)
        this->data[i] /= div;
}

void framedata_divide(framedata *this, framedata *div)
{
    for (int i = 0; i < this->cols*this->rows; i++)
        this->data[i] /=  div->data[i];
}

void framedata_free(framedata *frame)
{
    if (!frame)
        return;

    int status;
    if (frame->data)
        free(frame->data);

    if (frame->fptr)
        fits_close_file(frame->fptr, &status);

    free(frame);
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
