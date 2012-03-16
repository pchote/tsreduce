/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#include "framedata.h"
#include "helpers.h"

framedata framedata_new(const char *filename, framedata_type dtype)
{
    framedata this;
	int status = 0;
    if (fits_open_image(&this._fptr, filename, READONLY, &status))
    {
        char fitserr[128];
        while (fits_read_errmsg(fitserr))
            fprintf(stderr, "%s\n", fitserr);

        die("fits_open_image failed with error %d; %s", status, filename);
    }
    // Query the image size
    fits_read_key(this._fptr, TINT, "NAXIS1", &this.cols, NULL, &status);
    fits_read_key(this._fptr, TINT, "NAXIS2", &this.rows, NULL, &status);
    if (status)
        die("querying NAXIS failed");
    
    long fpixel[2] = {1,1}; // Read the entire image
    
    this.dtype = dtype;
    if (dtype == FRAMEDATA_INT)
    {
        this.dbl_data = NULL;
        this.data = (int *)malloc(this.cols*this.rows*sizeof(int));
        if (this.data == NULL)
            die("malloc failed");
        
        if (fits_read_pix(this._fptr, TINT, fpixel, this.cols*this.rows, 0, this.data, NULL, &status))
            die("fits_read_pix failed");
    }
    
    else if (dtype == FRAMEDATA_DBL)
    {
        this.data = NULL;
        this.dbl_data = (double *)malloc(this.cols*this.rows*sizeof(double));
        if (this.dbl_data == NULL)
            die("malloc failed");
        
        if (fits_read_pix(this._fptr, TDOUBLE, fpixel, this.cols*this.rows, 0, this.dbl_data, NULL, &status))
            die("fits_read_pix failed");
    }

    // Load image regions
    // TODO: have acquisition software save regions into a header key
    this.regions.has_overscan = (this.cols != this.rows);
    int *ir = this.regions.image_region;
    int *br = this.regions.bias_region;
    if (this.regions.has_overscan)
    {
        ir[0] = 0; ir[1] = 512; ir[2] = 0; ir[3] = 512;
        br[0] = 525; br[1] = 535; br[2] = 5; br[3] = 508;
    }
    else
    {
        ir[0] = 0; ir[1] = this.cols; ir[2] = 0; ir[3] = this.rows;
        br[0] = br[1] = br[2] = br[3] = 0;
    }

    this.regions.image_px = (ir[1] - ir[0])*(ir[3] - ir[2]);
    this.regions.bias_px = (br[1] - br[0])*(br[3] - br[2]);
    return this;
}

int framedata_get_header_int(framedata *this, const char *key)
{
    int ret, status = 0;
    if (fits_read_key(this->_fptr, TINT, key, &ret, NULL, &status))
        die("framedata_get_header_int failed");
    return ret;
}

int framedata_has_header_string(framedata *this, const char *key)
{
    int status = 0;
    char buf[128];
    fits_read_key(this->_fptr, TSTRING, key, &buf, NULL, &status);
    return status != KEY_NO_EXIST;
}

void framedata_get_header_string(framedata *this, const char *key, char *ret)
{
    int status = 0;
    if (fits_read_key(this->_fptr, TSTRING, key, ret, NULL, &status))
        die("framedata_get_header_string failed");
}

void framedata_subtract(framedata *this, framedata *other)
{
    if (this->cols != other->cols || this->rows != other->rows)
        die("Attempting to subtract frame with different size");
    
    if (this->dtype == FRAMEDATA_INT)
        for (int i = 0; i < this->cols*this->rows; i++)
            this->data[i] -= other->data[i];
    else if (this->dtype == FRAMEDATA_DBL)
        for (int i = 0; i < this->cols*this->rows; i++)
            this->dbl_data[i] -= other->dbl_data[i];
    else 
        die("Unknown datatype");
}

void framedata_add(framedata *this, framedata *other)
{
    if (this->cols != other->cols || this->rows != other->rows)
        die("Attempting to add frame with different size");
    
    if (this->dtype == FRAMEDATA_INT)
        for (int i = 0; i < this->cols*this->rows; i++)
            this->data[i] += other->data[i];
    else if (this->dtype == FRAMEDATA_DBL)
        for (int i = 0; i < this->cols*this->rows; i++)
            this->dbl_data[i] += other->dbl_data[i];
    else 
        die("Unknown datatype");
}

void framedata_multiply(framedata *this, int mul)
{
    if (this->dtype == FRAMEDATA_INT)
        for (int i = 0; i < this->cols*this->rows; i++)
            this->data[i] *= mul;
    else if (this->dtype == FRAMEDATA_DBL)
        for (int i = 0; i < this->cols*this->rows; i++)
            this->dbl_data[i] *= mul;
    else 
        die("Unknown datatype");
}

void framedata_divide_const(framedata *this, int div)
{
    if (this->dtype == FRAMEDATA_INT)
        for (int i = 0; i < this->cols*this->rows; i++)
            this->data[i] /= div;
    else if (this->dtype == FRAMEDATA_DBL)
        for (int i = 0; i < this->cols*this->rows; i++)
            this->dbl_data[i] /= div;
    else 
        die("Unknown datatype");
}

void framedata_divide(framedata *this, framedata *div)
{
    if (this->dtype == FRAMEDATA_INT)
        for (int i = 0; i < this->cols*this->rows; i++)
            this->data[i] /= div->data[i];
    else if (this->dtype == FRAMEDATA_DBL)
        for (int i = 0; i < this->cols*this->rows; i++)
            this->dbl_data[i] /=  div->dbl_data[i];
    else 
        die("Unknown datatype");
}

void framedata_free(framedata this)
{
    int status;
    if (this.data)
        free(this.data);
    if (this.dbl_data)
        free(this.dbl_data);
    fits_close_file(this._fptr, &status);
}

// Convenience function for calculating the mean signal in a sub-region of a frame
// Assumes that the frame type is double, and that the region is inside the frame
double mean_in_region(framedata *frame, int rgn[4])
{
    int num_px = (rgn[1] - rgn[0])*(rgn[3] - rgn[2]);
    double mean = 0;
    for (int j = rgn[2]; j < rgn[3]; j++)
        for (int i = rgn[0]; i < rgn[1]; i++)
            mean += frame->dbl_data[j*frame->cols + i]/num_px;
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
        frame->dbl_data[i] -= mean_bias;
}