/* Copyright 2010-2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE. */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <fitsio.h>

#include "framedata.h"
#include "reduction.h"

typedef enum 
{
    ADD,
    SUBTRACT,
    AVERAGE,
    MULTIPLY,
    DIVIDE,
} Mode;

// Subtracts the dark count, then normalizes the frame to average to unity
// (not that we actually care about the total photometric counts)
void normalize_flat(framedata *flat, void *data)
{
    framedata *dark = (framedata *)data;
    
    if (dark->dtype != FRAMEDATA_DBL || flat->dtype != FRAMEDATA_DBL)
        error("normalize_flat frames must be type DBL");
    
    if (dark->rows != flat->rows || dark->cols != flat->rows)
        error("normalize_flat frames must have same size");
    
    int flatexp = framedata_get_header_int(flat, "EXPTIME");
    int darkexp = framedata_get_header_int(dark, "EXPTIME");
    int n = flat->rows*flat->cols;
    
    // Calculate mean
    double mean = 0;
    for (int i = 0; i < n; i++)
    {
        flat->dbl_data[i] -= flatexp*1.0/darkexp*dark->dbl_data[i];
        mean += flat->dbl_data[i];
    }
    mean /= n;

    // Calculate standard deviation
    double std = 0;
    for (int i = 0; i < n; i++)
        std += (flat->dbl_data[i] - mean)*(flat->dbl_data[i] - mean);
    std = sqrt(std/n);
    
    // Recalculate the mean, excluding outliers at 3 sigma
    double mean_new = 0;
    for (int i = 0; i < n; i++)
        if (fabs(flat->dbl_data[i] - mean) < 3*std)
            mean_new += flat->dbl_data[i];
    mean_new /= n;

    // Normalize flat to unity counts
    for (int i = 0; i < n; i++)
        flat->dbl_data[i] /= mean_new;
}

#define MAX_FRAMES 100
// Create a flat field frame from the frames listed by the command `flatcmd',
// rejecting `minmax' highest and lowest pixels, with the dark frame `masterdark'
// saved to the file named `outname'
int create_flat(const char *flatcmd, int minmax, const char *masterdark, const char *outname)
{
    char **frames = (char **)malloc(MAX_FRAMES*sizeof(char*));
    for (int i = 0; i < MAX_FRAMES; i++)
        frames[i] = (char *)malloc(PATH_MAX*sizeof(char));
    
    int numflats = get_matching_files(flatcmd, frames, PATH_MAX, MAX_FRAMES);
    if (numflats < 2*minmax)
        error("Insufficient frames. %d found, %d will be discarded", numflats, 2*minmax);

    double *flat = (double *)malloc(512*512*sizeof(double));
    
    framedata dark = framedata_new(masterdark, FRAMEDATA_DBL);
        
    // Load the flat frames, discarding the 5 outermost pixels for each
    load_reject_minmax( (const char **)frames, numflats, dark.rows, dark.cols, minmax, minmax, flat, normalize_flat, (void *)&dark);
    framedata_free(dark);
    
    // Create a new fits file
    fitsfile *out;
    int status = 0;
    char outbuf[2048];
    sprintf(outbuf, "!%s", outname);
    fits_create_file(&out, outbuf, &status);
    
    /* Create the primary array image (16-bit short integer pixels */
	long size[2] = { 512, 512 };
	fits_create_img(out, DOUBLE_IMG, 2, size, &status);
    
    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, 512*512, flat, &status))
        error("fits_write_img failed with status %d", status);
    
    fits_close_file(out, &status);
    free(flat);
    
    for (int i = 0; i < MAX_FRAMES; i++)
        free(frames[i]);
    free(frames);

    return 0;
}

void do_nothing(framedata *flat, void *data) {}
int create_dark(const char *darkcmd, int minmax, const char *outname)
{
    char **frames = (char **)malloc(MAX_FRAMES*sizeof(char*));
    for (int i = 0; i < MAX_FRAMES; i++)
        frames[i] = (char *)malloc(PATH_MAX*sizeof(char));
    
    int numdarks = get_matching_files(darkcmd, frames, PATH_MAX, MAX_FRAMES);
    if (numdarks < 2*minmax)
        error("Insufficient frames. %d found, at least %d are required", numdarks, 2*minmax);
    
    framedata base = framedata_new(frames[0], FRAMEDATA_INT);
    int exptime = framedata_get_header_int(&base, "EXPTIME");
    double *dark = (double *)malloc(base.rows*base.cols*sizeof(double));
    

    // Load the flat frames, discarding the 5 outermost pixels for each
    load_reject_minmax( (const char **)frames, numdarks, base.rows, base.cols, minmax, minmax, dark, do_nothing, NULL);
    
    // Create a new fits file
    fitsfile *out;
    int status = 0;
    char outbuf[2048];
    sprintf(outbuf, "!%s", outname);
    
    fits_create_file(&out, outbuf, &status);
    
    /* Create the primary array image (16-bit short integer pixels */
	long size[2] = { 512, 512 };
	fits_create_img(out, DOUBLE_IMG, 2, size, &status);
    
    fits_update_key(out, TINT, "EXPTIME", &exptime, "Actual integration time (sec)", &status);
    
    // Write the frame data to the image
    if (fits_write_img(out, TDOUBLE, 1, base.rows*base.cols, dark, &status))
        error("fits_write_img failed with status %d", status);

    fits_close_file(out, &status);
    free(dark);
    
    framedata_free(base);

    for (int i = 0; i < MAX_FRAMES; i++)
        free(frames[i]);
    free(frames);
    
    return 0;
}

int main( int argc, char *argv[] )
{
    // First arg gives the mode (add / avg)
    // Second arg gives the file to take the metadata from
    // Last arg gives the file to save to
    
    if (argc > 3)
    {                
        // `frametool create-flat "/bin/ls dome-*.fits.gz" 5 master-dark.fits.gz master-dome.fits.gz`
        if (argc == 6 && strncmp(argv[1], "create-flat", 11) == 0)
        {
            create_flat(argv[2], atoi(argv[3]), argv[4], argv[5]);
            return 0;
        }
        
        // `frametool create-dark "/bin/ls dark-*.fits.gz" 5 master-dark.fits.gz`
        if (argc == 5 && strncmp(argv[1], "create-dark", 11) == 0)
        {
            create_dark(argv[2], atoi(argv[3]), argv[4]);
            return 0;
        }
        // `frametool reduce rawframe.fits.gz master-dark.fits.gz master-flat.fits.gz reduced.fits.gz`
        if (argc == 6 && strncmp(argv[1], "reduce", 6) == 0)
        {
            framedata base = framedata_new(argv[2], FRAMEDATA_DBL);
            framedata dark = framedata_new(argv[3], FRAMEDATA_DBL);
            framedata flat = framedata_new(argv[4], FRAMEDATA_DBL);
            // Subtract dark counts
            if (dark.dbl_data != NULL)
                framedata_subtract(&base, &dark);
            
            // Flat field image
            if (flat.dbl_data != NULL)
                framedata_divide(&base, &flat);
            
            // Create a new fits file           
            fitsfile *out;
            int status = 0;
            char outbuf[2048];
            sprintf(outbuf, "!%s", argv[5]);            
            fits_create_file(&out, outbuf, &status);
            
            
            for (int i = 0; i < base.rows*base.cols; i++)
                printf("%f\n",base.dbl_data[i]);
            
            /* Create the primary array image (16-bit short integer pixels */
            long size[2] = { 512, 512 };
            fits_create_img(out, DOUBLE_IMG, 2, size, &status);
                        
            // Write the frame data to the image
            if (fits_write_img(out, TDOUBLE, 1, base.rows*base.cols, base.dbl_data, &status))
                error("fits_write_img failed with status %d", status);
            
            fits_close_file(out, &status);

            
            
            printf("Saved to %s\n", argv[argc-1]);

            return 0;
        }
        
        Mode mode;
        if (strncmp(argv[1], "add",3) == 0)
            mode = ADD;
        else if (strncmp(argv[1], "subtract",7) == 0)
            mode = SUBTRACT;
        else if(strncmp(argv[1], "avg",3) == 0)
            mode = AVERAGE;
        else if(strncmp(argv[1], "multiply",8) == 0)
            mode = MULTIPLY; 
        else if(strncmp(argv[1], "divide",6) == 0)
            mode = DIVIDE; 
        else
            error("Invalid mode");
        
        printf("Opening file %s\n", argv[2]);
        framedata base = framedata_new(argv[2], FRAMEDATA_INT);
        
        
        if (mode == DIVIDE)
        {
            framedata_divide_const(&base, atoi(argv[3]));
        }
        else if (mode == MULTIPLY)
        {
            framedata_multiply(&base, atoi(argv[3]));
        }
        else // add, subtract, average
        {
            for (int i = 3; i < argc-1; i++)
            {
                printf("Adding file %s\n", argv[i]);
                framedata other = framedata_new(argv[i], FRAMEDATA_INT);
                if (mode == SUBTRACT)
                    framedata_subtract(&base, &other);
                else
                    framedata_add(&base, &other);
                framedata_free(other);
            }
            
            if (mode == AVERAGE)
            {
                for (int i = 0; i < base.cols*base.rows; i++)
                    base.data[i] /= argc - 3;
            }
        }
        
        fitsfile *out;
    	int status = 0;
    	
    	// Create a new fits file
        char outbuf[2048];
        sprintf(outbuf, "!%s(%s)", argv[argc-1], argv[2]);
    	fits_create_file(&out, outbuf, &status);

        // Write the frame data to the image
    	if (fits_write_img(out, TINT, 1, base.cols*base.rows, base.data, &status))
    	{
            printf("status: %d\n", status);
    	    error("fits_write_img failed");
    	}    
    	
    	fits_close_file(out, &status);
        printf("Saved to %s\n", argv[argc-1]);

        framedata_free(base);
    }
    else
        error("Invalid args");
    return 0;
}
