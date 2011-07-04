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
#include <libgen.h>

#include "framedata.h"
#include "reduction.h"

double2 process_target(target r, framedata *frame, double exptime)
{
    double2 ret = {0,0};
    double2 bg, xy;
    target last;
    int n = 0;
    // Iterate until we move less than 2px, or reach 10 iterations
    do
    {
        last = r;
        bg = calculate_background(r, frame);
        xy = center_aperture(r, bg, frame);
        printf("%d: (%f,%f) -> (%f,%f) [%f,%f]\n", n, r.x, r.y, xy.x, xy.y, bg.x, bg.y);

        r.x = (int)xy.x;
        r.y = (int)xy.y;
        
        if (r.x - r.r < 0 || r.x + r.r >= frame->cols || r.y - r.r < 0 || r.y + r.r >= frame->rows)
        {
            fprintf(stderr, "Aperture outside chip - skipping\n");
            return ret;
        }
    } while (++n < 10 && (xy.x-last.x)*(xy.x-last.x) + (xy.y-last.y)*(xy.y-last.y) >= 4);
    
    if (n >= 10)
    {
        printf("Aperture centering did not converge - skipping\n");
        return ret;
    }

    ret.y = bg.x*M_PI*r.r*r.r / exptime;
    ret.x = integrate_aperture(xy,r.r, frame) / exptime - ret.y;

    return ret;
}

#define MAX_OBS 10000
int main( int argc, char *argv[] )
{
    if (argc > 1)
    {
        target targets[10];
        int numtargets = 0;
        
        // basename and dirname may modify underlying memory
        char dirbuf[PATH_MAX];
        char namebuf[PATH_MAX];
        strncpy(dirbuf, argv[1], sizeof(dirbuf));
        strncpy(namebuf, argv[1], sizeof(namebuf));
        
        // change into the data file dir so FrameDir: is correct
        chdir(dirname(dirbuf));
        FILE *data = fopen(basename(namebuf), "r+");
        if (data == NULL)
            error("Error opening data file");
        
        record records[MAX_OBS];
        int numrecords = 0;
        
        char pattern[128];
        pattern[0] = '\0';
        
        time_t start_time = 0;
        
        // Processed dark template
        framedata dark;
        dark.data = NULL;
        dark.dbl_data = NULL;

        framedata flat;
        flat.dbl_data = NULL;
        
        char linebuf[1024];
        // Read any existing config and data
        while (fgets(linebuf, sizeof(linebuf)-1, data) != NULL)
        {            
            if (!strncmp(linebuf,"# FramePattern:", 15))
                sscanf(linebuf, "# FramePattern: %s\n", pattern);
            else if (!strncmp(linebuf,"# FrameDir:", 11))
            {
                char datadir[PATH_MAX];
                sscanf(linebuf, "# FrameDir: %s\n", datadir);
                chdir(datadir);
            }
            else if (!strncmp(linebuf,"# DarkTemplate:", 15))
            {
                char darkbuf[128];
                sscanf(linebuf, "# DarkTemplate: %s\n", darkbuf);
                
                dark = framedata_new(darkbuf, FRAMEDATA_DBL);
            }
            else if (!strncmp(linebuf,"# FlatTemplate:", 15))
            {
                char flatbuf[128];
                sscanf(linebuf, "# FlatTemplate: %s\n", flatbuf);
                
                flat = framedata_new(flatbuf, FRAMEDATA_DBL);
            }
            else if (!strncmp(linebuf,"# Target:", 9))
            {
                 sscanf(linebuf, "# Target: (%lf, %lf, %lf, %lf, %lf)\n",
                    &targets[numtargets].x,
                    &targets[numtargets].y,
                    &targets[numtargets].r,
                    &targets[numtargets].s1,
                    &targets[numtargets].s2
                 );
                numtargets++;
            }
            else if (!strncmp(linebuf,"# ReferenceTime:", 16))
            {
                struct tm t;
                strptime(linebuf, "# ReferenceTime: %Y-%m-%d %H:%M:%S\n", &t);
                start_time = timegm(&t);
            }
            else if (linebuf[0] == '#')
                continue;

            sscanf(linebuf, "%lf %lf %lf %lf %lf %lf %lf %lf %s\n",
                &records[numrecords].time,
                &records[numrecords].star[0],
                &records[numrecords].sky[0],
                &records[numrecords].star[1],
                &records[numrecords].sky[1],
                &records[numrecords].star[2],
                &records[numrecords].sky[2],
                &records[numrecords].ratio,
                records[numrecords].filename
            );
            numrecords++;
        }
                
        // Iterate through the list of files matching the filepattern
        char filenamebuf[PATH_MAX+8];
        sprintf(filenamebuf, "/bin/ls %s", pattern);
        FILE *ls = popen(filenamebuf, "r");
        if (ls == NULL)
            error("failed to list directory");

        while (fgets(filenamebuf, sizeof(filenamebuf)-1, ls) != NULL)
        {
            // Strip the newline character from the end of the filename
            filenamebuf[strlen(filenamebuf)-1] = '\0';
            char *filename = filenamebuf;

            // Check whether the frame has been processed
            int processed = FALSE;
            for (int i = 0; i < numrecords; i++)
                if (strcmp(filename, records[i].filename) == 0)
                {
                    processed = TRUE;
                    break;
                }
            printf("%s: processed %d\n", filename, processed);
            
            if (processed)
                continue;
                     
            framedata frame = framedata_new(filename, FRAMEDATA_DBL);
            
            int exptime = framedata_get_header_int(&frame, "EXPTIME");

            // Calculate time at the middle of the exposure relative to ReferenceTime
            double midtime = 0;
            if (framedata_has_header_string(&frame, "UTC-BEG"))
            {
                char datebuf[128], timebuf[128], datetimebuf[257];
                struct tm t;
                framedata_get_header_string(&frame, "UTC-DATE", datebuf);
                framedata_get_header_string(&frame, "UTC-BEG", timebuf);
                
                sprintf(datetimebuf, "%s %s", datebuf, timebuf);
                strptime(datetimebuf, "%Y-%m-%d %H:%M:%S", &t);
                time_t frame_time = timegm(&t);
                double start = difftime(frame_time, start_time);
                
                framedata_get_header_string(&frame, "UTC-END", timebuf);
                sprintf(datetimebuf, "%s %s", datebuf, timebuf);
                strptime(datetimebuf, "%Y-%m-%d %H:%M:%S", &t);
                frame_time = timegm(&t);
                double end = difftime(frame_time, start_time);
                
                midtime = (start + end)/2;
            }
            // Handle oldstyle frames
            else if (framedata_has_header_string(&frame, "GPSTIME"))
            {
                char datebuf[128];
                struct tm t;
                framedata_get_header_string(&frame, "GPSTIME", datebuf);
                strptime(datebuf, "%Y-%m-%d %H:%M:%S", &t);
                time_t frame_time = timegm(&t);

                midtime = difftime(frame_time, start_time) + exptime/2;
            }
            else
            {
                fprintf(stderr, "%s: no valid time header found - skipping\n", filename);
                continue;
            }
            
            // Subtract dark counts
            if (dark.dbl_data != NULL)
                framedata_subtract(&frame, &dark);
            
            // Flat field image
            if (flat.dbl_data != NULL)
                framedata_divide(&frame, &flat);
            
            // Process targets
            fprintf(data, "%f ", midtime);
            double2 target = process_target(targets[0], &frame, exptime);
            fprintf(data, "%f %f ", target.x, target.y);
            double comparison = 0;
            for (int i = 1; i < numtargets; i++)
            {
                double2 ret = process_target(targets[i], &frame, exptime);
                comparison += ret.x;
                fprintf(data, "%f %f ", ret.x, ret.y);
            }
            
            fprintf(data, "%f %s\n", target.x / comparison, filename);
            
            framedata_free(frame);
        }
        if (dark.data != NULL)
            framedata_free(dark);
        
        pclose(ls);
        fclose(data);
    }
    return 0;
}
