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

#include "aperture.h"
#include "helpers.h"
#include "framedata.h"

// Find the center of the star within the inner circle
//   Takes the search circle and imagedata
//   Returns x,y coordinates for the star center
double2 center_aperture(target reg, double2 bg2, framedata *frame)
{
    // Round to the nearest pixel
    int x = (int)reg.x;
    int y = (int)reg.y;
    int r = (int)reg.r;
    double bg = bg2.x;
    double std = bg2.y;
    
    // Calculate x and y marginals (sum the image into 1d lines in x and y)
    double total = 0;
    double *xm = (double *)malloc(2*r*sizeof(double));
    if (xm == NULL)
        error("malloc failed");
    for (int i = 0; i < 2*r; i++)
        xm[i] = 0;

    double *ym = (double *)malloc(2*r*sizeof(double));
    if (ym == NULL)
        error("malloc failed");
    
    for (int j = 0; j < 2*r; j++)
        ym[j] = 0;
        
    for (int j = 0; j < 2*r; j++)
        for (int i = 0; i < 2*r; i++)
        {            
            // Ignore points outside the circle
            if ((i-r)*(i-r) + (j-r)*(j-r) > r*r)
                continue;

            double px = frame->dbl_data[frame->cols*(y+j-r) + (x+i-r)] - bg;
            if (fabs(px) < 3*std)
                continue;

            xm[i] += px;
            ym[j] += px;
            total += px;
        }

    // Calculate x and y moments
    double xc = 0;
    double yc = 0;
    for (int i = 0; i < 2*r; i++)
        xc += (double)i*xm[i] / total;
    
    for (int j = 0; j < 2*r; j++)
        yc += (double)j*ym[j] / total;
    
    free(xm);
    free(ym);
    
    double2 ret = {xc + x - r,yc + y - r};
    return ret;
}


// Iterate center_aperture to converge on the best aperture position
double2 converge_aperture(target r, framedata *frame)
{
    double2 error = {-1,-1};
    double2 bg, xy;
    target last;

    int n = 0;
    double move = 0;
    double2 pos[20];
    // Iterate until we move less than 1/16 of a pixel or reach 20 iterations
    do
    {
        // Most apertures converge within 3 iterations.
        // If we haven't converged by 20, we are probably in a cycle
        if (n == 20)
        {
            xy.x = xy.y = 0;
            for (int i = 15; i < 20; i++)
            {
                xy.x += pos[i].x;
                xy.y += pos[i].y;
            }
            xy.x /= 5;
            xy.y /= 5;
            printf("Converge failed - using average of last 5 iterations: (%f,%f)\n", xy.x, xy.y);
            break;
        }

        if (r.x - r.r < 0 || r.x + r.r >= frame->cols || r.y - r.r < 0 || r.y + r.r >= frame->rows)
        {
            fprintf(stderr, "Aperture outside chip - skipping\n");
            return error;
        }

        last = r;
        bg = calculate_background(r, frame);
        xy = center_aperture(r, bg, frame);
        pos[n] = xy;

        printf("%d: (%f,%f) -> (%f,%f) [%f,%f]\n", n, r.x, r.y, xy.x, xy.y, bg.x, bg.y);
        r.x = xy.x;
        r.y = xy.y;

        move = (xy.x-last.x)*(xy.x-last.x) + (xy.y-last.y)*(xy.y-last.y);
        n++;
    } while (move >= 0.00390625);

    return xy;
}


// Calculate the mode intensity and standard deviation within an annulus
double2 calculate_background(target r, framedata *frame)
{
    int minx = (int)fmax(floor(r.x - r.s2), 0);
    int maxx = (int)fmin(ceil(r.x + r.s2), frame->cols-1);
    int miny = (int)fmax(floor(r.y - r.s2), 0);
    int maxy = (int)fmin(ceil(r.y + r.s2), frame->rows-1);
    
    // Copy pixels into a flat list that can be sorted
    // Allocate enough space to store the entire target region, but only copy pixels
    // within the annulus.
    double *data = (double *)malloc((maxy - miny + 1)*(maxx - minx + 1)*sizeof(double));
    int n = 0;
    for (int j = miny; j <= maxy; j++)
        for (int i = minx; i <= maxx; i++)
        {
            double d2 = (r.x-i)*(r.x-i) + (r.y-j)*(r.y-j);
            if (d2 > r.s1*r.s1 && d2 < r.s2*r.s2)
                data[n++] = frame->dbl_data[frame->cols*j + i];
        }

    // Sort data into ascending order
    qsort(data, n, sizeof(double), compare_double);

    // Calculate initial mean value
    double mean = 0;
    for (int i = 0; i < n; i++)
        mean += data[i];
    mean /= n;

    // Calculate initial stddev
    double std = 0;
    for (int i = 0; i < n; i++)
        std += (data[i] - mean)*(data[i] - mean);
    std = sqrt(std/n);

    int oldN = n;
    // Discard pixels brighter than mean + 10*stddev
    while (data[n-1] > mean + 10*std)
    {
        printf("discarding bright sky pixel %f\n", data[n-1]);
        n--;
    }

    // Calculate improved mean and median
    if (n != oldN)
    {
        mean = std = 0;
        for (int i = 0; i < n; i++)
            mean += data[i];
        mean /= n;

        for (int i = 0; i < n; i++)
            std += (data[i] - mean)*(data[i] - mean);
        std = sqrt(std/n);
    }
    
    // Calculate median
    double median = (data[n/2] + data[n/2+1])/2;

    //printf("mean: %f median: %f std: %f\n",mean,median,std);
    
    free(data);
    double2 ret = {3*mean - 2*median, std};
    return ret;
}

// Finds the intersection point of the line defined by p1 and p2 (both x,y)
// with the circle c (x,y,r).
//   Assumes that there is only one intersection (one point inside, one outside)
//   Returns (x,y) of the intersection
// See logbook 07/02/11 for calculation workthrough
static double2 line_circle_intersection(double x, double y, double r, double2 p0, double2 p1)
{
    // Line from p1 to p2
    double2 dp = {p1.x - p0.x, p1.y - p0.y};
    
    // Line from c to p1
    double2 dc = {p0.x - x, p0.y - y};
    
    // Polynomial coefficients
    double a = dp.x*dp.x + dp.y*dp.y;
    double b = 2*(dc.x*dp.x + dc.y*dp.y);
    double c = dc.x*dc.x + dc.y*dc.y - r*r;

    // Solve for line parameter x.
    double d = sqrt(b*b - 4*a*c);
    double x1 = (-b + d)/(2*a);
    double x2 = (-b - d)/(2*a);
    
    // The solution we want will be 0<=x<=1
    double sol = (x1 >= 0 && x1 <= 1) ? x1 : x2;
    
    double2 ret = {p0.x + sol*dp.x, p0.y + sol*dp.y};
    return ret;
}

// Calculate the area inside a chord, defined by p1,p2 on the edge of a circle radius r
static double chord_area(double2 p1, double2 p2, double r)
{
    // b is 0.5*the length of the chord defined by p1 and p2
    double b = sqrt((p2.x-p1.x)*(p2.x-p1.x) + (p2.y-p1.y)*(p2.y-p1.y))/2;
    return r*r*asin(b/r) - b*sqrt(r*r-b*b);
}

// Calculate the area of a polygon defined by a list of points
//   Returns the area
static double polygon_area(double2 v[], int nv)
{
    double a = 0;
    int n = 0;
    for (int i = 0; i < nv; i++)
    {
        n = i == 0 ? nv-1 : i-1;
        a += v[n].x*v[i].y - v[i].x*v[n].y;
    }
    return fabs(a/2);
}


// Convenience function to get a corner, with indices that go outside the indexed range
static double2 corners[4] = {
    {0,0},
    {0,1},
    {1,1},
    {1,0}
};

static double2 c(int k)
{
    while (k < 0) k += 4;
    while (k > 3) k -= 4;
    return corners[k];
}

// Calculate the intesection between the unit pixel (with TL corner at the origin)
// and the aperture defined by x,y,r.
//   Returns a number between 0 and 1 specifying the intersecting area
// Origin is defined as the bottom-left of the pixel
static double pixel_aperture_intesection(double x, double y, double r)
{
    // We don't yet handle the extra cases for r < 1
    if (r < 1)
        return 0;

    // Determine which pixels are inside and outside the aperture
    int hit[4];
    int numhit = 0;
    for (int k = 0; k < 4; k++)
        if ((corners[k].x - x)*(corners[k].x - x) + (corners[k].y - y)*(corners[k].y - y) <= r*r)
            hit[numhit++] = k;
    
    switch (numhit)
    {
        // If 0 edges are inside the aperture, but aperture is inside pixel then the aperture is completely inside
        case 0:
        {
            // Check that aperture doesn't intersect with each edge of the pixel in turn
            // If so, the intersecting area is a sector
            // See whiteboard photo from 13/7/11 for working

            // Aperture must be within r distance from the pixel for it to intersect
            if (x < -r || y < -r || x > r + 1 || y > r + 1)
                return 0;

            // Aperture is centered inside the pixel, but intersects no edges.
            // It must be fully contained inside the pixel
            // We don't yet handle this case, so return 0
            if (x >= 0 && y >= 0 && x <= 1 && y <= 1)
                return 0;

            // Aperture cannot intersect the pixel without having one of the
            // vertices inside if it lies in one of these corner regions
            if ((x < 0 && y < 0) || (x < 0 && y > 1) || (x > 1 && y > 1) || (x > 1 && y < 0))
                return 0;

            // Remaining aperture positions may potentially intersect the pixel twice, without
            // containing a vertex. Check each in turn.
            // TODO: Check each edge for intersections

            return 0;
        }
        case 4: return 1;
        case 1:
        {
            // Intersection points
            double2 pts[3] =
            {
                line_circle_intersection(x,y,r, c(hit[0] - 1), c(hit[0])),
                c(hit[0]),
                line_circle_intersection(x,y,r, c(hit[0]), c(hit[0] + 1))
            };

            // Area is triangle + chord
            return polygon_area(pts, 3) + chord_area(pts[0], pts[2], r);    
        }
        break;
        case 2:
        {
            // Find the first inside the aperture
            int first = (hit[1] - hit[0] == 3) ? hit[1] : hit[0];
            
            // Intersection points
            double2 pts[4] =
            {
                line_circle_intersection(x,y,r, c(first-1), c(first)),
                c(first),
                c(first+1),
                line_circle_intersection(x,y,r, c(first+1), c(first+2))
            };
            
            // Area is a quadralateral + chord
            return polygon_area(pts,4) + chord_area(pts[0], pts[3], r);
        }
        break;
        case 3:
        {
            int outside = 3;
            for (int k = 0; k < numhit; k++)
                if (hit[k] != k)
                {
                    outside = k;
                    break;
                }
    
            // Intersection points
            double2 pts[3] =
            {
                line_circle_intersection(x,y,r, c(outside-1), c(outside)),
                c(outside),
                line_circle_intersection(x,y,r, c(outside), c(outside+1))
            };
            
            // Area is square - triangle + chord
            return 1 - polygon_area(pts,3) + chord_area(pts[0], pts[2], r);
        }
        break;
    }
    return 0;
}

static double pixel_value(framedata *frame, int i, int j)
{
    return frame->dbl_data[i + frame->cols*j];
}

// Integrates the flux within the specified aperture, 
// accounting for partially covered pixels.
//   Takes the aperture (x,y,r) and the image data (2d numpy array)
//   Returns the contained flux (including background)
double integrate_aperture(double2 xy, double r, framedata *frame)
{
    double total = 0;
    int bx = floor(xy.x), by = floor(xy.y), br = ceil(r) + 1;
    for (int i = bx-br; i < bx+br; i++)
        for (int j = by-br; j < by+br; j++)
        {
            // Split pixel into 4 subpixels using the QUADPX algorithm
            // TODO: work through the details
            // TODO: test different values of H
            double H = 1; // hardness parameter
            double PN = pixel_value(frame, i, j+1);
            double PS = pixel_value(frame, i, j-1);
            double PE = pixel_value(frame, i+1, j);
            double PW = pixel_value(frame, i-1, j);
            double PC = pixel_value(frame, i, j);
            double S = PN + PS + PE + PW;
            double A = (1 - H)*PC/4.0;
            double B = (H/S)*(PC/2.0);
            total += ((PE+PN)*B + A)*pixel_aperture_intesection(2*(xy.x-i)+1, 2*(xy.y-j)+1, 2*r);
            total += ((PN+PW)*B + A)*pixel_aperture_intesection(2*(xy.x-i), 2*(xy.y-j)+1, 2*r);
            total += ((PW+PS)*B + A)*pixel_aperture_intesection(2*(xy.x-i), 2*(xy.y-j), 2*r);
            total += ((PS+PE)*B + A)*pixel_aperture_intesection(2*(xy.x-i)+1, 2*(xy.y-j), 2*r);
        }

    return total;
}


