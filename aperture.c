/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "aperture.h"
#include "helpers.h"
#include "framedata.h"

extern int verbosity;

// Optimize the centering of a circular aperture (radius r, initial position xy)
// over a star in frame
// Input position is modified in place, or an error code returned 
static int center_aperture_inner(double2 *xy, double rd, double sky_intensity, double sky_std_dev, framedata *frame)
{
    // Square within the frame to sample
    uint16_t x = floor(xy->x);
    uint16_t y = floor(xy->y);
    uint16_t r = ceil(rd) + 1;

    // Calculate x and y marginals (sum the image into 1d lines in x and y)
    double *xm = calloc(2*r, sizeof(double));
    double *ym = calloc(2*r, sizeof(double));
    if (!xm || !ym)
        return error("\tmalloc failed");

    double total = 0;
    for (uint16_t j = 0; j < 2*r; j++)
        for (uint16_t i = 0; i < 2*r; i++)
        {            
            // Ignore pixels outside the circle
            if ((i-r)*(i-r) + (j-r)*(j-r) > r*r)
                continue;

            double px = frame->data[frame->cols*(y+j-r) + (x+i-r)] - sky_intensity;
            if (fabs(px) < 3*sky_std_dev)
                continue;

            xm[i] += px;
            ym[j] += px;
            total += px;
        }

    if (total == 0)
        return error("\tNo signal inside aperture");

    // Calculate x and y moments
    double xc = 0, yc = 0;
    for (size_t i = 0; i < 2*r; i++)
    {
        xc += i*xm[i] / total;
        yc += i*ym[i] / total;
    }

    free(xm);
    free(ym);

    xy->x = xc + x - r;
    xy->y = yc + y - r;
    return 0;
}


// Iterate center_aperture to converge on the best aperture position
int center_aperture(aperture r, framedata *frame, double2 *center)
{
    // Allow up to 20 iterations to center
    double2 pos[20];
    double move;
    uint8_t n = 0;

    // Set initial center from the given aperture
    pos[0].x = r.x;
    pos[0].y = r.y;

    // Iterate until we move less than 1/16 of a pixel or reach 20 iterations
    do
    {
        n++;

        // Most apertures converge within 3 iterations.
        // If we haven't converged by 20, we are probably in a cycle
        if (n == 19)
        {
            pos[n].x = pos[n].y = 0;
            for (uint8_t i = 14; i < 19; i++)
            {
                pos[n].x += pos[i].x;
                pos[n].y += pos[i].y;
            }
            pos[n].x /= 5;
            pos[n].y /= 5;
            if (verbosity >= 1)
                error("\tConverge failed - using average of last 5 iterations: (%f,%f)", pos[n].x, pos[n].y);
            break;
        }

        // Set initial position from the previous iteration
        pos[n] = pos[n-1];

        // Check that the position is valid
        int16_t tx = floor(pos[n].x);
        int16_t ty = floor(pos[n].y);
        int16_t tr = ceil(r.r) + 1;
        if (tx < tr || tx + tr >= frame->cols ||
            ty < tr || ty + tr >= frame->rows)
            return (verbosity >= 1) ? error("\tAperture outside chip - skipping %f %f", pos[n].x, pos[n].y) : 1;

        // Calculate new background
        double sky_intensity, sky_std_dev;
        if (calculate_background(r, frame, &sky_intensity, &sky_std_dev))
            return error("\tcalculate_background failed");

        // Iterate centering
        if (center_aperture_inner(&pos[n], r.r, sky_intensity, sky_std_dev, frame))
            return error("\tcenter_aperture_inner failed");

        // Calculate shift
        move = (pos[n].x - pos[n-1].x)*(pos[n].x - pos[n-1].x) + (pos[n].y - pos[n-1].y)*(pos[n].y - pos[n-1].y);
    } while (move >= 0.00390625);

    *center = pos[n];
    return 0;
}


// Calculate the mode intensity and standard deviation within an annulus
// Pixels brighter than 10 x standard deviation are discarded
// This provides a robust method for calculating the background sky intensity and uncertainty
int calculate_background(aperture r, framedata *frame, double *sky_mode, double *sky_std_dev)
{
    uint16_t minx = (uint16_t)fmax(floor(r.x - r.s2), 0);
    uint16_t maxx = (uint16_t)fmin(ceil(r.x + r.s2), frame->cols-1);
    uint16_t miny = (uint16_t)fmax(floor(r.y - r.s2), 0);
    uint16_t maxy = (uint16_t)fmin(ceil(r.y + r.s2), frame->rows-1);

    // Copy pixels into a flat list that can be sorted
    // Allocate enough space to store the entire aperture region, but only copy pixels
    // within the annulus.
    double *data = calloc((maxy - miny + 1)*(maxx - minx + 1), sizeof(double));
    if (!data)
        return error("malloc failed");

    size_t n = 0;
    for (uint16_t j = miny; j <= maxy; j++)
        for (uint16_t i = minx; i <= maxx; i++)
        {
            double d2 = (r.x - i)*(r.x - i) + (r.y - j)*(r.y - j);
            if (d2 > r.s1*r.s1 && d2 < r.s2*r.s2)
                data[n++] = frame->data[frame->cols*j + i];
        }

    // Sort data into ascending order
    qsort(data, n, sizeof(double), compare_double);

    // Calculate initial mean value
    double mean = 0;
    for (size_t i = 0; i < n; i++)
        mean += data[i];
    mean /= n;

    // Calculate initial stddev
    double std = 0;
    for (size_t i = 0; i < n; i++)
        std += (data[i] - mean)*(data[i] - mean);
    std = sqrt(std/n);

    // Discard pixels brighter than mean + 10*stddev
    size_t filtered_n = n;
    while (data[filtered_n - 1] > mean + 10*std)
        filtered_n--;

    if (verbosity >= 1 && filtered_n != n)
        printf("\tdiscarding %zu bright sky pixels\n", n - filtered_n);

    // Recalculate mean and median ignoring discarded pixels
    if (n != filtered_n)
    {
        mean = std = 0;
        for (size_t i = 0; i < filtered_n; i++)
            mean += data[i];
        mean /= filtered_n;

        for (size_t i = 0; i < filtered_n; i++)
            std += (data[i] - mean)*(data[i] - mean);
        std = sqrt(std/filtered_n);
    }

    // Set return values and clean up
    double median = (data[filtered_n/2] + data[filtered_n/2+1])/2;
    if (sky_mode)
        *sky_mode = 3*mean - 2*median;

    if (sky_std_dev)
        *sky_std_dev = std;

    free(data);
    return 0;
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
static double polygon_area(double2 v[], size_t nv)
{
    double a = 0;
    size_t n = 0;
    for (size_t i = 0; i < nv; i++)
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

static double2 c(int8_t k)
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
    uint8_t hit[4];
    uint8_t numhit = 0;
    for (uint8_t k = 0; k < 4; k++)
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
                line_circle_intersection(x, y, r, c(hit[0] - 1), c(hit[0])),
                c(hit[0]),
                line_circle_intersection(x, y, r, c(hit[0]), c(hit[0] + 1))
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
                line_circle_intersection(x, y, r, c(first - 1), c(first)),
                c(first),
                c(first+1),
                line_circle_intersection(x, y, r, c(first + 1), c(first + 2))
            };
            
            // Area is a quadralateral + chord
            return polygon_area(pts, 4) + chord_area(pts[0], pts[3], r);
        }
        break;
        case 3:
        {
            int outside = 3;
            for (uint8_t k = 0; k < numhit; k++)
                if (hit[k] != k)
                {
                    outside = k;
                    break;
                }
    
            // Intersection points
            double2 pts[3] =
            {
                line_circle_intersection(x, y, r, c(outside - 1), c(outside)),
                c(outside),
                line_circle_intersection(x, y, r, c(outside), c(outside + 1))
            };
            
            // Area is square - triangle + chord
            return 1 - polygon_area(pts, 3) + chord_area(pts[0], pts[2], r);
        }
        break;
    }
    return 0;
}

// Integrates the flux within the specified aperture, 
// accounting for partially covered pixels.
//   Takes the aperture (x,y,r) and the image data
//   Returns the contained flux (including background)
double integrate_aperture(double2 xy, double r, framedata *frame)
{
    double total = 0;
    uint16_t bx = floor(xy.x);
    uint16_t by = floor(xy.y);
    uint16_t br = ceil(r) + 1;

    for (int32_t i = bx - br; i < bx + br; i++)
        for (int32_t j = by - br; j < by + br; j++)
            total += pixel_aperture_intesection(xy.x - i, xy.y - j, r)*frame->data[i + frame->cols*j];

    return total;
}

void integrate_aperture_and_noise(double2 xy, double r, framedata *frame, framedata *dark, double readnoise, double gain, double *signal, double *noise)
{
    *signal = 0;
    *noise = 0;

    uint16_t bx = floor(xy.x);
    uint16_t by = floor(xy.y);
    uint16_t br = ceil(r) + 1;

    for (int32_t i = bx - br; i < bx + br; i++)
        for (int32_t j = by - br; j < by + br; j++)
        {
            double area = pixel_aperture_intesection(xy.x - i, xy.y - j, r);
            double flux = frame->data[i + frame->cols*j];
            double darkflux = dark ? dark->data[i + frame->cols*j] : 0;

            *signal += area*flux;
            *noise += area*(readnoise*readnoise + (flux + darkflux)/gain);
        }

    // Convert variance to standard deviation
    *noise = sqrt(*noise);
}

double estimate_fwhm(framedata *frame, double2 xy, double bg, double max_radius)
{
    double center_profile = frame->data[frame->cols*((uint16_t)xy.y) + (uint16_t)xy.x] - bg;
    double last_intensity = 0;
    double last_profile = center_profile;

    double fwhm = 0;
    for (uint16_t radius = 1; radius <= (uint16_t)max_radius; radius++)
    {
        double intensity = integrate_aperture(xy, radius, frame) - bg*M_PI*radius*radius;
        double profile = (intensity - last_intensity) / (M_PI*(2*radius - 1));

        if (profile < center_profile/2)
        {
            double last_radius = radius - 1;
            fwhm = 2*(last_radius + (radius - last_radius)*(center_profile/2 - last_profile)/(profile - last_profile));
            break;
        }

        last_intensity = intensity;
        last_profile = profile;
    }

    return fwhm;
}
