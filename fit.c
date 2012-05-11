/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <stdlib.h>
#include <math.h>
#include "helpers.h"

// Calculates the reduced row echelon form using Gauss-Jordan elimination of a m x n matrix A.
// Returns non-zero on error
static int rref(double *A, int m, int n)
{
    double eps = 1e-10; // Limit for singular values

    // Don't reduce into a diagonal form directly; pick the column with the largest
    // value as the pivot for each row and keep track of the row for each column
    // so we can quickly reorder at the end of the process.
    // This reduces numerical problems while avoiding the need to swap columns.
    int *prows = (int *)malloc(m*sizeof(int));
    if (prows == NULL)
        return error("malloc failed");

    // Loop over diagonal index / column
    for (int i = 0; i < m; i++)
    {
        // Pick the largest value in the remaining (leftmost square) rows as the next pivot
        int prow = i, pcol = i;
        double pval = -1;
        for (int j = i; j < m; j++)
            for (int k = 0; k < m; k++)
                if (fabs(A[j*n + k]) >= pval)
                {
                    pval = fabs(A[j*n + k]);
                    prow = j;
                    pcol = k;
                }

        // Swap pivot into the current row
        if (prow != i)
        {
            for (int j = 0; j < n; j++)
            {
                double t = A[i*n + j];
                A[i*n + j] = A[prow*n + j];
                A[prow*n + j] = t;
            }
            prow = i;
        }

        // Store the pivot row for later
        prows[pcol] = prow;

        // Check that we don't have a singular matrix
        if (fabs(A[prow*n + pcol]) < eps)
        {
            free(prows);
            return error("Matrix is singular");
        }

        // Normalize row by pivot
        for (int j = 0; j < n; j++)
        {
            if (j == pcol) continue;
            A[prow*n + j] /= A[prow*n + pcol];
        }
        A[prow*n + pcol] = 1.0;

        // Reduce column
        for (int l = 0; l < m; l++)
        {
            if (l == prow)
                continue;

            for (int j = 0; j < n; j++)
            {
                if (j == pcol) continue;
                A[l*n + j] -= A[l*n + pcol]*A[prow*n + j];
            }
            A[l*n + pcol] = 0;
        }
    }

    // Swap rows to make A the identity
    for (int i = 0; i < m; i++)
    {
        // Don't swap with self or rows that have been completed
        if (prows[i] <= i)
            continue;

        int k = prows[i];
        for (int j = 0; j < n; j++)
        {
            double t = A[i*n + j];
            A[i*n + j] = A[k*n + j];
            A[k*n + j] = t;
        }

        // Update remaining rows with new indices
        prows[i] = i;
        for (int j = i+1; j < n; j++)
            if (prows[j] == i)
            {
                prows[j] = k;
                break;
            }
    }

    free(prows);
    return 0;
}

/*
 * Calculate a linear least-squares fit of numParams basis functions defined by evaluate_basis to (x,y)
 * e (if non-NULL) specifies the estimated standard deviation in y, and is used to weight the fit
 */
static int fit(double *x, double *y, double *e, int c, double *params, int numParams,
               void (*evaluate_basis)(double x, double *basis, int n, void *user), void *user)
{
    int n = numParams;
    int m = numParams + 1;
    double *A = (double *)malloc(m*n*sizeof(double));
    if (A == NULL)
        return error("malloc failed");

    // Initialize to zero
    for (int i = 0; i < m*n; i++)
        A[i] = 0;

    double *basis = (double *)malloc(n*sizeof(double));
    for (int i = 0; i < c; i++)
    {
        // Evaluate basis functions at x[i]
        evaluate_basis(x[i], basis, n, user);

        // Estimated variance in y
        double var = e ? e[i]*e[i] : 1.0;

        for (int j = 0; j < n; j++)
        {
            // Calculate Ajk contribution from i
            for (int k = 0; k < n; k++)
                A[j*m + k] += basis[j]*basis[k]/var;

            // Calculate b_j contribution from i
            A[j*m + n] += y[i]*basis[j]/var;
        }
    }

    // Solve for the coefficients
    if (rref(A, n, m))
    {
        free(A);
        return error("fit failed");
    }

    // Copy coeffs to output
    for (int i = 0; i < n; i++)
        params[i] = A[i*m + n];

    free(A);
    return 0;
}

static void polynomial_fit(double x, double *basis, int n, void *user)
{
    basis[0] = 1;
    for (int j = 1; j < n; j++)
        basis[j] = basis[j-1]*x;
}

static void sinusoidal_fit(double t, double *basis, int n, void *freqs)
{
    double *f = (double *)freqs;
    for (int j = 0; j < n/2; j++)
    {
        double phase = 2*M_PI*f[j]*t;
        basis[2*j] = cos(phase);
        basis[2*j+1] = sin(phase);
    }
}

// Calculate a polynomial fit to the given x,y data
int fit_polynomial(double *x, double *y, double *e, int c, double *coeffs, int degree)
{
    return fit(x, y, e, c, coeffs, degree + 1, polynomial_fit, NULL);
}

/*
 * Takes a list of frequencies 0..N-1
 * Returns a list of amplitudes 0..2*N-1; alternating between cos and sin for each freq in freqs
 */
int fit_sinusoids(double *x, double *y, double *e, int c, double *freqs, int numFreqs, double *amplitudes)
{
    return fit(x, y, e, c, amplitudes, 2*numFreqs, sinusoidal_fit, freqs);
}