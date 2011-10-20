/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <stdlib.h>
#include <math.h>

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

    // Loop over diagonal index / column
    for (int i = 0; i < m; i++)
    {
        // Pick the largest value in the remaining (leftmost square) rows as the next pivot
        int prow, pcol;
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

        // Store the pivot column for later
        prows[pcol] = prow;

        // Check that we don't have a singular matrix
        if (fabs(A[prow*n + pcol]) < eps)
        {
            free(prows);
            return 1;
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
        if (prows[i] == i)
            continue;

        int k = prows[i];
        for (int j = 0; j < n; j++)
        {
            double t = A[i*n + j];
            A[i*n + j] = A[k*n + j];
            A[k*n + j] = t;
        }

        prows[i] = prows[k];
        prows[k] = k;
    }

    free(prows);
    return 0;
}