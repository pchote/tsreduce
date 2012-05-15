/* Copyright 2010-2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <xpa.h>
#include <unistd.h>
#include <dirent.h>
#include <regex.h>
#include <math.h>
#include "helpers.h"

// Helper function to free a 2d char array allocated using malloc etc.
void free_2d_array(char **array, int len)
{
    for (int i = 0; i < len; i++)
        free(array[i]);
    free(array);
}

/*
 * Cast an array of doubles to an array of floats reusing the same memory
 */
float *cast_double_array_to_float(double *d_ptr, size_t count)
{
    float *f_ptr = (float *)d_ptr;
    for (size_t i = 0; i < count; i++)
        f_ptr[i] = d_ptr[i];
    return f_ptr;
}

/*
 * versionsort isn't defined on all platforms (e.g osx), so define our own copy
 * Taken from GNU strverscmp.c, licenced under GPLv2 or later
 */
#define ISDIGIT(c) ((unsigned int) (c) - '0' <= 9)
#define S_N    0x0
#define S_I    0x4
#define S_F    0x8
#define S_Z    0xC
#define CMP    2
#define LEN    3
int tsreduce_versionsort(const void *a, const void *b)
{
    const unsigned char *p1 = (const unsigned char *) (*(const struct dirent **) a)->d_name;
    const unsigned char *p2 = (const unsigned char *) (*(const struct dirent **) b)->d_name;
    unsigned char c1, c2;
    int state;
    int diff;

    /* Symbol(s)    0       [1-9]   others  (padding)
     Transition   (10) 0  (01) d  (00) x  (11) -   */
    static const unsigned int next_state[] =
    {
        /* state    x    d    0    - */
        /* S_N */  S_N, S_I, S_Z, S_N,
        /* S_I */  S_N, S_I, S_I, S_I,
        /* S_F */  S_N, S_F, S_F, S_F,
        /* S_Z */  S_N, S_F, S_Z, S_Z
    };

    static const int result_type[] =
    {
        /* state   x/x  x/d  x/0  x/-  d/x  d/d  d/0  d/-
         0/x  0/d  0/0  0/-  -/x  -/d  -/0  -/- */

        /* S_N */  CMP, CMP, CMP, CMP, CMP, LEN, CMP, CMP,
        CMP, CMP, CMP, CMP, CMP, CMP, CMP, CMP,
        /* S_I */  CMP, -1,  -1,  CMP,  1,  LEN, LEN, CMP,
        1,  LEN, LEN, CMP, CMP, CMP, CMP, CMP,
        /* S_F */  CMP, CMP, CMP, CMP, CMP, LEN, CMP, CMP,
        CMP, CMP, CMP, CMP, CMP, CMP, CMP, CMP,
        /* S_Z */  CMP,  1,   1,  CMP, -1,  CMP, CMP, CMP,
        -1,  CMP, CMP, CMP
    };

    if (p1 == p2)
        return 0;

    c1 = *p1++;
    c2 = *p2++;
    /* Hint: '0' is a digit too.  */
    state = S_N | ((c1 == '0') + (ISDIGIT (c1) != 0));

    while ((diff = c1 - c2) == 0 && c1 != '\0')
    {
        state = next_state[state];
        c1 = *p1++;
        c2 = *p2++;
        state |= (c1 == '0') + (ISDIGIT (c1) != 0);
    }

    state = result_type[state << 2 | ((c2 == '0') + (ISDIGIT (c2) != 0))];

    switch (state)
    {
        case CMP:
            return diff;

        case LEN:
            while (ISDIGIT (*p1++))
                if (!ISDIGIT (*p2++))
                    return 1;

            return ISDIGIT (*p2) ? -1 : diff;

        default:
            return state;
    }
}

// Find files that match the given regex.
// outList is set to point to an array (which must be later freed by the caller) of filenames
// Returns the number of files matched or negative on error
int get_matching_files(const char *pattern, char ***outList)
{
    // Compile the pattern into a regex
    regex_t regex;
    int regerr = 0;
    if ((regerr = regcomp(&regex, pattern, REG_EXTENDED | REG_NOSUB)))
    {
        char errbuf[1024];
        regerror(regerr, &regex, errbuf, 1024);
        return -error("Error compiling `%s` into a regular expression: %s", pattern, errbuf);
    }

    // Find the matching files
    struct dirent **files;
    int numFiles = scandir(".", &files, 0, tsreduce_versionsort);
    
    char **ret = (char **)malloc(numFiles*sizeof(char*));
    if (ret == NULL)
    {
        regfree(&regex);
        free(files);
        return -error("Memory allocation for %d filenames failed.", numFiles);
    }

    int numMatched = 0;
    for (int i = 0; i < numFiles; i++)
    {
        if (!regexec(&regex, files[i]->d_name, 0, NULL, 0))
        {
            ret[numMatched] = strdup(files[i]->d_name);
            if (ret[numMatched] == NULL)
            {
                regfree(&regex);
                for (int j = i; j < numFiles; j++)
                    free(files[j]);
                free(files);
                return -error("Memory allocation failed.");
            }
            numMatched++;
        }
        free(files[i]);
    }
    free(files);
    regfree(&regex);

    *outList = ret;
    return numMatched;
}

// Find the filename of the first file that (alphabetically) matches the given regex pattern
// Result is copied into filenamebuf. Returns 1 on success, 0 on failure.
int get_first_matching_file(char *pattern, char *filenamebuf, int buflen)
{
    // Compile the pattern into a regex
    regex_t regex;
    int regerr = 0;
    if ((regerr = regcomp(&regex, pattern, REG_EXTENDED | REG_NOSUB)))
    {
        char errbuf[1024];
        regerror(regerr, &regex, errbuf, 1024);
        return error("Error compiling `%s` into a regular expression: %s", pattern, errbuf);
    }

    // Find the first matching file
    struct dirent **matched;
    int numMatched = scandir(".", &matched, 0, alphasort);
    int found = 0;
    for (int i = 0; i < numMatched; i++)
    {
        if (!found && !regexec(&regex, matched[i]->d_name, 0, NULL, 0))
        {
            found = 1;
            strncpy(filenamebuf, matched[i]->d_name, buflen);
            filenamebuf[NAME_MAX-1] = '\0';
        }
        free(matched[i]);
    }
    free(matched);
    regfree(&regex);
    return found;
}

// qsort() doubles in ascending order
int compare_double(const void *a, const void *b)
{
    const double *da = (const double *)a;
    const double *db = (const double *)b;

    return (*da > *db) - (*da < *db);
}

// Prints an vararg error to stderr then returns 1
int error(const char * format, ...)
{
	va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	fprintf(stderr, "\n");
    return 1;
}

// Prints a vararg error and exits
void die(const char * format, ...)
{
    va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	fprintf(stderr, "\n");
    exit(1);
}

static int ds9_available(char *title)
{
    char *names[1];
    char *errs[1];
    int valid = XPAAccess(NULL, title, NULL, NULL, names, errs, 1);
    if (errs[0] != NULL)
    {
        valid = 0;
        free(errs[0]);
    }
    if (names[0]) free(names[0]);

    return valid;
}

int init_ds9(char *title)
{
    if (!ds9_available(title))
    {
        char buf[128];
        snprintf(buf, 128, "ds9 -title %s&", title);
        system(buf);

        int wait = 0;
        while (!ds9_available(title))
        {
            if (wait++ == 10)
                return 0;

            printf("Waiting for ds9... %d\n", wait);
            sleep(1);
        }
    }
    return 1;
}

// Send a command to ds9 (with optional data) and ignore any response
int tell_ds9(char *title, char *command, void *data, int dataSize)
{
    char *names[1];
    char *errs[1];
    int valid = XPASet(NULL, title, command, NULL, data, dataSize, names, errs, 1);
    if (errs[0] != NULL)
    {
        valid = 0;
        printf("Error: %s", errs[0]);
        free(errs[0]);
    }
    if (names[0]) free(names[0]);
    return !valid;
}

// Send a command to ds9 and return its response.
// outbuf is a pointer to a string, which must be later freed by the caller.
int ask_ds9(char *title, char *command, char **outbuf)
{
    char *names[1];
    char *errs[1];
    char *ret[1];
    int len[1];

    int valid = XPAGet(NULL, title, command, NULL, ret, len, names, errs, 1);

    if( errs[0] != NULL )
    {
        valid = 0;
        printf("Error: %s", errs[0]);
        free(errs[0]);
    }
    else
        *outbuf = ret[0]; // The caller can free this later

    if (names[0]) free(names[0]);
    return !valid;
}

// Calculate the amplitude spectrum of the signal defined by numData points in (time, data)
// in the frequency range (fmin, fmax)
// numOut results are stored in (outFreq, outPower)
void calculate_amplitude_spectrum(double fmin, double fmax, double *t, double *data, int numData, double *outFreq, double *outAmpl, int numOut)
{
    double df = (fmax-fmin)/numOut;
    for (int j = 0; j < numOut; j++)
    {
        double real = 0;
        double imag = 0;
        outFreq[j] = fmin + j*df;

        for (int i = 0; i < numData; i++)
        {
            double phase = -outFreq[j]*2*M_PI*(t[i]-t[0]);
            real += data[i]*cos(phase)/numData;
            imag += data[i]*sin(phase)/numData;
        }

        outAmpl[j] = 2*sqrt(real*real + imag*imag);
    }
}

void calculate_amplitude_spectrum_float(float fmin, float fmax, float *t, float *data, int numData, float *outFreq, float *outAmpl, int numOut)
{
    double df = (fmax-fmin)/numOut;
    for (int j = 0; j < numOut; j++)
    {
        double real = 0;
        double imag = 0;
        outFreq[j] = fmin + j*df;

        for (int i = 0; i < numData; i++)
        {
            double phase = -outFreq[j]*2*M_PI*(t[i]-t[0]);
            real += data[i]*cos(phase)/numData;
            imag += data[i]*sin(phase)/numData;
        }

        outAmpl[j] = 2*sqrt(real*real + imag*imag);
    }
}

