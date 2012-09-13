/* Copyright 2010-2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <dirent.h>
#include <regex.h>
#include <math.h>
#include "helpers.h"

#if (defined _WIN32 || defined _WIN64)

#include <windows.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/stat.h>

// scandir() for win32
// this tool should make life easier for people writing for both unix and wintel
// written by Tom Torfs, 2002/10/31
// donated to the public domain; use this code for anything you like as long as
// it is understood there are absolutely *NO* warranties of any kind, even implied
int scandir(const char *dirname,
            struct dirent ***namelist,
            int (*select)(const struct dirent *),
            int (*compar)(const struct dirent **,
                          const struct dirent **))
{
	WIN32_FIND_DATA wfd;
	HANDLE hf;
	struct dirent **plist, **newlist;
	struct dirent d;
	int numentries = 0;
	int allocentries = 255;
	int i;
	char path[FILENAME_MAX];

	i = strlen(dirname);

	if (i > sizeof path - 5)
		return -1;

	strcpy(path, dirname);
	if (i>0 && dirname[i-1]!='\\' && dirname[i-1]!='/')
		strcat(path, "\\");
	strcat(path, "*.*");

	hf = FindFirstFile(path, &wfd);
	if (hf == INVALID_HANDLE_VALUE)
		return -1;

	plist = malloc(sizeof *plist * allocentries);
	if (plist==NULL)
	{
		FindClose(hf);
		return -1;
	}

	do
	{
		if (numentries==allocentries)
		{
			allocentries *= 2;
			newlist = realloc(plist, sizeof *plist * allocentries);
			if (newlist==NULL)
			{
				for (i=0; i<numentries; i++)
                    free(plist[i]);
				free(plist);
				FindClose(hf);
				return -1;
			}
			plist = newlist;
		}

		strncpy(d.d_name, wfd.cFileName, sizeof d.d_name);
		d.d_ino = 0;
		d.d_namlen = strlen(wfd.cFileName);
		d.d_reclen = sizeof d;

		if (select==NULL || select(&d))
		{
			plist[numentries] = malloc(sizeof d);
			if (plist[numentries]==NULL)
			{
				for (i=0; i<numentries; i++)
                    free(plist[i]);
				free(plist);
				FindClose(hf);
				return -1;
			};
			memcpy(plist[numentries], &d, sizeof d);
			numentries++;
		}
	}
	while (FindNextFile(hf, &wfd));

	FindClose(hf);

	if (numentries==0)
	{
		free(plist);
		*namelist = NULL;
	}
	else
	{
		newlist = realloc(plist, sizeof *plist * numentries);
		if (newlist!=NULL)
            plist = newlist;

		if (compar!=NULL)
			qsort(plist, numentries, sizeof *plist, compar);

		*namelist = plist;
	}
	return numentries;
}

char *strndup(const char *s, size_t max)
{
#ifdef _WIN64
    size_t len = strnlen(s, max);
#else
    size_t len = strlen(s);
    if (len > max)
        len = max;
#endif
    char *ret = malloc(len + 1);
    if (ret == NULL)
        return NULL;

    memcpy(ret, s, len);
    ret[len] = '\0';
    return ret;
}

/*
 * realpath() Win32 implementation, supports non standard glibc extension
 * This file has no copyright assigned and is placed in the Public Domain.
 * Written by Nach M. S. September 8, 2005
 */

char *realpath(const char *path, char resolved_path[PATH_MAX])
{
    char *return_path = 0;

    if (path) //Else EINVAL
    {
        if (resolved_path)
        {
            return_path = resolved_path;
        }
        else
        {
            //Non standard extension that glibc uses
            return_path = malloc(PATH_MAX);
        }

        if (return_path) //Else EINVAL
        {
            //This is a Win32 API function similar to what realpath() is supposed to do
            size_t size = GetFullPathNameA(path, PATH_MAX, return_path, 0);

            //GetFullPathNameA() returns a size larger than buffer if buffer is too small
            if (size > PATH_MAX)
            {
                if (return_path != resolved_path) //Malloc'd buffer - Unstandard extension retry
                {
                    size_t new_size;

                    free(return_path);
                    return_path = malloc(size);

                    if (return_path)
                    {
                        new_size = GetFullPathNameA(path, size, return_path, 0); //Try again

                        if (new_size > size) //If it's still too large, we have a problem, don't try again
                        {
                            free(return_path);
                            return_path = 0;
                            errno = ENAMETOOLONG;
                        }
                        else
                        {
                            size = new_size;
                        }
                    }
                    else
                    {
                        //I wasn't sure what to return here, but the standard does say to return EINVAL
                        //if resolved_path is null, and in this case we couldn't malloc large enough buffer
                        errno = EINVAL;
                    }
                }
                else //resolved_path buffer isn't big enough
                {
                    return_path = 0;
                    errno = ENAMETOOLONG;
                }
            }

            //GetFullPathNameA() returns 0 if some path resolve problem occured
            if (!size)
            {
                if (return_path != resolved_path) //Malloc'd buffer
                {
                    free(return_path);
                }

                return_path = 0;

                //Convert MS errors into standard errors
                switch (GetLastError())
                {
                    case ERROR_FILE_NOT_FOUND:
                        errno = ENOENT;
                        break;

                    case ERROR_PATH_NOT_FOUND: case ERROR_INVALID_DRIVE:
                        errno = ENOTDIR;
                        break;

                    case ERROR_ACCESS_DENIED:
                        errno = EACCES;
                        break;

                    default: //Unknown Error
                        errno = EIO;
                        break;
                }
            }

            //If we get to here with a valid return_path, we're still doing good
            if (return_path)
            {
                struct stat stat_buffer;

                //Make sure path exists, stat() returns 0 on success
                if (stat(return_path, &stat_buffer))
                {
                    if (return_path != resolved_path)
                    {
                        free(return_path);
                    }

                    return_path = 0;
                    //stat() will set the correct errno for us
                }
                //else we succeeded!
            }
        }
        else
        {
            errno = EINVAL;
        }
    }
    else
    {
        errno = EINVAL;
    }

    return return_path;
}

#else
#include <xpa.h>
#endif

// Helper function to free a 2d char array allocated using malloc etc.
void free_2d_array(char **array, int len)
{
    for (int i = 0; i < len; i++)
        free(array[i]);
    free(array);
}


// Cast an array of doubles to an array of floats reusing the same memory
float *cast_double_array_to_float(double *d_ptr, size_t count)
{
    float *f_ptr = (float *)d_ptr;
    for (size_t i = 0; i < count; i++)
        f_ptr[i] = d_ptr[i];
    return f_ptr;
}


// versionsort isn't defined on all platforms (e.g osx), so define our own copy
// Taken from GNU strverscmp.c, licenced under GPLv2 or later
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

    //  Symbol(s)    0       [1-9]   others  (padding)
    // Transition   (10) 0  (01) d  (00) x  (11) -
    static const unsigned int next_state[] =
    {
        // state    x    d    0    - 
        /* S_N */  S_N, S_I, S_Z, S_N,
        /* S_I */  S_N, S_I, S_I, S_I,
        /* S_F */  S_N, S_F, S_F, S_F,
        /* S_Z */  S_N, S_F, S_Z, S_Z
    };

    static const int result_type[] =
    {
        //  state   x/x  x/d  x/0  x/-  d/x  d/d  d/0  d/-
        // 0/x  0/d  0/0  0/-  -/x  -/d  -/0  -/- 

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
    // Hint: '0' is a digit too.
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
    // TODO: Use native filehandling functionality under windows
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

// Find the filename of the first file that matches the given regex pattern
// Returns an allocated string, or NULL on failure
char *get_first_matching_file(char *pattern)
{
    char **frame_paths;
    int num_frames = get_matching_files(pattern, &frame_paths);
    if (num_frames <= 0)
        return NULL;

    char *match = strdup(frame_paths[0]);
    free_2d_array(frame_paths, num_frames);
    return match;
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
#if (defined _WIN32 || defined _WIN64)
    return 0;
#else
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
#endif
}

int init_ds9(char *title)
{
#if (defined _WIN32 || defined _WIN64)
    return 0;
#else
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
#endif
}

// Send a command to ds9 (with optional data) and ignore any response
int tell_ds9(char *title, char *command, void *data, int dataSize)
{
#if (defined _WIN32 || defined _WIN64)
    return error("ds9 communication not implemented");
#else
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
#endif
}

// Send a command to ds9 and return its response.
// outbuf is a pointer to a string, which must be later freed by the caller.
int ask_ds9(char *title, char *command, char **outbuf)
{
#if (defined _WIN32 || defined _WIN64)
    return error("ds9 communication not implemented");
#else
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
#endif
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

