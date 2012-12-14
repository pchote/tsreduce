/* Copyright 2010-2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE. */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <dirent.h>
#include <regex.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#ifdef USE_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif

#include "helpers.h"
#include "astro_convert.h"

#if (defined _WIN32)
#include <windows.h>
#include <errno.h>
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

#endif

// Cross platform equivalent of realpath()
char *canonicalize_path(const char *path)
{
#if (defined _WIN32)
    char path_buf[MAX_PATH], *ptr;
    GetFullPathName(path, MAX_PATH, path_buf, &ptr);

    // Replace all '\' in path cd ../tsreducewith '/'
    char *i;
    while ((i = strstr(path_buf, "\\")))
        i[0] = '/';
#else
    char path_buf[PATH_MAX];

    // Expand tilde to home dir
    if (path[0] == '~')
    {
        char *home = getenv("HOME");
        char temp_buf[PATH_MAX];
        snprintf(temp_buf, PATH_MAX, "%s%s", home, &path[1]);
        realpath(temp_buf, path_buf);
    }
    else
        realpath(path, path_buf);
#endif
    return strdup(path_buf);
}

// Cross platform equivalent of gmtime_r()
static void ts_gmtime(ts_time in, struct tm *out)
{
#ifdef _WIN64
    _gmtime_s(out, &in.time);
#elif defined _WIN32
    *out = *localtime(&in.time);
#else
    gmtime_r(&in.time, out);
#endif
}

// Cross platform equivalent of timegm()
static ts_time ts_timegm(struct tm *t, uint16_t ms)
{
    ts_time tt;
#ifdef _WIN64
    tt.time = _mkgmtime(t);
#elif defined _WIN32
    tt.time = mktime(t);
#else
    tt.time = timegm(t);
#endif
    tt.ms = ms;
    return tt;
}

ts_time parse_time(const char *string)
{
    struct tm t;
    unsigned int ms = 0;
    sscanf(string, "%d-%d-%d %d:%d:%d.%u", &t.tm_year, &t.tm_mon, &t.tm_mday, &t.tm_hour, &t.tm_min, &t.tm_sec, &ms);
    t.tm_year -= 1900;
    t.tm_mon -= 1;

    return ts_timegm(&t, ms);
}

ts_time parse_date_time(const char *date, const char *time)
{
    size_t datetime_len = strlen(date) + strlen(time) + 2;

    char *datetime = malloc(datetime_len*sizeof(char));
    snprintf(datetime, datetime_len, "%s %s", date, time);
    ts_time ret = parse_time(datetime);
    free(datetime);

    return ret;
}

void serialize_time(ts_time t, char buf[24])
{
    struct tm time;
    ts_gmtime(t, &time);
    strftime(buf, 24, "%Y-%m-%d %H:%M:%S", &time);
    if (t.ms)
        snprintf(buf+19, 5,".%03d", t.ms);
}

double ts_time_to_utc_hour(ts_time t)
{
    struct tm st;
    ts_gmtime(t, &st);
    return st.tm_hour + st.tm_min / 60.0 + st.tm_sec / 3600.0 + t.ms / 3600000.0;
}

double ts_difftime(ts_time a, ts_time b)
{
    double ret = difftime(a.time, b.time);
    ret += (a.ms - b.ms)/1000.0;
    return ret;
}

double ts_time_to_bjd(ts_time t, double ra, double dec, double epoch)
{
    struct tm tt;
    ts_gmtime(t, &tt);

    // Convert from UT to TT
    tt.tm_sec += utcttoffset(t.time);
    double jd = tmtojd(&tt) + t.ms/86400000.0;

    // Calculate BJD
    double2 reference_coords = precess((double2){ra, dec}, epoch, tmtoyear(&tt));
    return jdtobjd(jd, reference_coords);
}

// Helper function to free a 2d char array allocated using malloc etc.
void free_2d_array(char **array, size_t len)
{
    for (size_t i = 0; i < len; i++)
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
#ifndef ISDIGIT
#define ISDIGIT(c) ((unsigned int) (c) - '0' <= 9)
#endif

#define S_N    0x0
#define S_I    0x4
#define S_F    0x8
#define S_Z    0xC
#define CMP    2
#define LEN    3
int ts_versionsort(const struct dirent **a, const struct dirent **b)
{
    const char *p1 = (*a)->d_name;
    const char *p2 = (*b)->d_name;
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
size_t get_matching_files(const char *pattern, char ***outList)
{
    // Compile the pattern into a regex
    regex_t regex;
    int regerr = 0;
    if ((regerr = regcomp(&regex, pattern, REG_EXTENDED | REG_NOSUB)))
    {
        char errbuf[1024];
        regerror(regerr, &regex, errbuf, 1024);
        fprintf(stderr, "Error compiling `%s` into a regular expression: %s", pattern, errbuf);
        return 0;
    }

    // Find the matching files
    // TODO: Use native filehandling functionality under windows
    struct dirent **files;
    int numFiles = scandir(".", &files, 0, ts_versionsort);
    
    if (numFiles == 0)
    {
        *outList = NULL;
        return 0;
    }

    char **ret = (char **)malloc(numFiles*sizeof(char*));
    if (ret == NULL)
    {
        regfree(&regex);
        free(files);
        fprintf(stderr, "Memory allocation for %d filenames failed.", numFiles);
        return 0;
    }

    size_t numMatched = 0;
    for (size_t i = 0; i < numFiles; i++)
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
                fprintf(stderr, "Memory allocation failed.");
                return 0;
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
    size_t num_frames = get_matching_files(pattern, &frame_paths);
    if (num_frames == 0)
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

int compare_float(const void *a, const void *b)
{
    const float *da = (const float *)a;
    const float *db = (const float *)b;

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

// Sleep for ms milliseconds
void millisleep(int ms)
{
#if (defined _WIN32 || defined _WIN64)
    Sleep(ms);
#else
    nanosleep(&(struct timespec){ms / 1000, (ms % 1000)*1e6}, NULL);
#endif
}

// Run an external process and return the error code.
// stdout is a pointer to a string which is allocated to hold the output
int ts_exec_read(const char *cmd, char **output_ptr)
{
    *output_ptr = NULL;
    FILE *process = popen(cmd, "r");
    if (!process)
        return error("Error invoking read process: %s", cmd);

    size_t output_size = 512;
    size_t output_step = 512;
    char *output = calloc(output_size, sizeof(char));
    if (!output)
        return error("Allocation error");

    char buffer[128];
    while (!feof(process))
    {
        if (fgets(buffer, 128, process) != NULL)
        {
            size_t current_length = strlen(output);
            size_t append_length = strlen(buffer);
            size_t new_length = current_length + append_length;

            // Increase buffer size
            if (new_length >= output_size)
            {
                while (output_size < new_length)
                    output_size += output_step;
                output = realloc(output, output_size*sizeof(char));
                if (!output)
                    return error("Allocation error");

                for (size_t i = current_length; i < output_size; i++)
                    output[i] = '\0';
            }
            strcat(output, buffer);
        }
    }

    *output_ptr = output;
    return pclose(process);
}

// Runs an external process and pipes the given data to stdin
int ts_exec_write(const char *cmd, const void *restrict data, size_t size)
{
#if (defined _WIN32)
    const char *mode = "wb";
#else
    const char *mode = "w";
#endif
    FILE *process = popen(cmd, mode);
    if (!process)
        return error("Error invoking write process: %s", cmd);

    if (size > 0)
        fwrite(data, size, 1, process);

    return pclose(process);
}

int init_ds9()
{
    // xpa exit code = 1 if ds9 is available
#if (defined _WIN32)
    const char *test_command = "xpaaccess tsreduce > nul";
#else
    const char *test_command = "xpaaccess tsreduce > /dev/null";
#endif
    bool available = ts_exec_write(test_command, NULL, 0);

    if (!available)
    {
#if (defined _WIN32)
        const char *ds9_command = "start /b ds9 -title tsreduce";
#else
        const char *ds9_command = "ds9 -title tsreduce&";
#endif
        printf("Starting ds9...\n");
        ts_exec_write(ds9_command, NULL, 0);

        do
        {
            printf("Waiting...\n");
            millisleep(1000);
            available = ts_exec_write(test_command, NULL, 0);
        } while (!available);
    }
    return 0;
}

char *prompt_user_input(char *message, char *fallback)
{
    char *prompt = NULL;
    if (fallback)
    {
        prompt = malloc(strlen(message) + strlen(fallback) + 6);
        sprintf(prompt, "%s [%s]: ", message, fallback);
    }
    else
    {
        prompt = malloc(strlen(message) + 3);
        sprintf(prompt, "%s: ", message);
    }

#ifdef USE_READLINE
    // Prevent trailing space being added after completion
	rl_completion_suppress_append = 1;
    rl_bind_key('\t', rl_complete);

    char *input = readline(prompt);

    // Encountered EOF
    if (!input)
    {
        printf("\n");
        return strdup(strdup(fallback ? fallback : ""));
    }
#else
    char inputbuf[1024];
    printf("%s", prompt);
    fgets(inputbuf, 1024, stdin);

    // Strip trailing newline
    size_t len = strlen(inputbuf);
    if (len >= 1)
        inputbuf[len - 1] = '\0';

    char *input = strdup(inputbuf);
#endif

    // Empty string: Return fallback
    if (strlen(input) == 0)
    {
        free(input);
        return strdup(strdup(fallback ? fallback : ""));
    }

#ifdef USE_READLINE
    add_history(input);
#endif
    return input;
}

// Calculate the amplitude spectrum between (freq_min, freq_max)
// for the signal time/mma (count points)
// length results are stored in (freq, ampl)
void calculate_amplitude_spectrum(double *time, double *mma, size_t count,
                                  double freq_min, double freq_max,
                                  double *freq, double *ampl, size_t length)
{
    double df = (freq_max - freq_min)/length;
    for (size_t j = 0; j < length; j++)
    {
        double real = 0;
        double imag = 0;
        freq[j] = freq_min + j*df;

        for (size_t i = 0; i < count; i++)
        {
            double phase = -freq[j]*2*M_PI*(time[i] - time[0]);
            real += mma[i]*cos(phase)/count;
            imag += mma[i]*sin(phase)/count;
        }

        ampl[j] = 2*sqrt(real*real + imag*imag);
    }
}

bool region_contains(uint16_t r[4], size_t x, size_t y)
{
    return x >= r[0] && x < r[1] && y >= r[2] && y < r[3];
}

// Convenience function for calculating the mean signal in a sub-region of a frame
// Assumes that the frame type is double, and that the region is inside the frame
double region_mean(uint16_t r[4], double *data, size_t stride)
{
    int num_px = (r[1] - r[0])*(r[3] - r[2]);
    double mean = 0;
    for (uint16_t j = r[2]; j < r[3]; j++)
        for (uint16_t i = r[0]; i < r[1]; i++)
            mean += data[j*stride + i]/num_px;
    return mean;
}
