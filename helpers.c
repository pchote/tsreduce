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
#include "helpers.h"

// Write the first numFiles filenames that match the system command cmd into
// the array files (allocated by the caller).
int get_matching_files(const char *cmd, char **files, int numFiles)
{    
    FILE *p = popen(cmd, "r");
    if (p == NULL)
        return error("failed to run command `%s`", cmd);
    
    int n = 0;
    char *buf = (char *)malloc(PATH_MAX*sizeof(char));
    while (fgets(buf, PATH_MAX-1, p) != NULL)
    {
        if (n == numFiles)
        {
            fprintf(stderr, "Warning: Exceeded file limit of %d frames. Extra frames will be discarded", numFiles);
            break;
        }
        
        int len = strlen(buf);
        if (len <= 0)
            continue;
        // Copy then strip the trailing newline
        strncpy(files[n], buf, len);
        files[n][len-1] = '\0';
        n++;
    }
    free(buf);
    pclose(p);
    return n;
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
    return valid;
}

int ask_ds9(char *title, char *command, char *outbuf, int outlen)
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
    {
        if (outlen < len[0])
        {
            valid = 0;
            printf("Error: message buffer too small for returned ds9 output");
        }

        strncpy(outbuf, ret[0], outlen);
        outbuf[outlen - 1] = '\0';
        free(ret[0]);
    }
    return valid;
}
