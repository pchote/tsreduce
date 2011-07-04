/* Copyright 2010-2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

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