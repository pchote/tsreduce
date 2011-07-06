/*
* Copyright 2010, 2011 Paul Chote
* This file is part of Puoko-nui, which is free software. It is made available
* to you under the terms of version 3 of the GNU General Public License, as
* published by the Free Software Foundation. For more information, see LICENSE.
*/

#include <stdarg.h>

#ifndef HELPERS_H
#define HELPERS_H

int get_matching_files(const char *cmd, char **files, int numFiles);
int error(const char * format, ...);
void die(const char * format, ...);
int init_ds9(char *);
int command_ds9(char *title, char *command, void *data, int dataSize);
#endif