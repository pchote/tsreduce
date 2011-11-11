/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */


#ifndef REDUCTION_H
#define REDUCTION_H

int create_flat(const char *pattern, int minmax, const char *masterdark, const char *outname);
int create_dark(const char *pattern, int minmax, const char *outname);


int reduce_single_frame(char *framePath, char *darkPath, char *flatPath, char *outPath);
int update_reduction(char *dataPath);
int create_reduction_file(char *framePath, char *framePattern, char *darkTemplate, char *flatTemplate, char *filePath);


#endif
