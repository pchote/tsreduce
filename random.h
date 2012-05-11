/*
 * Copyright 2010, 2011 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <stdint.h>

#ifndef RANDOM_H
#define RANDOM_H

typedef struct
{
    uint32_t state[624];
    int left;
    uint32_t *next;
} random_generator;

random_generator *random_create(uint32_t seed);
void random_free(random_generator *g);
uint32_t random_uint32(random_generator *g);

#endif
