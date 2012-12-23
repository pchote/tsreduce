/*
 * Copyright 2010, 2011, 2012 Paul Chote
 * This file is part of Puoko-nui, which is free software. It is made available
 * to you under the terms of version 3 of the GNU General Public License, as
 * published by the Free Software Foundation. For more information, see LICENSE.
 */

#ifndef RANDOM_H
#define RANDOM_H

#include <stdint.h>

typedef struct random_generator random_generator;

random_generator *random_create(uint32_t seed);
void random_free(random_generator *g);
uint32_t random_uint32(random_generator *g);
uint32_t random_uint32_max(random_generator *g, uint32_t n);
double random_double(random_generator *g);
double random_normal(random_generator *g, double mu, double sigma);
void random_shuffle_double_array(random_generator *g, double *a, size_t n);

#endif
