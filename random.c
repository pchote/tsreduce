/* Based on 2002 improvement of MT19937
 * by Takuji Nishimura and Makoto Matsumoto and made
 * avalable under the terms of the BSD license. 
 */

#include <stdlib.h>
#include "random.h"

// Period parameters 
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   // constant vector a
#define UMASK 0x80000000UL // most significant w-r bits
#define LMASK 0x7fffffffUL // least significant r bits
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

random_generator *random_create(uint32_t seed)
{
    random_generator *ret = (random_generator *)malloc(sizeof(random_generator));

    ret->state[0] = seed;
    for (int j = 1; j < N; j++)
        ret->state[j] = (1812433253UL * (ret->state[j-1] ^ (ret->state[j-1] >> 30)) + j); 

    ret->left = 1;
    return ret;
}

void random_free(random_generator *g)
{
    free(g);
}

static void random_next_state(random_generator *g)
{
    uint32_t *p = g->state;

    g->left = N;
    g->next = g->state;

    for (int j = N - M + 1; --j; p++) 
        *p = p[M] ^ TWIST(p[0], p[1]);

    for (int j = M; --j; p++) 
        *p = p[M-N] ^ TWIST(p[0], p[1]);

    *p = p[M-N] ^ TWIST(p[0], g->state[0]);
}

// Get a random number between [0,2^32-1)
uint32_t random_uint32(random_generator *g)
{    
    if (--g->left == 0)
        random_next_state(g);

    uint32_t y = *(g->next)++;

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

// Get a uniformly distributed random number between [0,n)
uint32_t random_uint32_max(random_generator *g, uint32_t n)
{
    // Avoid modulo bias by taking only the range that maps equally into [0,n)
    uint32_t max = UINT32_MAX - UINT32_MAX % n;
    uint32_t r;
    do
    {
        r = random_uint32(g);
    } while (r >= max);

    return (uint32_t)(r % n);
}

// Shuffle an array of doubles in place using the Fisher–Yates algorithm
void shuffle_double_array(double *a, size_t n, random_generator *g)
{
    for (size_t i = n - 1; i > 0; i--)
    {
        size_t j = random_uint32_max(g, i + 1);
        double t = a[j];
        a[j] = a[i];
        a[i] = t;
    }
}
