#ifndef DEF_HASH
#define DEF_HASH

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

#define rand_nb (0.5 * (sqrt(5)-1))
#define floor_rand_nb floor(rand_nb * INT16_MAX)

inline uint16_t hash1(uint8_t* kmer){
     return ((((uint16_t)kmer[0]) << 2) | (kmer[1] >> 6));
}

inline uint16_t hash2_bis(uint8_t* kmer){
    return ((((uint16_t)kmer[1]) << 2) | (kmer[2] >> 6));
}

inline uint16_t dbj2(uint8_t c1, uint8_t c2, int size_bf){
    unsigned long h = 177573 + c1;
    h = ((h << 5) + h) + c2;

    return ((uint16_t)h)%size_bf;
}

inline uint16_t sdbm(uint8_t c1, uint8_t c2, int size_bf){
    unsigned long h = c1;
    h = c2 + (h << 6) + (h << 16) - h;

    return ((uint16_t)h)%size_bf;
}

inline uint16_t KnutDivision(uint32_t substring_prefix, int size_bf){
    return (substring_prefix*(substring_prefix+3))%size_bf;
}

inline uint16_t MutiplicationMethod(uint32_t substring_prefix, int size_bf){
    return (((uint16_t)(substring_prefix*floor_rand_nb)) >> 5)%size_bf;
}

#endif
