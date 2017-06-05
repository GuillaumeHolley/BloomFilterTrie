#include "useful_macros.h"

extern inline void time_spent(struct timeval *start_time, struct timeval *end_time, struct timeval *resulting_time);
extern inline off_t fsize(const char *filename);

int uint32_t_cmp(const void *a, const void *b){

    if (*((uint32_t*)a) > *((uint32_t*)b)) return 1;
    if (*((uint32_t*)a) < *((uint32_t*)b))  return -1;

    return 0;
}
