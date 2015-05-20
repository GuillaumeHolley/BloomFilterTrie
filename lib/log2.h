#ifndef DEF_LOG2
#define DEF_LOG2

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#include "./../lib/useful_macros.h"

static const int MultiplyDeBruijnBitPosition2[32] =
{
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};

int get_log2_upper_power2(uint32_t pos_rounded_power_of_2);
int get_nb_bytes_power2(uint32_t pos);
int round_up_next_highest_power2(uint32_t pos);

#endif
