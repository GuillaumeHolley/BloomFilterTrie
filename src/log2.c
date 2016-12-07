#include "log2.h"

const int MultiplyDeBruijnBitPosition_power2[32] =
{
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};

extern inline int get_log2_upper_power2(uint32_t pos_rounded_power_of_2);
extern inline uint32_t round_up_next_highest_power2(uint32_t pos);
extern inline int get_nb_bytes_power2_comp(uint32_t pos);
extern inline int get_nb_bytes_power2_annot(uint32_t pos);
extern inline int get_nb_bytes_power2_annot_bis(uint32_t pos, uint32_t pow2_pos);

