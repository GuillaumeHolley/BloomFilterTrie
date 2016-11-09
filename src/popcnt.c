#include "popcnt.h"

extern inline uint8_t reverse_word_8(uint8_t v);
extern inline int popcnt_8(uint8_t v);
extern inline int popcnt_8_par(const uint8_t* v, int start, int end);
