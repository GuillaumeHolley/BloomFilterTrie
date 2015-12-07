#include "./../lib/log2.h"

int get_log2_upper_power2(uint32_t pos_rounded_power_of_2){
    return MultiplyDeBruijnBitPosition_power2[(uint32_t)(pos_rounded_power_of_2 * 0x077CB531U) >> 27];
}

int round_up_next_highest_power2(uint32_t pos){

    pos--;
    pos |= pos >> 1;
    pos |= pos >> 2;
    pos |= pos >> 4;
    pos |= pos >> 8;
    pos |= pos >> 16;
    pos++;

    return pos;
}

int get_nb_bytes_power2_comp(uint32_t pos){

    int log2_upper_pos = get_log2_upper_power2(pos);

    if (log2_upper_pos > 6) return CEIL(log2_upper_pos-6, 7) + 1;
    else return 1;
}

int get_nb_bytes_power2_annot(uint32_t pos){

    int pow2_pos = round_up_next_highest_power2(pos);

    if (pow2_pos == pos) return CEIL(get_log2_upper_power2(pow2_pos)+1, 6);
    else return CEIL(get_log2_upper_power2(pow2_pos), 6);
}
