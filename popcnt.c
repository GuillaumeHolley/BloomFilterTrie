#include "popcnt.h"

uint8_t reverse_word_8(uint8_t v) {
  return rev_LSB[v & 0xf] | rev_MSB[v >> 4];
}

int popcnt_8(uint8_t v) {
    return tbl[v & 0xf] + tbl[v >> 4];
}

int popcnt_8_par(uint8_t* v, int start, int end) {
    int hamming_weight = 0;
    int i=0;

    for (i=start; i<end; i++){
        hamming_weight += tbl[v[i] & 0xf] + tbl[v[i] >> 4];
    }

    return hamming_weight;
}

/*int popcount_wikipedia_3(uint64_t *buf, int n) {
  int cnt=0;
  uint64_t x;
  do {
    x = *buf;
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits
    cnt += (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24)+...
    buf++;
  } while (--n);
  return cnt;
}*/
