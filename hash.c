#include "hash.h"


uint16_t hash2(uint8_t kmer0_1, uint8_t kmer0_2, uint8_t kmer1_1, uint8_t kmer1_2, uint8_t kmer2){
    uint16_t hash_v = A_OR_G[0][kmer0_1] | A_OR_G[1][kmer0_2] | A_OR_G[2][kmer1_1] | A_OR_G[3][kmer1_2] | A_OR_G[4][kmer2];

    if ((count_G[kmer0_1]+count_G[kmer0_2]+count_G[kmer1_1]+count_G[kmer1_2]+count_G[kmer2])
        > (count_A[kmer0_1]+count_A[kmer0_2]+count_A[kmer1_1]+count_A[kmer1_2]+count_A[kmer2])){
            hash_v |= 1;
    }

    return hash_v;
}

uint16_t hash3(uint8_t kmer0_1, uint8_t kmer0_2, uint8_t kmer1_1, uint8_t kmer1_2, uint8_t kmer2){
    uint16_t hash_v = A_OR_C[0][kmer0_1] | A_OR_C[1][kmer0_2] | A_OR_C[2][kmer1_1] | A_OR_C[3][kmer1_2] | A_OR_C[4][kmer2];

    if ((count_C[kmer0_1]+count_C[kmer0_2]+count_C[kmer1_1]+count_C[kmer1_2]+count_C[kmer2])
        > (count_A[kmer0_1]+count_A[kmer0_2]+count_A[kmer1_1]+count_A[kmer1_2]+count_A[kmer2])){
            hash_v |= 1;
    }

    return hash_v;
}

uint16_t hash1(uint8_t* kmer){
     return ((((uint16_t)kmer[1]) << 2) | (kmer[2] >> 6));
}

uint16_t hash2_bis(uint8_t kmer0_1, uint8_t kmer0_2, uint8_t kmer1_1, uint8_t kmer1_2, uint8_t kmer2){
    uint16_t hash_v = A_OR_G[0][kmer0_1] | A_OR_G[1][kmer0_2] | A_OR_G[2][kmer1_1] | A_OR_G[3][kmer1_2] | A_OR_G[4][kmer2];
    return (hash_v >> 1);
}

uint16_t hash3_bis(uint8_t kmer0_1, uint8_t kmer0_2, uint8_t kmer1_1, uint8_t kmer1_2, uint8_t kmer2){
    uint16_t hash_v = A_OR_C[0][kmer0_1] | A_OR_C[1][kmer0_2] | A_OR_C[2][kmer1_1] | A_OR_C[3][kmer1_2] | A_OR_C[4][kmer2];
    return (hash_v >> 1);
}

uint16_t hash1_bis(uint8_t* kmer){
     return ((((uint16_t)kmer[1]) << 1) | (kmer[2] >> 7));
}

uint16_t hash4_bis(uint8_t* kmer){
     return ((((uint16_t)kmer[2]) << 1) | (kmer[3] >> 7));
}

uint16_t hash5_bis(uint8_t* kmer){
     return (((((uint16_t)kmer[1]) << 4) | (kmer[2] >> 4))) & 0x1ff;
}
