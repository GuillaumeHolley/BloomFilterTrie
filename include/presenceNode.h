#ifndef DEF_PRESENCE_NODE
#define DEF_PRESENCE_NODE

/* ===================================================================================================================================
*  INCLUDES
*  ===================================================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "useful_macros.h"
#include "Node.h"
#include "CC.h"
#include "write_to_disk.h"
#include "extract_kmers.h"

/* ===================================================================================================================================
*  INLINE FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

inline int suffix_matching(CC* cc, int position, uint8_t mask, uint8_t suffix, uint8_t size_suffix_bit);

inline int suffix_matching(CC* cc, int position, uint8_t mask, uint8_t suffix, uint8_t size_suffix_bit){

    ASSERT_NULL_PTR(cc, "suffix_matching()")

    if (size_suffix_bit==8){
        if ((cc->filter3[position] & mask) == suffix) return 1;
    }
    else if (IS_ODD(position)){
        if (((cc->filter3[position/2] >> 4) & mask) == suffix) return 1;
    }
    else if (((cc->filter3[position/2] & 15) & mask) == suffix) return 1;

    return 0;
}

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

void presenceNeighborsRight(Node* node, BFT_Root* root, uint8_t*  kmer, int size_kmer, resultPresence* res);

void presenceNeighborsLeft(Node* node, BFT_Root* root, uint8_t*  kmer, int size_kmer, resultPresence* res);

void presenceKmer(Node* node, BFT_Root* root, uint8_t* kmer, int size_kmer, int posCC_start_search,
                  int verify_UC_or_not, resultPresence* res);

resultPresence* isKmerPresent(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer);

void findCluster(CC* cc, int pos_filter2, int* cpt_node_return, int* pos_extra_filter3,
                 int* hamming_weight_0, resultPresence* res, info_per_level*  info_per_lvl);

int* get_bf_presence_per_cc(BFT_Root* root);

bool prefix_matching_comp(BFT_Root* bft, uint8_t* prefix_comp, int length_prefix, size_t (*f)(BFT_kmer*, BFT_Root*, va_list), ...);

size_t v_prefix_matching(Node* node, BFT_Root* root, int lvl_node,
                         uint8_t* shifted_prefix, uint8_t* prefix, int size_prefix_shifted, int size_prefix,
                         BFT_kmer* bft_kmer, size_t (*f)(BFT_kmer*, BFT_Root*, va_list), va_list args);

size_t v_prefix_matching_custom(Node* node, BFT_Root* root, int lvl_node,
                                uint8_t* shifted_prefix, uint8_t* prefix, int size_prefix_shifted, int size_prefix,
                                BFT_kmer* bft_kmer, resultPresence** junction_prefix, size_t (*f)(BFT_kmer*, BFT_Root*, va_list), va_list args);

int memcmp_bits(uint8_t* s1, uint8_t* s2, int length);

#endif
