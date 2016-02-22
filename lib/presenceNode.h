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

#include "./../lib/useful_macros.h"
#include "./../lib/Node.h"
#include "./../lib/CC.h"
#include "./../lib/write_to_disk.h"

/* ===================================================================================================================================
*  INLINE FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

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

void presenceNeighborsRight(Node*  node, BFT_Root* root, uint8_t*  kmer, int size_kmer, resultPresence* res);
void presenceNeighborsLeft(Node*  node, BFT_Root* root, uint8_t*  kmer, int size_kmer, resultPresence* res);
void presenceKmer(Node*  node, BFT_Root* root, uint8_t*  kmer, int size_kmer, int posCC_start_search,
                  int verify_UC_or_not, resultPresence* res);

void presenceSeed(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  seed, int size_seed, int size_kmer,
                  int posCC_start_search, int id_genome, int search_for_id_genome, int stop_if_2idgenomes_found,
                  resultPresenceSeed* res, annotation_inform* ann_inf, annotation_array_elem* annot_sorted, Duo* d1);

resultPresence* isKmerPresent(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer);

void isKmerPresent_bis(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer, resultPresence* res);

void isSeedPresent(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  seed, int size_seed, int size_kmer,
                   int id_genome, int search_for_id_genome, int stop_if_2idgenomes_found, annotation_inform* ann_inf,
                   annotation_array_elem* annot_sorted, Duo* d1);

void countKmers_Node(Node*  node, int lvl_node, int size_kmer, int id_genome,
                     int search_for_id_genome, int stop_if_2idgenomes_found, info_per_level*  info_per_lvl,
                     annotation_inform* ann_inf, annotation_array_elem* annot_sorted, Duo* d1);
void countKmers_CC(CC*  cc, int lvl_cont, int size_kmer, int id_genome, int search_for_id_genome,
                   int stop_if_2idgenomes_found, info_per_level*  info_per_lvl,
                   annotation_inform* ann_inf, annotation_array_elem* annot_sorted, Duo* d1);

void findCluster(CC* cc, int pos_filter2, int* cpt_node_return, int* pos_extra_filter3,
                 int* hamming_weight_0, resultPresence* res, info_per_level*  info_per_lvl);

int* get_bf_presence_per_cc(BFT_Root* root);

#endif
