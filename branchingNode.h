#ifndef DEF_BRANCHING_NODE
#define DEF_BRANCHING_NODE

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
#include "CC.h"
#include "presenceNode.h"
#include "insertNode.h"

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

int getBranchingNode(Node* n, Node* root, Node* tree_branching_nodes, uint8_t* kmer, int size_kmer, int bucket, int pos_in_bucket,
                     int size_kmer_root, ptrs_on_func* restrict func_on_types, uint16_t** skip_node_root, annotation_inform* ann_inf, resultPresence* res);

int isBranchingLeft(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types, uint16_t** skip_node_root);
int isBranchingRight(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types, uint16_t** skip_node_root);
resultPresence* getRightNeighbors(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types, uint16_t** skip_node_root);
resultPresence* getLeftNeighbors(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types, uint16_t** skip_node_root);
void extractSuffix(uint8_t* kmer_tmp, int size_kmer, int size_kmer_array, int shifting_prefix, int it_bucket, uint8_t* suffix_start, ptrs_on_func* restrict func_on_types);

#endif
