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

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

void presenceNeighborsRight(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types, resultPresence* res, uint16_t** skip_node_root);
void presenceNeighborsLeft(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types, resultPresence* res, uint16_t** skip_node_root);
void presenceKmer(Node* restrict node, uint8_t* restrict kmer, int size_kmer, int posCC_start_search, ptrs_on_func* restrict func_on_types, resultPresence* res);
resultPresence* isKmerPresent(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types);

#endif
