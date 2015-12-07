#ifndef DEF_INSERT_NODE
#define DEF_INSERT_NODE

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
#include "./../lib/presenceNode.h"
#include "./../lib/CC.h"
#include "./../lib/retrieveAnnotation.h"

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

void insertKmer_Node(Node* restrict node, Root* restrict root, int lvl_node, uint8_t* restrict suffix, int size_suffix,
                     uint8_t* kmer, uint32_t id_genome, int size_id_genome, int pos_CC_start_search,
                     annotation_inform* ann_inf, resultPresence* res, annotation_array_elem* annot_sorted);

Node* insertKmer_Node_special(Root* restrict root, resultPresence* restrict pres, int lvl_cont, uint8_t* restrict suffix, int size_suffix,
                              uint8_t* restrict kmer, uint32_t id_genome, int size_id_genome, annotation_inform* ann_inf,
                              annotation_array_elem* annot_sorted);

#endif
