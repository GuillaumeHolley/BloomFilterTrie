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

#include "useful_macros.h"
#include "presenceNode.h"
#include "CC.h"
#include "retrieveAnnotation.h"

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

void insertKmer_Node(Node* restrict node, Node* restrict root, uint8_t* restrict suffix, int size_suffix, uint8_t* kmer, int size_kmer, int id_genome,
                     ptrs_on_func* restrict func_on_types, annotation_inform* ann_inf, resultPresence* res, annotation_array_elem* annot_sorted);

Node* insertKmer_Node_special(Node* root, resultPresence* restrict pres, uint8_t* restrict suffix, int size_suffix, uint8_t* restrict kmer, int size_kmer, int id_genome,
                              ptrs_on_func* restrict func_on_types, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

#endif
