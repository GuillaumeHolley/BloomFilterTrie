#ifndef DEF_DELETE_COLORS_NODE
#define DEF_DELETE_COLORS_NODE

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
#include "./../lib/branchingNode.h"
#include "./../lib/CC.h"
#include "./../lib/retrieveAnnotation.h"
#include "./../lib/marking.h"

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

int deleteColors_simplePath(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array, int shifting_suffix, int id_genome,
                            uint16_t** skip_node_root, ptrs_on_func* restrict func_on_types, annotation_inform* ann_inf, resultPresence* res, annotation_array_elem* annot_sorted);

int deleteColors_from_branchingNodes(Node* n, Node* root, uint8_t* kmer, int size_kmer, int bucket, int pos_in_bucket, int size_kmer_root, int id_genome,
                                     ptrs_on_func* restrict func_on_types, uint16_t** skip_node_root, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

int resize_annotation_Node(Node* n, int size_kmer, ptrs_on_func* restrict func_on_types);
int resize_annotation_CC(CC* cc, int size_kmer, ptrs_on_func* restrict func_on_types);
int resize_annotation_UC(UC* uc, int size_substring, int nb_children);

#endif
