#pragma once

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

/*int deleteColors_simplePath(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array, int shifting_suffix, uint32_t id_genome,
                            uint16_t** skip_node_root, info_per_level*  info_per_lvl, annotation_inform* ann_inf, resultPresence* res, annotation_array_elem* annot_sorted);

int deleteColors_from_branchingNodes(Node* n, Node* root, uint8_t* kmer, int size_kmer, int bucket, int pos_in_bucket, int size_kmer_root, uint32_t id_genome,
                                     info_per_level*  info_per_lvl, uint16_t** skip_node_root, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

int resize_annotation_Node(Node* n, int size_kmer, info_per_level*  info_per_lvl);
int resize_annotation_CC(CC* cc, int size_kmer, info_per_level*  info_per_lvl);
int resize_annotation_UC(UC* uc, int size_substring, int nb_children);*/
