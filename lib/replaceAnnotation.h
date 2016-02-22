#pragma once

/* ===================================================================================================================================
*  INCLUDES AND DEFINES
*  ===================================================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include "./../lib/useful_macros.h"
#include "./../lib/Node.h"
#include "./../lib/CC.h"
#include "./../lib/UC.h"

const uint8_t MASK_POWER_8[8];

/* ===================================================================================================================================
*  INLINE FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

void load_annotation_from_Node(Node*  node, int lvl_node, int size_kmer, int longest_annot, info_per_level*  info_per_lvl,
                               Pvoid_t* PJArray, annotation_array_elem* annot_sorted, annotation_inform* ann_inf);
void load_annotation_from_CC(CC*  cc, int lvl_cc, int size_kmer, int longest_annot, info_per_level*  info_per_lvl,
                             Pvoid_t* PJArray, annotation_array_elem* annot_sorted, annotation_inform* ann_inf);
void load_annotation_from_UC(UC* uc, int size_substring, int nb_children, int longest_annot, Pvoid_t* PJArray,
                             annotation_array_elem* annot_sorted, annotation_inform* ann_inf);

int compress_annotation_from_Node(Node*  node, int lvl_node, int size_kmer, info_per_level*  info_per_lvl,
                                  Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted, annotation_inform* ann_inf);
int compress_annotation_from_CC(CC*  cc, int lvl_cc, int size_kmer, info_per_level*  info_per_lvl,
                                Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted, annotation_inform* ann_inf);
int compress_annotation_from_UC(UC* uc, int size_substring, int nb_children, Pvoid_t* PJArray,
                                annotation_array_elem* old_annot_sorted, annotation_inform* ann_inf);
