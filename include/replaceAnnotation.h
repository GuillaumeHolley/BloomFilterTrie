#pragma once

/* ===================================================================================================================================
*  INCLUDES AND DEFINES
*  ===================================================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include "useful_macros.h"
#include "Node.h"
#include "CC.h"
#include "UC.h"

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

int compress_annotation_from_UC(UC* uc, int size_substring, int nb_children,
                                Pvoid_t* new_annots, annotation_array_elem* old_annots,
                                annotation_inform* ann_inf, bool insert_comp_annot);

