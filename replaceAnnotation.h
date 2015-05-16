#ifndef DEF_REPLACE_ANNOTATION
#define DEF_REPLACE_ANNOTATION

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

void load_annotation_from_Node(Node* restrict node, int size_kmer, ptrs_on_func* restrict func_on_types, Pvoid_t* PJArray, annotation_array_elem* annot_sorted);
void load_annotation_from_CC(CC* restrict cc, int size_kmer, ptrs_on_func* restrict func_on_types, Pvoid_t* PJArray, annotation_array_elem* annot_sorted);
void load_annotation_from_UC(UC* uc, int size_substring, int nb_children, Pvoid_t* PJArray, annotation_array_elem* annot_sorted);

int compress_annotation_from_Node(Node* restrict node, int size_kmer, ptrs_on_func* restrict func_on_types, Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted);
int compress_annotation_from_CC(CC* restrict cc, int size_kmer, ptrs_on_func* restrict func_on_types, Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted);
int compress_annotation_from_UC(UC* uc, int size_substring, int nb_children, Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted);

void uncompress_annotation_from_Node(Node* restrict node, int size_kmer, ptrs_on_func* restrict func_on_types, annotation_array_elem* annot_sorted);
void uncompress_annotation_from_CC(CC* restrict cc, int size_kmer, ptrs_on_func* restrict func_on_types, annotation_array_elem* annot_sorted);
void uncompress_annotation_from_UC(UC* uc, int size_substring, int nb_children, annotation_array_elem* annot_sorted);

#endif
