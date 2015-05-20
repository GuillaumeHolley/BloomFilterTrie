#ifndef DEF_RETRIEVE_ANNOTATION
#define DEF_RETRIEVE_ANNOTATION

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#include "./../lib/Node.h"
#include "./../lib/UC_annotation.h"
#include "./../lib/CC.h"
#include "./../lib/branchingNode.h"

annotation_array_elem* retrieve_annotation_left(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array, int id_genome_avoid,
                        uint16_t** skip_node_root, ptrs_on_func* restrict func_on_types, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

annotation_array_elem* retrieve_annotation_right(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array, int shifting_suffix,
                        int id_genome_avoid, uint16_t** skip_node_root, ptrs_on_func* restrict func_on_types, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

annotation_array_elem* retrieve_annotation(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array, int shifting_suffix,
                        int id_genome_avoid, uint16_t** skip_node_root, ptrs_on_func* restrict func_on_types, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

void modify_annotations(Node* root, UC* uc, int size_substring, int nb_children, int position, int id_genome, uint8_t* kmer, int size_kmer,
                        ptrs_on_func* restrict func_on_types, annotation_inform* ann_inf, annotation_array_elem* annot_sorted, int special_case);

#endif
