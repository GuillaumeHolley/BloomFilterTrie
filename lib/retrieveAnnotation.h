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

/*annotation_array_elem* retrieve_annotation_left(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array, uint32_t id_genome_avoid,
                        uint16_t** skip_node_root, info_per_level* restrict info_per_lvl, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

annotation_array_elem* retrieve_annotation_right(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array, int shifting_suffix,
                        uint32_t id_genome_avoid, uint16_t** skip_node_root, info_per_level* restrict info_per_lvl, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

annotation_array_elem* retrieve_annotation(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array, int shifting_suffix,
                        uint32_t id_genome_avoid, uint16_t** skip_node_root, info_per_level* restrict info_per_lvl, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);
*/

void modify_annotations(Root* root, UC* uc, int size_substring, int nb_children, int position, uint32_t id_genome,
                        int size_id_genome, uint8_t* kmer, annotation_inform* ann_inf, annotation_array_elem* annot_sorted,
                        int special_case);

#endif
