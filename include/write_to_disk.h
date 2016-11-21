#pragma once

/* ===================================================================================================================================
*  INCLUDES AND DEFINES
*  ===================================================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "CC.h"
#include "bft.h"

#define CC_VERTEX 1
#define UC_VERTEX 2
#define UC_PRESENT 4

void write_kmers_2disk(BFT_Root* root, char* filename, bool compressed_output);

void write_BFT_Root(BFT_Root*  root, char* filename, bool write_root_only);
void write_Node(Node*  node, BFT_Root* root, int lvl_node, FILE* file, int size_kmer);
void write_UC(UC* uc, BFT_Root* root, FILE* file, int size_substring, uint16_t nb_children, bool compressed_header);
void write_CC(CC*  cc, BFT_Root* root, int lvl_cc, FILE* file, int size_kmer);

BFT_Root* read_BFT_Root(char* filename);
BFT_Root* read_BFT_Root_offset(char* filename, long int offset_read);

void read_Node(Node*  node, BFT_Root* root, int lvl_node, FILE* file, int size_kmer);
void read_UC(UC* uc, BFT_Root* root, FILE* file, int size_substring, uint16_t nb_children);
void read_CC(CC*  cc, BFT_Root* root, int lvl_cc, FILE* file, int size_kmer);

void append_BFT_2_comp_annots(BFT_Root* root, char* filename, bool write_root_only);

void read_BFT_replace_comp_annots(BFT_Root* root, char* filename_bft, char* filename_new_comp_colors, Pvoid_t* new_annots,
                                  int len_longest_annot);
void read_BFT_replace_comp_annots_bis(BFT_Root* root, char* filename_bft, char* filename_new_comp_colors, Pvoid_t* new_annots,
                                  int len_longest_annot, bool comp_annot_only);
void read_Node_replace_comp_annots(Node* node, BFT_Root* root, int lvl_node, FILE* file, int size_kmer, Pvoid_t* new_annots,
                                   bool comp_annot_only);
void read_CC_replace_comp_annots(CC* cc, BFT_Root* root, int lvl_cc, FILE* file, int size_kmer, Pvoid_t* new_annots,
                                 bool comp_annot_only);

void write_annotation_array_elem(char* filename, annotation_array_elem* annot_sorted, int size_array);
void read_annotation_array_elem(char* filename, annotation_array_elem** annot_sorted, int* size_array);

