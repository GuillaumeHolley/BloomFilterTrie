#pragma once

/* ===================================================================================================================================
*  INCLUDES AND DEFINES
*  ===================================================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "./../lib/CC.h"
#include "./../lib/UC_annotation.h"
#include "./../lib/interface.h"

#define CC_VERTEX 1
#define UC_VERTEX 2
#define UC_PRESENT 4

void write_kmers_2disk(BFT_Root* root, char* filename, bool compressed_output);

void write_BFT_Root_sparse(BFT_Root*  root, char* filename);
void write_Node_sparse(Node*  node, BFT_Root* root, int lvl_node, FILE* file, int size_kmer);
void write_UC_sparse(UC* uc, BFT_Root* root, FILE* file, int size_substring, uint16_t nb_children);
void write_CC_sparse(CC*  cc, BFT_Root* root, int lvl_cc, FILE* file, int size_kmer);

BFT_Root* read_BFT_Root_sparse(char* filename);
void read_Node_sparse(Node*  node, BFT_Root* root, int lvl_node, FILE* file, int size_kmer);
void read_UC_sparse(UC* uc, FILE* file, int size_substring, uint16_t nb_children);
void read_CC_sparse(CC*  cc, BFT_Root* root, int lvl_cc, FILE* file, int size_kmer);
