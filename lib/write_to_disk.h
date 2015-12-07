#ifndef DEF_WRITE2DISK
#define DEF_WRITE2DISK

/* ===================================================================================================================================
*  INCLUDES AND DEFINES
*  ===================================================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "./../lib/CC.h"
#include "./../lib/UC_annotation.h"

#define CC_VERTEX 1
#define UC_VERTEX 2
#define UC_PRESENT 4

/*void write_Root(Root* restrict root, char* filename);
void write_Node(Node* restrict node, Root* root, int lvl_node, FILE* file, int size_kmer);
void write_CC(CC* restrict cc, Root* root, int lvl_cc, FILE* file, int size_kmer);
void write_UC(UC* uc, FILE* file, int size_substring, uint16_t nb_children);*/

void write_Root_sparse(Root* restrict root, char* filename);
void write_Node_sparse(Node* restrict node, Root* root, int lvl_node, FILE* file, int size_kmer);
void write_UC_sparse(UC* uc, Root* root, FILE* file, int size_substring, uint16_t nb_children);
void write_CC_sparse(CC* restrict cc, Root* root, int lvl_cc, FILE* file, int size_kmer);

/*Root* read_Root(char* filename);
void read_Node(Node* restrict node, Root* root, int lvl_node, FILE* file, int size_kmer);
void read_UC(UC* uc, FILE* file, int size_substring, uint16_t nb_children);
void read_CC(CC* restrict cc, Root* root, int lvl_cc, FILE* file, int size_kmer);*/

Root* read_Root_sparse(char* filename);
void read_Node_sparse(Node* restrict node, Root* root, int lvl_node, FILE* file, int size_kmer);
void read_UC_sparse(UC* uc, FILE* file, int size_substring, uint16_t nb_children);
void read_CC_sparse(CC* restrict cc, Root* root, int lvl_cc, FILE* file, int size_kmer);

#endif
