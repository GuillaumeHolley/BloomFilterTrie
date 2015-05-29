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

void write_Root(Root* restrict root, char* filename, ptrs_on_func* restrict func_on_types);
void write_Node(Node* restrict node, FILE* file, int size_kmer, ptrs_on_func* restrict func_on_types);
void write_UC(UC* uc, FILE* file, int size_substring, uint16_t nb_children);
void write_CC(CC* restrict cc, FILE* file, int size_kmer, ptrs_on_func* restrict func_on_types);

Root* read_Root(char* filename);
void read_Node(Node* restrict node, FILE* file, int size_kmer, ptrs_on_func* restrict func_on_types);
void read_UC(UC* uc, FILE* file, int size_substring, uint16_t nb_children);
void read_CC(CC* restrict cc, FILE* file, int size_kmer, ptrs_on_func* restrict func_on_types);

#endif
