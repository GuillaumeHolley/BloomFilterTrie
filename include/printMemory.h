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

#define FREE 0 //Debugging purpose: boolean indicating if the memory has to be free when iterating over the tree

/* ===================================================================================================================================
*  STRUCTURES DECLARATION
*  ===================================================================================================================================
*/

//memory_Used is a structure used to iterate on the tree and have information about a tree structure (Node, CC or UC) and all its children
typedef struct{
    double memory; //memory used, includes pointer, does not care about internal/external fragmentation or allocator class size rounding
    double nb_node_visited; //Number of nodes it is possible to visit from this tree structure (itself included)
    double nb_CC_visited; //Number of CCs it is possible to visit from this tree structure (itself included)
    double nb_UCptr_visited; //Number of UCs it is possible to visit from this tree structure (itself included)
    double nb_kmers_in_UCptr; //Number of suffixes stored in children's field and UCs from this tree structure (itself included)
    double nb_pointers_used; //Number of pointers used from this tree structure
    double size_biggest_annot; //Size of the biggest annotation in bytes from this tree structure

    double kmers_in_CC63; //Number of suffixes stored in nodes where suffixes in UCs are of length 63
    double nb_CC63; //Number of CC stored in nodes where suffixes in UCs are of length 63
    double kmers_in_CC54; //Number of suffixes stored in nodes where suffixes in UCs are of length 54
    double nb_CC54; //Number of CC stored in nodes where suffixes in UCs are of length 54
    double kmers_in_CC45; //Number of suffixes stored in nodes where suffixes in UCs are of length 45
    double nb_CC45; //Number of CC stored in nodes where suffixes in UCs are of length 45
    double kmers_in_CC36; //Number of suffixes stored in nodes where suffixes in UCs are of length 36
    double nb_CC36; //Number of CC stored in nodes where suffixes in UCs are of length 36
    double kmers_in_CC27; //Number of suffixes stored in nodes where suffixes in UCs are of length 27
    double nb_CC27; //Number of CC stored in nodes where suffixes in UCs are of length 27
    double kmers_in_CC18; //Number of suffixes stored in nodes where suffixes in UCs are of length 18
    double nb_CC18; //Number of CC stored in nodes where suffixes in UCs are of length 18
    double kmers_in_CC9; //Number of suffixes stored in nodes where suffixes in UCs are of length 9
    double nb_CC9; //Number of CC stored in nodes where suffixes in UCs are of length 9
} memory_Used;

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

memory_Used* create_memory_Used();
void add_memory_Used(memory_Used*  mem1, memory_Used*  mem2);
memory_Used* printMemoryUsedFromNode(Node*  node, int lvl_node, int size_kmer, info_per_level* info_per_lvl);
memory_Used* printMemoryUsedFrom_CC(CC*  cc, int lvl_cc, int size_kmer, info_per_level* info_per_lvl);
memory_Used* printMemoryUsedFrom_UC(UC*  uc);
void printMemory(memory_Used* mem);
