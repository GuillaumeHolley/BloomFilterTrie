/** \file include/snippets.h
 * Code snippets using a BFT.
 * The purpose of this file is to give examples of how to use the functions of the BFT API (bft.h).
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <stdarg.h>
#include <stdbool.h>

#include "bft.h"

#include "list.h"
#include "extract_kmers.h"

/** Flag for marked vertices indicating the vertex has not been visited before.*/
#ifndef V_NOT_VISITED
#define V_NOT_VISITED 0
#endif

/** Flag for marked vertices indicating the vertex has been visited before.*/
#ifndef V_VISITED
#define V_VISITED 1
#endif

/** @name K-mer extraction functions
 *  These functions extract k-mers stored in a BFT to disk.
 *  The last argument of extract_pangenome_kmers_to_disk() is of type BFT_func_ptr,
 *  you can use extract_core_kmers(), extract_dispensable_kmers() or extract_singleton_kmers()
 *  to extract core, dispensable or singleton k-mers.
 */

///@{
size_t extract_core_kmers(BFT_kmer* kmer, BFT* graph, va_list args);
size_t extract_dispensable_kmers(BFT_kmer* kmer, BFT* graph, va_list args);
size_t extract_singleton_kmers(BFT_kmer* kmer, BFT* graph, va_list args);
void extract_pangenome_kmers_to_disk(BFT* graph, char* filename_output, BFT_func_ptr f);
///@}

/** @name Path extraction functions
 *  These functions extract simple (non branching) paths of k-mers stored in a BFT to disk.
 */

///@{
size_t extract_simple_paths(BFT_kmer* kmer, BFT* graph, va_list args);
void extract_simple_paths_to_disk(BFT* graph, char* filename_output);
size_t extract_core_simple_paths(BFT_kmer* kmer, BFT* graph, va_list args);
void extract_simple_core_paths_to_disk(BFT* graph, double core_ratio, char* filename_output);
///@}

/** @name Graph traversal functions
 *  These functions iterate over a colored de Bruijn graph stored as a BFT.
 */

///@{
size_t BFS(BFT_kmer* kmer, BFT* graph, va_list args);
size_t BFS_subgraph(BFT_kmer* kmer, BFT* graph, va_list args);
size_t DFS(BFT_kmer* kmer, BFT* graph, va_list args);
size_t DFS_subgraph(BFT_kmer* kmer, BFT* graph, va_list args);
bool is_in_subgraph(BFT_kmer* kmer, BFT* graph, int nb_id_genomes, const va_list args);
void cdbg_traversal(BFT* graph, BFT_func_ptr f, ...);
size_t nb_connected_components(BFT_kmer* kmer, BFT* graph, va_list args);
void get_nb_connected_component(BFT* graph, ...);
///@}
