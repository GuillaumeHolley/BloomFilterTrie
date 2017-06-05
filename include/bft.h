/** \file include/bft.h
 * Interface containing all functions to use a BFT.
 * Code snippets using this interface are provided in snippets.h.
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <stdarg.h>
#include <stdbool.h>

#include "intersection.h"
#include "insertNode.h"
#include "branchingNode.h"
#include "fasta.h"
#include "marking.h"
#include "replaceAnnotation.h"
#include "write_to_disk.h"
#include "extract_kmers.h"
#include "printMemory.h"
#include "file_io.h"
#include "useful_macros.h"

/** @copydoc BFT_Root */
typedef BFT_Root BFT;

/** Annotation associated with a BFT_kmer.
 *  A BFT_annotation contains the set of colors associated with a k-mer of a BFT_Root.
 */
typedef struct{
    uint8_t* annot;
    uint8_t* annot_ext;
    uint8_t* annot_cplx;

    int size_annot;
    int size_annot_cplx;

    uint8_t from_BFT;
} BFT_annotation;

/** Pointer on function used by iterate_over_kmers() and v_iterate_over_kmers().
    Such a function (user written) is called on every k-mer of a BFT.
*   @param bft_kmer is a k-mer from a BFT.
*   @param bft is the BFT from which bft_kmer is from.
*   @param args contains all additional parameters given to iterate_over_kmers() / v_iterate_over_kmers().
*   @return a size_t type object. It can be use to return an unsigned integer or a pointer.
*/
typedef size_t (*BFT_func_ptr) (BFT_kmer* bft_kmer, BFT* bft, va_list args);

inline uint8_t intersection_annots(const uint8_t a, const uint8_t b);
inline uint8_t union_annots(const uint8_t a, const uint8_t b);
inline uint8_t sym_difference_annots(const uint8_t a, const uint8_t b);

/** @name Graph functions
 *  These functions manipulate a colored de Bruijn graph stored in a BFT.
 */

///@{
BFT* create_cdbg(int k, int treshold_compression);
void free_cdbg(BFT* bft);
///@}

/** @name Insertion functions
 *  These functions insert genomes in a colored de Bruijn graph stored in a BFT.
 */

///@{
void insert_genomes_from_files(int nb_files, char** paths, BFT* bft, char* prefix_bft_filename);
void insert_kmers_new_genome(int nb_kmers, char** kmers, char* genome_name, BFT* bft);
void insert_kmers_last_genome(int nb_kmers, char** kmers, BFT* bft);
///@}

/** @name K-mer functions
 *  These functions manipulate k-mers.
 */

///@{
BFT_kmer* create_kmer(const char* kmer, int k);
BFT_kmer* create_empty_kmer();
void free_BFT_kmer(BFT_kmer* bft_kmer, int nb_bft_kmer);
void free_BFT_kmer_content(BFT_kmer* bft_kmer, int nb_bft_kmer);
void extract_kmers_to_disk(BFT* bft, char* filename_output, bool compressed_output);
size_t write_kmer_ascii_to_disk(BFT_kmer* bft_kmer, BFT* bft, va_list args);
size_t write_kmer_comp_to_disk(BFT_kmer* bft_kmer, BFT* bft, va_list args);
///@}

/** @name Annotation functions
 *  These functions manipulate annotations (color sets).
 */

///@{
BFT_annotation* create_BFT_annotation();
void free_BFT_annotation(BFT_annotation* bft_annot);
BFT_annotation* get_annotation(BFT_kmer* bft_kmer);
bool presence_genome(uint32_t id_genome, BFT_annotation* bft_annot, BFT* bft);

inline uint8_t intersection_annots(const uint8_t a, const uint8_t b){
    return a & b;
}

inline uint8_t union_annots(const uint8_t a, const uint8_t b){
    return a | b;
}

inline uint8_t sym_difference_annots(const uint8_t a, const uint8_t b){
    return a ^ b;
}

BFT_annotation* intersection_annotations(BFT* bft, uint32_t nb_annotations, ... );
BFT_annotation* union_annotations(BFT* bft, uint32_t nb_annotations, ... );
BFT_annotation* sym_difference_annotations(BFT* bft, uint32_t nb_annotations, ... );
uint32_t* get_list_id_genomes(BFT_annotation* bft_annot, BFT* bft);
uint32_t get_count_id_genomes(BFT_annotation* bft_annot, BFT* bft);
uint32_t* intersection_list_id_genomes(uint32_t* list_a, uint32_t* list_b);
///@}

/** @name Query functions
 *  These functions query for k-mers or sequences.
 */

///@{
BFT_kmer* get_kmer(const char* kmer, BFT* bft);
bool is_kmer_in_cdbg(BFT_kmer* bft_kmer);
uint32_t* query_sequence(BFT* bft, char* sequence, double threshold, bool canonical_search);
///@}

/** @name Pattern matching functions
 *  These functions provide pattern matching functionalities over the k-mers or paths of a colored de Bruijn graph stored as a BFT.
 */

///@{
bool prefix_matching(BFT* bft, char* prefix, BFT_func_ptr f, ...);
///@}

/** @name Marking functions
 *  These functions allow to mark k-mers of a colored de Bruijn graph with flags.
 */

///@{
void set_marking(BFT* bft);
void unset_marking(BFT* bft);
void set_flag_kmer(uint8_t flag, BFT_kmer* bft_kmer, BFT* bft);
uint8_t get_flag_kmer(BFT_kmer* bft_kmer, BFT* bft);
///@}

/** @name Traversal functions
 *  These functions allow to traverse a colored de Bruijn graph stored as a BFT.
 */

///@{
void set_neighbors_traversal(BFT* bft);
void unset_neighbors_traversal(BFT* bft);
BFT_kmer* get_neighbors(BFT_kmer* bft_kmer, BFT* bft);
BFT_kmer* get_predecessors(BFT_kmer* bft_kmer, BFT* bft);
BFT_kmer* get_successors(BFT_kmer* bft_kmer, BFT* bft);
///@}

/** @name Iteration functions
 *  These functions iterate over the k-mers of a colored de Bruijn graph stored as a BFT.
 */

///@{
void iterate_over_kmers(BFT* bft, BFT_func_ptr f, ... );
void v_iterate_over_kmers(BFT* bft, BFT_func_ptr f, va_list args);
///@}

/** @name Disk I/O functions
 *  These functions write and load a BFT from disk.
 */

///@{
void write_BFT(BFT* bft, char* filename, bool compress_annotations);
BFT* load_BFT(char* filename);
///@}

BFT* create_cdbg_from_bft_kmers(BFT_kmer** bft_kmers, uint32_t nb_bft_kmers, BFT* bft, bool add_colors);
void add_id_genomes(BFT_kmer* bft_kmer, BFT_annotation* bft_annot, BFT* bft, uint32_t* list_id_genomes);
