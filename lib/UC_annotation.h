#ifndef DEF_UC_ANNOTATION
#define DEF_UC_ANNOTATION

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#include <Judy.h>
#include "./../lib/useful_macros.h"

#define SIZE_CELL 8
#define KMER_LENGTH_MAX 63

#define SIZE_BYTE_EXT_ANNOT 3
#define SIZE_BYTE_CPLX_N 2

#define NB_MAX_ID_GENOMES 2000
#define NB_CELL_GENOMES CEIL(NB_MAX_ID_GENOMES+2,SIZE_CELL)

/* ===================================================================================================================================
*  STRUCTURES DECLARATION
*  ===================================================================================================================================
*/

static const uint8_t MASK_POWER_8[8] = {1, 2, 4, 8, 16, 32, 64, 128};

static const char* presence_genomes[16] = {"0,0,0,0,", "1,0,0,0,", "0,1,0,0,", "1,1,0,0,",
                                    "0,0,1,0,", "1,0,1,0,", "0,1,1,0,", "1,1,1,0,",
                                    "0,0,0,1,", "1,0,0,1,", "0,1,0,1,", "1,1,0,1,",
                                    "0,0,1,1,", "1,0,1,1,", "0,1,1,1,", "1,1,1,1,"};

typedef struct {
    //Represent the genomes IDs stored in the annotation, doesn't matter if they are present or not
    uint16_t id_stored[NB_MAX_ID_GENOMES];
    int nb_id_stored;

    uint8_t annotation[NB_CELL_GENOMES];

    //Last ID added to the current annotation before annotation
    int last_added;

    //Number of genomes present and not present.
    //Number of genomes present and not present for which the ID is more than 6 bits
    int nb_genome_pres;
    int nb_ext_pres;

    //Current mode of the annation, new optimal mode and the required size
    int current_mode;
    int min_mode;
    int min_size;
} annotation_inform;

typedef struct {
    int64_t last_index;
    uint8_t* annot_array;
    uint8_t size_annot;
} annotation_array_elem;

typedef struct{
    uint8_t* suffixes;
    uint8_t size_annot;
    uint8_t size_annot_cplx_nodes;
    uint16_t nb_extended_annot;
    uint16_t nb_cplx_nodes;
    uint16_t nb_children;
} UC;

typedef struct{
    uint8_t* suffixes;
    uint8_t size_annot;
    uint8_t size_annot_cplx_nodes;
    uint16_t nb_extended_annot;
    uint16_t nb_cplx_nodes;
} __attribute__ ((__packed__)) UC_9;



inline void reinit_annotation_inform(annotation_inform* ann_inf){
    ann_inf->nb_id_stored = 0;
    ann_inf->last_added = 0;
    ann_inf->nb_genome_pres = 0;
    ann_inf->nb_ext_pres = 0;
    ann_inf->current_mode = 0;
    ann_inf->min_mode = 0;
    ann_inf->min_size = 0;
    memset(ann_inf->annotation, 0, NB_CELL_GENOMES*sizeof(uint8_t));
    return;
}

inline uint8_t* min_size_per_annot_cplx(UC* uc, int nb_substrings, int size_substring){

    if (nb_substrings > 0){
        ASSERT_NULL_PTR(uc, "min_size_per_annot_cplx()")
        ASSERT_NULL_PTR(uc->suffixes, "min_size_per_annot_cplx()")
    }

    int i, z, base;
    int pos = 0;
    int size_line = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;

    uint8_t* sizes = calloc(nb_substrings, sizeof(uint8_t));
    ASSERT_NULL_PTR(sizes, "min_size_per_annot_cplx()")

    uint8_t* annot_cplx = &(uc->suffixes[nb_substrings * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]);

    for (i = 0; i < uc->nb_cplx_nodes; i++){

        z = uc->size_annot_cplx_nodes - 1;
        base = i * size_line;
        pos += ((((uint16_t)annot_cplx[base]) << SIZE_CELL) | ((uint16_t)annot_cplx[base + 1]));
        base += SIZE_BYTE_CPLX_N;

        while ((z >= 0) && (annot_cplx[base + z] == 0)) z--;
        sizes[i] = z+1;
    }

    return sizes;
}

inline uint8_t* min_size_per_annot_cplx_sub(uint8_t* annot_cplx, int nb_cplx, int size_cplx, int pos_start, int pos_end){

    if (nb_cplx > 0) ASSERT_NULL_PTR(annot_cplx, "min_size_per_annot_cplx()")
    if (pos_start > pos_end) ERROR("min_size_per_annot_cplx_sub(): pos_start > pos_end")

    int i, z, base;
    int size_line = SIZE_BYTE_CPLX_N + size_cplx;

    uint16_t pos = 0;

    uint8_t* sizes = calloc(pos_end - pos_start + 1, sizeof(uint8_t));
    ASSERT_NULL_PTR(sizes, "min_size_per_annot_cplx_sub()")

    for (i = 0; i < nb_cplx; i++){

        z = size_cplx - 1;
        base = i * size_line;
        pos += ((((uint16_t)annot_cplx[base]) << SIZE_CELL) | ((uint16_t)annot_cplx[base + 1]));

        if (pos >= pos_start){

            if (pos >= pos_end) goto OUT_LOOP;

            base += SIZE_BYTE_CPLX_N;

            while ((z >= 0) && (annot_cplx[base + z] == 0)) z--;
            sizes[i-pos_start] = z+1;
        }
    }

    OUT_LOOP: return sizes;
}

inline int max_size_annot_cplx_sub(uint8_t* annot_cplx, int nb_cplx, int size_cplx, int pos_start, int pos_end){

    if (nb_cplx > 0) ASSERT_NULL_PTR(annot_cplx, "max_size_annot_cplx_sub()")
    if (pos_start > pos_end) ERROR("max_size_annot_cplx_sub(): pos_start > pos_end")

    int i, z, base;
    int max_size = 0;
    int size_line = SIZE_BYTE_CPLX_N + size_cplx;

    uint16_t pos = 0;

    for (i = 0; i < nb_cplx * size_line; i += size_line){

        pos += ((((uint16_t)annot_cplx[i]) << SIZE_CELL) | ((uint16_t)annot_cplx[i + 1]));

        if (pos >= pos_start){

            if (pos >= pos_end) return max_size;

            base = i + SIZE_BYTE_CPLX_N;
            z = size_cplx - 1;

            while ((z >= 0) && (annot_cplx[base + z] == 0)) z--;
            max_size = MAX(max_size, z+1);
        }
    }

    return max_size;
}

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION UC.h
*  ===================================================================================================================================
*/

void insertKmer_UC(UC* restrict uc, uint8_t* restrict kmer, int id_genome, int size_kmer, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION annotation.h
*  ===================================================================================================================================
*/

void shift_extended_annot(UC* uc, int size_substring, int nb_substring, int pos_insert);

void insert_extend_annot(UC* uc, int size_substring, int nb_substring, int pos_insert, uint8_t annot, int shift_or_not);

void delete_extend_annots(UC* uc, int size_substring, int nb_substring, int pos_sub_start, int pos_sub_end, int delete_sub, int delete_ext_sub_array, int realloc_table);

uint8_t* get_extend_annot(UC* uc, int size_substring, int nb_substring, int pos_substring);
uint8_t** get_extend_annots(UC* uc, int size_substring, int nb_substring, int pos_substring_begin, int pos_substring_end);

uint8_t* realloc_annotation(UC* uc, int size_substring, int nb_substring, uint8_t new_size_annotation, int new_insertion, int pos_insert_extend);
void recopy_back_annot_extend(UC* uc, int size_substring, int nb_substring);
void create_annot_extended(UC* uc, int size_substring, int nb_substring);

int is_genome_present(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot, int size_annot, uint8_t* annot_sup,
                      int size_annot_sup, uint16_t id_genome);

void compute_best_mode(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot, int size_annot, uint8_t* annot_sup,
                       int size_annot_sup, uint16_t id_genome2insert);

void modify_mode_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup,
                            uint16_t id_genome2insert);

annotation_array_elem* intersection_annotations(uint8_t* annot1, int size_annot1, uint8_t* annot_sup1, int size_annot_sup1, uint8_t* annot2,
                            int size_annot2, uint8_t* annot_sup2, int size_annot_sup2, int id_genome_max, annotation_array_elem* annot_sorted);

void printAnnotation_CSV(FILE* file_output, uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup,
                         int id_genome_max, annotation_array_elem* annot_sorted);

annotation_array_elem* sort_annotations(Pvoid_t* PJArray, int* size_array);

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION annotation_special_nodes.h
*  ===================================================================================================================================
*/

void shift_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_insert);
void insert_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_insert, uint8_t* annot, int size_annot, int shift_or_not);
void delete_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_sub_start, int pos_sub_end, int delete_sub, int delete_ext_sub_array, int realloc_table);
uint8_t* get_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_substring);
uint8_t** get_annots_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_substring_begin, int pos_substring_end);
void create_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring);
void increase_size_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int new_size, int to_realloc);
void decrease_size_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int new_size);

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION UC_annotation.h
*  ===================================================================================================================================
*/

int get_annotation(UC* uc, uint8_t** annot, uint8_t** annot_ext, uint8_t** annot_cplx, int* size_annot, int* size_annot_cplx, int size_substring,
                   int nb_substring, int position);

#endif
