#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#include "useful_macros.h"
//#include "UC_annotation.h"
#include "annotation.h"

typedef struct{
    uint8_t* suffixes;
    UC_SIZE_ANNOT_T size_annot;
    UC_SIZE_ANNOT_CPLX_T size_annot_cplx_nodes;
    uint16_t nb_extended_annot;
    uint16_t nb_cplx_nodes;
    uint16_t nb_children;
} __attribute__ ((__packed__)) UC;

inline UC* createUC();
inline void initiateUC(UC* uc);
inline void initializeUC(UC* uc);
inline void resetUC(UC* uc);
inline void freeUC(UC* uc);

inline UC_SIZE_ANNOT_CPLX_T * min_size_per_annot_cplx(UC* uc, int nb_substrings, int size_substring);
inline UC_SIZE_ANNOT_CPLX_T * min_size_per_annot_cplx_sub(UC* uc, int nb_substrings, int size_substring, int pos_start, int pos_end);
inline int max_size_annot_cplx_sub(uint8_t* annot_cplx, int nb_cplx, int size_cplx, int pos_start, int pos_end);

// Pure UC inline functions

inline UC* createUC(){

    UC* uc = calloc(1,sizeof(UC));
    ASSERT_NULL_PTR(uc,"createUC()")

    uc->suffixes = NULL;

    return uc;
}

inline void initiateUC(UC* uc){

    ASSERT_NULL_PTR(uc,"initiateUC()")

    uc->nb_children &= 1;
    uc->suffixes = NULL;

    return;
}

inline void initializeUC(UC* uc){

    ASSERT_NULL_PTR(uc,"initializeUC()")

    uc->size_annot = 0;
    uc->size_annot_cplx_nodes = 0;
    uc->nb_extended_annot = 0;
    uc->nb_cplx_nodes = 0;
    uc->nb_children = 0;
    uc->suffixes = NULL;

    return;
}

inline void resetUC(UC* uc){

    ASSERT_NULL_PTR(uc,"resetUC()")

    uc->size_annot = 0;
    uc->size_annot_cplx_nodes = 0;
    uc->nb_extended_annot = 0;
    uc->nb_cplx_nodes = 0;

    uc->nb_children &= 1;

    uc->suffixes = NULL;

    return;
}

inline void freeUC(UC* uc){
    free(uc->suffixes);
    free(uc);
    return;
}

// UC complexe annotation management inline functions

inline UC_SIZE_ANNOT_CPLX_T * min_size_per_annot_cplx(UC* uc, int nb_substrings, int size_substring){

    if (nb_substrings > 0){
        ASSERT_NULL_PTR(uc, "min_size_per_annot_cplx()")
        ASSERT_NULL_PTR(uc->suffixes, "min_size_per_annot_cplx()")
    }

    int i, z, base;
    int pos = 0;
    int size_line = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;

    UC_SIZE_ANNOT_CPLX_T *sizes = calloc(nb_substrings, sizeof( UC_SIZE_ANNOT_CPLX_T ));
    ASSERT_NULL_PTR(sizes, "min_size_per_annot_cplx()")

    uint8_t* annot_cplx = &(uc->suffixes[nb_substrings * (size_substring + uc->size_annot)
                            + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]);

    for (i = 0; i < uc->nb_cplx_nodes; i++){

        z = uc->size_annot_cplx_nodes - 1;
        base = i * size_line;
        pos += ((((uint16_t)annot_cplx[base]) << SIZE_BITS_UINT_8T) | ((uint16_t)annot_cplx[base + 1]));
        base += SIZE_BYTE_CPLX_N;

        while ((z >= 0) && (annot_cplx[base + z] == 0)) z--;
        sizes[i] = z+1;
    }

    return sizes;
}

inline UC_SIZE_ANNOT_CPLX_T * min_size_per_annot_cplx_sub(UC* uc, int nb_substrings, int size_substring, int pos_start, int pos_end){

    if (nb_substrings > 0){
        ASSERT_NULL_PTR(uc, "min_size_per_annot_cplx_sub()")
        ASSERT_NULL_PTR(uc->suffixes, "min_size_per_annot_cplx_sub()")
    }

    if (pos_start > pos_end) ERROR("min_size_per_annot_cplx_sub(): pos_start > pos_end")

    int i, z, base;
    int size_line = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;

    uint16_t pos = 0;

    UC_SIZE_ANNOT_CPLX_T *sizes = calloc(pos_end - pos_start + 1, sizeof( UC_SIZE_ANNOT_CPLX_T ));
    ASSERT_NULL_PTR(sizes, "min_size_per_annot_cplx_sub()")

    uint8_t* annot_cplx = &(uc->suffixes[nb_substrings * (size_substring + uc->size_annot)
                            + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]);

    for (i = 0; i < uc->nb_cplx_nodes; i++){

        z = uc->size_annot_cplx_nodes - 1;
        base = i * size_line;
        pos += ((((uint16_t)annot_cplx[base]) << SIZE_BITS_UINT_8T) | ((uint16_t)annot_cplx[base + 1]));

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

        pos += ((((uint16_t)annot_cplx[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)annot_cplx[i + 1]));

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

// Pure UC functions

void insertKmer_UC(UC*  uc, uint8_t*  kmer, uint32_t id_genome,
                   int size_id_genome, int size_kmer, int pos_insertion, annotation_inform* ann_inf,
                   annotation_array_elem* annot_sorted);

int binary_search_UC(const UC*, int pos_start, int pos_end, const uint8_t* suf, int size_suf_byte,
                     uint8_t mask_for_last_byte);

int binary_search_UC_array(const uint8_t* uc_array, int size_annot, int pos_start, int pos_end,
                           const uint8_t* suf, int size_suf_byte, uint8_t mask_for_last_byte);

                           int get_annot(UC* uc, uint8_t** annot, uint8_t** annot_ext, uint8_t** annot_cplx, int* size_annot,
                   int* size_annot_cplx, int size_substring, int nb_substring, int position);

// UC annotation management functions

int get_annot(UC* uc, uint8_t** annot, uint8_t** annot_ext, uint8_t** annot_cplx,
                   int* size_annot, int* size_annot_cplx, int size_substring, int nb_substring, int position);

void get_annots(UC* uc, uint8_t*** annots, uint8_t*** annots_ext, uint8_t*** annots_cplx,
                   int** size_annots, int** size_annots_cplx, int size_substring, int nb_substring,
                   int position_start, int position_end);

// UC extended annotation management functions

void shift_extended_annot(UC* uc, int size_substring, int nb_substring, int pos_insert);

void insert_extend_annot(UC* uc, int size_substring, int nb_substring, int pos_insert,
                         uint8_t annot, int shift_or_not);

void delete_extend_annots(UC* uc, int size_substring, int nb_substring, int pos_sub_start,
                          int pos_sub_end, int delete_sub, int delete_ext_sub_array,
                          int realloc_table);

uint8_t* get_extend_annot(UC* uc, int size_substring, int nb_substring, int pos_substring);
uint8_t** get_extend_annots(UC* uc, int size_substring, int nb_substring, int pos_substring_begin,
                            int pos_substring_end);

uint8_t* realloc_annotation(UC* uc, int size_substring, int nb_substring, int new_size_annotation,
                            int new_insertion, int pos_insert_extend);
void recopy_back_annot_extend(UC* uc, int size_substring, int nb_substring);
void create_annot_extended(UC* uc, int size_substring, int nb_substring);

// UC complex annotation management functions

void shift_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_insert);
void insert_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_insert, uint8_t* annot,
                             int size_annot, int shift_or_not);
void delete_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_sub_start, int pos_sub_end,
                             int delete_sub, int delete_ext_sub_array, int realloc_table);
uint8_t* get_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_substring);
uint8_t** get_annots_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_substring_begin,
                                int pos_substring_end);
//void create_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring);
void increase_size_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int new_size, int to_realloc);
void decrease_size_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int new_size);

void create_annot_cplx_nodes_marked(UC* uc, int size_substring, int nb_substring); // For compression test
