#pragma once

#define _XOPEN_SOURCE 500

#define UC_SIZE_ANNOT_T int32_t
#define UC_SIZE_ANNOT_CPLX_T UC_SIZE_ANNOT_T

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>
#include <limits.h>
#include <string.h>
#include <math.h>

#include <Judy.h>

#include "default_param.h"
#include "useful_macros.h"
#include "log2.h"
#include "quicksort.h"

extern const uint8_t MASK_POWER_8[8];

typedef struct {
    //Represent the genomes IDs stored in the annotation, doesn't matter if they are present or not
    uint32_t* id_stored;
    uint32_t* size_id_stored;

    uint8_t* annotation;

    int nb_id_stored;

    //Last ID added to the current annotation before annotation
    int last_added;

    //Current mode of the annotation, new optimal mode and the required size
    int min_size;
    int min_mode;
    int current_mode;
    int size_annot;
    int comp_annot;

    uint8_t disabled_flags;

} annotation_inform;

typedef struct {
    int64_t last_index;
    uint8_t* annot_array;
    int size_annot;
} annotation_array_elem;

// Function declarations

int is_genome_present(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot,
                      int size_annot, uint8_t* annot_sup, int size_annot_sup, uint32_t id_genome);

int is_genome_present_from_end_annot(annotation_inform* ann_inf, annotation_array_elem* annot_sorted,
                                     uint8_t* annot, int size_annot, uint8_t* annot_sup,
                                     int size_annot_sup, uint32_t id_genome);

int get_last_genome_inserted(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot,
                      int size_annot, uint8_t* annot_sup, int size_annot_sup, uint32_t* last_id_genome);

int comp_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot, uint8_t* annot_sup,
                     int size_annot_sup);
int decomp_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot, uint8_t* annot_sup,
                       int size_annot_sup, int get_sizes);

void compute_best_mode(annotation_inform* ann_inf, annotation_array_elem* annot_sorted,
                       uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup,
                       uint32_t id_genome2insert, int size_id_genome);

void modify_mode_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot, uint8_t* annot_sup,
                            int size_annot_sup, uint32_t id_genome2insert, int size_id_genome);

annotation_array_elem* cmp_annots(uint8_t* annot1, int size_annot1, uint8_t* annot_sup1, int size_annot_sup1,
                                  uint8_t* annot2, int size_annot2,uint8_t* annot_sup2, int size_annot_sup2,
                                  uint32_t id_genome_max, uint8_t (*f)(const uint8_t, const uint8_t),
                                  annotation_array_elem* annot_sorted);

void printAnnotation_CSV(FILE* file_output, uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup,
                         uint32_t id_genome_max, annotation_array_elem* annot_sorted);

annotation_array_elem* sort_annotations(Pvoid_t* PJArray, int* size_array, uint32_t longest_annot);

void sort_annotations2(char* filename_annot_array_elem, Pvoid_t* JArray_annot, annotation_array_elem** root_comp_set_colors,
                       int* length_root_comp_set_colors, uint32_t longest_annot);

void sort_annotations3(Pvoid_t* JArray_annot, uint32_t longest_annot);

void replace_annots_comp(annotation_array_elem* comp_colors, Pvoid_t* JArray_annot, char* filename_new_comp_colors, uint32_t longest_annot);

void write_partial_comp_set_colors(char* filename_annot_array_elem, Pvoid_t* JArray_annot, uint32_t longest_annot);

void get_id_genomes_from_annot(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot,
                               int size_annot, uint8_t* annot_sup, int size_annot_sup);

int get_count_id_genomes_from_annot(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot,
                                    int size_annot, uint8_t* annot_sup, int size_annot_sup);

//Inline function declarations

inline annotation_inform* create_annotation_inform(int nb_id_genomes, bool disable_flag_0);
inline void reinit_annotation_inform(annotation_inform* ann_inf);
inline void free_annotation_inform(annotation_inform* ann_inf);
inline int size_annot_sub(uint8_t* annot, int size_substring, int size_annot);
inline int modify_annot_bis(uint8_t** current_annot, uint8_t* annot_sup, int* it_annot, int* size_current_annot,
                            uint32_t id_genome, int size_id_genome, uint8_t flag, uint8_t flag_ext);
inline UC_SIZE_ANNOT_T *min_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot);
inline int max_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot);
inline uint8_t* extract_from_annotation_array_elem(annotation_array_elem* annot_sorted, uint32_t position,
                                                   int* size_annot);
inline int getSize_from_annotation_array_elem(annotation_array_elem* annot_sorted, uint32_t position);
inline double getTotalSize_annotation_array_elem(annotation_array_elem* annot_sorted, int size_array);
inline int getMaxSize_annotation_array_elem(annotation_array_elem* annot_sorted);
inline void free_annotation_array_elem(annotation_array_elem** annot_sorted, int* size_array);

// Inline functions

inline annotation_inform* create_annotation_inform(int nb_id_genomes, bool disable_flag_0){

    annotation_inform* ann_inf = calloc(1, sizeof(annotation_inform));
    ASSERT_NULL_PTR(ann_inf,"create_annotation_inform()")

    if (nb_id_genomes <= 0){
        ann_inf->id_stored = malloc(NB_MAX_ID_GENOMES * sizeof(uint32_t));
        ASSERT_NULL_PTR(ann_inf->id_stored, "create_annotation_inform()");

        ann_inf->size_id_stored = malloc(NB_MAX_ID_GENOMES * sizeof(uint32_t));
        ASSERT_NULL_PTR(ann_inf->size_id_stored, "create_annotation_inform()");

        ann_inf->annotation = calloc(SIZE_MAX_BYTE_ANNOT, sizeof(uint8_t));
        ASSERT_NULL_PTR(ann_inf->annotation, "create_annotation_inform()");
    }
    else{
        ann_inf->id_stored = malloc(nb_id_genomes * sizeof(uint32_t));
        ASSERT_NULL_PTR(ann_inf->id_stored, "create_annotation_inform()");

        ann_inf->size_id_stored = malloc(nb_id_genomes * sizeof(uint32_t));
        ASSERT_NULL_PTR(ann_inf->size_id_stored, "create_annotation_inform()");

        ann_inf->annotation = calloc(CEIL(nb_id_genomes + 2, SIZE_BITS_UINT_8T), sizeof(uint8_t));
        ASSERT_NULL_PTR(ann_inf->annotation, "create_annotation_inform()");
    }

    if (disable_flag_0) ann_inf->disabled_flags = 1;

    return ann_inf;
}

inline void reinit_annotation_inform(annotation_inform* ann_inf){

    memset(ann_inf->annotation, 0, ann_inf->size_annot * sizeof(uint8_t));

    ann_inf->nb_id_stored = 0;
    ann_inf->last_added = 0;
    ann_inf->current_mode = 0;
    ann_inf->min_mode = 0;
    ann_inf->min_size = 0;
    ann_inf->size_annot = 0;
    ann_inf->comp_annot = 0;

    return;
}

inline void free_annotation_inform(annotation_inform* ann_inf){

    ASSERT_NULL_PTR(ann_inf,"free_annotation_inform()");

    free(ann_inf->id_stored);
    free(ann_inf->size_id_stored);
    free(ann_inf->annotation);
    free(ann_inf);

    return;
}

inline int size_annot_sub(uint8_t* annot, int size_substring, int size_annot){

    if ((annot == NULL) || (size_annot == 0)) return 0;

    uint8_t* z = annot + size_substring + size_annot - 1;
    while ((z >= annot) && (*z == 0)) z--;

    return z - annot + 1;
}

inline int modify_annot_bis(uint8_t** current_annot, uint8_t* annot_sup, int* it_annot, int* size_current_annot,
                            uint32_t id_genome, int size_id_genome, uint8_t flag, uint8_t flag_ext){

    ASSERT_NULL_PTR(current_annot, "modify_annot_bis()")
    ASSERT_NULL_PTR(it_annot, "modify_annot_bis()")

    if (id_genome < 0x40){

        (*current_annot)[*it_annot] = (id_genome << 2) | flag;
        (*it_annot)++;

        if (*it_annot >= *size_current_annot){
            *it_annot = 0;
            *size_current_annot = 1;
            (*current_annot) = annot_sup;
        }

        return 1;
    }
    else{

        if (size_id_genome <= 0) size_id_genome = get_nb_bytes_power2_annot(id_genome);

        int size_id_genome_cpy = size_id_genome;
        size_id_genome = (size_id_genome-1) * 6 - 2;

        uint8_t* curr_ann_tmp = *current_annot;
        int it_tmp = *it_annot;
        int size_tmp = *size_current_annot;

        curr_ann_tmp[it_tmp] = ((id_genome >> size_id_genome) & 0xfc) | flag;
        it_tmp++;
        size_id_genome -= 6;

        while (size_id_genome > 0){
            curr_ann_tmp[it_tmp] = ((id_genome >> size_id_genome) & 0xfc) | flag_ext;
            it_tmp++;
            size_id_genome -= 6;
        }

        if (it_tmp >= size_tmp){
            it_tmp = 0;
            size_tmp = 1;
            curr_ann_tmp = annot_sup;
        }

        if (size_id_genome <= 0){

            curr_ann_tmp[it_tmp] = (id_genome << 2) | flag_ext;
            it_tmp++;

            if (it_tmp >= size_tmp){
                it_tmp = 0;
                size_tmp = 1;
                curr_ann_tmp = annot_sup;
            }
        }

        *current_annot = curr_ann_tmp;
        *it_annot = it_tmp;
        *size_current_annot = size_tmp;

        return size_id_genome_cpy;
    }
}

inline UC_SIZE_ANNOT_T *min_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot){

    if (nb_substrings == 0) return 0;
    else ASSERT_NULL_PTR(annot, "min_size_per_sub()")

    UC_SIZE_ANNOT_T *sizes = calloc(nb_substrings, sizeof( UC_SIZE_ANNOT_T ));
    ASSERT_NULL_PTR(sizes, "min_size_per_sub()")

    if (size_annot == 0) return sizes;

    int i, k;
    int size_line = size_substring+size_annot;

    uint8_t* z;

    for (i=size_substring, k=0; i < nb_substrings * size_line; i += size_line, k++){

        z = annot + i + size_annot - 1;
        while ((z >= annot + i) && (*z == 0)) z--;
        sizes[k] = z - annot - i + 1;
    }

    return sizes;
}

inline int max_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot){

    if (nb_substrings == 0) return 0;
    else ASSERT_NULL_PTR(annot, "max_size_per_sub()")

    if (size_annot == 0) return 0;

    int size_line = size_substring + size_annot;
    int i = size_substring;
    int max_size = -1;

    uint8_t* z;

    for (; i < nb_substrings * size_line; i += size_line){

        z = annot + i + size_annot - 1;
        while ((z >= annot + i) && (*z == 0)) z--;

        max_size = MAX(max_size, z - annot - i + 1);
        if (max_size == size_annot) return max_size;
    }

    return max_size;
}

inline uint8_t* extract_from_annotation_array_elem(annotation_array_elem* annot_sorted, uint32_t position,
                                                   int* size_annot){

    ASSERT_NULL_PTR(annot_sorted, "extract_from_annotation_array_elem() 1")
    ASSERT_NULL_PTR(size_annot, "extract_from_annotation_array_elem() 2")

    int it_annot_sorted = 0;

    while(position > annot_sorted[it_annot_sorted].last_index) it_annot_sorted++;

    *size_annot = annot_sorted[it_annot_sorted].size_annot;

    if (it_annot_sorted == 0) return &(annot_sorted[0].annot_array[position * (*size_annot)]);
    else return &(annot_sorted[it_annot_sorted].annot_array[(position-annot_sorted[it_annot_sorted-1].last_index-1) * (*size_annot)]);
}

inline int getSize_from_annotation_array_elem(annotation_array_elem* annot_sorted, uint32_t position){

    ASSERT_NULL_PTR(annot_sorted, "getSize_from_annotation_array_elem()")

    int it_annot_sorted = 0;

    while(position > annot_sorted[it_annot_sorted].last_index) it_annot_sorted++;

    return annot_sorted[it_annot_sorted].size_annot;
}

inline double getTotalSize_annotation_array_elem(annotation_array_elem* annot_sorted, int size_array){

    double size_annot_array_elem = 0;

    if (annot_sorted != NULL){

        int64_t prev_pos = 0;

        for (int i=0; i<size_array; i++){
            if (annot_sorted[i].annot_array != NULL){
                size_annot_array_elem += (annot_sorted[i].last_index - prev_pos + 1) * annot_sorted[i].size_annot;
                prev_pos = annot_sorted[i].last_index;
            }
        }
    }

    return size_annot_array_elem;
}

inline int getMaxSize_annotation_array_elem(annotation_array_elem* annot_sorted){
    if (annot_sorted != NULL) return annot_sorted[0].size_annot;
    return 0;
}

inline void free_annotation_array_elem(annotation_array_elem** annot_sorted, int* size_array){

    if (*annot_sorted != NULL){

        for (int i = 0; i < *size_array; i++){
            if ((*annot_sorted)[i].annot_array != NULL) free((*annot_sorted)[i].annot_array);
        }

        free(*annot_sorted);
        *annot_sorted = NULL;

        *size_array = 0;
    }

    return;
}
