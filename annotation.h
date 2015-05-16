#ifndef DEF_ANNOTATION
#define DEF_ANNOTATION

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <string.h>
#include <math.h>

#include "useful_macros.h"
#include "popcnt.h"
#include "UC_annotation.h"
#include "log2.h"

/* ===================================================================================================================================
*  INLINE FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

inline int modify_annot_bis(uint8_t** current_annot, uint8_t* annot_sup, int* it_annot, int* size_current_annot, uint16_t id_genome, uint8_t flag, uint8_t flag_ext){

    ASSERT_NULL_PTR(current_annot, "modify_annot_bis()")
    ASSERT_NULL_PTR(it_annot, "modify_annot_bis()")

    if (id_genome < MASK_POWER_8[6]){

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

        (*current_annot)[*it_annot] = ((id_genome & 0xffc0) >> 4) | flag;
        (*it_annot)++;

        if (*it_annot >= *size_current_annot){
            *it_annot = 0;
            *size_current_annot = 1;
            (*current_annot) = annot_sup;
        }

        (*current_annot)[*it_annot] = ((id_genome & 0x3f) << 2) | flag_ext;
        (*it_annot)++;

        if (*it_annot >= *size_current_annot){
            *it_annot = 0;
            *size_current_annot = 1;
            (*current_annot) = annot_sup;
        }

        return 2;
    }
}

inline uint8_t* min_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot){

    if (nb_substrings == 0) return 0;
    else ASSERT_NULL_PTR(annot, "min_size_per_sub()")

    uint8_t* sizes = calloc(nb_substrings, sizeof(uint8_t));
    ASSERT_NULL_PTR(sizes, "min_size_per_sub()")

    if (size_annot == 0) return sizes;

    int i, z, base;
    int size_line = size_substring+size_annot;

    for (i=0; i<nb_substrings; i++){

        z = size_annot-1;
        base = i*size_line+size_substring;

        while ((z >= 0) && (annot[base+z] == 0)) z--;
        sizes[i] = z+1;
    }

    return sizes;
}

inline int max_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot){

    if (nb_substrings == 0) return 0;
    else ASSERT_NULL_PTR(annot, "max_size_per_sub()")

    if (size_annot == 0) return 0;

    int i, z;
    int size_line = size_substring+size_annot;
    int max_size = -1;

    for (i=size_substring; i<nb_substrings * size_line; i+=size_line){

        z = size_annot-1;
        while ((z >= 0) && (annot[i+z] == 0)) z--;

        max_size = MAX(max_size,z+1);
    }

    return max_size;
}

inline int size_annot_sub(uint8_t* annot, int size_substring, int size_annot){

    if (annot == NULL) return 0;
    if (size_annot == 0) return 0;

    int z = size_annot-1;
    while ((z >= 0) && (annot[size_substring+z] == 0)) z--;

    return z+1;
}

inline uint8_t* extract_from_annotation_array_elem(annotation_array_elem* annot_sorted, uint32_t position, int* size_annot){

    ASSERT_NULL_PTR(annot_sorted, "extract_from_annotation_array_elem()")
    ASSERT_NULL_PTR(size_annot, "extract_from_annotation_array_elem()")

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

inline void free_annotation_array_elem(annotation_array_elem* annot_sorted, int size_array){
    if (annot_sorted != NULL){
        int i = 0;
        for (i=0; i<size_array; i++){
            if (annot_sorted[i].annot_array != NULL) free(annot_sorted[i].annot_array);
        }
        free(annot_sorted);
    }
}

#endif
