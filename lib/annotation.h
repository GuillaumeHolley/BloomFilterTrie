#ifndef DEF_ANNOTATION
#define DEF_ANNOTATION

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <string.h>
#include <math.h>

#include "./../lib/useful_macros.h"
#include "./../lib/UC_annotation.h"
#include "./../lib/log2.h"

/* ===================================================================================================================================
*  INLINE FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

const uint8_t MASK_POWER_8[8];

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

            curr_ann_tmp[it_tmp] = ((id_genome >> size_id_genome) & 0xfc) | flag_ext;
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

inline int size_annot_sub(uint8_t* annot, int size_substring, int size_annot){

    if ((annot == NULL) || (size_annot == 0)) return 0;

    uint8_t* z = annot + size_substring + size_annot - 1;
    while ((z >= annot) && (*z == 0)) z--;

    return z - annot + 1;
}

inline uint8_t* extract_from_annotation_array_elem(annotation_array_elem* annot_sorted, uint32_t position, int* size_annot){

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

inline void free_annotation_array_elem(annotation_array_elem* annot_sorted, int size_array){

    if (annot_sorted != NULL){

        int i = 0;
        for (i=0; i<size_array; i++){
            if (annot_sorted[i].annot_array != NULL) free(annot_sorted[i].annot_array);
        }

        free(annot_sorted);
    }
}

inline double getTotalSize_annotation_array_elem(annotation_array_elem* annot_sorted, int size_array){

    double size_annot_array_elem = 0;

    if (annot_sorted != NULL){

        int64_t old_pos = 0;

        for (int i=0; i<size_array; i++){
            if (annot_sorted[i].annot_array != NULL){
                size_annot_array_elem += (annot_sorted[i].last_index - old_pos + 1) * annot_sorted[i].size_annot;
                old_pos = annot_sorted[i].last_index;
            }
        }
    }

    return size_annot_array_elem;
}

inline int getMaxSize_annotation_array_elem(annotation_array_elem* annot_sorted){
    if (annot_sorted != NULL) return annot_sorted[0].size_annot;
    else return 0;
}

#endif
