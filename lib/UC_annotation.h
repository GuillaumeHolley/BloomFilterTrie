#ifndef DEF_UC_ANNOTATION
#define DEF_UC_ANNOTATION

#define UC_SIZE_ANNOT_T int32_t
#define UC_SIZE_ANNOT_CPLX_T UC_SIZE_ANNOT_T

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#include <Judy.h>

#include "./../lib/default_param.h"
#include "./../lib/useful_macros.h"

/* ===================================================================================================================================
*  MACRO DECLARATION
*  ===================================================================================================================================
*/

#define FREE_5PTRS(_p1,_p2,_p3,_p4,_p5) \
     if ( (_p1) != NULL){               \
         free( (_p1) );                 \
         (_p1) = NULL;                  \
     }                                  \
     if ( (_p2) != NULL){               \
         free( (_p2) );                 \
         (_p2) = NULL;                  \
     }                                  \
     if ( (_p3) != NULL){               \
         free( (_p3) );                 \
         (_p3) = NULL;                  \
     }                                  \
     if ( (_p4) != NULL){               \
         free( (_p4) );                 \
         (_p4) = NULL;                  \
     }                                  \
     if ( (_p5) != NULL){               \
         free( (_p5) );                 \
         (_p5) = NULL;                  \
     }

/* ===================================================================================================================================
*  STRUCTURES DECLARATION
*  ===================================================================================================================================
*/

/*static const char* presence_genomes[16] = {"0,0,0,0,", "1,0,0,0,", "0,1,0,0,", "1,1,0,0,",
                                    "0,0,1,0,", "1,0,1,0,", "0,1,1,0,", "1,1,1,0,",
                                    "0,0,0,1,", "1,0,0,1,", "0,1,0,1,", "1,1,0,1,",
                                    "0,0,1,1,", "1,0,1,1,", "0,1,1,1,", "1,1,1,1,"};*/

/*typedef struct {
    //Represent the genomes IDs stored in the annotation, doesn't matter if they are present or not
    uint32_t id_stored[NB_MAX_ID_GENOMES];
    uint32_t size_id_stored[NB_MAX_ID_GENOMES];
    int nb_id_stored;

    //Last ID added to the current annotation before annotation
    int last_added;

    //Current mode of the annation, new optimal mode and the required size
    int min_size;
    int min_mode;
    int current_mode;

    uint8_t annotation[SIZE_MAX_BYTE_ANNOT];
    int size_annot;

} annotation_inform;*/

typedef struct {
    //Represent the genomes IDs stored in the annotation, doesn't matter if they are present or not
    uint32_t* id_stored;
    uint32_t* size_id_stored;

    uint8_t* annotation;

    int nb_id_stored;

    //Last ID added to the current annotation before annotation
    int last_added;

    //Current mode of the annation, new optimal mode and the required size
    int min_size;
    int min_mode;
    int current_mode;
    int size_annot;
    int comp_annot;

} annotation_inform;

typedef struct {
    int64_t last_index;
    uint8_t* annot_array;
    int size_annot;
} annotation_array_elem;

typedef struct{
    uint8_t* suffixes;
    UC_SIZE_ANNOT_T size_annot;
    UC_SIZE_ANNOT_CPLX_T size_annot_cplx_nodes;
    uint16_t nb_extended_annot;
    uint16_t nb_cplx_nodes;
    uint16_t nb_children;
} __attribute__ ((__packed__)) UC;

inline annotation_inform* create_annotation_inform(int nb_id_genomes){

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

        ann_inf->annotation = calloc(CEIL(nb_id_genomes+2,SIZE_BITS_UINT_8T), sizeof(uint8_t));
        ASSERT_NULL_PTR(ann_inf->annotation, "create_annotation_inform()");
    }

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

inline UC_SIZE_ANNOT_CPLX_T *min_size_per_annot_cplx_sub(UC* uc, int nb_substrings, int size_substring, int pos_start, int pos_end){

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

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION UC.h
*  ===================================================================================================================================
*/

void insertKmer_UC(UC* restrict uc, uint8_t* restrict kmer, uint32_t id_genome,
                   int size_id_genome, int size_kmer, int pos_insertion, annotation_inform* ann_inf,
                   annotation_array_elem* annot_sorted);

int binary_search_UC(const UC*, int pos_start, int pos_end, const uint8_t* suf, int size_suf_byte,
                     uint8_t mask_for_last_byte);

int binary_search_UC_array(const uint8_t* uc_array, int size_annot, int pos_start, int pos_end,
                           const uint8_t* suf, int size_suf_byte, uint8_t mask_for_last_byte);

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION annotation.h
*  ===================================================================================================================================
*/

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

int is_genome_present(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot,
                      int size_annot, uint8_t* annot_sup, int size_annot_sup, uint32_t id_genome);

int is_genome_present_from_end_annot(annotation_inform* ann_inf, annotation_array_elem* annot_sorted,
                                     uint8_t* annot, int size_annot, uint8_t* annot_sup,
                                     int size_annot_sup, uint32_t id_genome);

int get_last_genome_inserted(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot,
                      int size_annot, uint8_t* annot_sup, int size_annot_sup, int* last_id_genome);

int comp_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot, uint8_t* annot_sup,
                     int size_annot_sup);
int decomp_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot, uint8_t* annot_sup,
                       int size_annot_sup, int get_sizes);

void compute_best_mode(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot,
                       int size_annot, uint8_t* annot_sup, int size_annot_sup, uint32_t id_genome2insert,
                       int size_id_genome);

void modify_mode_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot, uint8_t* annot_sup,
                            int size_annot_sup, uint32_t id_genome2insert, int size_id_genome);

annotation_array_elem* intersection_annotations(uint8_t* annot1, int size_annot1, uint8_t* annot_sup1,
                                                int size_annot_sup1, uint8_t* annot2, int size_annot2,
                                                uint8_t* annot_sup2, int size_annot_sup2, uint32_t id_genome_max,
                                                annotation_array_elem* annot_sorted);

void printAnnotation_CSV(FILE* file_output, uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup,
                         uint32_t id_genome_max, annotation_array_elem* annot_sorted);

annotation_array_elem* sort_annotations(Pvoid_t* PJArray, int* size_array, uint32_t longest_annot);

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION annotation_special_nodes.h
*  ===================================================================================================================================
*/

void shift_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_insert);
void insert_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_insert, uint8_t* annot,
                             int size_annot, int shift_or_not);
void delete_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_sub_start, int pos_sub_end,
                             int delete_sub, int delete_ext_sub_array, int realloc_table);
uint8_t* get_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_substring);
uint8_t** get_annots_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_substring_begin,
                                int pos_substring_end);
void create_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring);
void increase_size_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int new_size, int to_realloc);
void decrease_size_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int new_size);

void create_annot_cplx_nodes_marked(UC* uc, int size_substring, int nb_substring); // For compression test

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION UC_annotation.h
*  ===================================================================================================================================
*/

int get_annotation(UC* uc, uint8_t** annot, uint8_t** annot_ext, uint8_t** annot_cplx, int* size_annot,
                   int* size_annot_cplx, int size_substring, int nb_substring, int position);

void get_annotations(UC* uc, uint8_t*** annots, uint8_t*** annots_ext, uint8_t*** annots_cplx,
                   int** size_annots, int** size_annots_cplx, int size_substring, int nb_substring,
                   int position_start, int position_end);

#endif
