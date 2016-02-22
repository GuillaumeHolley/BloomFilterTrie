#include "./../lib/UC_annotation.h"

extern int size_annot_sub(uint8_t* annot, int size_substring, int size_annot);

int get_annot(UC* uc, uint8_t** annot, uint8_t** annot_ext, uint8_t** annot_cplx,
                   int* size_annot, int* size_annot_cplx, int size_substring, int nb_substring, int position){

    if (position >= nb_substring) ERROR("get_annot(): position >= nb_substring")

    ASSERT_NULL_PTR(uc,"get_annot()")
    ASSERT_NULL_PTR(uc->suffixes,"get_annot()")

    uint8_t* ext_annot = NULL;
    uint8_t* cplx_annot = NULL;

    if (uc->size_annot == 0){ // No standard annotation

        *annot = NULL;
        *size_annot = 0;

        if ((ext_annot = get_extend_annot(uc, size_substring, nb_substring, position)) == NULL){ // No extended annotation

            *annot_ext = NULL;

            if ((cplx_annot = get_annot_cplx_nodes(uc, size_substring, nb_substring, position)) == NULL){ // No complex annotation

                *annot_cplx = NULL;
                *size_annot_cplx = 0;

                return 0;
            }
            else {
                *annot_cplx = cplx_annot;
                *size_annot_cplx = size_annot_sub(*annot_cplx, 0, uc->size_annot_cplx_nodes);
            }
        }
        else{
            *annot_ext = ext_annot;
            *size_annot_cplx = 0;
        }
    }
    else{

        *annot = &(uc->suffixes[position * (size_substring + uc->size_annot) + size_substring]);

        if ((ext_annot = get_extend_annot(uc, size_substring, nb_substring, position)) == NULL){// No extended annotation

            *annot_ext = NULL;

            if ((cplx_annot = get_annot_cplx_nodes(uc, size_substring, nb_substring, position)) == NULL){ // No complex annotation

                *annot_cplx = NULL;
                *size_annot_cplx = 0;
                *size_annot = size_annot_sub(*annot, 0, uc->size_annot);

                if (*size_annot == 0) return 0;
            }
            else {
                *annot_cplx = cplx_annot;
                *size_annot_cplx = size_annot_sub(*annot_cplx, 0, uc->size_annot_cplx_nodes);
                *size_annot = 0;
            }
        }
        else{
            *annot_ext = ext_annot;
            *annot_cplx = NULL;
            *size_annot = uc->size_annot;
            *size_annot_cplx = 0;
        }
    }

    return 1;
}

void get_annots(UC* uc, uint8_t*** annots, uint8_t*** annots_ext, uint8_t*** annots_cplx,
                   int** size_annots, int** size_annots_cplx, int size_substring, int nb_substring,
                   int position_start, int position_end){

    if (nb_substring <= 0) return;

    if ((position_start >= nb_substring) || (position_start < 0)){
        ERROR("get_annots(): position_start >= nb_substring or < 0")
    }

    if ((position_end >= nb_substring) || (position_end < 0)){
        ERROR("get_annots(): position_end >= nb_substring or < 0")
    }

    ASSERT_NULL_PTR(uc,"get_annot()")
    ASSERT_NULL_PTR(uc->suffixes,"get_annot()")

    int size_line;

    int i = 0;
    int nb_suffixes = position_end - position_start + 1;

    *annots = malloc(nb_suffixes * sizeof(uint8_t*));
    ASSERT_NULL_PTR(*annots, "get_annots()")

    *size_annots_cplx = calloc(nb_suffixes, sizeof(int));
    ASSERT_NULL_PTR(*size_annots_cplx, "get_annots()")

    *size_annots = calloc(nb_suffixes, sizeof(int));
    ASSERT_NULL_PTR(*size_annots, "get_annots()")

    *annots_ext = NULL;
    *annots_cplx = NULL;

    if (uc->size_annot == 0){

        for (i = 0; i < nb_suffixes; i++) (*annots)[i] = NULL;

        if ((*annots_ext = get_extend_annots(uc, size_substring, nb_substring,
                                             position_start, position_end)) == NULL){

            if ((*annots_cplx = get_annots_cplx_nodes(uc, size_substring, nb_substring,
                                                     position_start, position_end)) != NULL){
                free(*size_annots_cplx);
                *size_annots_cplx = min_size_per_annot_cplx_sub(uc, size_substring, nb_substring,
                                                                position_start, position_end);
            }
        }
    }
    else{

        size_line = size_substring + uc->size_annot;

        for (i = position_start; i <= position_end; i++)
            (*annots)[i - position_start] = &(uc->suffixes[i * size_line + size_substring]);

        if ((*annots_ext = get_extend_annots(uc, size_substring, nb_substring,
                                             position_start, position_end)) == NULL){

            for (i = 0; i < nb_suffixes; i++)
                (*size_annots)[i] = size_annot_sub((*annots)[i], 0, uc->size_annot);

            if ((*annots_cplx = get_annots_cplx_nodes(uc, size_substring, nb_substring,
                                                     position_start, position_end)) != NULL){

                for (i = 0; i < nb_suffixes; i++)
                    (*size_annots_cplx)[i] = size_annot_sub((*annots_cplx)[i], 0, uc->size_annot_cplx_nodes);
            }
        }
        else{
            for (i = 0; i < nb_suffixes; i++){
                if ((*annots_ext)[i] != NULL) (*size_annots)[i] = uc->size_annot;
                else (*size_annots)[i] = size_annot_sub((*annots)[i], 0, uc->size_annot);
            }
        }
    }

    return;
}
