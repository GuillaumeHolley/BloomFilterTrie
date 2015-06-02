#include "./../lib/UC_annotation.h"

extern int size_annot_sub(uint8_t* annot, int size_substring, int size_annot);

int get_annotation(UC* uc, uint8_t** annot, uint8_t** annot_ext, uint8_t** annot_cplx,
                   int* size_annot, int* size_annot_cplx, int size_substring, int nb_substring, int position){

    if (position >= nb_substring) ERROR("get_annotation(): position >= nb_substring")

    ASSERT_NULL_PTR(uc,"get_annotation()")
    ASSERT_NULL_PTR(uc->suffixes,"get_annotation()")

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
