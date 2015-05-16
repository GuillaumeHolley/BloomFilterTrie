#include "marking.h"

void create_marking_Node_4states(Node* n, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(n, "create_marking_Node_4states()")

    if (n->CC_array != NULL){

        int i = -1;

        do {
            i++;
            create_marking_CC_4states(&(((CC*)n->CC_array)[i]), size_kmer, func_on_types);
        }
        while ((((CC*)n->CC_array)[i].type & 0x1) == 0);
    }

    if (n->UC_array.suffixes != NULL) create_marking_UC_4states(&(n->UC_array), size_kmer, func_on_types);

    return;
}

void create_marking_CC_4states(CC* cc, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(cc,"create_marking_CC_4states()")

    UC* uc;

    int tot;
    int sup_bytes;
    int size_substring;

    int i = 0;
    int nb_skp = CEIL(cc->nb_elem, NB_CHILDREN_PER_SKP);

    if (size_kmer == SIZE_SEED){

        sup_bytes = CEIL(NB_CHILDREN_PER_SKP * 2, SIZE_CELL);

        for (i=0; i<nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            tot = NB_CHILDREN_PER_SKP * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT +
                    uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

            uc->suffixes = realloc(uc->suffixes, (tot + sup_bytes) * sizeof(uint8_t));
            ASSERT_NULL_PTR(uc->suffixes, "create_marking_CC_4states")

            memset(&(uc->suffixes[tot]), 0, sup_bytes * sizeof(uint8_t));
        }
    }
    else{

        size_substring = func_on_types[(size_kmer/SIZE_SEED)-1].size_kmer_in_bytes_minus_1;

        for (i=0; i<nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            tot = uc->nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT +
                    uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

            sup_bytes = CEIL(uc->nb_children * 2, SIZE_CELL);

            uc->suffixes = realloc(uc->suffixes, (tot + sup_bytes) * sizeof(uint8_t));
            ASSERT_NULL_PTR(uc->suffixes, "create_marking_CC_4states")

            memset(&(uc->suffixes[tot]), 0, sup_bytes * sizeof(uint8_t));
        }
    }

    for (i=0; i<cc->nb_Node_children; i++) create_marking_Node_4states(&(cc->children_Node_container[i]), size_kmer-SIZE_SEED, func_on_types);

    return;
}

void create_marking_UC_4states(UC* uc, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(uc, "create_marking_UC_4states()")

    int nb_elt = uc->nb_children >> 1;

    int tot = nb_elt * (func_on_types[(size_kmer/SIZE_SEED)-1].size_kmer_in_bytes + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT +
            uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

    int sup_bytes = CEIL(nb_elt * 2, SIZE_CELL);

    uc->suffixes = realloc(uc->suffixes, (tot + sup_bytes) * sizeof(uint8_t));
    ASSERT_NULL_PTR(uc->suffixes, "create_marking_UC_4states")

    memset(&(uc->suffixes[tot]), 0, sup_bytes * sizeof(uint8_t));

    return;
}

void delete_marking_Node_4states(Node* n, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(n, "delete_marking_Node_4states()")

    if (n->CC_array != NULL){

        int i = -1;

        do {
            i++;
            delete_marking_CC_4states(&(((CC*)n->CC_array)[i]), size_kmer, func_on_types);
        }
        while ((((CC*)n->CC_array)[i].type & 0x1) == 0);
    }

    if (n->UC_array.suffixes != NULL)
        delete_marking_UC_4states(&(n->UC_array), func_on_types[(size_kmer/SIZE_SEED)-1].size_kmer_in_bytes, n->UC_array.nb_children >> 1);

    return;
}

void delete_marking_CC_4states(CC* cc, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(cc,"delete_marking_CC_4states()")

    UC* uc;

    int i = 0;
    int nb_skp = CEIL(cc->nb_elem, NB_CHILDREN_PER_SKP);

    if (size_kmer == SIZE_SEED){

        for (i=0; i<nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            uc->suffixes = realloc(uc->suffixes, (NB_CHILDREN_PER_SKP * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT +
                    uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
            ASSERT_NULL_PTR(uc->suffixes, "delete_marking_CC_4states")
        }
    }
    else{

        int size_substring = func_on_types[(size_kmer/SIZE_SEED)-1].size_kmer_in_bytes_minus_1;

        for (i=0; i<nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            uc->suffixes = realloc(uc->suffixes, (uc->nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT +
                    uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
            ASSERT_NULL_PTR(uc->suffixes, "delete_marking_CC_4states")
        }
    }

    for (i=0; i<cc->nb_Node_children; i++) delete_marking_Node_4states(&(cc->children_Node_container[i]), size_kmer-SIZE_SEED, func_on_types);

    return;
}

void delete_marking_UC_4states(UC* uc, int size_substring, int nb_children){

    ASSERT_NULL_PTR(uc, "delete_marking_UC_4states()")
    ASSERT_NULL_PTR(uc->suffixes, "delete_marking_UC_4states()")

    uc->suffixes = realloc(uc->suffixes, (nb_children * (size_substring + uc->size_annot) +
                uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
    ASSERT_NULL_PTR(uc->suffixes, "delete_marking_UC_4states")

    return;
}

void mark_UC_4states(UC* uc, int size_substring, int nb_children, int position, uint8_t flag){

    ASSERT_NULL_PTR(uc, "mark_UC_4states()")
    ASSERT_NULL_PTR(uc->suffixes, "mark_UC_4states()")

    uc->suffixes[nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT +
                uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes) + position/4] |= flag << ((position % 4) * 2);

    return;
}

uint8_t get_mark_UC_4states(UC* uc, int size_substring, int nb_children, int position){

    ASSERT_NULL_PTR(uc, "get_mark_UC_4states()")
    ASSERT_NULL_PTR(uc->suffixes, "get_mark_UC_4states()")

    return (uc->suffixes[nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes) + position/4] >> ((position % 4) * 2)) & 0x3;
}

uint8_t get_mark_UC_4states_bis(UC* uc, int position, int nb_bytes_before_marking){

    ASSERT_NULL_PTR(uc, "get_mark_UC_4states_bis()")
    ASSERT_NULL_PTR(uc->suffixes, "get_mark_UC_4states_bis()")

    return (uc->suffixes[nb_bytes_before_marking + position/4] >> ((position % 4) * 2)) & 0x3;
}

int count_flag_mark_UC_4states(UC* uc, int size_substring, int nb_children){

    ASSERT_NULL_PTR(uc, "count_flag_mark_UC_4states()")
    ASSERT_NULL_PTR(uc->suffixes, "count_flag_mark_UC_4states()")

    int i = 0;
    int count = 0;

    uint8_t* marking = &(uc->suffixes[nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)]);

    for (i=0; i < CEIL(nb_children * 2, SIZE_CELL); i++) count += popcnt_8(marking[i] & MASK_COUNT_FLAG_1);

    return count;
}
