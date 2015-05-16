#include "replaceAnnotation.h"

void load_annotation_from_Node(Node* restrict node, int size_kmer, ptrs_on_func* restrict func_on_types, Pvoid_t* PJArray, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(node,"load_annotation_from_Node()")
    ASSERT_NULL_PTR(func_on_types,"load_annotation_from_Node()")

    int i = -1;

    if ((CC*)node->CC_array != NULL){
        do {
            i++;
            load_annotation_from_CC(&(((CC*)node->CC_array)[i]), size_kmer, func_on_types, PJArray, annot_sorted);
        }
        while ((((CC*)node->CC_array)[i].type & 0x1) == 0);
    }

    if (node->UC_array.suffixes != NULL)
        load_annotation_from_UC(&(node->UC_array), func_on_types[(size_kmer/SIZE_SEED)-1].size_kmer_in_bytes, node->UC_array.nb_children >> 1, PJArray, annot_sorted);
}

void load_annotation_from_CC(CC* restrict cc, int size_kmer, ptrs_on_func* restrict func_on_types, Pvoid_t* PJArray, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(cc,"load_annotation_from_CC()")

    UC* uc;

    int i;

    int level = (size_kmer/SIZE_SEED)-1;
    int nb_skp = CEIL(cc->nb_elem, NB_CHILDREN_PER_SKP);

    if (size_kmer != SIZE_SEED){

        for (i = 0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);
            load_annotation_from_UC(uc, func_on_types[level].size_kmer_in_bytes_minus_1, uc->nb_children, PJArray, annot_sorted);
        }
    }
    else{
        for (i = 0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            if (i != nb_skp-1) load_annotation_from_UC(uc, 0, NB_CHILDREN_PER_SKP, PJArray, annot_sorted);
            else load_annotation_from_UC(uc, 0, cc->nb_elem - i * NB_CHILDREN_PER_SKP, PJArray, annot_sorted);
        }
    }

    for (i = 0; i < cc->nb_Node_children; i++)
        load_annotation_from_Node(&(cc->children_Node_container[i]), size_kmer-SIZE_SEED, func_on_types, PJArray, annot_sorted);

    return;
}

void load_annotation_from_UC(UC* uc, int size_substring, int nb_children, Pvoid_t* PJArray, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(uc, "load_annotation_from_UC()")

    if (uc->suffixes != NULL){

        PWord_t PValue;

        int i = 0, j = 0;
        int size_line = size_substring + uc->size_annot;

        int size_checksum;
        int size_tmp;
        uint32_t position;

        uint8_t** extended_annots = get_extend_annots(uc, size_substring, nb_children, 0, nb_children-1);

        uint8_t* min_sizes = min_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        uint8_t* annot = calloc(NB_CELL_GENOMES + CEIL(NB_CELL_GENOMES, SIZE_CELL), sizeof(uint8_t)); //annot is always on max NB_CELL_GENOMES bytes + checksum
        ASSERT_NULL_PTR(annot, "load_annotation_from_UC()")

        for (i=0; i<nb_children; i++){

            if ((extended_annots != NULL) && (extended_annots[i] != NULL) && (extended_annots[i][0] != 0)){
                min_sizes[i] = uc->size_annot;
                memcpy(&(annot[1]), &(uc->suffixes[i * size_line + size_substring]), min_sizes[i] * sizeof(uint8_t));
                annot[min_sizes[i] + 1] = extended_annots[i][0];
                min_sizes[i]++;
            }
            else memcpy(&(annot[1]), &(uc->suffixes[i * size_line + size_substring]), min_sizes[i] * sizeof(uint8_t));

            if ((annot[1] & 0x3) == 3){

                position = 0;
                j = 1;

                while ((j < min_sizes[i]+1) && (annot[j] != 0)){
                    if (j == 1) position = annot[1] >> 2;
                    else position |= ((uint32_t)(annot[j] >> 1)) << (6+(j-2)*7);
                    j++;
                }

                uint8_t* annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &size_tmp);
                memcpy(&(annot[1]), annot_tmp, size_tmp*sizeof(uint8_t));
                annot[0] = size_tmp;
                min_sizes[i] = size_tmp;
            }
            else annot[0] = min_sizes[i];

            size_checksum = CEIL(min_sizes[i], SIZE_CELL-1);

            memset(&(annot[annot[0] + 1]), 1, size_checksum * sizeof(uint8_t));

            for (j=1; j <= min_sizes[i]; j++){
                if (annot[j] == 0){
                    annot[j] = 254;
                    annot[min_sizes[i] + 1 + (j-1)/(SIZE_CELL-1)] |= MASK_POWER_8[(j-1) % (SIZE_CELL-1) + 1];
                }
            }

            JSLI(PValue, *PJArray, annot);
            if (PValue == PJERR) ERROR("Out of memory -- exit\n")

            *PValue += min_sizes[i];

            memset(annot, 0, (min_sizes[i] + size_checksum + 2) * sizeof(uint8_t));
        }

        free(annot);
        free(min_sizes);
        if (extended_annots != NULL) free(extended_annots);
    }
}

int compress_annotation_from_Node(Node* restrict node, int size_kmer, ptrs_on_func* restrict func_on_types, Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted){

    ASSERT_NULL_PTR(node,"compress_annotation_from_Node()")
    ASSERT_NULL_PTR(func_on_types,"compress_annotation_from_Node()")

    int i = -1;
    int tot_sizes_annot = 0;

    if (node->CC_array != NULL){
        do {
            i++;
            tot_sizes_annot += compress_annotation_from_CC(&(((CC*)node->CC_array)[i]), size_kmer, func_on_types, PJArray, old_annot_sorted);
        }
        while ((((CC*)node->CC_array)[i].type & 0x1) == 0);
    }

    if (node->UC_array.suffixes != NULL)
        tot_sizes_annot += compress_annotation_from_UC(&(node->UC_array), func_on_types[(size_kmer/SIZE_SEED)-1].size_kmer_in_bytes,
                                                      node->UC_array.nb_children >> 1, PJArray, old_annot_sorted);

    return tot_sizes_annot;
}

int compress_annotation_from_CC(CC* restrict cc, int size_kmer, ptrs_on_func* restrict func_on_types, Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted){

    ASSERT_NULL_PTR(cc,"compress_annotation_from_CC()")

    int i;

    int level = (size_kmer/SIZE_SEED)-1;
    int tot_sizes_annot = 0;
    int nb_skp = CEIL(cc->nb_elem,NB_CHILDREN_PER_SKP);

    UC* uc;

    if (size_kmer != SIZE_SEED){
        for (i=0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);
            tot_sizes_annot += compress_annotation_from_UC(uc, func_on_types[level].size_kmer_in_bytes_minus_1, uc->nb_children, PJArray, old_annot_sorted);
        }
    }
    else {
        for (i=0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            if (i != nb_skp-1) tot_sizes_annot += compress_annotation_from_UC(uc, 0, NB_CHILDREN_PER_SKP, PJArray, old_annot_sorted);
            else tot_sizes_annot += compress_annotation_from_UC(uc, 0, cc->nb_elem - i * NB_CHILDREN_PER_SKP, PJArray, old_annot_sorted);
        }
    }

    for (i = 0; i < cc->nb_Node_children; i++)
        tot_sizes_annot += compress_annotation_from_Node(&(cc->children_Node_container[i]), size_kmer-SIZE_SEED, func_on_types, PJArray, old_annot_sorted);

    return tot_sizes_annot;
}

int compress_annotation_from_UC(UC* uc, int size_substring, int nb_children, Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted){

    ASSERT_NULL_PTR(uc, "compress_annotation_from_UC()")

    if (uc->suffixes != NULL){

        PWord_t PValue;

        int i = 0, j = 0, z = 0;
        int size_max = 0;
        int size_line = size_substring + uc->size_annot;

        int new_size_line;
        int size_tmp;
        int size_checksum;
        int size_line_annot;
        int max_size_annot;

        uint32_t id = 0;

        uint8_t** extended_annots = get_extend_annots(uc, size_substring, nb_children, 0, nb_children-1);
        uint8_t* min_sizes = min_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        uint8_t* annot;
        uint8_t* annot_tmp;
        uint8_t* annotations;
        uint8_t* new_tab_sub;

        uint32_t* positions = calloc(nb_children, sizeof(uint32_t));
        ASSERT_NULL_PTR(positions,"compress_annotation_from_UC()")

        if (extended_annots != NULL) max_size_annot = uc->size_annot + 1;
        else max_size_annot = max_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        for (z=0; z<nb_children; z++){

            annot = &(uc->suffixes[z * size_line + size_substring]);

            if ((annot[0] & 0x3) == 3){

                j = 0;

                while ((j < min_sizes[z]) && (annot[j] != 0)){
                    if (j == 0) positions[z] = annot[0] >> 2;
                    else positions[z] |= ((uint32_t)(annot[j] >> 1)) << (6+(j-1)*7);
                    j++;
                }

                if ((extended_annots != NULL) && (extended_annots[z] != NULL) && (extended_annots[z][0] != 0))
                    positions[z] |= ((uint32_t)(extended_annots[z][0] >> 1)) << (6+(j-1)*7);

                annot_tmp = extract_from_annotation_array_elem(old_annot_sorted, positions[z], &size_tmp);
                max_size_annot = MAX(max_size_annot, size_tmp);
            }
        }

        size_line_annot = max_size_annot + CEIL(max_size_annot, SIZE_CELL-1) + 2;

        annotations = calloc(nb_children * size_line_annot, sizeof(uint8_t));
        ASSERT_NULL_PTR(annotations,"compress_annotation_from_UC()")

        for (z=0; z<nb_children; z++){

            annot = &(annotations[z * size_line_annot]);

            if ((extended_annots != NULL) && (extended_annots[z] != NULL) && (extended_annots[z][0] != 0)){

                min_sizes[z] = uc->size_annot;
                memcpy(&(annot[1]), &(uc->suffixes[z * size_line + size_substring]), min_sizes[z] * sizeof(uint8_t));
                min_sizes[z]++;
                annot[min_sizes[z]] = extended_annots[z][0];
            }
            else memcpy(&(annot[1]), &(uc->suffixes[z * size_line + size_substring]), min_sizes[z] * sizeof(uint8_t));

            if ((annot[1] & 0x3) == 3){
                annot_tmp = extract_from_annotation_array_elem(old_annot_sorted, positions[z], &size_tmp);
                memcpy(&(annot[1]), annot_tmp, size_tmp * sizeof(uint8_t));
                min_sizes[z] = size_tmp;
                annot[0] = size_tmp;
            }
            else annot[0] = min_sizes[z];

            size_checksum = CEIL(min_sizes[z], SIZE_CELL-1);

            memset(&(annot[annot[0] + 1]), 1, size_checksum * sizeof(uint8_t));

            for (j=1; j <= min_sizes[z]; j++){
                if (annot[j] == 0){
                    annot[j] = 254;
                    annot[min_sizes[z] + 1 + (j-1)/(SIZE_CELL-1)] |= MASK_POWER_8[(j-1) % (SIZE_CELL-1) + 1];
                }
            }

            JSLG(PValue, *PJArray, annot);
            if (PValue == PJERR) ERROR("Out of memory -- exit")
            if (PValue == NULL) ERROR("compress_annotation_from_UC() 1");

            if (*PValue == UINT32_MAX) size_max = MAX(size_max, min_sizes[z]);
            else id = MAX(id, *PValue);
        }

        size_max = MAX(size_max, get_nb_bytes_power2(round_up_next_highest_power2((int)id+1)));

        new_size_line = size_substring + size_max;
        uc->size_annot = size_max;

        new_tab_sub = calloc(nb_children * new_size_line, sizeof(uint8_t));
        ASSERT_NULL_PTR(new_tab_sub,"compress_annotation_from_UC()")

        for (z = 0; z < nb_children; z++){

            memcpy(&(new_tab_sub[z * new_size_line]), &(uc->suffixes[z * size_line]), size_substring * sizeof(uint8_t));

            annot = &(annotations[z * size_line_annot]);

            JSLG(PValue, *PJArray, annot);
            if (PValue == PJERR) ERROR("Out of memory -- exit\n")
            if (PValue == NULL) ERROR("compress_annotation_from_UC() 2");

            if (*PValue != UINT32_MAX){

                id = *PValue;

                new_tab_sub[z * new_size_line + size_substring] = ((id & 0xff) << 2) | 0x3;

                for (i=1; i<get_nb_bytes_power2(round_up_next_highest_power2(id+1)); i++)
                    new_tab_sub[z * new_size_line + size_substring + i] = (((id >> (i*6 + (i-1))) & 0xff) << 1) | 0x1;
            }
            else{
                for (j=1; j<=min_sizes[z]; j++)
                    if ((annot[j] == 254) && ((annot[min_sizes[z] + 1 + (j-1)/(SIZE_CELL-1)] & MASK_POWER_8[(j-1) % (SIZE_CELL-1) + 1]) != 0)) annot[j] = 0;

                memcpy(&(new_tab_sub[z * new_size_line + size_substring]), &(annot[1]), min_sizes[z] * sizeof(uint8_t));
            }
        }

        free(annotations);
        free(positions);
        free(min_sizes);
        if (extended_annots != NULL) free(extended_annots);

        free(uc->suffixes);
        uc->suffixes = new_tab_sub;
        uc->nb_extended_annot = 0;

        create_annot_extended(uc, size_substring, nb_children);

        return nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;
    }

    return 0;
}

void uncompress_annotation_from_Node(Node* restrict node, int size_kmer, ptrs_on_func* restrict func_on_types, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(node,"uncompress_annotation_from_Node()")
    ASSERT_NULL_PTR(func_on_types,"uncompress_annotation_from_Node()")

    int i = -1;

    if (node->CC_array != NULL){
        do {
            i++;
            uncompress_annotation_from_CC(&(((CC*)node->CC_array)[i]), size_kmer, func_on_types, annot_sorted);
        }
        while ((((CC*)node->CC_array)[i].type & 0x1) == 0);
    }

    if (node->UC_array.suffixes != NULL)
        uncompress_annotation_from_UC(&(node->UC_array), func_on_types[(size_kmer/SIZE_SEED)-1].size_kmer_in_bytes, node->UC_array.nb_children >> 1, annot_sorted);

    return;
}

void uncompress_annotation_from_CC(CC* restrict cc, int size_kmer, ptrs_on_func* restrict func_on_types, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(cc,"uncompress_annotation_from_CC()")

    int i;

    int level = (size_kmer/SIZE_SEED)-1;
    int nb_skp = CEIL(cc->nb_elem,NB_CHILDREN_PER_SKP);

    UC* uc;

    if (size_kmer != SIZE_SEED){

        for (i=0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);
            uncompress_annotation_from_UC(uc, func_on_types[level].size_kmer_in_bytes_minus_1, uc->nb_children, annot_sorted);
        }
    }
    else {
        for (i=0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            if (i != nb_skp-1) uncompress_annotation_from_UC(uc, 0, NB_CHILDREN_PER_SKP, annot_sorted);
            else uncompress_annotation_from_UC(uc, 0, cc->nb_elem - i * NB_CHILDREN_PER_SKP, annot_sorted);
        }
    }

    for (i = 0; i < cc->nb_Node_children; i++)
        uncompress_annotation_from_Node(&(cc->children_Node_container[i]), size_kmer-SIZE_SEED, func_on_types, annot_sorted);

    return;
}

void uncompress_annotation_from_UC(UC* uc, int size_substring, int nb_children, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(uc, "uncompress_annotation_from_UC()")

    if (uc->suffixes != NULL){

        int j = 0, z = 0;
        int size_line = size_substring + uc->size_annot;

        int new_size_line;
        int size_tmp;
        int max_size_annot;

        uint8_t** extended_annots = get_extend_annots(uc, size_substring, nb_children, 0, nb_children-1);
        uint8_t* min_sizes = min_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        uint8_t* annot;
        uint8_t* annot_tmp;
        uint8_t* new_tab_sub;

        uint32_t* positions = calloc(nb_children, sizeof(uint32_t));
        ASSERT_NULL_PTR(positions,"uncompress_annotation_from_UC()")

        if (extended_annots != NULL) max_size_annot = uc->size_annot + 1;
        else max_size_annot = max_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        for (z=0; z<nb_children; z++){

            annot = &(uc->suffixes[z * size_line + size_substring]);

            if ((annot[0] & 0x3) == 3){

                j = 0;

                while ((j < min_sizes[z]) && (annot[j] != 0)){
                    if (j == 0) positions[z] = annot[0] >> 2;
                    else positions[z] |= ((uint32_t)(annot[j] >> 1)) << (6+(j-1)*7);
                    j++;
                }

                if ((extended_annots != NULL) && (extended_annots[z] != NULL) && (extended_annots[z][0] != 0))
                    positions[z] |= ((uint32_t)(extended_annots[z][0] >> 1)) << (6+(j-1)*7);

                annot_tmp = extract_from_annotation_array_elem(annot_sorted, positions[z], &size_tmp);
                max_size_annot = MAX(max_size_annot, size_tmp);
            }
        }

        new_size_line = size_substring + max_size_annot;

        new_tab_sub = calloc(nb_children * new_size_line, sizeof(uint8_t));
        ASSERT_NULL_PTR(new_tab_sub,"uncompress_annotation_from_UC()")

        for (z=0; z<nb_children; z++){

            memcpy(&(new_tab_sub[z * new_size_line]), &(uc->suffixes[z * size_line]), size_substring * sizeof(uint8_t));

            annot = &(new_tab_sub[z * new_size_line + size_substring]);

            if ((extended_annots != NULL) && (extended_annots[z] != NULL) && (extended_annots[z][0] != 0)){

                memcpy(annot, &(uc->suffixes[z * size_line + size_substring]), uc->size_annot * sizeof(uint8_t));
                annot[uc->size_annot] = extended_annots[z][0];
            }
            else memcpy(annot, &(uc->suffixes[z * size_line + size_substring]), min_sizes[z] * sizeof(uint8_t));

            if ((annot[0] & 0x3) == 3){
                annot_tmp = extract_from_annotation_array_elem(annot_sorted, positions[z], &size_tmp);
                memcpy(annot, annot_tmp, size_tmp * sizeof(uint8_t));
            }
        }

        free(positions);
        free(min_sizes);
        if (extended_annots != NULL) free(extended_annots);

        free(uc->suffixes);
        uc->suffixes = new_tab_sub;
        uc->size_annot = max_size_annot;
        uc->nb_extended_annot = 0;

        create_annot_extended(uc, size_substring, nb_children);
    }

    return;
}
