#include "replaceAnnotation.h"

void load_annotation_from_Node(Node* node, int lvl_node, int size_kmer, int longest_annot, info_per_level* info_per_lvl,
                               Pvoid_t* PJArray, annotation_array_elem* annot_sorted, annotation_inform* ann_inf){

    ASSERT_NULL_PTR(node,"load_annotation_from_Node()")
    ASSERT_NULL_PTR(info_per_lvl,"load_annotation_from_Node()")

    int i = -1;

    if ((CC*)node->CC_array != NULL){

        do {

            i++;

            load_annotation_from_CC(&(((CC*)node->CC_array)[i]), lvl_node, size_kmer, longest_annot, info_per_lvl,
                                    PJArray, annot_sorted, ann_inf);
        }
        while (IS_EVEN(((CC*)node->CC_array)[i].type));
    }

    if (node->UC_array.suffixes != NULL){

        load_annotation_from_UC(&(node->UC_array), info_per_lvl[lvl_node].size_kmer_in_bytes,
                                node->UC_array.nb_children >> 1, longest_annot, PJArray, annot_sorted, ann_inf);
    }

    return;
}

void load_annotation_from_CC(CC*  cc, int lvl_cc, int size_kmer, int longest_annot, info_per_level*  info_per_lvl,
                             Pvoid_t* PJArray, annotation_array_elem* annot_sorted, annotation_inform* ann_inf){

    ASSERT_NULL_PTR(cc,"load_annotation_from_CC()")

    UC* uc;

    int i;

    int nb_skp = CEIL(cc->nb_elem, info_per_lvl[lvl_cc].nb_ucs_skp);

    if (size_kmer != NB_CHAR_SUF_PREF){

        for (i = 0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            load_annotation_from_UC(uc, info_per_lvl[lvl_cc].size_kmer_in_bytes_minus_1, uc->nb_children,
                                    longest_annot, PJArray, annot_sorted, ann_inf);

            free(uc->suffixes);
            uc->suffixes = NULL;
        }
    }
    else{
        for (i = 0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            if (i != nb_skp-1){

                load_annotation_from_UC(uc, 0, info_per_lvl[lvl_cc].nb_ucs_skp, longest_annot,
                                        PJArray, annot_sorted, ann_inf);
            }
            else{

                load_annotation_from_UC(uc, 0, cc->nb_elem - i * info_per_lvl[lvl_cc].nb_ucs_skp,
                                        longest_annot, PJArray, annot_sorted, ann_inf);
            }

            free(uc->suffixes);
            uc->suffixes = NULL;
        }
    }

    for (i = 0; i < cc->nb_Node_children; i++){

        load_annotation_from_Node(&(cc->children_Node_container[i]), lvl_cc-1, size_kmer-NB_CHAR_SUF_PREF,
                                  longest_annot, info_per_lvl, PJArray, annot_sorted, ann_inf);

        freeNode(&(cc->children_Node_container[i]), lvl_cc-1, info_per_lvl);
    }

    free(cc->children_Node_container);
    cc->children_Node_container = NULL;

    return;
}

void load_annotation_from_UC(UC* uc, int size_substring, int nb_children, int longest_annot,
                             Pvoid_t* PJArray, annotation_array_elem* annot_sorted, annotation_inform* ann_inf){

    ASSERT_NULL_PTR(uc, "load_annotation_from_UC()")

    if (uc->suffixes != NULL){

        PWord_t PValue;

        int i = 0, j = 0;
        int bits_per_byte_checksum = SIZE_BITS_UINT_8T - 1;
        int size_line = size_substring + uc->size_annot;

        int shift_mode_3;
        int size_checksum;

        uint32_t position;

        uint64_t size_decomp;

        uint8_t* annot_tmp;
        uint8_t* annot_start;

        uint8_t* annot = calloc(longest_annot + CEIL(longest_annot, bits_per_byte_checksum) + 4, sizeof(uint8_t));
        ASSERT_NULL_PTR(annot, "load_annotation_from_UC()")

        uint8_t** extended_annots = get_extend_annots(uc, size_substring, nb_children, 0, nb_children-1);

        UC_SIZE_ANNOT_T real_size;

        UC_SIZE_ANNOT_T *min_sizes = min_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        for (i=0; i<nb_children; i++){

            if ((extended_annots != NULL) && (extended_annots[i] != NULL) && (extended_annots[i][0] != 0)){
                memcpy(&annot[3], &(uc->suffixes[i * size_line + size_substring]), uc->size_annot * sizeof(uint8_t));
                annot[uc->size_annot + 3] = extended_annots[i][0];
                min_sizes[i] = uc->size_annot + 1;
            }
            else memcpy(&annot[3], &(uc->suffixes[i * size_line + size_substring]), min_sizes[i] * sizeof(uint8_t));

            if ((annot[3] & 0x3) == 3){

                j = 4;
                shift_mode_3 = 5;
                position = annot[3] >> 2;

                while ((j < min_sizes[i]+3) && annot[j]){
                    position |= ((uint32_t)(annot[j] & 0xfe)) << shift_mode_3;
                    shift_mode_3 += 7;
                    j++;
                }

                real_size = min_sizes[i];
                annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &min_sizes[i]);

                j = decomp_annotation(ann_inf, annot_tmp, min_sizes[i], NULL, 0, 1);
                size_decomp = (j != 0 ? j : min_sizes[i]);

                reinit_annotation_inform(ann_inf);
            }
            else {
                size_decomp = min_sizes[i];
                j = comp_annotation(ann_inf, &annot[3], min_sizes[i], NULL, 0);
                if (j) min_sizes[i] = j;
                real_size = min_sizes[i];
            }

            annot[0] = (min_sizes[i] >> 8) | 0x3;
            annot[1] = (min_sizes[i] >> 2) | 0x3;
            annot[2] = (min_sizes[i] << 4) | 0x3;

            size_checksum = CEIL(real_size, bits_per_byte_checksum);

            memset(&annot[real_size + 3], 1, size_checksum * sizeof(uint8_t));

            j = 0;
            annot_start = annot + 3;
            annot_tmp = annot_start + real_size;

            while (annot_start < annot_tmp){

                if (*annot_start == 0x0){
                    *annot_start = 0xfe;
                    annot_tmp[j/bits_per_byte_checksum] |= MASK_POWER_8[j % bits_per_byte_checksum + 1];
                    j++;
                }
                else if (*annot_start == 0xfe) j++;

                annot_start++;
            }

            annot_tmp[CEIL(j, bits_per_byte_checksum)] = '\0';

            JSLI(PValue, *PJArray, annot);
            if (PValue == PJERR) ERROR("Out of memory -- exit\n")

            #if defined (_WORDx64)
                if (*PValue == 0) *PValue = (size_decomp << 32) | size_decomp;
                else *PValue += size_decomp << 32;
            #else
                if (*PValue == NULL){
                    *PValue = malloc(sizeof(uint64_t));
                    ASSERT_NULL_PTR(*PValue, "load_annotation_from_UC()")

                    **((uint64_t**)PValue) = (size_decomp << 32) | size_decomp;
                }
                else **((uint64_t**)PValue) += size_decomp << 32;
            #endif

            memset(annot, 0, (real_size + size_checksum + 4) * sizeof(uint8_t));
        }

        free(annot);
        free(min_sizes);
        if (extended_annots != NULL) free(extended_annots);
    }
}

int compress_annotation_from_UC(UC* uc, int size_substring, int nb_children,
                                Pvoid_t* new_annots, annotation_array_elem* old_annots,
                                annotation_inform* ann_inf, bool insert_comp_annot){

    ASSERT_NULL_PTR(uc, "compress_annotation_from_UC() 5")

    if (uc->suffixes != NULL){

        PWord_t PValue;

        int i = 0, j = 0, k = 0, z = 0;
        int size_max = 0;
        int size_max_annot_cmp = 0;
        int longest_annot = 0;
        int bits_per_byte_checksum = SIZE_BITS_UINT_8T - 1;
        int size_line = size_substring + uc->size_annot;

        int shift_mode_3;
        int new_size_line;
        int size_tmp;
        int size_checksum;
        int size_line_annot;
        int max_size_annot;
        int size_current_annot;

        UC_SIZE_ANNOT_T real_size;

        UC_SIZE_ANNOT_T *min_sizes = min_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        uint8_t* annot;
        uint8_t* annot_tmp;
        uint8_t* annot_start;
        uint8_t* annotations;
        uint8_t* new_tab_sub;
        uint8_t* new_tab_sub_tmp;

        uint8_t** extended_annots = get_extend_annots(uc, size_substring, nb_children, 0, nb_children-1);

        uint32_t id = 0;

        uint32_t* positions = calloc(nb_children, sizeof(uint32_t));
        ASSERT_NULL_PTR(positions,"compress_annotation_from_UC() 4")

        if (extended_annots != NULL) max_size_annot = uc->size_annot + 1;
        else max_size_annot = max_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        for (z = 0; z < nb_children; z++){

            annot = &(uc->suffixes[z * size_line + size_substring]);

            if ((annot[0] & 0x3) == 3){

                j = 1;
                shift_mode_3 = 5;
                positions[z] = annot[0] >> 2;

                while ((j < min_sizes[z]) && annot[j]){
                    positions[z] |= ((uint32_t)(annot[j] & 0xfe)) << shift_mode_3;
                    shift_mode_3 += 7;
                    j++;
                }

                if ((extended_annots != NULL) && (extended_annots[z] != NULL) && (extended_annots[z][0]))
                    positions[z] |= ((uint32_t)(extended_annots[z][0] & 0xfe)) << shift_mode_3;

                annot_tmp = extract_from_annotation_array_elem(old_annots, positions[z], &size_tmp);

                max_size_annot = MAX(max_size_annot, size_tmp);

                if (insert_comp_annot) size_max_annot_cmp = MAX(size_max_annot_cmp, size_tmp);
                else {

                    size_current_annot = decomp_annotation(ann_inf, annot_tmp, size_tmp, NULL, 0, 1);

                    if (size_current_annot) size_max_annot_cmp = MAX(size_max_annot_cmp, size_current_annot);
                    else size_max_annot_cmp = MAX(size_max_annot_cmp, size_tmp);

                    reinit_annotation_inform(ann_inf);
                }
            }
        }

        size_line_annot = max_size_annot + CEIL(max_size_annot, bits_per_byte_checksum) + 4;

        annotations = calloc(nb_children * size_line_annot, sizeof(uint8_t));
        ASSERT_NULL_PTR(annotations,"compress_annotation_from_UC() 3")

        for (z = 0; z < nb_children; z++){

            annot = &(annotations[z * size_line_annot]);

            if ((extended_annots != NULL) && (extended_annots[z] != NULL) && extended_annots[z][0]){

                memcpy(&annot[3], &(uc->suffixes[z * size_line + size_substring]), uc->size_annot * sizeof(uint8_t));
                annot[uc->size_annot+3] = extended_annots[z][0];
                min_sizes[z] = uc->size_annot + 1;
            }
            else memcpy(&annot[3], &(uc->suffixes[z * size_line + size_substring]), min_sizes[z] * sizeof(uint8_t));

            if ((annot[3] & 0x3) == 3){
                real_size = min_sizes[z];
                extract_from_annotation_array_elem(old_annots, positions[z], &min_sizes[z]);
            }
            else {
                size_max_annot_cmp = MAX(size_max_annot_cmp, min_sizes[z]);
                j = comp_annotation(ann_inf, &annot[3], min_sizes[z], NULL, 0);
                if (j) min_sizes[z] = j;
                real_size = min_sizes[z];
            }

            annot[0] = (min_sizes[z] >> 8) | 0x3;
            annot[1] = (min_sizes[z] >> 2) | 0x3;
            annot[2] = (min_sizes[z] << 4) | 0x3;

            size_checksum = CEIL(real_size, bits_per_byte_checksum);

            memset(&annot[real_size + 3], 1, size_checksum * sizeof(uint8_t));

            j = 0;
            annot_start = annot + 3;
            annot_tmp = annot_start + real_size;

            while (annot_start < annot_tmp){

                if (*annot_start == 0x0){
                    *annot_start = 0xfe;
                    annot_tmp[j/bits_per_byte_checksum] |= MASK_POWER_8[j % bits_per_byte_checksum + 1];
                    j++;
                }
                else if (*annot_start == 0xfe) j++;

                annot_start++;
            }

            annot_tmp[CEIL(j, bits_per_byte_checksum)] = '\0';

            JSLG(PValue, *new_annots, annot);
            if (PValue == PJERR) ERROR("Out of memory -- exit")

            #if defined (_WORDx64)
                if (PValue == NULL) size_max = MAX(size_max, min_sizes[z]);
                else id = MAX(id, *PValue);
            #else
                if (PValue == NULL) size_max = MAX(size_max, min_sizes[z]);
                else id = MAX(id, **((uint64_t**)PValue));
            #endif
        }

        size_max = MAX(MAX(size_max, get_nb_bytes_power2_comp(round_up_next_highest_power2((int)id+1))), size_max_annot_cmp);

        new_size_line = size_substring + size_max;

        new_tab_sub = calloc(nb_children * new_size_line, sizeof(uint8_t));
        ASSERT_NULL_PTR(new_tab_sub,"compress_annotation_from_UC() 2")

        for (z = 0; z < nb_children; z++){

            memcpy(&new_tab_sub[z * new_size_line], &(uc->suffixes[z * size_line]), size_substring * sizeof(uint8_t));

            new_tab_sub_tmp = &new_tab_sub[z * new_size_line + size_substring];
            annot = &annotations[z * size_line_annot];

            JSLG(PValue, *new_annots, annot);
            if (PValue == PJERR) ERROR("Out of memory -- exit\n")

            if (PValue != NULL){

                #if defined (_WORDx86)
                    PValue = *PValue;
                #endif

                id = *PValue;
                size_current_annot = get_nb_bytes_power2_comp(round_up_next_highest_power2(id+1));
                longest_annot = MAX(longest_annot, size_current_annot);

                new_tab_sub_tmp[0] = ((id & 0xff) << 2) | 0x3;

                for (i = 1, k = 0; i < size_current_annot; i++, k++) new_tab_sub_tmp[i] = (((id >> (i * 6 + k)) & 0xff) << 1) | 0x1;
            }
            else if ((annot[3] & 0x3) == 3){

                if (insert_comp_annot){

                    annot_tmp = extract_from_annotation_array_elem(old_annots, positions[z], &size_tmp);
                    memcpy(new_tab_sub_tmp, annot_tmp, size_tmp * sizeof(uint8_t));
                    longest_annot = MAX(longest_annot, size_tmp);
                }
                else{

                    annot_tmp = extract_from_annotation_array_elem(old_annots, positions[z], &size_tmp);
                    longest_annot = MAX(longest_annot, decomp_annotation(ann_inf, annot_tmp, size_tmp, NULL, 0, 1));

                    uint8_t flag1 = annot_tmp[0] & 0x3;

                    if (flag1){

                        int it_annot = 0, size_current_annot = size_max;
                        uint8_t flag2 = (flag1 == 1 ? 2 : 1);

                        for (int y = 0; y < ann_inf->nb_id_stored; y++){

                            modify_annot_bis(&new_tab_sub_tmp, NULL, &it_annot, &size_current_annot,
                                             ann_inf->id_stored[y], ann_inf->size_id_stored[y], flag1, flag2);
                        }
                    }
                    else{
                        memcpy(new_tab_sub_tmp, annot_tmp, size_tmp * sizeof(uint8_t));
                        longest_annot = MAX(longest_annot, size_tmp);
                    }

                    reinit_annotation_inform(ann_inf);
                }
            }
            else if (insert_comp_annot){

                j = 0;
                annot_start = annot + 3;
                annot_tmp = annot_start + min_sizes[z];

                while (annot_start < annot_tmp){

                    if (*annot_start == 0xfe){
                        if (annot_tmp[j/bits_per_byte_checksum] & MASK_POWER_8[j % bits_per_byte_checksum + 1]) *annot_start = 0x0;
                        j++;
                    }

                    annot_start++;
                }

                memcpy(new_tab_sub_tmp, &annot[3], min_sizes[z] * sizeof(uint8_t));
                longest_annot = MAX(longest_annot, min_sizes[z]);
            }
            else {

                if ((extended_annots != NULL) && (extended_annots[z] != NULL) && extended_annots[z][0]){

                    memcpy(new_tab_sub_tmp, &(uc->suffixes[z * size_line + size_substring]), uc->size_annot * sizeof(uint8_t));
                    new_tab_sub_tmp[uc->size_annot] = extended_annots[z][0];
                    longest_annot = MAX(longest_annot, uc->size_annot + 1);
                }
                else {

                    real_size = size_annot_sub(&(uc->suffixes[z * size_line + size_substring]), 0, uc->size_annot);
                    memcpy(new_tab_sub_tmp, &(uc->suffixes[z * size_line + size_substring]), real_size * sizeof(uint8_t));
                    longest_annot = MAX(longest_annot, real_size);
                }
            }
        }

        free(annotations);
        free(positions);
        free(min_sizes);
        free(extended_annots);
        free(uc->suffixes);

        uc->suffixes = new_tab_sub;
        uc->size_annot = size_max;
        uc->nb_extended_annot = 0;

        if (longest_annot < uc->size_annot) realloc_annotation(uc, size_substring, nb_children, longest_annot, 0, 0);
        create_annot_extended(uc, size_substring, nb_children);

        return nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;
    }

    return 0;
}

/*void load_annotation_from_Node(Node*  node, int lvl_node, int size_kmer, int longest_annot, info_per_level*  info_per_lvl,
                               Pvoid_t* PJArray, annotation_array_elem* annot_sorted, annotation_inform* ann_inf){

    ASSERT_NULL_PTR(node,"load_annotation_from_Node()")
    ASSERT_NULL_PTR(info_per_lvl,"load_annotation_from_Node()")

    int i = -1;

    if ((CC*)node->CC_array != NULL){
        do {
            i++;
            load_annotation_from_CC(&(((CC*)node->CC_array)[i]), lvl_node, size_kmer, longest_annot, info_per_lvl,
                                    PJArray, annot_sorted, ann_inf);
        }
        while (IS_EVEN(((CC*)node->CC_array)[i].type));
    }

    if (node->UC_array.suffixes != NULL){
        load_annotation_from_UC(&(node->UC_array), info_per_lvl[lvl_node].size_kmer_in_bytes,
                                node->UC_array.nb_children >> 1, longest_annot, PJArray, annot_sorted, ann_inf);
    }
}

void load_annotation_from_CC(CC*  cc, int lvl_cc, int size_kmer, int longest_annot, info_per_level*  info_per_lvl,
                             Pvoid_t* PJArray, annotation_array_elem* annot_sorted, annotation_inform* ann_inf){

    ASSERT_NULL_PTR(cc,"load_annotation_from_CC()")

    UC* uc;

    int i;

    int nb_skp = CEIL(cc->nb_elem, info_per_lvl[lvl_cc].nb_ucs_skp);

    if (size_kmer != NB_CHAR_SUF_PREF){

        for (i = 0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);
            load_annotation_from_UC(uc, info_per_lvl[lvl_cc].size_kmer_in_bytes_minus_1, uc->nb_children,
                                    longest_annot, PJArray, annot_sorted, ann_inf);
        }
    }
    else{
        for (i = 0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            if (i != nb_skp-1){
                load_annotation_from_UC(uc, 0, info_per_lvl[lvl_cc].nb_ucs_skp, longest_annot,
                                        PJArray, annot_sorted, ann_inf);
            }
            else{
                load_annotation_from_UC(uc, 0, cc->nb_elem - i * info_per_lvl[lvl_cc].nb_ucs_skp, longest_annot,
                                        PJArray, annot_sorted, ann_inf);
            }
        }
    }

    for (i = 0; i < cc->nb_Node_children; i++){
        load_annotation_from_Node(&(cc->children_Node_container[i]), lvl_cc-1, size_kmer-NB_CHAR_SUF_PREF,
                                  longest_annot, info_per_lvl, PJArray, annot_sorted, ann_inf);
    }

    return;
}

void load_annotation_from_UC(UC* uc, int size_substring, int nb_children, int longest_annot,
                             Pvoid_t* PJArray, annotation_array_elem* annot_sorted, annotation_inform* ann_inf){

    ASSERT_NULL_PTR(uc, "load_annotation_from_UC()")

    if (uc->suffixes != NULL){

        PWord_t PValue;

        uint64_t size_decomp;

        int i = 0, j = 0;
        int bits_per_byte_checksum = SIZE_BITS_UINT_8T - 1;
        int size_line = size_substring + uc->size_annot;

        int shift_mode_3;
        int size_checksum;
        int size_tmp;

        uint32_t position;

        uint8_t* annot_tmp;
        uint8_t* annot_start;
        uint8_t** extended_annots = get_extend_annots(uc, size_substring, nb_children, 0, nb_children-1);

        UC_SIZE_ANNOT_T *min_sizes = min_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        uint8_t* annot = calloc(longest_annot + CEIL(longest_annot, bits_per_byte_checksum) + 4, sizeof(uint8_t));
        ASSERT_NULL_PTR(annot, "load_annotation_from_UC()")

        for (i=0; i<nb_children; i++){

            if ((extended_annots != NULL) && (extended_annots[i] != NULL) && (extended_annots[i][0] != 0)){
                memcpy(&(annot[3]), &(uc->suffixes[i * size_line + size_substring]), uc->size_annot * sizeof(uint8_t));
                annot[uc->size_annot + 3] = extended_annots[i][0];
                min_sizes[i] = uc->size_annot + 1;
            }
            else memcpy(&(annot[3]), &(uc->suffixes[i * size_line + size_substring]), min_sizes[i] * sizeof(uint8_t));

            if ((annot[3] & 0x3) == 3){

                j = 4;
                shift_mode_3 = 6;
                position = annot[3] >> 2;

                while ((j < min_sizes[i]+3) && (annot[j] != 0)){
                    position |= ((uint32_t)(annot[j] >> 1)) << shift_mode_3;
                    shift_mode_3 += 7;
                    j++;
                }

                annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &size_tmp);
                memcpy(&(annot[3]), annot_tmp, size_tmp*sizeof(uint8_t));
                min_sizes[i] = size_tmp;

                if ((j = decomp_annotation(ann_inf, annot_tmp, size_tmp, NULL, 0, 1)) != 0) size_decomp = j;
                else size_decomp = size_tmp;

                reinit_annotation_inform(ann_inf);
            }
            else {
                size_decomp = min_sizes[i];
                j = comp_annotation(ann_inf, &(annot[3]), min_sizes[i], NULL, 0);
                if (j) min_sizes[i] = j;
            }

            annot[0] = (min_sizes[i] >> 8) | 0x3;
            annot[1] = (min_sizes[i] >> 2) | 0x3;
            annot[2] = (min_sizes[i] << 4) | 0x3;

            size_checksum = CEIL(min_sizes[i], bits_per_byte_checksum);

            memset(&(annot[min_sizes[i] + 3]), 1, size_checksum * sizeof(uint8_t));

            //j = -1;
            //annot_start = annot + 3;
            //annot_tmp = annot_start - 1;
            //while ((j < min_sizes[i]) && ((annot_tmp = memchr(annot_tmp + 1, 0, min_sizes[i] - j - 1)) != NULL)){
            //    *annot_tmp = 254;
            //    j = annot_tmp - annot_start;
            //    annot_start[min_sizes[i] + j/bits_per_byte_checksum] |= MASK_POWER_8[j % bits_per_byte_checksum + 1];
            //}

            j = 0;
            annot_start = annot + 3;
            annot_tmp = annot_start + min_sizes[i];

            while (annot_start < annot_tmp){

                if (*annot_start == 0){
                    *annot_start = 254;
                    annot_tmp[j/bits_per_byte_checksum] |= MASK_POWER_8[j % bits_per_byte_checksum + 1];
                    j++;
                }
                else if (*annot_start == 254) j++;

                annot_start++;
            }

            annot_tmp[CEIL(j, bits_per_byte_checksum)] = '\0';

            JSLI(PValue, *PJArray, annot);
            if (PValue == PJERR) ERROR("Out of memory -- exit\n")

            #if defined (_WORDx64)
                if (*PValue == 0) *PValue = (size_decomp << 32) | size_decomp;
                else *PValue += size_decomp << 32;
            #else
                if (*PValue == NULL){
                    *PValue = malloc(sizeof(uint64_t));
                    ASSERT_NULL_PTR(*PValue, "load_annotation_from_UC()")

                    **((uint64_t**)PValue) = (size_decomp << 32) | size_decomp;
                }
                else **((uint64_t**)PValue) += size_decomp << 32;
            #endif

            memset(annot, 0, (min_sizes[i] + size_checksum + 4) * sizeof(uint8_t));
        }

        free(annot);
        free(min_sizes);
        if (extended_annots != NULL) free(extended_annots);
    }
}

int compress_annotation_from_Node(Node*  node, int lvl_node, int size_kmer, info_per_level*  info_per_lvl,
                                  Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted, annotation_inform* ann_inf){

    ASSERT_NULL_PTR(node,"compress_annotation_from_Node()")
    ASSERT_NULL_PTR(info_per_lvl,"compress_annotation_from_Node()")

    int i = -1;
    int tot_sizes_annot = 0;

    if (node->CC_array != NULL){

        do {
            i++;
            tot_sizes_annot += compress_annotation_from_CC(&(((CC*)node->CC_array)[i]), lvl_node, size_kmer,
                                                           info_per_lvl, PJArray, old_annot_sorted, ann_inf);
        }
        while (IS_EVEN(((CC*)node->CC_array)[i].type));
    }

    if (node->UC_array.suffixes != NULL){

        tot_sizes_annot += compress_annotation_from_UC(&(node->UC_array), info_per_lvl[lvl_node].size_kmer_in_bytes,
                                                      node->UC_array.nb_children >> 1, PJArray, old_annot_sorted, ann_inf);
    }

    return tot_sizes_annot;
}

int compress_annotation_from_CC(CC*  cc, int lvl_cc, int size_kmer, info_per_level*  info_per_lvl,
                                Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted, annotation_inform* ann_inf){

    ASSERT_NULL_PTR(cc,"compress_annotation_from_CC()")

    int i;

    int tot_sizes_annot = 0;
    int nb_skp = CEIL(cc->nb_elem, info_per_lvl[lvl_cc].nb_ucs_skp);

    UC* uc;

    if (size_kmer != NB_CHAR_SUF_PREF){

        for (i=0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);
            tot_sizes_annot += compress_annotation_from_UC(uc, info_per_lvl[lvl_cc].size_kmer_in_bytes_minus_1,
                                                           uc->nb_children, PJArray, old_annot_sorted, ann_inf);
        }
    }
    else {
        for (i=0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);

            if (i != nb_skp-1){

                tot_sizes_annot += compress_annotation_from_UC(uc, 0, info_per_lvl[lvl_cc].nb_ucs_skp,
                                                               PJArray, old_annot_sorted, ann_inf);
            }
            else{
                tot_sizes_annot += compress_annotation_from_UC(uc, 0, cc->nb_elem - i * info_per_lvl[lvl_cc].nb_ucs_skp,
                                                               PJArray, old_annot_sorted, ann_inf);
            }
        }
    }

    for (i = 0; i < cc->nb_Node_children; i++){
        tot_sizes_annot += compress_annotation_from_Node(&(cc->children_Node_container[i]), lvl_cc-1, size_kmer-NB_CHAR_SUF_PREF,
                                                         info_per_lvl, PJArray, old_annot_sorted, ann_inf);
    }

    return tot_sizes_annot;
}

int compress_annotation_from_UC(UC* uc, int size_substring, int nb_children, Pvoid_t* PJArray, annotation_array_elem* old_annot_sorted,
                                annotation_inform* ann_inf){

    ASSERT_NULL_PTR(uc, "compress_annotation_from_UC() 5")

    if (uc->suffixes != NULL){

        PWord_t PValue;

        int i = 0, j = 0, k = 0, z = 0;
        int size_max = 0;
        int bits_per_byte_checksum = SIZE_BITS_UINT_8T-1;
        int size_line = size_substring + uc->size_annot;

        int shift_mode_3;
        int new_size_line;
        int size_tmp;
        int size_checksum;
        int size_line_annot;
        int max_size_annot;

        uint32_t id = 0;

        uint8_t** extended_annots = get_extend_annots(uc, size_substring, nb_children, 0, nb_children-1);
        UC_SIZE_ANNOT_T *min_sizes = min_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        uint8_t* annot;
        uint8_t* annot_tmp;
        uint8_t* annot_start;
        uint8_t* annotations;
        uint8_t* new_tab_sub;
        uint8_t* new_tab_sub_tmp;

        uint8_t* current_annot;
        uint8_t flag1, flag2;
        int it_annot, size_current_annot;
        int size_max_annot_cmp = 0;
        int longest_annot = 0;

        uint32_t* positions = calloc(nb_children, sizeof(uint32_t));
        ASSERT_NULL_PTR(positions,"compress_annotation_from_UC() 4")

        if (extended_annots != NULL) max_size_annot = uc->size_annot + 1;
        else max_size_annot = max_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);

        for (z=0; z<nb_children; z++){

            annot = &(uc->suffixes[z * size_line + size_substring]);

            if ((annot[0] & 0x3) == 3){

                j = 1;
                shift_mode_3 = 6;
                positions[z] = annot[0] >> 2;

                while ((j < min_sizes[z]) && (annot[j] != 0)){
                    positions[z] |= ((uint32_t)(annot[j] >> 1)) << shift_mode_3;
                    shift_mode_3 += 7;
                    j++;
                }

                if ((extended_annots != NULL) && (extended_annots[z] != NULL) && (extended_annots[z][0] != 0))
                    positions[z] |= ((uint32_t)(extended_annots[z][0] >> 1)) << shift_mode_3;

                annot_tmp = extract_from_annotation_array_elem(old_annot_sorted, positions[z], &size_tmp);
                max_size_annot = MAX(max_size_annot, size_tmp);

                if ((j = decomp_annotation(ann_inf, annot_tmp, size_tmp, NULL, 0, 1)) != 0) size_tmp = j;

                size_max_annot_cmp = MAX(size_max_annot_cmp, size_tmp);
                reinit_annotation_inform(ann_inf);
            }
        }

        size_line_annot = max_size_annot + CEIL(max_size_annot, bits_per_byte_checksum) + 4;

        annotations = calloc(nb_children * size_line_annot, sizeof(uint8_t));
        ASSERT_NULL_PTR(annotations,"compress_annotation_from_UC() 3")

        for (z=0; z<nb_children; z++){

            annot = &(annotations[z * size_line_annot]);

            if ((extended_annots != NULL) && (extended_annots[z] != NULL) && (extended_annots[z][0] != 0)){

                memcpy(&(annot[3]), &(uc->suffixes[z * size_line + size_substring]), uc->size_annot * sizeof(uint8_t));
                annot[uc->size_annot+3] = extended_annots[z][0];
                min_sizes[z] = uc->size_annot + 1;
            }
            else memcpy(&(annot[3]), &(uc->suffixes[z * size_line + size_substring]), min_sizes[z] * sizeof(uint8_t));

            if ((annot[3] & 0x3) == 3){
                annot_tmp = extract_from_annotation_array_elem(old_annot_sorted, positions[z], &size_tmp);
                memcpy(&(annot[3]), annot_tmp, size_tmp * sizeof(uint8_t));
                min_sizes[z] = size_tmp;
            }
            else {
                size_max_annot_cmp = MAX(size_max_annot_cmp, min_sizes[z]);
                j = comp_annotation(ann_inf, &(annot[3]), min_sizes[z], NULL, 0);
                if (j != 0) min_sizes[z] = j;
            }

            annot[0] = (min_sizes[z] >> 8) | 0x3;
            annot[1] = (min_sizes[z] >> 2) | 0x3;
            annot[2] = (min_sizes[z] << 4) | 0x3;

            size_checksum = CEIL(min_sizes[z], bits_per_byte_checksum);

            memset(&annot[min_sizes[z] + 3], 1, size_checksum * sizeof(uint8_t));

            //j = -1;
            //annot_start = annot + 3;
            //annot_tmp = annot_start - 1;
            //while (((annot_tmp = memchr(annot_tmp + 1, 0, min_sizes[z] - j - 1)) != NULL)){
            //    *annot_tmp = 254;
            //    j = annot_tmp - annot_start;
            //    annot_start[min_sizes[z] + j/bits_per_byte_checksum] |= MASK_POWER_8[j % bits_per_byte_checksum + 1];
            //}

            j = 0;
            annot_start = annot + 3;
            annot_tmp = annot_start + min_sizes[z];

            while (annot_start < annot_tmp){

                if (*annot_start == 0){
                    *annot_start = 254;
                    annot_tmp[j/bits_per_byte_checksum] |= MASK_POWER_8[j % bits_per_byte_checksum + 1];
                    j++;
                }
                else if (*annot_start == 254) j++;

                annot_start++;
            }

            annot_tmp[CEIL(j, bits_per_byte_checksum)] = '\0';

            JSLG(PValue, *PJArray, annot);
            if (PValue == PJERR) ERROR("Out of memory -- exit")
            //if (PValue == NULL) ERROR("compress_annotation_from_UC() 1");
            //
            //#if defined (_WORDx64)
            //    if (*PValue == UINT32_MAX) size_max = MAX(size_max, min_sizes[z]);
            //    else id = MAX(id, *PValue);
            //#else
            //    if (**((uint64_t**)PValue) == UINT32_MAX) size_max = MAX(size_max, min_sizes[z]);
            //    else id = MAX(id, **((uint64_t**)PValue));
            //#endif

            #if defined (_WORDx64)
                if (PValue == NULL) size_max = MAX(size_max, min_sizes[z]);
                else id = MAX(id, *PValue);
            #else
                if (PValue == NULL) size_max = MAX(size_max, min_sizes[z]);
                else id = MAX(id, **((uint64_t**)PValue));
            #endif
        }

        size_max = MAX(MAX(size_max, get_nb_bytes_power2_comp(round_up_next_highest_power2((int)id+1))), size_max_annot_cmp);

        new_size_line = size_substring + size_max;
        uc->size_annot = size_max;

        new_tab_sub = calloc(nb_children * new_size_line, sizeof(uint8_t));
        ASSERT_NULL_PTR(new_tab_sub,"compress_annotation_from_UC() 2")

        for (z = 0; z < nb_children; z++){

            memcpy(&(new_tab_sub[z * new_size_line]), &(uc->suffixes[z * size_line]), size_substring * sizeof(uint8_t));

            annot = &(annotations[z * size_line_annot]);

            JSLG(PValue, *PJArray, annot);
            if (PValue == PJERR) ERROR("Out of memory -- exit\n")
            //if (PValue == NULL) ERROR("compress_annotation_from_UC() 6");

            //#if defined (_WORDx86)
            //    PValue = *PValue;
            //#endif
            //
            //if (*PValue != UINT32_MAX){

            new_tab_sub_tmp = &new_tab_sub[z * new_size_line + size_substring];

            if (PValue != NULL){

                #if defined (_WORDx86)
                    PValue = *PValue;
                #endif

                id = *PValue;
                size_current_annot = get_nb_bytes_power2_comp(round_up_next_highest_power2(id+1));
                longest_annot = MAX(longest_annot, size_current_annot);

                new_tab_sub_tmp[0] = ((id & 0xff) << 2) | 0x3;

                for (i=1, k=0; i < size_current_annot; i++, k++) new_tab_sub_tmp[i] = (((id >> (i*6 + k)) & 0xff) << 1) | 0x1;
            }
            else{

                //j = -1;
                //annot_start = annot + 3;
                //annot_tmp = annot_start - 1;
                //while (((annot_tmp = memchr(annot_tmp + 1, 254, min_sizes[z] - j - 1)) != NULL)){
                //    j = annot_tmp - annot_start;
                //    if (annot_start[min_sizes[z] + j/bits_per_byte_checksum] & MASK_POWER_8[j % bits_per_byte_checksum + 1]) *annot_tmp = 0;
                //}

                j = 0;
                annot_start = annot + 3;
                annot_tmp = annot_start + min_sizes[z];

                while (annot_start < annot_tmp){

                    if (*annot_start == 254){
                        if (annot_tmp[j/bits_per_byte_checksum] & MASK_POWER_8[j % bits_per_byte_checksum + 1]) *annot_start = 0;
                        j++;
                    }

                    annot_start++;
                }

                //decomp_annotation(ann_inf, annot_start, min_sizes[z], NULL, 0, 1);

                //if (ann_inf->comp_annot <= 0){
                    memcpy(new_tab_sub_tmp, &annot[3], min_sizes[z] * sizeof(uint8_t));
                    longest_annot = MAX(longest_annot, min_sizes[z]);
                //}
                //else {
                //
                //    current_annot = &(new_tab_sub[z * new_size_line + size_substring]);
                //    it_annot = 0;
                //    size_max_annot_cmp = 0;
                //    size_current_annot = size_max;
                //
                //    if ((annot_start[0] & 0x3) == 1){
                //        flag1 = 1;
                //        flag2 = 2;
                //    }
                //    else{
                //        flag1 = 2;
                //        flag2 = 1;
                //    }
                //
                //    for (i = 0; i < ann_inf->nb_id_stored; i++){
                //        size_max_annot_cmp += modify_annot_bis(&current_annot, NULL, &it_annot, &size_current_annot,
                //                                               ann_inf->id_stored[i], ann_inf->size_id_stored[i], flag1, flag2);
                //    }
                //
                //    longest_annot = MAX(longest_annot, size_max_annot_cmp);
                //}
                //
                //reinit_annotation_inform(ann_inf);
            }
        }

        free(annotations);
        free(positions);
        free(min_sizes);
        if (extended_annots != NULL) free(extended_annots);

        free(uc->suffixes);
        uc->suffixes = new_tab_sub;
        uc->nb_extended_annot = 0;

        if (longest_annot < uc->size_annot) realloc_annotation(uc, size_substring, nb_children, longest_annot, 0, 0);
        create_annot_extended(uc, size_substring, nb_children);

        return nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;
    }

    return 0;
}*/
