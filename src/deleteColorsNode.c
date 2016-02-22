#include "./../lib/deleteColorsNode.h"

/*int deleteColors_simplePath(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array, int shifting_suffix, int id_genome,
                            uint16_t** skip_node_root, info_per_level*  info_per_lvl, annotation_inform* ann_inf, resultPresence* res, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(root,"deleteColors_simplePath()")
    ASSERT_NULL_PTR(kmer_start,"deleteColors_simplePath()")
    ASSERT_NULL_PTR(kmer_start_tmp,"deleteColors_simplePath()")

    UC_SIZE_ANNOT_T last_id_genome;

    int size_annot_res_tmp_right, annot_present;
    int size_annot_res = 0;
    int count = 0;
    int size_annot = 0;
    int size_annot_cplx = 0;
    int i, j;
    int tot, ceil, inc_cplx;

    uint8_t q, z;
    uint8_t flag = 0xff, flag_tmp;

    uint8_t* annot_res = NULL;
    uint8_t* annot_res_tmp_right = NULL;

    uint8_t* annot = NULL;
    uint8_t* annot_ext = NULL;
    uint8_t* annot_cplx = NULL;

    uint8_t kmer_retrieve[size_kmer_array];

    resultPresence* res_tmp_right = NULL;

    annotation_array_elem* annot_array_elem;

    UC* uc;
    UC* uc_tmp_right;

    for (z=0; z<4; z++){ //for each possible right neighbor

        if (res[z].link_child != NULL){ //res[z] is the current neighbor

            if (flag == 0xff){

                kmer_start[size_kmer_array-1] |= z << shifting_suffix;

                uc = (UC*)res[z].container;

                annot_present = get_annot(uc, &annot, &annot_ext, &annot_cplx, &size_annot, &size_annot_cplx,
                                               res[z].posFilter2, res[z].posFilter3, res[z].pos_sub_bucket);

                size_annot_res = 0;

                if (annot_present != 0){

                    annot_res = malloc(MAX(size_annot+1, size_annot_cplx)*sizeof(uint8_t));
                    ASSERT_NULL_PTR(annot_res,"deleteColors_simplePath()")

                    if (size_annot != 0){
                        memcpy(annot_res, annot, size_annot * sizeof(uint8_t));
                        size_annot_res = size_annot;
                    }

                    if ((annot_ext != NULL) && (annot_ext[0] != 0)){
                        memcpy(&(annot_res[size_annot]), annot_ext, sizeof(uint8_t));
                        size_annot_res++;
                    }

                    if (size_annot_cplx != 0){
                        memcpy(annot_res, annot_cplx, size_annot_cplx * sizeof(uint8_t));
                        size_annot_res = size_annot_cplx;
                    }
                }

                if (size_annot_res != 0){

                    if (is_genome_present(ann_inf, annot_sorted, annot_res, size_annot_res, NULL, 0, id_genome) == 0) goto IT_BRANCHING_NODES;
                    else flag = get_mark_UC_4states(uc, res[z].posFilter2, res[z].posFilter3, res[z].pos_sub_bucket);
                }
                else{
                    memcpy(kmer_retrieve, kmer_start, size_kmer_array * sizeof(uint8_t));

                    if (id_genome != 0){
                        annot_array_elem = retrieve_annotation(root, kmer_start, kmer_retrieve, size_kmer_root, size_kmer_array, shifting_suffix, id_genome,
                                                                skip_node_root, info_per_lvl, ann_inf, annot_sorted);
                    }
                    else{
                        annot_array_elem = retrieve_annotation(root, kmer_start, kmer_retrieve, size_kmer_root, size_kmer_array, shifting_suffix, -1,
                                                                skip_node_root, info_per_lvl, ann_inf, annot_sorted);
                    }

                    for (i = annot_array_elem->size_annot - 1; i >= 0; i--){
                        if (annot_array_elem->annot_array[i] != 0){
                            for (j = 7; j >= 0; j--){
                                if ((annot_array_elem->annot_array[i] & MASK_POWER_8[j]) != 0){
                                    last_id_genome = i * SIZE_BITS_UINT_8T + j - 2;
                                    annot_array_elem->annot_array[i] &= ~(MASK_POWER_8[j]);
                                    goto ANNOTATE;
                                }
                            }
                        }
                    }

                    ANNOTATE: compute_best_mode(ann_inf, annot_sorted, annot_array_elem->annot_array, annot_array_elem->size_annot, NULL, 0, last_id_genome);

                    free(annot_res);
                    annot_res = malloc(ann_inf->min_size * sizeof(uint8_t));
                    ASSERT_NULL_PTR(annot_res,"deleteColors_simplePath()")

                    modify_mode_annotation(ann_inf, annot_res, ann_inf->min_size, NULL, 0, last_id_genome);

                    tot = res[z].posFilter3 * (res[z].posFilter2 + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT +
                            uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

                    ceil = CEIL(res[z].posFilter3 * 2, SIZE_BITS_UINT_8T);

                    if (ann_inf->min_size > uc->size_annot_cplx_nodes){

                        inc_cplx = uc->nb_cplx_nodes * (ann_inf->min_size - uc->size_annot_cplx_nodes);

                        uc->suffixes = realloc(uc->suffixes, (tot + inc_cplx + SIZE_BYTE_CPLX_N + ann_inf->min_size + ceil) * sizeof(uint8_t));
                        ASSERT_NULL_PTR(uc->suffixes,"deleteColors_simplePath()")

                        memmove(&(uc->suffixes[tot + inc_cplx + SIZE_BYTE_CPLX_N + ann_inf->min_size]), &(uc->suffixes[tot]), ceil * sizeof(uint8_t));

                        increase_size_annot_cplx_nodes(uc, res[z].posFilter2, res[z].posFilter3, ann_inf->min_size, 0);
                    }
                    else{
                        uc->suffixes = realloc(uc->suffixes, (tot + SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes + ceil) * sizeof(uint8_t));
                        ASSERT_NULL_PTR(uc->suffixes,"deleteColors_simplePath()")

                        memmove(&(uc->suffixes[tot + SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes]), &(uc->suffixes[tot]), ceil * sizeof(uint8_t));
                    }

                    insert_annot_cplx_nodes(uc, res[z].posFilter2, res[z].posFilter3, res[z].pos_sub_bucket, annot_res, ann_inf->min_size, 1);

                    free(annot_array_elem->annot_array);
                    free(annot_array_elem);

                    reinit_annotation_inform(ann_inf, annot_array_elem->size_annot);

                    goto IT_BRANCHING_NODES;
                }
            }

            IT_SPL_NODES: if (flag < 2){ //if the neighbor is not branching

                res_tmp_right = getRightNeighbors(root, kmer_start, size_kmer_root, info_per_lvl, skip_node_root);

                for (q=0; q<4; q++){

                    if (res_tmp_right[q].link_child != NULL){

                        kmer_start[size_kmer_array-1] |= q << shifting_suffix;

                        uc_tmp_right = (UC*)res_tmp_right[q].container;

                        flag_tmp = get_mark_UC_4states(uc_tmp_right, res_tmp_right[q].posFilter2, res_tmp_right[q].posFilter3, res_tmp_right[q].pos_sub_bucket);

                        if (flag_tmp < 2){ //if the right neighbor of the right neighbor is not branching

                            annot_present = get_annot(uc_tmp_right, &annot, &annot_ext, &annot_cplx, &size_annot, &size_annot_cplx,
                               res_tmp_right[q].posFilter2, res_tmp_right[q].posFilter3, res_tmp_right[q].pos_sub_bucket);

                            size_annot_res_tmp_right = 0;

                            if (annot_present != 0){

                                annot_res_tmp_right = malloc(MAX(size_annot+1, size_annot_cplx)*sizeof(uint8_t));
                                ASSERT_NULL_PTR(annot_res_tmp_right,"getBranchingNode()")

                                if (size_annot != 0){
                                    memcpy(annot_res_tmp_right, annot, size_annot * sizeof(uint8_t));
                                    size_annot_res_tmp_right = size_annot;
                                }

                                if ((annot_ext != NULL) && (annot_ext[0] != 0)){
                                    memcpy(&(annot_res_tmp_right[size_annot]), annot_ext, sizeof(uint8_t));
                                    size_annot_res_tmp_right++;
                                }

                                if (size_annot_cplx != 0){
                                    memcpy(annot_res_tmp_right, annot_cplx, size_annot_cplx * sizeof(uint8_t));
                                    size_annot_res_tmp_right = size_annot_cplx;
                                }
                            }

                            if (size_annot_res_tmp_right == 0){
                                memcpy(kmer_retrieve, kmer_start, size_kmer_array * sizeof(uint8_t));

                                if (id_genome != 0){
                                    annot_array_elem = retrieve_annotation(root, kmer_start, kmer_retrieve, size_kmer_root, size_kmer_array, shifting_suffix, id_genome,
                                                                            skip_node_root, info_per_lvl, ann_inf, annot_sorted);
                                }
                                else{
                                    annot_array_elem = retrieve_annotation(root, kmer_start, kmer_retrieve, size_kmer_root, size_kmer_array, shifting_suffix, -1,
                                                                            skip_node_root, info_per_lvl, ann_inf, annot_sorted);
                                }

                                for (i = annot_array_elem->size_annot - 1; i >= 0; i--){
                                    if (annot_array_elem->annot_array[i] != 0){
                                        for (j = 7; j >= 0; j--){
                                            if ((annot_array_elem->annot_array[i] & MASK_POWER_8[j]) != 0){
                                                last_id_genome = i * SIZE_BITS_UINT_8T + j - 2;
                                                annot_array_elem->annot_array[i] &= ~(MASK_POWER_8[j]);
                                                goto ANNOTATE2;
                                            }
                                        }
                                    }
                                }

                                ANNOTATE2: compute_best_mode(ann_inf, annot_sorted, annot_array_elem->annot_array, annot_array_elem->size_annot, NULL, 0, last_id_genome);

                                free(annot_res);
                                annot_res = malloc(ann_inf->min_size * sizeof(uint8_t));
                                ASSERT_NULL_PTR(annot_res,"getBranchingNode()")

                                modify_mode_annotation(ann_inf, annot_res, ann_inf->min_size, NULL, 0, last_id_genome);

                                tot = res_tmp_right[q].posFilter3 * (res_tmp_right[q].posFilter2 + uc_tmp_right->size_annot) + uc_tmp_right->nb_extended_annot *
                                        SIZE_BYTE_EXT_ANNOT + uc_tmp_right->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc_tmp_right->size_annot_cplx_nodes);

                                ceil = CEIL(res_tmp_right[q].posFilter3 * 2, SIZE_BITS_UINT_8T);

                                if (ann_inf->min_size > uc_tmp_right->size_annot_cplx_nodes){

                                    inc_cplx = uc_tmp_right->nb_cplx_nodes * (ann_inf->min_size - uc_tmp_right->size_annot_cplx_nodes);

                                    uc_tmp_right->suffixes = realloc(uc_tmp_right->suffixes, (tot + inc_cplx + SIZE_BYTE_CPLX_N + ann_inf->min_size + ceil) * sizeof(uint8_t));
                                    ASSERT_NULL_PTR(uc_tmp_right->suffixes,"deleteColors_simplePath()")

                                    memmove(&(uc_tmp_right->suffixes[tot + inc_cplx + SIZE_BYTE_CPLX_N + ann_inf->min_size]),
                                            &(uc_tmp_right->suffixes[tot]),
                                            ceil * sizeof(uint8_t));

                                    increase_size_annot_cplx_nodes(uc_tmp_right, res_tmp_right[q].posFilter2, res_tmp_right[q].posFilter3, ann_inf->min_size, 0);

                                    reinit_annotation_inform(ann_inf, annot_array_elem->size_annot);
                                }
                                else{
                                    uc_tmp_right->suffixes = realloc(uc_tmp_right->suffixes, (tot + SIZE_BYTE_CPLX_N + uc_tmp_right->size_annot_cplx_nodes +
                                                                     ceil) * sizeof(uint8_t));
                                    ASSERT_NULL_PTR(uc_tmp_right->suffixes,"deleteColors_simplePath()")

                                    memmove(&(uc_tmp_right->suffixes[tot + SIZE_BYTE_CPLX_N + uc_tmp_right->size_annot_cplx_nodes]),
                                            &(uc_tmp_right->suffixes[tot]),
                                            ceil * sizeof(uint8_t));
                                }

                                insert_annot_cplx_nodes(uc_tmp_right, res_tmp_right[q].posFilter2, res_tmp_right[q].posFilter3, res_tmp_right[q].pos_sub_bucket,
                                                        annot_res, ann_inf->min_size, 1);

                                kmer_start[size_kmer_array-1] &= ~(((uint8_t)0x3) << shifting_suffix);

                                free(annot_array_elem->annot_array);
                                free(annot_array_elem);

                                uc_tmp_right = NULL;

                                reinit_annotation_inform(ann_inf, annot_array_elem->size_annot);

                                goto WRONG_NODE;
                            }
                            else if (is_genome_present(ann_inf, annot_sorted, annot_res_tmp_right, size_annot_res_tmp_right, NULL, 0, id_genome) == 0){

                                kmer_start[size_kmer_array-1] &= ~(((uint8_t)0x3) << shifting_suffix);

                                uc_tmp_right = NULL;

                                reinit_annotation_inform(ann_inf, annot_array_elem->size_annot);

                                goto WRONG_NODE;
                            }

                            if (size_annot_res == size_annot_res_tmp_right){ //If sizes are different, annotations have to be different

                                if (memcmp(annot_res, annot_res_tmp_right, size_annot_res*sizeof(uint8_t)) == 0){ //If annotations are the same

                                    //Erase the annotation
                                    mark_UC_4states(uc, res[z].posFilter2, res[z].posFilter3, res[z].pos_sub_bucket, 1);

                                    count++;
                                }
                            }

                            //Current annotation becomes the one of the right neighbor
                            memcpy(&(res[z]), &(res_tmp_right[q]), sizeof(resultPresence));

                            flag = flag_tmp;

                            uc = uc_tmp_right;

                            size_annot_res = size_annot_res_tmp_right;

                            memcpy(annot_res, annot_res_tmp_right, size_annot_res);

                            free(annot_res_tmp_right);
                            free(res_tmp_right);

                            annot_res_tmp_right = NULL;
                            res_tmp_right = NULL;

                            uc_tmp_right = NULL;

                            goto IT_SPL_NODES;
                        }
                        else goto IT_BRANCHING_NODES;
                    }

                    WRONG_NODE: continue;
                }
            }
            else goto IT_BRANCHING_NODES;
        }

        IT_BRANCHING_NODES: flag = 0xff;

        if (res_tmp_right != NULL){
            free(res_tmp_right);
            res_tmp_right = NULL;
        }

        if (annot_res != NULL){
            free(annot_res);
            annot_res = NULL;
        }

        if (annot_res_tmp_right != NULL){
            free(annot_res_tmp_right);
            annot_res_tmp_right = NULL;
        }

        memcpy(kmer_start, kmer_start_tmp, size_kmer_array*sizeof(uint8_t));
    }

    initialize_resultPresence(res);

    return count;
}

int deleteColors_from_branchingNodes(Node* n, Node* root, uint8_t* kmer, int size_kmer, int bucket, int pos_in_bucket, int size_kmer_root, uint32_t id_genome,
                                     info_per_level*  info_per_lvl, uint16_t** skip_node_root, annotation_inform* ann_inf, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(n,"deleteColors_from_branchingNodes()")
    ASSERT_NULL_PTR(root,"deleteColors_from_branchingNodes()")
    ASSERT_NULL_PTR(kmer,"deleteColors_from_branchingNodes()")
    ASSERT_NULL_PTR(ann_inf,"deleteColors_from_branchingNodes()")
    ASSERT_NULL_PTR(skip_node_root,"deleteColors_from_branchingNodes()")
    ASSERT_NULL_PTR(info_per_lvl,"deleteColors_from_branchingNodes()")

    CC* cc;
    UC* uc;
    resultPresence* res;

    int i = -1, j = 0, k = 0, count = 0, level = (size_kmer/SIZE_BITS_UINT_8T)-1, size_kmer_array = CEIL(size_kmer_root*2,SIZE_BITS_UINT_8T);

    uint16_t size_bf, nb_elt, it_filter2;

    uint64_t new_substring;

    int it_filter3, first_bit, it_bucket, last_shift, last_it_children_bucket, nb_cell_children, shifting_UC, tot;
    int it_children_pos_bucket, it_children_bucket, it_node, it_substring, size_line, size_new_substring, size_new_substring_bytes;

    int shifting_suffix = SIZE_BITS_UINT_8T - (size_kmer_array*SIZE_BITS_UINT_8T - (size_kmer_root-1)*2);
    int shifting1 = (NB_CHAR_SUF_PREF*2)-SIZE_BITS_UINT_8T+pos_in_bucket;
    int shifting2 = shifting1-SIZE_BITS_UINT_8T;
    int shifting3 = SIZE_BITS_UINT_8T-shifting2;

    uint8_t mask = ~(MASK_POWER_8[shifting3]-1);

    uint8_t s, p;

    uint8_t kmer_tmp[size_kmer_array];
    uint8_t kmer_start[size_kmer_array];
    uint8_t kmer_start_tmp[size_kmer_array];

    if ((CC*)n->CC_array != NULL){
        do {
            i++;
            cc = &(((CC*)n->CC_array)[i]);

            s = (cc->type >> 1) & 0x1f;
            p = NB_CHAR_SUF_PREF*2-s;

            it_filter2 = 0;
            it_filter3 = 0;
            it_bucket = 0;
            it_node = 0;
            it_children_pos_bucket = 0;
            it_children_bucket = 0;
            it_substring = 0;

            size_bf = cc->type >> 7;
            last_it_children_bucket = -1;

            if (s==8){
                if (info_per_lvl[level].level_min == 1){
                    for (it_filter2=0; it_filter2<MASK_POWER_16[p]; it_filter2++){
                        if ((cc->BF_filter2[size_bf+it_filter2/SIZE_BITS_UINT_8T] & (MASK_POWER_8[it_filter2%SIZE_BITS_UINT_8T])) != 0){

                            first_bit = 1;

                            while((it_filter3 < cc->nb_elem) && (((cc->extra_filter3[it_filter3/SIZE_BITS_UINT_8T] & MASK_POWER_8[it_filter3%SIZE_BITS_UINT_8T]) == 0) || (first_bit == 1))){

                                memcpy(kmer_tmp, kmer, size_kmer_array*sizeof(uint8_t));

                                new_substring = (it_filter2 << SIZE_BITS_UINT_8T) | cc->filter3[it_filter3];
                                new_substring = (new_substring >> 2) | ((new_substring & 0x3) << 16);
                                kmer_tmp[bucket] |= new_substring >> shifting1;
                                kmer_tmp[bucket+1] = new_substring >> shifting2;
                                kmer_tmp[bucket+2] = new_substring << shifting3;
                                it_bucket = bucket+2;
                                if (shifting3 == 0) it_bucket++;

                                if (size_kmer != NB_CHAR_SUF_PREF){

                                    if ((nb_elt = getNbElts(cc, it_filter3)) != 0){
                                        it_children_bucket = it_filter3/NB_UC_PER_SKP;

                                        if (it_children_bucket != last_it_children_bucket){

                                            it_children_pos_bucket = 0;
                                            uc = &(((UC*)cc->children)[it_children_bucket]);
                                            size_line = info_per_lvl[level].size_kmer_in_bytes_minus_1 + uc->size_annot;
                                            last_it_children_bucket = it_children_bucket;

                                            tot = uc->nb_children * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                                        + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);
                                        }

                                        for (j=it_children_pos_bucket*size_line; j<(it_children_pos_bucket+nb_elt)*size_line; j+=size_line){

                                            if (get_mark_UC_4states_bis(uc, j/size_line, tot) >= 2){

                                                extractSuffix(kmer_tmp, size_kmer, size_kmer_array, shifting3, it_bucket, &(uc->suffixes[j]), &(info_per_lvl[level]));

                                                memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));
                                                res = getRightNeighbors(root, kmer_start, size_kmer_root, info_per_lvl, skip_node_root);
                                                memcpy(kmer_start_tmp, kmer_start, size_kmer_array*sizeof(uint8_t));

                                                count += deleteColors_simplePath(root, kmer_start, kmer_start_tmp, size_kmer_root, size_kmer_array, shifting_suffix,
                                                                                 id_genome, skip_node_root, info_per_lvl, ann_inf, res, annot_sorted);

                                                free(res);

                                                for (k=0; k<bucket+3; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);
                                                if (shifting3 != 0) kmer_tmp[bucket+2] &= mask;
                                                memset(&(kmer_tmp[bucket+3]), 0, (size_kmer_array-bucket-3)*sizeof(uint8_t));
                                            }
                                        }

                                        it_children_pos_bucket += nb_elt;
                                    }
                                    else{
                                        count += deleteColors_from_branchingNodes(&(cc->children_Node_container[it_node]), root, kmer_tmp, size_kmer-NB_CHAR_SUF_PREF, it_bucket,
                                                                                  shifting2, size_kmer_root, id_genome, info_per_lvl, skip_node_root, ann_inf, annot_sorted);
                                        it_node++;
                                    }
                                }
                                else{
                                    if ((isBranchingRight(root, kmer_tmp, size_kmer_root, info_per_lvl, skip_node_root) > 1) ||
                                        (isBranchingLeft(root, kmer_tmp, size_kmer_root, info_per_lvl, skip_node_root) > 1)){
                                            //count++;
                                            memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));
                                        }
                                }

                                it_filter3++;
                                first_bit=0;
                            }
                        }
                    }
                }
                else{
                    int cpt_pv = 0;
                    int nb_skp = CEIL(cc->nb_elem,NB_UC_PER_SKP);

                    nb_cell_children = info_per_lvl[level].size_kmer_in_bytes_minus_1-1;

                    for (it_filter2=0; it_filter2<MASK_POWER_16[p]; it_filter2++){
                        if ((cc->BF_filter2[size_bf+it_filter2/SIZE_BITS_UINT_8T] & (MASK_POWER_8[it_filter2%SIZE_BITS_UINT_8T])) != 0){

                            first_bit = 1;

                            while (it_children_bucket < nb_skp){
                                uc = &(((UC*)cc->children)[it_children_bucket]);
                                size_line = info_per_lvl[level].size_kmer_in_bytes_minus_1 + uc->size_annot;

                                tot = uc->nb_children * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                            + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

                                while (cpt_pv < uc->nb_children){

                                    memcpy(kmer_tmp, kmer, size_kmer_array*sizeof(uint8_t));

                                    new_substring = (it_filter2 << SIZE_BITS_UINT_8T) | cc->filter3[it_filter3];
                                    new_substring = (new_substring >> 2) | ((new_substring & 0x3) << 16);
                                    kmer_tmp[bucket] |= new_substring >> shifting1;
                                    kmer_tmp[bucket+1] = new_substring >> shifting2;
                                    kmer_tmp[bucket+2] = new_substring << shifting3;
                                    it_bucket = bucket+2;
                                    if (shifting3 == 0) it_bucket++;

                                    if ((nb_elt = getNbElts(cc, it_filter3)) == 0){

                                        if (((cc->children_Node_container[it_node].UC_array.nb_children & 0x1) == 0) || (first_bit == 1)){

                                            first_bit=0;
                                            count += deleteColors_from_branchingNodes(&(cc->children_Node_container[it_node]), root, kmer_tmp, size_kmer-NB_CHAR_SUF_PREF, it_bucket,
                                                                                      shifting2, size_kmer_root, id_genome, info_per_lvl, skip_node_root, ann_inf, annot_sorted);
                                            it_node++;
                                        }
                                        else goto OUT_LOOP_S8;
                                    }
                                    else{
                                        if ((uc->suffixes[cpt_pv*size_line+nb_cell_children] < 0x80)  || (first_bit == 1)){

                                            first_bit=0;

                                            for (j=cpt_pv*size_line; j<(cpt_pv+nb_elt)*size_line; j+=size_line){

                                                if (get_mark_UC_4states_bis(uc, j/size_line, tot) >= 2){

                                                    extractSuffix(kmer_tmp, size_kmer, size_kmer_array, shifting3, it_bucket, &(uc->suffixes[j]), &(info_per_lvl[level]));
                                                    kmer_tmp[size_kmer_array-1] &= 0x7f;
                                                    memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));
                                                    res = getRightNeighbors(root, kmer_start, size_kmer_root, info_per_lvl, skip_node_root);
                                                    memcpy(kmer_start_tmp, kmer_start, size_kmer_array*sizeof(uint8_t));

                                                    count += deleteColors_simplePath(root, kmer_start, kmer_start_tmp, size_kmer_root, size_kmer_array, shifting_suffix,
                                                                                     id_genome, skip_node_root, info_per_lvl, ann_inf, res, annot_sorted);

                                                    free(res);

                                                    for (k=0; k<bucket+3; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);
                                                    if (shifting3 != 0) kmer_tmp[bucket+2] &= mask;
                                                    memset(&(kmer_tmp[bucket+3]), 0, (size_kmer_array-bucket-3)*sizeof(uint8_t));
                                                }
                                            }

                                            cpt_pv += nb_elt;
                                        }
                                        else goto OUT_LOOP_S8;
                                    }

                                    it_filter3++;
                                }

                                cpt_pv = 0;
                                it_children_bucket++;
                            }
                        }

                        OUT_LOOP_S8: continue;
                    }
                }
            }
            else {
                if (info_per_lvl[level].level_min == 1){
                    for (it_filter2=0; it_filter2<MASK_POWER_16[p]; it_filter2++){
                        if ((cc->BF_filter2[size_bf+it_filter2/SIZE_BITS_UINT_8T] & (MASK_POWER_8[it_filter2%SIZE_BITS_UINT_8T])) != 0){

                            first_bit = 1;

                            while((it_filter3 < cc->nb_elem) && (((cc->extra_filter3[it_filter3/SIZE_BITS_UINT_8T] & MASK_POWER_8[it_filter3%SIZE_BITS_UINT_8T]) == 0) || (first_bit == 1))){

                                memcpy(kmer_tmp, kmer, size_kmer_array*sizeof(uint8_t));

                                if (IS_ODD(it_filter3)) new_substring = (it_filter2 << 4) | (cc->filter3[it_filter3/2] >> 4);
                                else new_substring = (it_filter2 << 4) | (cc->filter3[it_filter3/2] & 0xf);

                                new_substring = (new_substring >> 2) | ((new_substring & 0x3) << 16);
                                kmer_tmp[bucket] |= new_substring >> shifting1;
                                kmer_tmp[bucket+1] = new_substring >> shifting2;
                                kmer_tmp[bucket+2] = new_substring << shifting3;
                                it_bucket = bucket+2;
                                if (shifting3 == 0) it_bucket++;

                                if (size_kmer != NB_CHAR_SUF_PREF){

                                    if ((nb_elt = getNbElts(cc, it_filter3)) != 0){
                                        it_children_bucket = it_filter3/NB_UC_PER_SKP;

                                        if (it_children_bucket != last_it_children_bucket){

                                            it_children_pos_bucket = 0;
                                            uc = &(((UC*)cc->children)[it_children_bucket]);

                                            size_line = info_per_lvl[level].size_kmer_in_bytes_minus_1 + uc->size_annot;
                                            last_it_children_bucket = it_children_bucket;

                                            tot = uc->nb_children * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                                        + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);
                                        }

                                        for (j=it_children_pos_bucket*size_line; j<(it_children_pos_bucket+nb_elt)*size_line; j+=size_line){

                                            if (get_mark_UC_4states_bis(uc, j/size_line, tot) >= 2){

                                                extractSuffix(kmer_tmp, size_kmer, size_kmer_array, shifting3, it_bucket, &(uc->suffixes[j]), &(info_per_lvl[level]));
                                                memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));
                                                res = getRightNeighbors(root, kmer_start, size_kmer_root, info_per_lvl, skip_node_root);
                                                memcpy(kmer_start_tmp, kmer_start, size_kmer_array*sizeof(uint8_t));

                                                count += deleteColors_simplePath(root, kmer_start, kmer_start_tmp, size_kmer_root, size_kmer_array, shifting_suffix,
                                                                                 id_genome, skip_node_root, info_per_lvl, ann_inf, res, annot_sorted);

                                                free(res);

                                                for (k=0; k<bucket+3; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);
                                                if (shifting3 != 0) kmer_tmp[bucket+2] &= mask;
                                                memset(&(kmer_tmp[bucket+3]), 0, (size_kmer_array-bucket-3)*sizeof(uint8_t));
                                            }
                                        }

                                        it_children_pos_bucket += nb_elt;
                                    }
                                    else{
                                        count += deleteColors_from_branchingNodes(&(cc->children_Node_container[it_node]), root, kmer_tmp, size_kmer-NB_CHAR_SUF_PREF, it_bucket,
                                                                                  shifting2, size_kmer_root, id_genome, info_per_lvl, skip_node_root, ann_inf, annot_sorted);
                                        it_node++;
                                    }
                                }
                                else{
                                    if ((isBranchingRight(root, kmer_tmp, size_kmer_root, info_per_lvl, skip_node_root) > 1) ||
                                        (isBranchingLeft(root, kmer_tmp, size_kmer_root, info_per_lvl, skip_node_root) > 1)){
                                            //count++;
                                            memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));
                                        }
                                }

                                it_filter3++;
                                first_bit=0;
                            }
                        }
                    }
                }
                else{
                    int cpt_pv = 0;
                    int nb_skp = CEIL(cc->nb_elem,NB_UC_PER_SKP);

                    nb_cell_children = info_per_lvl[level].size_kmer_in_bytes_minus_1-1;

                    for (it_filter2=0; it_filter2<MASK_POWER_16[p]; it_filter2++){
                        if ((cc->BF_filter2[size_bf+it_filter2/SIZE_BITS_UINT_8T] & (MASK_POWER_8[it_filter2%SIZE_BITS_UINT_8T])) != 0){

                            first_bit = 1;

                            while (it_children_bucket < nb_skp){
                                uc = &(((UC*)cc->children)[it_children_bucket]);
                                size_line = info_per_lvl[level].size_kmer_in_bytes_minus_1 + uc->size_annot;

                                tot = uc->nb_children * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                            + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

                                while (cpt_pv < uc->nb_children){

                                    memcpy(kmer_tmp, kmer, size_kmer_array*sizeof(uint8_t));

                                    if (IS_ODD(it_filter3)) new_substring = (it_filter2 << 4) | (cc->filter3[it_filter3/2] >> 4);
                                    else new_substring = (it_filter2 << 4) | (cc->filter3[it_filter3/2] & 0xf);

                                    new_substring = (new_substring >> 2) | ((new_substring & 0x3) << 16);
                                    kmer_tmp[bucket] |= new_substring >> shifting1;
                                    kmer_tmp[bucket+1] = new_substring >> shifting2;
                                    kmer_tmp[bucket+2] = new_substring << shifting3;
                                    it_bucket = bucket+2;
                                    if (shifting3 == 0) it_bucket++;

                                    if ((nb_elt = getNbElts(cc, it_filter3)) == 0){

                                        if (((cc->children_Node_container[it_node].UC_array.nb_children & 0x1) == 0) || (first_bit == 1)){

                                            first_bit=0;
                                            count += deleteColors_from_branchingNodes(&(cc->children_Node_container[it_node]), root, kmer_tmp, size_kmer-NB_CHAR_SUF_PREF, it_bucket,
                                                                                      shifting2, size_kmer_root, id_genome, info_per_lvl, skip_node_root, ann_inf, annot_sorted);
                                            it_node++;
                                        }
                                        else goto OUT_LOOP_S4;
                                    }
                                    else{
                                        if ((uc->suffixes[cpt_pv*size_line+nb_cell_children] < 0x80)  || (first_bit == 1)){

                                            first_bit=0;

                                            for (j=cpt_pv*size_line; j<(cpt_pv+nb_elt)*size_line; j+=size_line){

                                                if (get_mark_UC_4states_bis(uc, j/size_line, tot) >= 2){

                                                    extractSuffix(kmer_tmp, size_kmer, size_kmer_array, shifting3, it_bucket, &(uc->suffixes[j]), &(info_per_lvl[level]));
                                                    kmer_tmp[size_kmer_array-1] &= 0x7f;

                                                    memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));
                                                    res = getRightNeighbors(root, kmer_start, size_kmer_root, info_per_lvl, skip_node_root);
                                                    memcpy(kmer_start_tmp, kmer_start, size_kmer_array*sizeof(uint8_t));

                                                    count += deleteColors_simplePath(root, kmer_start, kmer_start_tmp, size_kmer_root, size_kmer_array, shifting_suffix,
                                                                                     id_genome, skip_node_root, info_per_lvl, ann_inf, res, annot_sorted);

                                                    free(res);

                                                    for (k=0; k<bucket+3; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);
                                                    if (shifting3 != 0) kmer_tmp[bucket+2] &= mask;
                                                    memset(&(kmer_tmp[bucket+3]), 0, (size_kmer_array-bucket-3)*sizeof(uint8_t));
                                                }
                                            }

                                            cpt_pv += nb_elt;
                                        }
                                        else goto OUT_LOOP_S4;
                                    }

                                    it_filter3++;
                                }

                                cpt_pv = 0;
                                it_children_bucket++;
                            }
                        }

                        OUT_LOOP_S4: continue;
                    }
                }
            }
        }
        while ((((CC*)n->CC_array)[i].type & 0x1) == 0);
    }

    if (n->UC_array.suffixes != NULL){

        size_line = info_per_lvl[level].size_kmer_in_bytes + n->UC_array.size_annot;
        nb_elt = n->UC_array.nb_children >> 1;

        tot = nb_elt * size_line + n->UC_array.nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                    + n->UC_array.nb_cplx_nodes * (SIZE_BYTE_CPLX_N + n->UC_array.size_annot_cplx_nodes);

        for (j=0; j<nb_elt * size_line; j += size_line){

            if (get_mark_UC_4states_bis(&(n->UC_array), j/size_line, tot) >= 2){

                it_substring = 0;
                it_bucket = bucket;

                memcpy(kmer_tmp, kmer, size_kmer_array*sizeof(uint8_t));

                while (it_substring < info_per_lvl[level].size_kmer_in_bytes){

                    it_substring += sizeof(uint64_t);
                    new_substring = 0;

                    if (it_substring > info_per_lvl[level].size_kmer_in_bytes){

                        size_new_substring = size_kmer*2-((it_substring-sizeof(uint64_t))*SIZE_BITS_UINT_8T);
                        size_new_substring_bytes = CEIL(size_new_substring, SIZE_BITS_UINT_8T);
                        for (k=0; k<size_new_substring_bytes; k++) new_substring = (new_substring << 8) | reverse_word_8(n->UC_array.suffixes[j+(it_substring-sizeof(uint64_t))+k]);
                        new_substring >>= info_per_lvl[level].size_kmer_in_bytes*SIZE_BITS_UINT_8T - size_new_substring;
                    }
                    else{

                        size_new_substring = sizeof(uint64_t)*SIZE_BITS_UINT_8T;
                        size_new_substring_bytes = sizeof(uint64_t);
                        for (k=0; k<size_new_substring_bytes; k++) new_substring = (new_substring << 8) | reverse_word_8(n->UC_array.suffixes[j+(it_substring-sizeof(uint64_t))+k]);
                    }

                    shifting_UC = SIZE_BITS_UINT_8T-pos_in_bucket;

                    for (k=it_bucket; k<it_bucket+size_new_substring_bytes; k++){

                        last_shift = size_new_substring-shifting_UC;

                        if (last_shift >= 0) kmer_tmp[k] |= new_substring >> last_shift;
                        else kmer_tmp[k] |= new_substring << abs(last_shift);

                        shifting_UC += SIZE_BITS_UINT_8T;
                    }

                    it_bucket+=size_new_substring_bytes;
                }

                for (k=0; k<size_kmer_array; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);

                memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));
                res = getRightNeighbors(root, kmer_start, size_kmer_root, info_per_lvl, skip_node_root);
                memcpy(kmer_start_tmp, kmer_start, size_kmer_array*sizeof(uint8_t));

                count += deleteColors_simplePath(root, kmer_start, kmer_start_tmp, size_kmer_root, size_kmer_array, shifting_suffix, id_genome,
                                                 skip_node_root, info_per_lvl, ann_inf, res, annot_sorted);

                free(res);
            }
        }
    }

    return count;
}

int resize_annotation_Node(Node* n, int size_kmer, info_per_level*  info_per_lvl){

    ASSERT_NULL_PTR(n,"resize_annotation_Node()")

    int count = 0;
    int level = (size_kmer/NB_CHAR_SUF_PREF)-1;

    if (n->CC_array != NULL){

        int i = -1;

        do {
            i++;
            count += resize_annotation_CC(&(((CC*)n->CC_array)[i]), size_kmer, info_per_lvl);
        }
        while ((((CC*)n->CC_array)[i].type & 0x1) == 0);
    }

    if (n->UC_array.suffixes != NULL) count += resize_annotation_UC(&(n->UC_array), info_per_lvl[level].size_kmer_in_bytes, n->UC_array.nb_children >> 1);

    return count;
}

int resize_annotation_CC(CC* cc, int size_kmer, info_per_level*  info_per_lvl){

    ASSERT_NULL_PTR(cc,"resize_annotation_CC()")

    int i = 0;
    int level = (size_kmer/NB_CHAR_SUF_PREF)-1;
    int count = 0;

    UC* uc;


    for (i=0; i< CEIL(cc->nb_elem, NB_UC_PER_SKP); i++){

        uc = &(((UC*)cc->children)[i]);

        if (size_kmer != NB_CHAR_SUF_PREF) count += resize_annotation_UC(uc, info_per_lvl[level].size_kmer_in_bytes_minus_1, uc->nb_children);
        else count += resize_annotation_UC(uc, 0, MIN(NB_UC_PER_SKP, cc->nb_elem - i *NB_UC_PER_SKP));
    }

    for (i=0; i<cc->nb_Node_children; i++) count += resize_annotation_Node(&(cc->children_Node_container[i]), size_kmer-NB_CHAR_SUF_PREF, info_per_lvl);

    return count;
}

int resize_annotation_UC(UC* uc, int size_substring, int nb_children){

    ASSERT_NULL_PTR(uc,"resize_annotation_UC()")

    int tot;
    int flag;
    int size;
    int size_line;
    int new_possible_size;
    int old_nb_extended_annot;

    int max = 0;
    int max_tmp = uc->size_annot_cplx_nodes;
    int count = 0;
    int nb_annot_suf = 0;

    int* new_order;

    uint8_t* annot;
    uint8_t* marked;
    uint8_t* min_size;
    uint8_t* sizes_annot;
    uint8_t* sizes_plus_positions;

    uint8_t** extended_annot;

    uint16_t j;
    uint16_t pos;

    if (uc->suffixes != NULL){

        size_line = size_substring + uc->size_annot;

        min_size = min_size_per_sub(uc->suffixes, nb_children, size_substring, uc->size_annot);
        extended_annot = get_extend_annots(uc, size_substring, nb_children, 0, nb_children-1);
        marked = calloc(nb_children, sizeof(uint8_t));
        sizes_plus_positions = calloc(nb_children * 3, sizeof(uint8_t));
        sizes_annot = calloc(SIZE_MAX_BYTE_ANNOT, sizeof(uint8_t));

        tot = nb_children * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
            + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

        for (j = 0; j < nb_children; j++){

            flag = get_mark_UC_4states_bis(uc, j, tot);

            if (IS_ODD(flag)){

                nb_annot_suf++;

                marked[j] = 1;

                sizes_plus_positions[j * 3 + 1] = j >> 8;
                sizes_plus_positions[j * 3 + 2] = j & 0xff;

                if ((extended_annot != NULL) && (extended_annot[j] != NULL) && (extended_annot[j][0] != 0)){
                    max_tmp = MAX(max_tmp, min_size[j]+1);
                    sizes_plus_positions[j * 3] = min_size[j] + 1;
                    sizes_annot[min_size[j] + 1]++;
                }
                else{
                    max_tmp = MAX(max_tmp, min_size[j]);
                    sizes_plus_positions[j * 3] = min_size[j];
                    sizes_annot[min_size[j]]++;
                }
            }
        }

        new_order = quicksort_init(sizes_plus_positions, 3, 0, nb_children-1);

        //Best possible size: no annotations with kmers, no extended annot, necessary annotations are cplx
        new_possible_size = nb_children * size_substring + (uc->nb_cplx_nodes + nb_annot_suf) * (SIZE_BYTE_CPLX_N + max_tmp);

        j = 0;
        pos = 0;

        //While the best possible size is bigger than the actual size, we iterate over the bigger annotations first
        while ((new_possible_size >= tot) && (j < nb_children)){

            pos = (((uint16_t)sizes_plus_positions[j * 3 + 1]) << SIZE_BITS_UINT_8T) | (((uint16_t)sizes_plus_positions[j * 3 + 2]) | 0xff);

            if (marked[pos] == 1){

                size = sizes_plus_positions[j];

                nb_annot_suf--;
                sizes_annot[size]--;

                new_possible_size += nb_children * MAX(size - max, 0) - SIZE_BYTE_CPLX_N - max_tmp;

                if (sizes_annot[size] == 0){

                    while ((sizes_annot[size] == 0) && (size > 0)) size--;
                    new_possible_size -= nb_annot_suf * (max_tmp - size);
                    max_tmp = size;
                }

                max = MAX(max, size);
                marked[pos] = 0;
            }

            j++;
        }

        if (new_possible_size < tot){

            count += nb_annot_suf;

            if (max_tmp > uc->size_annot_cplx_nodes){
                uc->suffixes = realloc(uc->suffixes, (nb_children * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                        + (uc->nb_cplx_nodes + nb_annot_suf) * (SIZE_BYTE_CPLX_N + max_tmp)) * sizeof(uint8_t));
                ASSERT_NULL_PTR(uc->suffixes, "resize_annotation_UC()")

                increase_size_annot_cplx_nodes(uc, size_substring, nb_children, max_tmp, 0);
            }
            else{
                uc->suffixes = realloc(uc->suffixes, (nb_children * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                        + (uc->nb_cplx_nodes + nb_annot_suf) * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
                ASSERT_NULL_PTR(uc->suffixes, "resize_annotation_UC()")
            }

            annot = malloc(max * sizeof(uint8_t));
            ASSERT_NULL_PTR(annot, "resize_annotation_UC()")

            for (j = 0; j < nb_children; j++){

                if (marked[j] == 1){

                    memcpy(annot, &(uc->suffixes[j * size_line + size_substring]), min_size[j] * sizeof(uint8_t));

                    if ((extended_annot != NULL) && (extended_annot[j] != NULL) && (extended_annot[j][0] != 0)){

                        annot[min_size[j]] = extended_annot[j][0];
                        insert_annot_cplx_nodes(uc, size_substring, nb_children, j, annot, min_size[j]+1, 1);
                    }
                    else insert_annot_cplx_nodes(uc, size_substring, nb_children, j, annot, min_size[j], 1);

                    memset(&(uc->suffixes[j * size_line + size_substring]), 0, uc->size_annot * sizeof(uint8_t));
                }

                memmove(&(uc->suffixes[j * (size_substring + max)]), &(uc->suffixes[j * size_line]), (size_substring + max) * sizeof(uint8_t));
            }

            memmove(&(uc->suffixes[nb_children * (size_substring + max)]),
                    &(uc->suffixes[nb_children * size_line]),
                    (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));

            uc->size_annot = max;
            size_line = size_substring + uc->size_annot;
            old_nb_extended_annot = uc->nb_extended_annot;

            for (j = 0; j < nb_children; j++){
                if ((marked[j] == 1) && (extended_annot != NULL) && (extended_annot[j] != NULL))
                    delete_extend_annots(uc, size_substring, nb_children, j, j, 0, 0, 0);
            }

            memmove(&(uc->suffixes[nb_children * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
                    &(uc->suffixes[nb_children * size_line + old_nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
                    uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes) * sizeof(uint8_t));

            uc->suffixes = realloc(uc->suffixes, (nb_children * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                    + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes) * sizeof(uint8_t)));

            free(annot);
        }
        else delete_marking_UC_4states(uc, size_substring, nb_children);

        if (extended_annot != NULL) free(extended_annot);
        free(min_size);
        free(marked);
        free(sizes_plus_positions);
        free(sizes_annot);
        free(new_order);
    }

    return count;
}*/
