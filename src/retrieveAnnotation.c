#include "retrieveAnnotation.h"

/*annotation_array_elem* retrieve_annotation_right(BFT_Root* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array,
                                                 int shifting_suffix, uint32_t id_genome_avoid, uint16_t** skip_node_root,
                                                 info_per_level*  info_per_lvl, annotation_inform* ann_inf, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(root,"retrieve_annotation_right()")
    ASSERT_NULL_PTR(kmer_start,"retrieve_annotation_right()")
    ASSERT_NULL_PTR(kmer_start_tmp,"retrieve_annotation_right()")

    int annot_present;
    int size_annot = 0;
    int size_annot_cplx = 0;
    int lvl_root = (root->k / NB_CHAR_SUF_PREF) - 1;

    uint8_t* kmer = NULL;
    uint8_t* annot_res = NULL;
    uint8_t* annot = NULL;
    uint8_t* annot_ext = NULL;
    uint8_t* annot_cplx = NULL;

    uint8_t z;

    annotation_array_elem* ann_arr_elem;

    resultPresence* res = getRightNeighbors(root, lvl_root, kmer_start, size_kmer_root, info_per_lvl, skip_node_root);
    //memcpy(kmer_start, kmer_start_tmp, size_kmer_array*sizeof(uint8_t));

    for (z=0; z<4; z++){ //for each possible right neighbor

        if (res[z].link_child != NULL){ //res[z] is the current neighbor

            kmer_start[size_kmer_array-1] |= z << shifting_suffix;

            annot_present = get_annot((UC*)res[z].container, &annot, &annot_ext, &annot_cplx, &size_annot, &size_annot_cplx,
                                           res[z].posFilter2, res[z].posFilter3, res[z].pos_sub_bucket);

            if (annot_present != 0){

                annot_res = malloc(MAX(size_annot+1, size_annot_cplx)*sizeof(uint8_t));
                ASSERT_NULL_PTR(annot_res,"retrieve_annotation_right()")

                if (size_annot != 0) memcpy(annot_res, annot, size_annot * sizeof(uint8_t));

                if ((annot_ext != NULL) && (annot_ext[0] != 0)){
                    memcpy(&(annot_res[size_annot]), annot_ext, sizeof(uint8_t));
                    size_annot++;
                }

                if (size_annot_cplx != 0){
                    memcpy(annot_res, annot_cplx, size_annot_cplx * sizeof(uint8_t));
                    size_annot = size_annot_cplx;
                }

                if ((id_genome_avoid != -1) && (is_genome_present(ann_inf, annot_sorted, annot_res, size_annot, NULL, 0, id_genome_avoid) == 1)){
                    free(annot_res);
                    kmer_start[size_kmer_array-1] &= ~(z << shifting_suffix);
                    continue;
                }
                else{
                    ann_arr_elem = malloc(sizeof(annotation_array_elem));
                    ASSERT_NULL_PTR(ann_arr_elem,"retrieve_annotation_right()")

                    ann_arr_elem->annot_array = annot_res;
                    ann_arr_elem->size_annot = size_annot;

                    free(res);

                    return ann_arr_elem;
                }
            }
            else{
                kmer = malloc(size_kmer_array * sizeof(uint8_t));
                ASSERT_NULL_PTR(kmer,"retrieve_annotation_right()")

                memcpy(kmer, kmer_start, size_kmer_array*sizeof(uint8_t));

                ann_arr_elem = retrieve_annotation_right(root, kmer, kmer_start, size_kmer_root, size_kmer_array, shifting_suffix, id_genome_avoid,
                                                skip_node_root, info_per_lvl, ann_inf, annot_sorted);

                free(kmer);
                free(res);

                return ann_arr_elem;
            }

            //memcpy(kmer_start, kmer_start_tmp, size_kmer_array*sizeof(uint8_t));
        }
        //else memcpy(kmer_start, kmer_start_tmp, size_kmer_array*sizeof(uint8_t));
    }

    free(res);

    return NULL;
}

annotation_array_elem* retrieve_annotation_left(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root,
                                                int size_kmer_array, uint32_t id_genome_avoid, uint16_t** skip_node_root,
                                                info_per_level*  info_per_lvl, annotation_inform* ann_inf, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(root,"retrieve_annotation_left()")
    ASSERT_NULL_PTR(kmer_start,"retrieve_annotation_left()")
    ASSERT_NULL_PTR(kmer_start_tmp,"retrieve_annotation_left()")

    int annot_present;
    int size_annot = 0;
    int size_annot_cplx = 0;

    uint8_t* kmer = NULL;
    uint8_t* annot_res = NULL;
    uint8_t* annot = NULL;
    uint8_t* annot_ext = NULL;
    uint8_t* annot_cplx = NULL;

    uint8_t z;

    annotation_array_elem* ann_arr_elem;

    resultPresence* res = getLeftNeighbors(root, kmer_start, size_kmer_root, info_per_lvl, skip_node_root);
    //memcpy(kmer_start, kmer_start_tmp, size_kmer_array*sizeof(uint8_t));

    for (z=0; z<4; z++){ //for each possible right neighbor

        if (res[z].link_child != NULL){ //res[z] is the current neighbor

            kmer_start[0] |= z;

            annot_present = get_annot((UC*)res[z].container, &annot, &annot_ext, &annot_cplx, &size_annot, &size_annot_cplx,
                                           res[z].posFilter2, res[z].posFilter3, res[z].pos_sub_bucket);

            if (annot_present != 0){

                annot_res = malloc(MAX(size_annot+1, size_annot_cplx)*sizeof(uint8_t));
                ASSERT_NULL_PTR(annot_res,"retrieve_annotation_left()")

                if (size_annot != 0) memcpy(annot_res, annot, size_annot * sizeof(uint8_t));

                if ((annot_ext != NULL) && (annot_ext[0] != 0)){
                    memcpy(&(annot_res[size_annot]), annot_ext, sizeof(uint8_t));
                    size_annot++;
                }

                if (size_annot_cplx != 0){
                    memcpy(annot_res, annot_cplx, size_annot_cplx * sizeof(uint8_t));
                    size_annot = size_annot_cplx;
                }

                if ((id_genome_avoid != -1) && (is_genome_present(ann_inf, annot_sorted, annot_res, size_annot, NULL, 0, id_genome_avoid) == 1)){
                    free(annot_res);
                    kmer_start[0] &= 0xfc;
                    continue;
                }
                else{
                    ann_arr_elem = malloc(sizeof(annotation_array_elem));
                    ASSERT_NULL_PTR(ann_arr_elem,"retrieve_annotation_right()")

                    ann_arr_elem->annot_array = annot_res;
                    ann_arr_elem->size_annot = size_annot;

                    free(res);

                    return ann_arr_elem;
                }
            }
            else{
                kmer = malloc(size_kmer_array * sizeof(uint8_t));
                ASSERT_NULL_PTR(kmer,"retrieve_annotation_left()")

                memcpy(kmer, kmer_start, size_kmer_array*sizeof(uint8_t));

                ann_arr_elem = retrieve_annotation_left(root, kmer, kmer_start, size_kmer_root, size_kmer_array, id_genome_avoid,
                                                skip_node_root, info_per_lvl, ann_inf, annot_sorted);

                free(kmer);
                free(res);

                return ann_arr_elem;
            }

            //memcpy(kmer_start, kmer_start_tmp, size_kmer_array*sizeof(uint8_t));
        }
        //else memcpy(kmer_start, kmer_start_tmp, size_kmer_array*sizeof(uint8_t));
    }

    free(res);

    ann_arr_elem = malloc(sizeof(annotation_array_elem));
    ASSERT_NULL_PTR(ann_arr_elem,"retrieve_annotation_right()")

    ann_arr_elem->annot_array = NULL;
    ann_arr_elem->size_annot = 0;

    return ann_arr_elem;
}

annotation_array_elem* retrieve_annotation(Node* root, uint8_t* kmer_start, uint8_t* kmer_start_tmp, int size_kmer_root, int size_kmer_array, int shifting_suffix,
                                           uint32_t id_genome_avoid, uint16_t** skip_node_root, info_per_level*  info_per_lvl, annotation_inform* ann_inf,
                                           annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(root,"retrieve_annotation()")
    ASSERT_NULL_PTR(kmer_start,"retrieve_annotation()")
    ASSERT_NULL_PTR(kmer_start_tmp,"retrieve_annotation()")

    annotation_array_elem* annot_left = NULL;
    annotation_array_elem* annot_right = NULL;
    annotation_array_elem* annot_left_right = NULL;

    annot_left = retrieve_annotation_left(root, kmer_start, kmer_start_tmp, size_kmer_root,
                                          size_kmer_array, id_genome_avoid, skip_node_root,
                                          info_per_lvl, ann_inf, annot_sorted);

    memcpy(kmer_start, kmer_start_tmp, size_kmer_array * sizeof(uint8_t));

    annot_right = retrieve_annotation_right(root, kmer_start, kmer_start_tmp, size_kmer_root,
                                            size_kmer_array, shifting_suffix, id_genome_avoid,
                                            skip_node_root, info_per_lvl, ann_inf, annot_sorted);

    memcpy(kmer_start, kmer_start_tmp, size_kmer_array * sizeof(uint8_t));

    annot_left_right = intersection_annotations(annot_left->annot_array, annot_left->size_annot, NULL,
                                                0, annot_right->annot_array, annot_right->size_annot,
                                                NULL, 0, id_genome_avoid-1, annot_sorted);

    if (annot_left->annot_array != NULL) free(annot_left->annot_array);
    if (annot_right->annot_array != NULL) free(annot_right->annot_array);
    free(annot_left);
    free(annot_right);

    return annot_left_right;
}*/

void modify_annotations(BFT_Root* root, UC* uc, int size_substring, int nb_children, int position, uint32_t id_genome,
                        int size_id_genome, uint8_t* kmer, int special_case){

    ASSERT_NULL_PTR(uc,"modify_annotations()")

    int size_annot;
    int size_annot_cplx;
    int res_get_annotation;

    uint8_t* annot = NULL;
    uint8_t* annot_ext = NULL;
    uint8_t* annot_cplx = NULL;

    annotation_array_elem* annot_sorted = root->comp_set_colors;
    annotation_inform* ann_inf = root->ann_inf;

    annotation_array_elem* annot_array_elem = NULL;

    res_get_annotation = get_annot(uc, &annot, &annot_ext, &annot_cplx, &size_annot,
                                        &size_annot_cplx, size_substring, nb_children, position);

    if (res_get_annotation == 1){

        if (size_annot_cplx != 0) compute_best_mode(ann_inf, annot_sorted, annot_cplx,
                                                    size_annot_cplx, NULL, 0, id_genome, size_id_genome);
        else compute_best_mode(ann_inf, annot_sorted, annot, size_annot, annot_ext, 1, id_genome, size_id_genome);
    }
    else{
        /*int size_kmer_array = CEIL(root->k*2, SIZE_BITS_UINT_8T);
        int shifting_suffix = SIZE_BITS_UINT_8T - (size_kmer_array * SIZE_BITS_UINT_8T - (root->k - 1) * 2);
        uint8_t kmer_copy[size_kmer_array];

        memcpy(kmer_copy, kmer, size_kmer_array * sizeof(uint8_t));

        annot_array_elem = retrieve_annotation(&(root->node), kmer_copy, kmer, root->k, size_kmer_array, shifting_suffix,
                                               id_genome, NULL, root->info_per_lvl, ann_inf, annot_sorted);

        annot = annot_array_elem->annot_array;
        size_annot = annot_array_elem->size_annot;

        compute_best_mode(ann_inf, annot_sorted, annot, size_annot, annot_ext, 1, id_genome, size_id_genome);*/

        ERROR("modify_annotations()")
    }

    if (ann_inf->last_added != id_genome){

        if (size_annot_cplx != 0){

            if (ann_inf->min_size > size_annot_cplx){
                increase_size_annot_cplx_nodes(uc, size_substring, nb_children, ann_inf->min_size, 1);
                annot_cplx = get_annot_cplx_nodes(uc, size_substring, nb_children, position);
            }

            modify_mode_annotation(ann_inf, annot_cplx, ann_inf->min_size, NULL, 0, id_genome, size_id_genome);
        }
        else {

            if (ann_inf->min_size > uc->size_annot){
                if ((special_case == 0) || (annot_ext == NULL) || (ann_inf->min_size > uc->size_annot+1)){
                    annot_ext = realloc_annotation(uc, size_substring, nb_children, ann_inf->min_size, 0, position);
                }
            }


            modify_mode_annotation(ann_inf, &(uc->suffixes[position * (size_substring + uc->size_annot) + size_substring]),
                                   uc->size_annot, annot_ext, 1, id_genome, size_id_genome);

            if ((annot_ext != NULL) && (annot_ext[0] == 0))
                delete_extend_annots(uc, size_substring, nb_children, position, position, 0, 0, 1);
        }
    }

    if (annot_array_elem != NULL){
        free(annot_array_elem->annot_array);
        free(annot_array_elem);
    }

    reinit_annotation_inform(ann_inf);

    return;
}
