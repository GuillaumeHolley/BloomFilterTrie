#include "./../lib/branchingNode.h"

/* ---------------------------------------------------------------------------------------------------------------
*  insertKmer_Node(node, kmer, size_kmer, id_genome, info_per_lvl, ann_inf)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert a suffix or modify a suffix annotation into a node of the tree.
*  ---------------------------------------------------------------------------------------------------------------
*  node: pointer on a Node structure to insert the suffix
*  kmer: pointer on the suffix to insert
*  size_kmer: size of the suffix to insert, in char.
*  id_genome: genome identity of the suffix to insert
*  info_per_lvl: ptr on info_per_level structure, contains information to manipulate CCs field CC->children_type
*  ann_inf: ptr on annotation_inform structure, used to make the transition between reading and modifying an annot
*  ---------------------------------------------------------------------------------------------------------------
*/
int isBranchingRight(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer){

    ASSERT_NULL_PTR(node,"isBranchingRight()")
    ASSERT_NULL_PTR(root,"isBranchingRight()")
    ASSERT_NULL_PTR(kmer,"isBranchingRight()")

    int i, k;
    int count = 0;
    int size_line;

    uint16_t nb_elt;

    uint8_t mask;

    CC* cc;
    UC* uc;

    info_per_level*  info_per_lvl = &(root->info_per_lvl[lvl_node]);

    __builtin_prefetch (info_per_lvl, 0, 0);

    resultPresence* res = malloc(4*sizeof(resultPresence));
    uint8_t* right_shifting = calloc(info_per_lvl->size_kmer_in_bytes, sizeof(uint8_t));

    ASSERT_NULL_PTR(res,"isBranchingRight()")
    ASSERT_NULL_PTR(right_shifting,"isBranchingRight()")

    if (info_per_lvl->root == 1){
        for (i=0; i < info_per_lvl->size_kmer_in_bytes; i++){
            right_shifting[i] = kmer[i] >> 2;
            if (i+1 < info_per_lvl->size_kmer_in_bytes) right_shifting[i] |= kmer[i+1] << 6;
        }
    }
    else memcpy(right_shifting, kmer, info_per_lvl->size_kmer_in_bytes*sizeof(uint8_t));

    presenceNeighborsRight(node, root, right_shifting, size_kmer, res);

    if (size_kmer == NB_CHAR_SUF_PREF){

        for (i=0; i < 4; i++)
            if (res[i].link_child != NULL) count++;

        goto FREE_STRUCT;
    }
    else{
        int j=0;
        int nb_cell = info_per_lvl->size_kmer_in_bytes;
        int nb_cell_to_delete = 2;

        if (size_kmer == 45) nb_cell_to_delete++;

        for (j=0; j<nb_cell-nb_cell_to_delete; j++){
            right_shifting[j] = right_shifting[j+2] >> 2;
            if (j+3 < nb_cell) right_shifting[j] |= right_shifting[j+3] << 6;
        }

        right_shifting[j-1] &= info_per_lvl->mask_shift_kmer;

        for (i=0; i < 4; i++){
            if (res[i].link_child != NULL){
                if (res[i].children_type_leaf == 0){
                    if (res[i].container_is_UC == 0){
                        count = isBranchingRight((Node*)res[i].link_child, root, lvl_node-1, right_shifting,
                                                 size_kmer-NB_CHAR_SUF_PREF);
                        goto FREE_STRUCT;
                    }
                    else count++;
                }
                else{
                    cc = (CC*)res[i].container;
                    uc = &(((UC*)cc->children)[res[i].bucket]);
                    nb_elt = getNbElts(cc, res[i].pos_extra_filter3, (cc->type >> 6) & 0x1);

                    nb_cell = info_per_lvl->size_kmer_in_bytes_minus_1;
                    size_line = nb_cell+uc->size_annot;

                    mask = root->info_per_lvl[lvl_node-1].mask_shift_kmer;
                    if (mask == 0xff) mask = 0;

                    for (k=res[i].pos_sub_bucket*size_line; k<(res[i].pos_sub_bucket+nb_elt)*size_line; k+=size_line){
                        if (memcmp(&(uc->suffixes[k]), right_shifting, (nb_cell-1)*sizeof(uint8_t)) == 0){
                            if ((uc->suffixes[k+nb_cell-1] & mask) == right_shifting[nb_cell-1]){
                                count++;
                                if (count == 4) goto FREE_STRUCT;
                            }
                        }
                    }
                }
            }
        }
    }

    FREE_STRUCT: free(res);
    free(right_shifting);

    return count;
}

resultPresence* getRightNeighbors(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer){

    ASSERT_NULL_PTR(node,"getRightNeighbors()")
    ASSERT_NULL_PTR(root,"getRightNeighbors()")
    ASSERT_NULL_PTR(kmer,"getRightNeighbors()")

    CC* cc;
    UC* uc;

    info_per_level*  info_per_lvl = &(root->info_per_lvl[lvl_node]);

    int i = 0;
    int count = 0;
    int shifting_suffix = SIZE_BITS_UINT_8T
                            - (info_per_lvl->size_kmer_in_bytes_minus_1*SIZE_BITS_UINT_8T
                            - (size_kmer-NB_CHAR_SUF_PREF)*2) - 2;
    int size_line;
    int bucket;
    int pos_sub_bucket;
    int k;

    uint16_t nb_elt;

    uint8_t mask;
    uint8_t nuc;

    __builtin_prefetch (info_per_lvl, 0, 0);

    resultPresence* res = malloc(4 * sizeof(resultPresence));
    uint8_t* right_shifting = calloc(info_per_lvl->size_kmer_in_bytes, sizeof(uint8_t));

    ASSERT_NULL_PTR(res,"getRightNeighbors()")
    ASSERT_NULL_PTR(right_shifting,"getRightNeighbors()")

    if (info_per_lvl->root == 1){
        for (i=0; i < info_per_lvl->size_kmer_in_bytes; i++){
            kmer[i] >>= 2;
            if (i+1 < info_per_lvl->size_kmer_in_bytes) kmer[i] |= kmer[i+1] << 6;
        }
    }

    memcpy(right_shifting, kmer, info_per_lvl->size_kmer_in_bytes*sizeof(uint8_t));

    presenceNeighborsRight(node, root, right_shifting, size_kmer, res);

    if (size_kmer == NB_CHAR_SUF_PREF){
        free(right_shifting);
        return res;
    }
    else{
        int j=0;
        int nb_cell = info_per_lvl->size_kmer_in_bytes;
        for (j=0; j<nb_cell-2; j++){
            right_shifting[j] = right_shifting[j+2] >> 2;
            if (j+3 < nb_cell) right_shifting[j] |= right_shifting[j+3] << 6;
        }
        right_shifting[j-1] &= info_per_lvl->mask_shift_kmer;

        for (i=0; i < 4; i++){
            if (res[i].link_child != NULL){
                if (res[i].children_type_leaf == 0){
                    if (res[i].container_is_UC == 0){
                        resultPresence* res_tmp = getRightNeighbors((Node*)res[i].link_child, root, lvl_node - 1, right_shifting,
                                                                    size_kmer-NB_CHAR_SUF_PREF);
                        free(res);
                        free(right_shifting);
                        return res_tmp;
                    }
                }
                else{
                    cc = (CC*)res[i].container;
                    uc = &(((UC*)cc->children)[res[i].bucket]);

                    bucket = res[i].bucket;
                    pos_sub_bucket = res[i].pos_sub_bucket;
                    nb_elt = getNbElts(res[i].container, res[i].pos_extra_filter3, (cc->type >> 6) & 0x1);

                    nb_cell = info_per_lvl->size_kmer_in_bytes_minus_1;
                    size_line = nb_cell+uc->size_annot;

                    mask = root->info_per_lvl[lvl_node-1].mask_shift_kmer;
                    if (mask == 0xff) mask = 0;

                    res[i].link_child = NULL;
                    res[i].container = NULL;

                    for (k=pos_sub_bucket*size_line; k<(pos_sub_bucket+nb_elt)*size_line; k+=size_line){ //Iterate over every suffix stored
                        if (memcmp(&(uc->suffixes[k]), right_shifting, (nb_cell-1)*sizeof(uint8_t)) == 0){
                            if ((uc->suffixes[k+nb_cell-1] & mask) == right_shifting[nb_cell-1]){

                                nuc = (uc->suffixes[k+nb_cell-1] >> shifting_suffix) & 0x3;
                                res[nuc].link_child = &(((UC*)cc->children)[bucket].suffixes[k]);
                                res[nuc].container = uc;
                                res[nuc].posFilter2 = nb_cell;
                                res[nuc].pos_sub_bucket = k/size_line;

                                res[nuc].posFilter3 = uc->nb_children;

                                count++;
                                if (count == 4) break;
                            }
                        }
                    }

                    free(right_shifting);
                    return res;
                }
            }
        }
    }

    free(right_shifting);
    return res;
}

/* ---------------------------------------------------------------------------------------------------------------
*  insertKmer_Node(node, kmer, size_kmer, id_genome, info_per_lvl, ann_inf)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert a suffix or modify a suffix annotation into a node of the tree.
*  ---------------------------------------------------------------------------------------------------------------
*  node: pointer on a Node structure to insert the suffix
*  kmer: pointer on the suffix to insert
*  size_kmer: size of the suffix to insert, in char.
*  id_genome: genome identity of the suffix to insert
*  info_per_lvl: ptr on info_per_level structure, contains information to manipulate CCs field CC->children_type
*  ann_inf: ptr on annotation_inform structure, used to make the transition between reading and modifying an annot
*  ---------------------------------------------------------------------------------------------------------------
*/
int isBranchingLeft(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer){

    ASSERT_NULL_PTR(node,"isBranchingLeft()")
    ASSERT_NULL_PTR(root,"isBranchingLeft()")
    ASSERT_NULL_PTR(kmer,"isBranchingLeft()")

    int i = 0;
    int count = 0;
    int k;
    int size_line;

    uint16_t nb_elt;

    CC* cc;
    UC* uc;

    info_per_level*  info_per_lvl = &(root->info_per_lvl[lvl_node]);

    __builtin_prefetch (info_per_lvl, 0, 0);

    resultPresence* res = malloc(4*sizeof(resultPresence));
    ASSERT_NULL_PTR(res,"isBranchingLeft()")

    uint8_t* left_shifting = calloc(info_per_lvl->size_kmer_in_bytes, sizeof(uint8_t));
    ASSERT_NULL_PTR(left_shifting,"isBranchingLeft()")

    if (info_per_lvl->root == 1){
        uint8_t shifting = ((uint8_t)0xff) >> (info_per_lvl->size_kmer_in_bytes*SIZE_BITS_UINT_8T - size_kmer*2);

        for (i = info_per_lvl->size_kmer_in_bytes-1; i >= 0; i--){
            left_shifting[i] = kmer[i] << 2;
            if (i > 0) left_shifting[i] |= kmer[i-1] >> 6;
        }

        left_shifting[info_per_lvl->size_kmer_in_bytes-1] &= shifting;
    }
    else memcpy(left_shifting, kmer, info_per_lvl->size_kmer_in_bytes*sizeof(uint8_t));

    presenceNeighborsLeft(node, root, left_shifting, size_kmer, res);

    if (size_kmer == NB_CHAR_SUF_PREF){
        for (i=0; i < 4; i++){
            if (res[i].link_child != NULL) count = 1;
        }
        goto FREE_STRUCT;
    }

    int j=0;
    int nb_cell = info_per_lvl->size_kmer_in_bytes;
    int nb_cell_to_delete = 2;

    if (size_kmer == 45) nb_cell_to_delete++;

    for (j=0; j<nb_cell-nb_cell_to_delete; j++){
        left_shifting[j] = left_shifting[j+2] >> 2;
        if (j+3 < nb_cell) left_shifting[j] |= left_shifting[j+3] << 6;
    }

    left_shifting[j-1] &= info_per_lvl->mask_shift_kmer;

    if (info_per_lvl->root == 1){

        for (i=0; i < 4; i++){

            if (res[i].link_child != NULL){

                if (res[i].children_type_leaf == 0){

                    if (res[i].container_is_UC != 0) count++;
                    else count += isBranchingLeft((Node*)res[i].link_child, root, lvl_node-1, left_shifting,
                                                  size_kmer-NB_CHAR_SUF_PREF);
                }
                else{
                    if (size_kmer == NB_CHAR_SUF_PREF) count++;
                    else {
                        cc = (CC*)res[i].container;
                        uc = &(((UC*)cc->children)[res[i].bucket]);

                        nb_cell = info_per_lvl->size_kmer_in_bytes_minus_1;

                        nb_elt = getNbElts(cc, res[i].pos_extra_filter3, (cc->type >> 6) & 0x1);
                        size_line = nb_cell+uc->size_annot;

                        k = binary_search_UC(uc, res[i].pos_sub_bucket, res[i].pos_sub_bucket + nb_elt - 1, left_shifting, nb_cell, 0xff);
                        if (memcmp(&(uc->suffixes[k * size_line]), left_shifting, nb_cell * sizeof(uint8_t)) == 0) count++;
                    }
                }
            }
        }
    }
    else if (info_per_lvl->level_min == 1){

        for (i=0; i < 4; i++){

            if (res[i].link_child != NULL){

                if (res[i].children_type_leaf == 0){

                    if (res[i].container_is_UC != 0) count = 1;
                    else count = isBranchingLeft((Node*)res[i].link_child, root, lvl_node - 1, left_shifting,
                                                 size_kmer-NB_CHAR_SUF_PREF);

                    break;
                }
                else{
                    if (size_kmer == NB_CHAR_SUF_PREF){
                        count = 1;
                        break;
                    }
                    else {
                        cc = (CC*)res[i].container;
                        uc = &(((UC*)cc->children)[res[i].bucket]);

                        nb_cell = info_per_lvl->size_kmer_in_bytes_minus_1;

                        nb_elt = getNbElts(cc, res[i].pos_extra_filter3, (cc->type >> 6) & 0x1);
                        size_line = nb_cell+uc->size_annot;

                        k = binary_search_UC(uc, res[i].pos_sub_bucket, res[i].pos_sub_bucket + nb_elt - 1, left_shifting, nb_cell, 0xff);

                        if (memcmp(&(uc->suffixes[k * size_line]), left_shifting, nb_cell * sizeof(uint8_t)) == 0){
                            count = 1;
                            goto FREE_STRUCT;
                        }
                    }
                }
            }
        }
    }
    else{
        for (i=0; i < 4; i++){

            if (res[i].link_child != NULL){

                if (res[i].children_type_leaf == 0){

                    if (res[i].container_is_UC != 0) count = 1;
                    else count = isBranchingLeft((Node*)res[i].link_child, root, lvl_node-1, left_shifting,
                                                 size_kmer-NB_CHAR_SUF_PREF);
                    break;
                }
                else{
                    if (size_kmer == NB_CHAR_SUF_PREF){
                        count = 1;
                        break;
                    }
                    else {
                        cc = (CC*)res[i].container;
                        uc = &(((UC*)cc->children)[res[i].bucket]);

                        nb_cell = info_per_lvl->size_kmer_in_bytes_minus_1;

                        nb_elt = getNbElts(cc, res[i].pos_extra_filter3, (cc->type >> 6) & 0x1);
                        size_line = nb_cell+uc->size_annot;

                        for (k=res[i].pos_sub_bucket*size_line; k<(res[i].pos_sub_bucket+nb_elt)*size_line; k+=size_line){

                            if (memcmp(&(uc->suffixes[k]), left_shifting, (nb_cell-1)*sizeof(uint8_t)) == 0){

                                if ((uc->suffixes[k+nb_cell-1] & 0x7f) == left_shifting[nb_cell-1]){
                                    count = 1;
                                    goto FREE_STRUCT;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    FREE_STRUCT: free(res);
    free(left_shifting);

    return count;
}

/* ---------------------------------------------------------------------------------------------------------------
*  insertKmer_Node(node, kmer, size_kmer, id_genome, info_per_lvl, ann_inf)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert a suffix or modify a suffix annotation into a node of the tree.
*  ---------------------------------------------------------------------------------------------------------------
*  node: pointer on a Node structure to insert the suffix
*  kmer: pointer on the suffix to insert
*  size_kmer: size of the suffix to insert, in char.
*  id_genome: genome identity of the suffix to insert
*  info_per_lvl: ptr on info_per_level structure, contains information to manipulate CCs field CC->children_type
*  ann_inf: ptr on annotation_inform structure, used to make the transition between reading and modifying an annot
*  ---------------------------------------------------------------------------------------------------------------
*/
resultPresence* getLeftNeighbors(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer){

    ASSERT_NULL_PTR(node,"getLeftNeighbors()")
    ASSERT_NULL_PTR(root,"getLeftNeighbors()")
    ASSERT_NULL_PTR(kmer,"getLeftNeighbors()")

    CC* cc;
    UC* uc;

    info_per_level*  info_per_lvl = &(root->info_per_lvl[lvl_node]);

    int i = 0, j = 0, k;
    int size_line;
    int nb_cell = info_per_lvl->size_kmer_in_bytes;

    uint16_t nb_elt;

    resultPresence* res_tmp = NULL;
    resultPresence* res = malloc(4*sizeof(resultPresence));
    ASSERT_NULL_PTR(res,"getLeftNeighbors()")

    __builtin_prefetch (info_per_lvl, 0, 0);

    uint8_t* left_shifting = calloc(info_per_lvl->size_kmer_in_bytes, sizeof(uint8_t));
    ASSERT_NULL_PTR(left_shifting,"getLeftNeighbors()")

    if (info_per_lvl->root == 1){

        uint8_t shifting = ((uint8_t)0xff) >> (info_per_lvl->size_kmer_in_bytes*SIZE_BITS_UINT_8T - size_kmer*2);

        for (i = info_per_lvl->size_kmer_in_bytes-1; i >= 0; i--){
            kmer[i] <<= 2;
            if (i > 0) kmer[i] |= kmer[i-1] >> 6;
        }

        kmer[info_per_lvl->size_kmer_in_bytes-1] &= shifting;
    }

    memcpy(left_shifting, kmer, info_per_lvl->size_kmer_in_bytes*sizeof(uint8_t));

    presenceNeighborsLeft(node, root, left_shifting, size_kmer, res);

    if (size_kmer == NB_CHAR_SUF_PREF){
        for (i=0; i < 4; i++){
            if (res[i].link_child != NULL){
                res_tmp = malloc(sizeof(resultPresence));
                memcpy(res_tmp, &(res[i]), sizeof(resultPresence));
                break;
            }
        }
        goto FREE_STRUCT;
    }

    for (j=0; j<nb_cell-2; j++){
        left_shifting[j] = left_shifting[j+2] >> 2;
        if (j+3 < nb_cell) left_shifting[j] |= left_shifting[j+3] << 6;
    }
    left_shifting[j-1] &= info_per_lvl->mask_shift_kmer;

    if (info_per_lvl->root == 1){
        for (i=0; i < 4; i++){
            if (res[i].link_child != NULL){
                if (res[i].children_type_leaf == 0){
                    if (res[i].container_is_UC == 0){
                        if ((res_tmp = getLeftNeighbors((Node*)res[i].link_child, root, lvl_node-1, left_shifting,
                                                        size_kmer-NB_CHAR_SUF_PREF)) != NULL){;
                            memcpy(&(res[i]), res_tmp, sizeof(resultPresence));
                            free(res_tmp);
                        }
                        else res[i].link_child = NULL;
                    }
                }
                else{

                    cc = (CC*)res[i].container;
                    uc = &(((UC*)cc->children)[res[i].bucket]);

                    nb_cell = info_per_lvl->size_kmer_in_bytes_minus_1;

                    nb_elt = getNbElts(cc, res[i].pos_extra_filter3, (cc->type >> 6) & 0x1);
                    size_line = nb_cell+uc->size_annot;

                    res[i].link_child = NULL;

                    k = binary_search_UC(uc, res[i].pos_sub_bucket, res[i].pos_sub_bucket + nb_elt - 1, left_shifting, nb_cell, 0xff);

                    if (memcmp(&(uc->suffixes[k * size_line]), left_shifting, nb_cell * sizeof(uint8_t)) == 0){
                        res[i].link_child = &(uc->suffixes[k * size_line]);
                        res[i].container = uc;
                        res[i].posFilter3 = uc->nb_children;
                        res[i].posFilter2 = nb_cell;
                        res[i].pos_sub_bucket = k;
                    }
                }
            }
        }

        free(left_shifting);
        return res;
    }
    else if (info_per_lvl->level_min == 1){

        for (i=0; i < 4; i++){

            if (res[i].link_child != NULL){

                if (res[i].children_type_leaf == 0){

                    if (res[i].container_is_UC == 0){
                        res_tmp = getLeftNeighbors((Node*)res[i].link_child, root, lvl_node - 1, left_shifting,
                                                   size_kmer-NB_CHAR_SUF_PREF);
                    }
                    else {
                        res_tmp = create_resultPresence();
                        memcpy(res_tmp, &(res[i]), sizeof(resultPresence));
                    }

                    break;
                }
                else{
                    res_tmp = create_resultPresence();

                    if (size_kmer == NB_CHAR_SUF_PREF){
                        memcpy(res_tmp, &(res[i]), sizeof(resultPresence));
                        break;
                    }
                    else {
                        cc = (CC*)res[i].container;
                        uc = &(((UC*)cc->children)[res[i].bucket]);

                        nb_cell = info_per_lvl->size_kmer_in_bytes_minus_1;

                        nb_elt = getNbElts(cc, res[i].pos_extra_filter3, (cc->type >> 6) & 0x1);
                        size_line = nb_cell+uc->size_annot;

                        k = binary_search_UC(uc, res[i].pos_sub_bucket, res[i].pos_sub_bucket + nb_elt - 1, left_shifting, nb_cell, 0xff);

                        if (memcmp(&(uc->suffixes[k * size_line]), left_shifting, nb_cell * sizeof(uint8_t)) == 0){
                                res_tmp->link_child = &(uc->suffixes[k * size_line]);
                                res_tmp->container = uc;
                                res_tmp->posFilter3 = uc->nb_children;
                                res_tmp->posFilter2 = nb_cell;
                                res_tmp->pos_sub_bucket = k;

                                goto FREE_STRUCT;
                        }
                    }
                }
            }
        }
    }
    else{
        for (i=0; i < 4; i++){

            if (res[i].link_child != NULL){

                if (res[i].children_type_leaf == 0){

                    if (res[i].container_is_UC == 0){
                        res_tmp = getLeftNeighbors((Node*)res[i].link_child, root, lvl_node - 1, left_shifting,
                                                   size_kmer-NB_CHAR_SUF_PREF);
                    }
                    else {
                        res_tmp = create_resultPresence();
                        memcpy(res_tmp, &(res[i]), sizeof(resultPresence));
                    }

                    break;
                }
                else{

                    res_tmp = create_resultPresence();

                    if (size_kmer == NB_CHAR_SUF_PREF){
                        memcpy(res_tmp, &(res[i]), sizeof(resultPresence));
                        break;
                    }
                    else {
                        cc = (CC*)res[i].container;
                        uc = &(((UC*)cc->children)[res[i].bucket]);

                        nb_cell = info_per_lvl->size_kmer_in_bytes_minus_1;

                        nb_elt = getNbElts(cc, res[i].pos_extra_filter3, (cc->type >> 6) & 0x1);
                        size_line = nb_cell+uc->size_annot;

                        for (k=res[i].pos_sub_bucket*size_line; k<(res[i].pos_sub_bucket+nb_elt)*size_line; k+=size_line){ //Iterate over every suffix stored

                            if (memcmp(&(uc->suffixes[k]), left_shifting, (nb_cell-1)*sizeof(uint8_t)) == 0){ //if there is a match

                                if ((uc->suffixes[k+nb_cell-1] & 0x7f) == left_shifting[nb_cell-1]){

                                    res_tmp->link_child = &(uc->suffixes[k]);
                                    res_tmp->container = uc;
                                    res_tmp->posFilter3 = uc->nb_children;
                                    res_tmp->posFilter2 = nb_cell;
                                    res_tmp->pos_sub_bucket = k/size_line;

                                    goto FREE_STRUCT;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    FREE_STRUCT: free(res);
    free(left_shifting);

    return res_tmp;
}

void extractSuffix(uint8_t* kmer_tmp, int size_kmer, int size_kmer_array, int shifting_prefix, int it_bucket,
                   uint8_t* suffix_start, int size_suffix_start_bytes, int shifting_suffix){

    ASSERT_NULL_PTR(kmer_tmp,"extractSuffix()")
    ASSERT_NULL_PTR(suffix_start,"extractSuffix()")

    int k, size_new_substring, size_new_substring_bytes, last_shift, shifting;
    uint64_t new_substring;

    int it_substring = 0;
    int shifting_suffix_cpy = shifting_suffix;

    while (it_substring < size_suffix_start_bytes){

        it_substring += sizeof(uint64_t);
        new_substring = 0;

        if (it_substring > size_suffix_start_bytes){

            size_new_substring = (size_kmer-NB_CHAR_SUF_PREF)*2-((it_substring-sizeof(uint64_t))*SIZE_BITS_UINT_8T);
            size_new_substring_bytes = CEIL(size_new_substring, SIZE_BITS_UINT_8T);

            for (k=0; k<size_new_substring_bytes; k++)
                new_substring = (new_substring << 8) | reverse_word_8(suffix_start[(it_substring-sizeof(uint64_t))+k]);

            new_substring >>= (size_suffix_start_bytes - (it_substring-sizeof(uint64_t))) * SIZE_BITS_UINT_8T - size_new_substring;
        }
        else{

            size_new_substring = sizeof(uint64_t)*SIZE_BITS_UINT_8T;
            size_new_substring_bytes = sizeof(uint64_t);

            for (k=0; k<size_new_substring_bytes; k++)
                new_substring = (new_substring << 8) | reverse_word_8(suffix_start[(it_substring-sizeof(uint64_t))+k]);
        }

        shifting = shifting_prefix;
        if (shifting == 0) shifting = 8;
        else size_new_substring_bytes++;

        if (shifting_suffix_cpy > 0){
            size_new_substring -= shifting_suffix_cpy;
            shifting_suffix_cpy = 0;
        }

        for (k=it_bucket; k<it_bucket+size_new_substring_bytes; k++){

            last_shift = size_new_substring - shifting;

            if (last_shift >= 0) kmer_tmp[k] |= new_substring >> last_shift;
            else kmer_tmp[k] |= new_substring << abs(last_shift);

            shifting += SIZE_BITS_UINT_8T;
        }

        if ((shifting%SIZE_BITS_UINT_8T != 0)) size_new_substring_bytes--;

        it_bucket += size_new_substring_bytes;
    }

    for (k=0; k<size_kmer_array; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);
}
