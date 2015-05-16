#include "presenceNode.h"

/* ---------------------------------------------------------------------------------------------------------------
*  presenceKmer(node, kmer, size_kmer, nb_CC_node, func_on_types)
*  ---------------------------------------------------------------------------------------------------------------
*  Search for a prefix into a node.
*  ---------------------------------------------------------------------------------------------------------------
*  node: pointer on a Node structure
*  kmer: pointer on the prefix to search
*  size_kmer: size of the prefix, in char.
*  nb_CC_node: number of CCs in this node
*  func_on_types: ptr on ptrs_on_func structure, contains information to manipulate CCs field CC->children_type
*  ---------------------------------------------------------------------------------------------------------------
*/
void presenceNeighborsLeft(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types, resultPresence* res, uint16_t** skip_node_root){

    if (size_kmer == 0) ERROR("presenceNeighborsLeft(): this case should not happen")
    ASSERT_NULL_PTR(node,"presenceNeighborsLeft()")
    ASSERT_NULL_PTR(kmer,"presenceNeighborsLeft()")

    //Allocate and initialize the structure resultPresence
    int last_count_node = -1;
    int cpt_node_tmp = -1;

    int i;

    uint8_t substring[3];
    uint8_t nuc2add;

    //Keep in resultPresence->substring only the prefix we are looking for
    if (func_on_types->root == 1) substring[0] = reverse_word_8(kmer[0]) & 0x3f;
    else substring[0] = reverse_word_8(kmer[0]);

    substring[1] = reverse_word_8(kmer[1]);
    substring[2] = reverse_word_8(kmer[2]) & 0xc0;

    //Compute hash-functions of the Bloom filter for the prefix
    uint16_t substring_prefix = (substring[0] << SIZE_CELL) | substring[1];
    uint16_t hash1_v = dbj2(substring[1], substring[0] & 0x3f, MODULO_HASH);
    uint16_t hash2_v = sdbm(substring[1], substring[0] & 0x3f, MODULO_HASH);

    CC* cc;
    UC* uc;
    UC* uc_tmp;

    for (i=0; i<4; i++) initialize_resultPresence(&(res[i]));

    //We first iterate on the CCs
    i = -1;
    if ((CC*)node->CC_array != NULL){
        do {
            i++;
            cc = &(((CC*)node->CC_array)[i]);

            //First look if the prefix is in the Bloom filter
            if ((cc->BF_filter2[hash1_v/SIZE_CELL] & MASK_POWER_8[hash1_v%SIZE_CELL]) == 0) continue;
            if ((cc->BF_filter2[hash2_v/SIZE_CELL] & MASK_POWER_8[hash2_v%SIZE_CELL]) == 0) continue;

            substring[0] = substring_prefix >> 6;
            substring[1] = (substring_prefix << 2) | (substring[2] >> 6);
            substring[2] = (substring_prefix >> 8) & 0xc0;

            //At this point, this prefix is said present by the BF, maybe it's a True Positive, maybe it's a False Positive
            //If it is a FP, the prefix will have to be inserted in this CC (FP recycling)
            uint16_t size_bf = cc->type >> 8;
            uint8_t s = (cc->type >> 2) & 0x3f; //length v of the prefixes stored in this CC (p_vs)
            uint8_t p = SIZE_SEED*2-s; //length u of the prefixes stored in this CC (p_us)

            uint8_t word_tmp;

            int size_filter2 = size_bf+(MASK_POWER_16[p]/SIZE_CELL); //Size BF + filter2 in bytes
            int size_filter2_n_skip = size_filter2; //Size BF + filter2 + SkipFilter2 (SkipFilter2 may not exist) in bytes
            int skip_filter2 = MASK_POWER_16[p]/0xf8;
            if (cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM) size_filter2_n_skip += skip_filter2;

            //Compute p_u, the first u char. of the prefix p we are looking for
            int posFilter2 = 0;
            if (p==10) posFilter2 = (((uint16_t)substring[0]) << 2) | (((uint16_t)substring[1]) >> 6);
            else posFilter2 = (((uint16_t)substring[0]) << 6) | (((uint16_t)substring[1]) >> 2);

            //If the prefix is present in filter2, need to compute the Hamming weight between position 0 and the position of p_u
            //in order to know how many p_u are lexicographically inferior or equal to the one we insert
            if ((cc->BF_filter2[size_bf+posFilter2/SIZE_CELL] & MASK_POWER_8[posFilter2%SIZE_CELL]) != 0){

                int k=size_bf+posFilter2/SIZE_CELL;
                int cnt = 0;
                int hamming_weight = 0;

                //If SkipFilter2 exists, we use it to accelerate the computation of the Hamming weight
                if (cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM){
                    int skip_posfilter2 = posFilter2/0xf8;

                    while ((cnt<skip_posfilter2) && (cnt < skip_filter2)){
                        hamming_weight += cc->BF_filter2[size_filter2+cnt];
                        cnt++;
                    }
                }

                //finish the computation of the Hamming weight
                hamming_weight += popcnt_8_par(cc->BF_filter2, size_bf+cnt*31, k);

                word_tmp = cc->BF_filter2[k];
                for (k=0; k<=posFilter2%8; k++)
                    if ((word_tmp & MASK_POWER_8[k]) != 0) hamming_weight++;

                posFilter2 = hamming_weight; //posFilter2 now contains this Hamming weight

                //At this point, we know that the p_v we are looking for should be situated in the third filter
                //at a slot (plus subsequent one) corresponding to the hamming_weight-th 1 in the extra_filter3
                int nb_skp = CEIL(cc->nb_elem, NB_CHILDREN_PER_SKP);
                int pos_extra_filter3 = INT_MAX;
                int hamming_weight_0 = 0;

                int j=0, m=0, sum=0;
                int nb_cell_3rdlist = CEIL(cc->nb_elem,SIZE_CELL);

                k=0;
                word_tmp = 0;
                hamming_weight = 0;

                //This loop tries to use SkipFilter3 (if it exists) to accelerate the Hamming weight comput.
                //in extra_filter3 by skipping counting 1s directly in extra_filter3
                while(m<cc->nb_elem/0xf8){
                    if ((sum = hamming_weight + cc->BF_filter2[size_filter2_n_skip+m]) < posFilter2){
                        hamming_weight = sum;
                        m++;
                    }
                    else break;
                }

                //Finish the Hamming weight computation directly in extra_filter3
                if (func_on_types->level_min == 1){
                    for (k=m*31; k<nb_cell_3rdlist; k++){

                        if ((sum = hamming_weight+popcnt_8(cc->extra_filter3[k])) >= posFilter2){
                            word_tmp = cc->extra_filter3[k];

                            int size_word = 7;
                            if (k == nb_cell_3rdlist-1) size_word = (cc->nb_elem-1)%8;

                            for (j=0; j<=size_word; j++){
                                if (((word_tmp >> j)&1) == 1){
                                    hamming_weight++;
                                    if (hamming_weight == posFilter2){
                                        pos_extra_filter3 = k*8+j;
                                        j++;
                                        break;
                                    }
                                }
                            }

                            if (k == nb_cell_3rdlist) k--;
                            goto MATCH;
                        }
                        else hamming_weight = sum;
                    }
                    //If the substring must be inserted in the last slot of the list
                    if (k == nb_cell_3rdlist) k--;

                    MATCH: if (pos_extra_filter3 == INT_MAX) pos_extra_filter3 = cc->nb_elem;

                    int a=0;
                    for (a=k; a<nb_cell_3rdlist; a++){
                        uint8_t word_tmp = cc->extra_filter3[a];
                        while (j<8){
                            if (((word_tmp >> j)&1) == 0){
                                if (a*8+j < cc->nb_elem) hamming_weight_0++;
                                else goto OUT_LOOP;
                            }
                            else goto OUT_LOOP;
                            j++;
                        }
                        j=0;
                    }
                }
                else{
                    int size_line_children;
                    int pos_children;
                    int pos_children_tmp;
                    int cpt_pv;
                    int cpt_node;
                    int nb_elem_in_pv;
                    int end, end_tmp;
                    int it, it_tmp;

                    int cpt_pv_tmp = 0;
                    int tmp_hamming_weight = hamming_weight;

                    uint8_t* children;

                    k = m*0xf8;
                    pos_children_tmp = k/NB_CHILDREN_PER_SKP;

                    (*func_on_types->count_Nodes_Children)(&(((CC*)node->CC_array)[i]), pos_children_tmp*NB_CHILDREN_PER_SKP, k, &cpt_pv_tmp, &cpt_node_tmp);
                    last_count_node = (cpt_node_tmp += (*func_on_types->count_nodes)(&(((CC*)node->CC_array)[i]), 0, pos_children_tmp*NB_CHILDREN_PER_SKP));

                    it = k - pos_children_tmp * NB_CHILDREN_PER_SKP;
                    it_tmp = it;
                    nb_elem_in_pv = 0;

                    cpt_pv = cpt_pv_tmp;
                    cpt_node = cpt_node_tmp;
                    pos_children = pos_children_tmp;

                    while (pos_children < nb_skp){

                        if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * NB_CHILDREN_PER_SKP;
                        else end = NB_CHILDREN_PER_SKP;

                        uc = &(((UC*)cc->children)[pos_children]);
                        size_line_children = func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot;
                        children = &(uc->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                        while (it < end){

                            it_tmp = it;

                            for (m=k; m<MIN(k+8, cc->nb_elem); m++){

                                if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, m)) == 0){
                                    tmp_hamming_weight += cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1;
                                    cpt_node++;
                                }
                                else{
                                    tmp_hamming_weight += children[cpt_pv*size_line_children] >> 7;
                                    cpt_pv += nb_elem_in_pv;
                                }

                                it++;
                            }

                            if (tmp_hamming_weight >= posFilter2){

                                if (pos_children_tmp == nb_skp - 1) end_tmp = cc->nb_elem - pos_children_tmp * NB_CHILDREN_PER_SKP;
                                else end_tmp = NB_CHILDREN_PER_SKP;

                                uc_tmp = &(((UC*)cc->children)[pos_children_tmp]);
                                size_line_children = func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot;
                                children = &(uc_tmp->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                                while (cpt_pv_tmp < uc_tmp->nb_children){

                                    if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, k)) == 0){
                                        hamming_weight += cc->children_Node_container[cpt_node_tmp].UC_array.nb_children & 0x1;
                                        cpt_node_tmp++;
                                    }
                                    else{
                                        hamming_weight += children[cpt_pv_tmp*size_line_children] >> 7;
                                        cpt_pv_tmp += nb_elem_in_pv;
                                    }

                                    if (hamming_weight == posFilter2){
                                        pos_extra_filter3 = k;
                                        goto MATCH2;
                                    }

                                    k++;
                                    it_tmp++;
                                }
                            }
                            else{
                                hamming_weight = tmp_hamming_weight;
                                cpt_pv_tmp = cpt_pv;
                                cpt_node_tmp = cpt_node;
                                pos_children_tmp = pos_children;
                                k = m;
                            }
                        }

                        it = 0;
                        it_tmp = 0;
                        cpt_pv = 0;
                        cpt_pv_tmp = cpt_pv;
                        pos_children++;
                        pos_children_tmp = pos_children;
                    }

                    it_tmp = it;

                    MATCH2: if (pos_extra_filter3 == INT_MAX) pos_extra_filter3 = cc->nb_elem;

                    if (pos_children_tmp < nb_skp){

                        if (pos_children_tmp == nb_skp - 1) end_tmp = cc->nb_elem - pos_children_tmp * NB_CHILDREN_PER_SKP;
                        else end_tmp = NB_CHILDREN_PER_SKP;

                        it_tmp++;
                        k++;

                        if (it_tmp >= end_tmp){
                            it = 0;
                            cpt_pv = 0;
                            pos_children = pos_children_tmp+1;
                        }
                        else{
                            it = it_tmp;
                            cpt_pv = cpt_pv_tmp;
                            pos_children = pos_children_tmp;
                        }

                        cpt_node = cpt_node_tmp;

                        while (pos_children < nb_skp){

                            uc = &(((UC*)cc->children)[pos_children]);
                            size_line_children = func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot;
                            children = &(uc->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                            if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * NB_CHILDREN_PER_SKP;
                            else end = NB_CHILDREN_PER_SKP;

                            while (it < end){

                                if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, k)) == 0){
                                    if ((cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1) == 0) hamming_weight_0++;
                                    else goto OUT_LOOP;
                                    cpt_node++;
                                }
                                else{
                                    if (children[cpt_pv*size_line_children] < 0x80) hamming_weight_0++;
                                    else goto OUT_LOOP;
                                    cpt_pv += nb_elem_in_pv;
                                }

                                k++;
                                it++;
                            }

                            it = 0;
                            cpt_pv = 0;
                            pos_children++;
                        }
                    }
                }

                //We now compare the p_vs between start position and end position in filter3 to see if
                //the prefix we look for is present
                OUT_LOOP:if (pos_extra_filter3 < cc->nb_elem){
                    uint8_t suffix;
                    uint8_t tmp;

                    if (s==8){ //if the length of p_vs is 8bits
                        suffix = (substring[1] << 2) | (substring[2] >> 6); //compute p_v

                        if (func_on_types->root == 1){
                            suffix &= 0xfc;
                            uint8_t next_suffix = suffix | 0x3;
                            int last_k = 0;
                            int last_node = 0;
                            int clust;

                            for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){ //iterate between start and end position
                                if((cc->filter3[k] & 0xfc) == suffix){ //if there is a match

                                    nuc2add = cc->filter3[k] & 0x3;
                                    res[nuc2add].bucket = k/NB_CHILDREN_PER_SKP;
                                    clust = k/SIZE_CLUST_SKIP_NODES;
                                    uc = &(((UC*)cc->children)[res[nuc2add].bucket]);

                                    if (size_kmer != SIZE_SEED){
                                        if ((*func_on_types->is_child)((void*)cc, k)){
                                            res[nuc2add].pos_sub_bucket = (*func_on_types->count_children)((void*)cc, res[nuc2add].bucket*NB_CHILDREN_PER_SKP, k);
                                            res[nuc2add].pos_extra_filter3 = k;
                                            res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket*(func_on_types->size_kmer_in_bytes_minus_1+uc->size_annot)]);
                                            res[nuc2add].container = &(((CC*)node->CC_array)[i]);
                                            res[nuc2add].children_type_leaf = 1;
                                        }
                                        else{
                                            //res->link_child = &(cc->children_Node_container[(*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), 0, k)]);
                                            if (clust-1 >= 0){
                                                if ((skip_node_root != NULL) && (k-last_k > k-clust*SIZE_CLUST_SKIP_NODES)){
                                                    last_node = skip_node_root[i][clust-1] +
                                                                    (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), clust*SIZE_CLUST_SKIP_NODES, k);
                                                }
                                                else last_node += (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), last_k, k);
                                            }
                                            else last_node += (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), last_k, k);

                                            res[nuc2add].link_child = &(cc->children_Node_container[last_node]);
                                            last_k = k;
                                        }
                                    }
                                    else{
                                        res[nuc2add].pos_extra_filter3 = k;
                                        res[nuc2add].pos_sub_bucket = k%NB_CHILDREN_PER_SKP;
                                        res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket*uc->size_annot]);
                                        res[nuc2add].container = uc;
                                        res[nuc2add].posFilter2 = 0;
                                        res[nuc2add].posFilter3 = NB_CHILDREN_PER_SKP;
                                        res[nuc2add].children_type_leaf = 1;
                                    }
                                }
                                else if (cc->filter3[k] > next_suffix) break;
                            }
                        }
                        else{
                            for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){ //iterate between start and end position
                                if(cc->filter3[k] == suffix){ //if there is a match
                                    nuc2add = cc->filter3[k] & 0x3;

                                    if (size_kmer != SIZE_SEED){
                                        if ((*func_on_types->is_child)((void*)cc, k) == 1){

                                            res[nuc2add].bucket = k/NB_CHILDREN_PER_SKP;
                                            uc = &(((UC*)cc->children)[res[nuc2add].bucket]);

                                            res[nuc2add].pos_extra_filter3 = k;
                                            res[nuc2add].pos_sub_bucket = (*func_on_types->count_children)((void*)cc, res[nuc2add].bucket*NB_CHILDREN_PER_SKP, k);
                                            res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket*(func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot)]);
                                            res[nuc2add].container = &(((CC*)node->CC_array)[i]);
                                            res[nuc2add].children_type_leaf = 1;
                                        }
                                        else{
                                            if (last_count_node == -1) res[nuc2add].link_child = &(cc->children_Node_container[(*func_on_types->count_nodes)((void*)cc, 0, k)]);
                                            else res[nuc2add].link_child = &(cc->children_Node_container[MAX(cpt_node_tmp-1, 0)
                                                                             + (*func_on_types->count_nodes)((void*)cc, pos_extra_filter3, k)]);
                                        }
                                    }
                                    else{
                                        res[nuc2add].bucket = k/NB_CHILDREN_PER_SKP;
                                        uc = &(((UC*)cc->children)[res[nuc2add].bucket]);

                                        res[nuc2add].pos_extra_filter3 = k;
                                        res[nuc2add].pos_sub_bucket = k%NB_CHILDREN_PER_SKP;
                                        res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket * uc->size_annot]);
                                        res[nuc2add].container = uc;
                                        res[nuc2add].posFilter2 = 0;
                                        res[nuc2add].posFilter3 = NB_CHILDREN_PER_SKP;
                                        res[nuc2add].children_type_leaf = 1;
                                    }

                                    return;
                                }
                                else if(cc->filter3[k] > suffix) break;
                            }
                        }
                    }
                    else {
                        suffix = ((substring[1] & 0x3) << 2) | (substring[2] >> 6);

                        if (func_on_types->root == 1){

                            suffix &= 0xfc;
                            uint8_t next_suffix = suffix | 0x3;
                            int last_k = 0;
                            int last_node = 0;
                            int clust;

                            for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){

                                if (IS_ODD(k)) tmp = cc->filter3[k/2] >> 4;
                                else tmp = cc->filter3[k/2] & 0xf;

                                if((tmp & 0xfc) == suffix){
                                    nuc2add = tmp & 0x3;
                                    res[nuc2add].bucket = k/NB_CHILDREN_PER_SKP;
                                    clust = k/SIZE_CLUST_SKIP_NODES;
                                    uc = &(((UC*)cc->children)[res[nuc2add].bucket]);

                                    if (size_kmer != SIZE_SEED){
                                        if ((*func_on_types->is_child)((void*)cc, k)){
                                            res[nuc2add].pos_sub_bucket = (*func_on_types->count_children)((void*)cc, res[nuc2add].bucket*NB_CHILDREN_PER_SKP, k);
                                            res[nuc2add].pos_extra_filter3 = k;
                                            res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket*(func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot)]);
                                            res[nuc2add].container = &(((CC*)node->CC_array)[i]);
                                            res[nuc2add].children_type_leaf = 1;
                                        }
                                        else{
                                            //res->link_child = &(cc->children_Node_container[(*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), 0, k)]);
                                            if (clust-1 >= 0){
                                                if ((skip_node_root != NULL) && (k-last_k > k-clust*SIZE_CLUST_SKIP_NODES)){
                                                    last_node = skip_node_root[i][clust-1] + (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), clust*SIZE_CLUST_SKIP_NODES, k);
                                                }
                                                else last_node += (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), last_k, k);
                                            }
                                            else last_node += (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), last_k, k);

                                            res[nuc2add].link_child = &(cc->children_Node_container[last_node]);
                                            last_k = k;
                                        }
                                    }
                                    else{
                                        res[nuc2add].pos_extra_filter3 = k;
                                        res[nuc2add].pos_sub_bucket = k%NB_CHILDREN_PER_SKP;
                                        res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket * uc->size_annot]);
                                        res[nuc2add].container = uc;
                                        res[nuc2add].posFilter2 = 0;
                                        res[nuc2add].posFilter3 = NB_CHILDREN_PER_SKP;
                                        res[nuc2add].children_type_leaf = 1;
                                    }
                                }
                                else if(tmp > next_suffix) break;
                            }
                        }
                        else{
                            for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){

                                if (IS_ODD(k)) tmp = cc->filter3[k/2] >> 4;
                                else tmp = cc->filter3[k/2] & 0xf;

                                if(tmp == suffix){
                                    nuc2add = tmp & 0x3;

                                    if (size_kmer != SIZE_SEED){
                                        if ((*func_on_types->is_child)((void*)cc, k)){

                                            res[nuc2add].bucket = k/NB_CHILDREN_PER_SKP;
                                            uc = &(((UC*)cc->children)[res[nuc2add].bucket]);

                                            res[nuc2add].pos_extra_filter3 = k;
                                            res[nuc2add].pos_sub_bucket = (*func_on_types->count_children)((void*)cc, (res[nuc2add].bucket)*NB_CHILDREN_PER_SKP, k);
                                            res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket*(func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot)]);
                                            res[nuc2add].container = &(((CC*)node->CC_array)[i]);
                                            res[nuc2add].children_type_leaf = 1;
                                        }
                                        else{
                                            if (last_count_node == -1) res[nuc2add].link_child = &(cc->children_Node_container[(*func_on_types->count_nodes)((void*)cc, 0, k)]);
                                            else res[nuc2add].link_child = &(cc->children_Node_container[MAX(cpt_node_tmp-1, 0)
                                                                             + (*func_on_types->count_nodes)((void*)cc, pos_extra_filter3, k)]);
                                        }
                                    }
                                    else{
                                        res[nuc2add].bucket = k/NB_CHILDREN_PER_SKP;
                                        uc = &(((UC*)cc->children)[res[nuc2add].bucket]);

                                        res[nuc2add].pos_extra_filter3 = k;
                                        res[nuc2add].pos_sub_bucket = k%NB_CHILDREN_PER_SKP;
                                        res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket * uc->size_annot]);
                                        res[nuc2add].container = uc;
                                        res[nuc2add].posFilter2 = 0;
                                        res[nuc2add].posFilter3 = NB_CHILDREN_PER_SKP;
                                        res[nuc2add].children_type_leaf = 1;
                                    }

                                    return;
                                }
                                else if(tmp > suffix) break;
                            }
                        }
                    }
                }

                return;
            }
            else return;
        }
        while ((((CC*)node->CC_array)[i].type & 0x1) == 0);
    }

    //If the prefix was not found in any CC Bloom filters, we search in the UC of the node
    if (node->UC_array.suffixes != NULL){
        int nb_cell = func_on_types->size_kmer_in_bytes;
        int size_line = nb_cell+node->UC_array.size_annot;

        int k;

        if (func_on_types->root == 1){
            i = 0;
            for (k=0; k<(node->UC_array.nb_children >> 1)*size_line; k+=size_line){ //Iterate over every suffix stored
                if (memcmp(&(node->UC_array.suffixes[k+1]), &(kmer[1]), (nb_cell-1)*sizeof(uint8_t)) == 0){
                    if ((node->UC_array.suffixes[k] & 0xfc) == (kmer[0] & 0xfc)){

                        nuc2add = node->UC_array.suffixes[k] & 0x3;
                        res[nuc2add].link_child = &(node->UC_array.suffixes[k]);
                        res[nuc2add].pos_sub_bucket = k/size_line;
                        res[nuc2add].container = &(node->UC_array);
                        res[nuc2add].posFilter2 = nb_cell;
                        res[nuc2add].posFilter3 = node->UC_array.nb_children >> 1;
                        res[nuc2add].container_is_UC = 1;
                        i++;
                        if (i == 4) return;
                    }
                }
            }
        }
        else{
            for (k=0; k<(node->UC_array.nb_children >> 1)*size_line; k+=size_line){ //Iterate over every suffix stored
                if (memcmp(&(node->UC_array.suffixes[k]), kmer, nb_cell) == 0){

                    nuc2add = node->UC_array.suffixes[k] & 0x3;
                    res[nuc2add].link_child = &(node->UC_array.suffixes[k]);
                    res[nuc2add].pos_sub_bucket = k/size_line;
                    res[nuc2add].container = &(node->UC_array);
                    res[nuc2add].posFilter2 = nb_cell;
                    res[nuc2add].posFilter3 = node->UC_array.nb_children >> 1;
                    res[nuc2add].container_is_UC = 1;
                    return;
                }
            }
        }
    }

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  presenceKmer(node, kmer, size_kmer, nb_CC_node, func_on_types)
*  ---------------------------------------------------------------------------------------------------------------
*  Search for a prefix into a node.
*  ---------------------------------------------------------------------------------------------------------------
*  node: pointer on a Node structure
*  kmer: pointer on the prefix to search
*  size_kmer: size of the prefix, in char.
*  nb_CC_node: number of CCs in this node
*  func_on_types: ptr on ptrs_on_func structure, contains information to manipulate CCs field CC->children_type
*  ---------------------------------------------------------------------------------------------------------------
*/
void presenceNeighborsRight(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types, resultPresence* res, uint16_t** skip_node_root){

    if (size_kmer == 0) ERROR("presenceNeighborsRight(): this case should not happen")
    ASSERT_NULL_PTR(node,"presenceNeighborsRight()")
    ASSERT_NULL_PTR(kmer,"presenceNeighborsRight()")

    //Allocate and initialize the structure resultPresence
    int last_count_node = -1;
    int cpt_node_tmp = -1;

    int i;

    uint8_t substring[3];
    uint8_t nuc2add;

    CC* cc;
    UC* uc;
    UC* uc_tmp;

    for (i=0; i<4; i++) initialize_resultPresence(&(res[i]));

    //Keep in resultPresence->substring only the prefix we are looking for
    if (size_kmer == SIZE_SEED){
        substring[0] = reverse_word_8(kmer[0]);
        substring[1] = reverse_word_8(kmer[1]) & 0xfc;
        substring[2] = reverse_word_8(kmer[2]) & 0xc0;
    }
    else{
        substring[0] = reverse_word_8(kmer[0]);
        substring[1] = reverse_word_8(kmer[1]);
        substring[2] = reverse_word_8(kmer[2]) & 0xc0;
    }

    //Compute hash-functions of the Bloom filter for the prefix
    uint16_t substring_prefix = (substring[0] << SIZE_CELL) | substring[1];
    uint16_t hash1_v = dbj2(substring[1], substring[0] & 0x3f, MODULO_HASH);
    uint16_t hash2_v = sdbm(substring[1], substring[0] & 0x3f, MODULO_HASH);

    //We first iterate on the CCs
    i = -1;

    if ((CC*)node->CC_array != NULL){
        do {
            i++;
            cc = &(((CC*)node->CC_array)[i]);

            //First look if the prefix is in the Bloom filter
            if ((cc->BF_filter2[hash1_v/SIZE_CELL] & MASK_POWER_8[hash1_v%SIZE_CELL]) == 0) continue;
            if ((cc->BF_filter2[hash2_v/SIZE_CELL] & MASK_POWER_8[hash2_v%SIZE_CELL]) == 0) continue;

            substring[0] = substring_prefix >> 6;
            substring[1] = (substring_prefix << 2) | (substring[2] >> 6);
            substring[2] = (substring_prefix >> 8) & 0xc0;

            //At this point, this prefix is said present by the BF, maybe it's a True Positive, maybe it's a False Positive
            //If it is a FP, the prefix will have to be inserted in this CC (FP recycling)
            uint16_t size_bf = cc->type >> 8;
            uint8_t s = (cc->type >> 2) & 0x3f; //length v of the prefixes stored in this CC (p_vs)
            uint8_t p = SIZE_SEED*2-s; //length u of the prefixes stored in this CC (p_us)

            uint8_t word_tmp;

            int size_filter2 = size_bf+(MASK_POWER_16[p]/SIZE_CELL); //Size BF + filter2 in bytes
            int size_filter2_n_skip = size_filter2; //Size BF + filter2 + SkipFilter2 (SkipFilter2 may not exist) in bytes
            int skip_filter2 = MASK_POWER_16[p]/0xf8;
            if (cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM) size_filter2_n_skip += skip_filter2;

            //Compute p_u, the first u char. of the prefix p we are looking for
            int posFilter2 = 0;
            if (p==10) posFilter2 = (((uint16_t)substring[0]) << 2) | (((uint16_t)substring[1]) >> 6);
            else posFilter2 = (((uint16_t)substring[0]) << 6) | (((uint16_t)substring[1]) >> 2);

            //If the prefix is present in filter2, need to compute the Hamming weight between position 0 and the position of p_u
            //in order to know how many p_u are lexicographically inferior or equal to the one we insert
            if ((cc->BF_filter2[size_bf+posFilter2/SIZE_CELL] & MASK_POWER_8[posFilter2%SIZE_CELL]) != 0){

                int k=size_bf+posFilter2/SIZE_CELL;
                int cnt = 0;
                int hamming_weight = 0;

                //If SkipFilter2 exists, we use it to accelerate the computation of the Hamming weight
                if (cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM){
                    int skip_posfilter2 = posFilter2/0xf8;

                    while ((cnt<skip_posfilter2) && (cnt < skip_filter2)){
                        hamming_weight += cc->BF_filter2[size_filter2+cnt];
                        cnt++;
                    }
                }

                //finish the computation of the Hamming weight
                hamming_weight += popcnt_8_par(cc->BF_filter2, size_bf+cnt*31, k);

                word_tmp = cc->BF_filter2[k];
                for (k=0; k<=posFilter2%8; k++)
                    if ((word_tmp & MASK_POWER_8[k]) != 0) hamming_weight++;

                posFilter2 = hamming_weight; //posFilter2 now contains this Hamming weight

                //At this point, we know that the p_v we are looking for should be situated in the third filter
                //at a slot (plus subsequent one) corresponding to the hamming_weight-th 1 in the extra_filter3
                int nb_skp = CEIL(cc->nb_elem, NB_CHILDREN_PER_SKP);
                int pos_extra_filter3 = INT_MAX;
                int hamming_weight_0 = 0;

                int j=0, m=0, sum=0;
                int nb_cell_3rdlist = CEIL(cc->nb_elem,SIZE_CELL);

                k=0;
                word_tmp = 0;
                hamming_weight = 0;

                //This loop tries to use SkipFilter3 (if it exists) to accelerate the Hamming weight comput.
                //in extra_filter3 by skipping counting 1s directly in extra_filter3
                while(m<cc->nb_elem/0xf8){
                    if ((sum = hamming_weight + cc->BF_filter2[size_filter2_n_skip+m]) < posFilter2){
                        hamming_weight = sum;
                        m++;
                    }
                    else break;
                }

                //Finish the Hamming weight computation directly in extra_filter3
                if (func_on_types->level_min == 1){
                    for (k=m*31; k<nb_cell_3rdlist; k++){

                        if ((sum = hamming_weight+popcnt_8(cc->extra_filter3[k])) >= posFilter2){
                            word_tmp = cc->extra_filter3[k];

                            int size_word = 7;
                            if (k == nb_cell_3rdlist-1) size_word = (cc->nb_elem-1)%8;

                            for (j=0; j<=size_word; j++){
                                if (((word_tmp >> j)&1) == 1){
                                    hamming_weight++;
                                    if (hamming_weight == posFilter2){
                                        pos_extra_filter3 = k*8+j;
                                        j++;
                                        break;
                                    }
                                }
                            }

                            if (k == nb_cell_3rdlist) k--;
                            goto MATCH;
                        }
                        else hamming_weight = sum;
                    }
                    //If the substring must be inserted in the last slot of the list
                    if (k == nb_cell_3rdlist) k--;

                    MATCH: if (pos_extra_filter3 == INT_MAX) pos_extra_filter3 = cc->nb_elem;

                    int a=0;
                    for (a=k; a<nb_cell_3rdlist; a++){
                        uint8_t word_tmp = cc->extra_filter3[a];
                        while (j<8){
                            if (((word_tmp >> j)&1) == 0){
                                if (a*8+j < cc->nb_elem) hamming_weight_0++;
                                else goto OUT_LOOP;
                            }
                            else goto OUT_LOOP;
                            j++;
                        }
                        j=0;
                    }
                }
                else{
                    int size_line_children;
                    int pos_children;
                    int pos_children_tmp;
                    int cpt_pv;
                    int cpt_node;
                    int nb_elem_in_pv;
                    int end, end_tmp;
                    int it, it_tmp;

                    int cpt_pv_tmp = 0;
                    int tmp_hamming_weight = hamming_weight;

                    uint8_t* children;

                    k = m*0xf8;
                    pos_children_tmp = k/NB_CHILDREN_PER_SKP;

                    (*func_on_types->count_Nodes_Children)(&(((CC*)node->CC_array)[i]), pos_children_tmp*NB_CHILDREN_PER_SKP, k, &cpt_pv_tmp, &cpt_node_tmp);
                    last_count_node = (cpt_node_tmp += (*func_on_types->count_nodes)(&(((CC*)node->CC_array)[i]), 0, pos_children_tmp*NB_CHILDREN_PER_SKP));

                    it = k - pos_children_tmp * NB_CHILDREN_PER_SKP;
                    it_tmp = it;
                    nb_elem_in_pv = 0;

                    cpt_pv = cpt_pv_tmp;
                    cpt_node = cpt_node_tmp;
                    pos_children = pos_children_tmp;

                    while (pos_children < nb_skp){

                        if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * NB_CHILDREN_PER_SKP;
                        else end = NB_CHILDREN_PER_SKP;

                        uc = &(((UC*)cc->children)[pos_children]);
                        size_line_children = func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot;
                        children = &(uc->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                        while (it < end){

                            it_tmp = it;

                            for (m=k; m<MIN(k+8, cc->nb_elem); m++){

                                if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, m)) == 0){
                                    tmp_hamming_weight += cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1;
                                    cpt_node++;
                                }
                                else{
                                    tmp_hamming_weight += children[cpt_pv*size_line_children] >> 7;
                                    cpt_pv += nb_elem_in_pv;
                                }

                                it++;
                            }

                            if (tmp_hamming_weight >= posFilter2){

                                if (pos_children_tmp == nb_skp - 1) end_tmp = cc->nb_elem - pos_children_tmp * NB_CHILDREN_PER_SKP;
                                else end_tmp = NB_CHILDREN_PER_SKP;

                                uc_tmp = &(((UC*)cc->children)[pos_children_tmp]);
                                size_line_children = func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot;
                                children = &(uc_tmp->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                                while (cpt_pv_tmp < uc_tmp->nb_children){

                                    if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, k)) == 0){
                                        hamming_weight += cc->children_Node_container[cpt_node_tmp].UC_array.nb_children & 0x1;
                                        cpt_node_tmp++;
                                    }
                                    else{
                                        hamming_weight += children[cpt_pv_tmp*size_line_children] >> 7;
                                        cpt_pv_tmp += nb_elem_in_pv;
                                    }

                                    if (hamming_weight == posFilter2){
                                        pos_extra_filter3 = k;
                                        goto MATCH2;
                                    }

                                    k++;
                                    it_tmp++;
                                }
                            }
                            else{
                                hamming_weight = tmp_hamming_weight;
                                cpt_pv_tmp = cpt_pv;
                                cpt_node_tmp = cpt_node;
                                pos_children_tmp = pos_children;
                                k = m;
                            }
                        }

                        it = 0;
                        it_tmp = 0;
                        cpt_pv = 0;
                        cpt_pv_tmp = cpt_pv;
                        pos_children++;
                        pos_children_tmp = pos_children;
                    }

                    it_tmp = it;

                    MATCH2: if (pos_extra_filter3 == INT_MAX) pos_extra_filter3 = cc->nb_elem;

                    if (pos_children_tmp < nb_skp){

                        if (pos_children_tmp == nb_skp - 1) end_tmp = cc->nb_elem - pos_children_tmp * NB_CHILDREN_PER_SKP;
                        else end_tmp = NB_CHILDREN_PER_SKP;

                        it_tmp++;
                        k++;

                        if (it_tmp >= end_tmp){
                            it = 0;
                            cpt_pv = 0;
                            pos_children = pos_children_tmp+1;
                        }
                        else{
                            it = it_tmp;
                            cpt_pv = cpt_pv_tmp;
                            pos_children = pos_children_tmp;
                        }

                        cpt_node = cpt_node_tmp;

                        while (pos_children < nb_skp){

                            uc = &(((UC*)cc->children)[pos_children]);
                            size_line_children = func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot;
                            children = &(uc->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                            if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * NB_CHILDREN_PER_SKP;
                            else end = NB_CHILDREN_PER_SKP;

                            while (it < end){

                                if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, k)) == 0){
                                    if ((cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1) == 0) hamming_weight_0++;
                                    else goto OUT_LOOP;
                                    cpt_node++;
                                }
                                else{
                                    if (children[cpt_pv*size_line_children] < 0x80) hamming_weight_0++;
                                    else goto OUT_LOOP;
                                    cpt_pv += nb_elem_in_pv;
                                }

                                k++;
                                it++;
                            }

                            it = 0;
                            cpt_pv = 0;
                            pos_children++;
                        }
                    }
                }

                //We now compare the p_vs between start position and end position in filter3 to see if
                //the prefix we look for is present
                OUT_LOOP:if (pos_extra_filter3 < cc->nb_elem){
                    uint8_t suffix;
                    uint8_t tmp;

                    if (s==8){ //if the length of p_vs is 8bits
                        suffix = (substring[1] << 2) | (substring[2] >> 6); //compute p_v

                        if (size_kmer == SIZE_SEED){

                            uint8_t next_suffix = suffix | 0xc;

                            for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){ //iterate between start and end position

                                if((cc->filter3[k] & 0xf3) == suffix){ //if there is a match

                                    nuc2add = (cc->filter3[k] >> 2) & 0x3;
                                    res[nuc2add].bucket = k/NB_CHILDREN_PER_SKP;
                                    uc = &(((UC*)cc->children)[res[nuc2add].bucket]);

                                    res[nuc2add].pos_sub_bucket = k%NB_CHILDREN_PER_SKP;
                                    res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket * uc->size_annot]);
                                    res[nuc2add].container = uc;
                                    res[nuc2add].posFilter2 = 0;
                                    res[nuc2add].posFilter3 = NB_CHILDREN_PER_SKP;
                                }
                                else if(cc->filter3[k] > next_suffix) break;
                            }
                        }
                        else{
                            for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){ //iterate between start and end position

                                if(cc->filter3[k] == suffix){ //if there is a match

                                    nuc2add = cc->filter3[k] & 0x3;

                                    if ((*func_on_types->is_child)((void*)cc, k) == 1){

                                        res[nuc2add].bucket = k/NB_CHILDREN_PER_SKP;
                                        uc = &(((UC*)cc->children)[res[nuc2add].bucket]);

                                        res[nuc2add].pos_extra_filter3 = k;
                                        res[nuc2add].pos_sub_bucket = (*func_on_types->count_children)((void*)cc, res[nuc2add].bucket*NB_CHILDREN_PER_SKP, k);
                                        res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket*(func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot)]);
                                        res[nuc2add].container = &(((CC*)node->CC_array)[i]);
                                        res[nuc2add].children_type_leaf = 1;
                                    }
                                    else{
                                        if ((skip_node_root != NULL) && (func_on_types->root == 1) && (last_count_node == -1)){

                                            int clust = k/SIZE_CLUST_SKIP_NODES;

                                            if (clust-1 >= 0){
                                                last_count_node = skip_node_root[i][clust-1] +
                                                                    (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), clust*SIZE_CLUST_SKIP_NODES, k);
                                            }
                                            else last_count_node = (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), 0, k);

                                            res[nuc2add].link_child = &(cc->children_Node_container[last_count_node]);
                                        }
                                        else{
                                            if (last_count_node == -1){
                                                res[nuc2add].link_child = &(cc->children_Node_container[(*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), 0, k)]);
                                            }
                                            else res[nuc2add].link_child = &(cc->children_Node_container[MAX(cpt_node_tmp-1, 0) +
                                                                             (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), pos_extra_filter3, k)]);

                                            if (cc->nb_Node_children == 0) ERROR("presenceKmer(s=8): the count of nodes in this CC is 0, but some nodes are detected\n");
                                        }
                                    }

                                    return;
                                }
                                else if(cc->filter3[k] > suffix) break;
                            }
                        }
                    }
                    else {
                        suffix = ((substring[1] & 0x3) << 2) | (substring[2] >> 6);
                        if (size_kmer == SIZE_SEED){

                            uint8_t next_suffix = suffix | 0xc;

                            for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){

                                if (IS_ODD(k)) tmp = cc->filter3[k/2] >> 4;
                                else tmp = cc->filter3[k/2] & 0xf;

                                if((tmp & 0x3) == suffix){

                                    nuc2add = (tmp >> 2) & 0x3;
                                    res[nuc2add].bucket = k/NB_CHILDREN_PER_SKP;
                                    uc = &(((UC*)cc->children)[res[nuc2add].bucket]);

                                    res[nuc2add].pos_sub_bucket = k%NB_CHILDREN_PER_SKP;
                                    res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket * uc->size_annot]);
                                    res[nuc2add].container = uc;
                                    res[nuc2add].posFilter2 = 0;
                                    res[nuc2add].posFilter3 = NB_CHILDREN_PER_SKP;
                                }
                                else if(tmp > next_suffix) break;
                            }
                        }
                        else{
                            for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){

                                if (IS_ODD(k)) tmp = cc->filter3[k/2] >> 4;
                                else tmp = cc->filter3[k/2] & 0xf;

                                if(tmp == suffix){

                                    nuc2add = tmp & 0x3;
                                    if ((*func_on_types->is_child)((void*)cc, k)){

                                        res[nuc2add].bucket = k/NB_CHILDREN_PER_SKP;
                                        uc = &(((UC*)cc->children)[res[nuc2add].bucket]);

                                        res[nuc2add].pos_extra_filter3 = k;
                                        res[nuc2add].pos_sub_bucket = (*func_on_types->count_children)((void*)cc, res[nuc2add].bucket*NB_CHILDREN_PER_SKP, k);
                                        res[nuc2add].link_child = &(uc->suffixes[res[nuc2add].pos_sub_bucket*(func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot)]);
                                        res[nuc2add].container = &(((CC*)node->CC_array)[i]);
                                        res[nuc2add].children_type_leaf = 1;
                                    }
                                    else{
                                        if ((skip_node_root != NULL) && (func_on_types->root == 1) && (last_count_node == -1)){

                                            int clust = k/SIZE_CLUST_SKIP_NODES;

                                            if (clust-1 >= 0){
                                                last_count_node = skip_node_root[i][clust-1] +
                                                                    (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), clust*SIZE_CLUST_SKIP_NODES, k);
                                            }
                                            else last_count_node = (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), 0, k);

                                            res[nuc2add].link_child = &(cc->children_Node_container[last_count_node]);
                                        }
                                        else{
                                            if (last_count_node == -1){
                                                res[nuc2add].link_child = &(cc->children_Node_container[(*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), 0, k)]);
                                            }
                                            else res[nuc2add].link_child = &(cc->children_Node_container[MAX(cpt_node_tmp-1, 0) +
                                                                             (*func_on_types->count_nodes)((void*)&(((CC*)node->CC_array)[i]), pos_extra_filter3, k)]);

                                            if (cc->nb_Node_children == 0) ERROR("presenceKmer(s=8): the count of nodes in this CC is 0, but some nodes are detected\n");
                                        }
                                    }

                                    return;
                                }
                                else if(tmp > suffix) break;
                            }
                        }
                    }
                }

                return;
            }
            else return;
        }
        while ((((CC*)node->CC_array)[i].type & 0x1) == 0);
    }

    //If the prefix was not found in any CC Bloom filters, we search in the UC of the node
    if (node->UC_array.suffixes != NULL){
        int nb_cell = func_on_types->size_kmer_in_bytes;
        int size_line = nb_cell+node->UC_array.size_annot;

        int k;
        uint8_t bits_left = (size_kmer*2)%SIZE_CELL;
        uint8_t mask2 = MASK_POWER_8[bits_left] - MASK_POWER_8[bits_left-2];
        uint8_t mask = func_on_types->mask_shift_kmer;

        if (mask == 0xff) mask = 0;

        i = 0;
        bits_left -= 2;

        for (k=0; k<(node->UC_array.nb_children >> 1)*size_line; k+=size_line){ //Iterate over every suffix stored

            if (memcmp(&(node->UC_array.suffixes[k]), kmer, (nb_cell-1)*sizeof(uint8_t)) == 0){

                if ((node->UC_array.suffixes[k+nb_cell-1] & mask) == kmer[nb_cell-1]){

                    nuc2add = (node->UC_array.suffixes[k+nb_cell-1] & mask2) >> bits_left;
                    res[nuc2add].link_child = &(node->UC_array.suffixes[k]);
                    res[nuc2add].container = &(node->UC_array);
                    res[nuc2add].posFilter2 = nb_cell;
                    res[nuc2add].posFilter3 = node->UC_array.nb_children >> 1;
                    res[nuc2add].pos_sub_bucket = k/size_line;
                    res[nuc2add].container_is_UC = 1;
                    i++;
                    if (i == 4) return;

                }
            }
        }
    }

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  presenceKmer(node, kmer, size_kmer, nb_CC_node, func_on_types)
*  ---------------------------------------------------------------------------------------------------------------
*  Search for a prefix into a node.
*  ---------------------------------------------------------------------------------------------------------------
*  node: pointer on a Node structure
*  kmer: pointer on the prefix to search
*  size_kmer: size of the prefix, in char.
*  nb_CC_node: number of CCs in this node
*  func_on_types: ptr on ptrs_on_func structure, contains information to manipulate CCs field CC->children_type
*  ---------------------------------------------------------------------------------------------------------------
*/
void presenceKmer(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types, resultPresence* res){

    if (size_kmer == 0) ERROR("presenceKmer(): this case should not happen")
    ASSERT_NULL_PTR(node,"presenceKmer()")
    ASSERT_NULL_PTR(kmer,"presenceKmer()")

    int i = -1;
    int cpt_node_tmp = -1;

    CC* cc;
    UC* uc;
    UC* uc_tmp;

    //Allocate and initialize the structure resultPresence
    initialize_resultPresence(res);

    //Keep in resultPresence->substring only the prefix we are looking for
    res->substring[0] = reverse_word_8(kmer[0]);
    res->substring[1] = reverse_word_8(kmer[1]);
    res->substring[2] = reverse_word_8(kmer[2]) & 0xc0;

    //Compute hash-functions of the Bloom filter for the prefix
    uint16_t substring_prefix = (res->substring[0] << SIZE_CELL) | res->substring[1];
    uint16_t hash1_v = dbj2(res->substring[1], res->substring[0] & 0x3f, MODULO_HASH);
    uint16_t hash2_v = sdbm(res->substring[1], res->substring[0] & 0x3f, MODULO_HASH);

    //We first iterate on the CCs
    if ((CC*)node->CC_array != NULL){
        do {
            i++;

            cc = &(((CC*)node->CC_array)[i]);

            //First look if the prefix is in the Bloom filter
            if ((cc->BF_filter2[hash1_v/SIZE_CELL] & MASK_POWER_8[hash1_v%SIZE_CELL]) == 0) continue;
            if ((cc->BF_filter2[hash2_v/SIZE_CELL] & MASK_POWER_8[hash2_v%SIZE_CELL]) == 0) continue;

            res->substring[0] = (substring_prefix >> 6) & 0xff;
            res->substring[1] = (substring_prefix << 2) | (res->substring[2] >> 6);
            res->substring[2] = (substring_prefix >> 8) & 0xc0;

            res->container = &(((CC*)node->CC_array)[i]);
            res->presBF = 1;

            //At this point, this prefix is said present by the BF, maybe it's a True Positive, maybe it's a False Positive
            //If it is a FP, the prefix will have to be inserted in this CC (FP recycling)
            uint16_t size_bf = cc->type >> 8;
            uint8_t s = (cc->type >> 2) & 0x3f; //length v of the prefixes stored in this CC (p_vs)
            uint8_t p = SIZE_SEED*2-s; //length u of the prefixes stored in this CC (p_us)

            uint8_t word_tmp;

            int size_filter2 = size_bf+(MASK_POWER_16[p]/SIZE_CELL); //Size BF + filter2 in bytes
            int size_filter2_n_skip = size_filter2; //Size BF + filter2 + SkipFilter2 (SkipFilter2 may not exist) in bytes
            int skip_filter2 = MASK_POWER_16[p]/0xf8;
            if (cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM) size_filter2_n_skip += skip_filter2;

            //Compute p_u, the first u char. of the prefix p we are looking for
            int posFilter2 = 0;
            if (p==10) posFilter2 = (((uint16_t)res->substring[0]) << 2) | (((uint16_t)res->substring[1]) >> 6);
            else posFilter2 = (((uint16_t)res->substring[0]) << 6) | (((uint16_t)res->substring[1]) >> 2);

            //If the prefix is present in filter2, need to compute the Hamming weight between position 0 and the position of p_u
            //in order to know how many p_u are lexicographically inferior or equal to the one we insert
            if ((cc->BF_filter2[size_bf+posFilter2/SIZE_CELL] & MASK_POWER_8[posFilter2%SIZE_CELL]) != 0){

                int k=size_bf+posFilter2/SIZE_CELL;
                int cnt = 0;
                int hamming_weight = 0;

                //If SkipFilter2 exists, we use it to accelerate the computation of the Hamming weight
                if (cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM){
                    int skip_posfilter2 = posFilter2/0xf8;

                    while ((cnt<skip_posfilter2) && (cnt < skip_filter2)){
                        hamming_weight += cc->BF_filter2[size_filter2+cnt];
                        cnt++;
                    }
                }

                //finish the computation of the Hamming weight
                hamming_weight += popcnt_8_par(cc->BF_filter2, size_bf+cnt*31, k);

                word_tmp = cc->BF_filter2[k];
                for (k=0; k<=posFilter2%8; k++)
                    if ((word_tmp & MASK_POWER_8[k]) != 0) hamming_weight++;

                res->presFilter2 = 1;
                res->posFilter2 = hamming_weight; //res->posFilter2 now contains this Hamming weight

                //At this point, we know that the p_v we are looking for should be situated in the third filter
                //at a slot (plus subsequent one) corresponding to the hamming_weight-th 1 in the extra_filter3
                int nb_skp = CEIL(cc->nb_elem, NB_CHILDREN_PER_SKP);
                int pos_extra_filter3 = INT_MAX;
                int hamming_weight_0 = 0;

                int j=0, m=0, sum=0;
                int nb_cell_3rdlist = CEIL(cc->nb_elem,SIZE_CELL);

                k=0;
                word_tmp = 0;
                hamming_weight = 0;

                //This loop tries to use SkipFilter3 (if it exists) to accelerate the Hamming weight comput.
                //in extra_filter3 by skipping counting 1s directly in extra_filter3
                while(m<cc->nb_elem/0xf8){
                    if ((sum = hamming_weight + cc->BF_filter2[size_filter2_n_skip+m]) < res->posFilter2){
                        hamming_weight = sum;
                        m++;
                    }
                    else break;
                }

                //Finish the Hamming weight computation directly in extra_filter3
                if (func_on_types->level_min == 1){
                    for (k=m*31; k<nb_cell_3rdlist; k++){

                        if ((sum = hamming_weight+popcnt_8(cc->extra_filter3[k])) >= res->posFilter2){
                            word_tmp = cc->extra_filter3[k];

                            int size_word = 7;
                            if (k == nb_cell_3rdlist-1) size_word = (cc->nb_elem-1)%8;

                            for (j=0; j<=size_word; j++){
                                if (((word_tmp >> j)&1) == 1){
                                    hamming_weight++;
                                    if (hamming_weight == res->posFilter2){
                                        pos_extra_filter3 = k*8+j;
                                        j++;
                                        break;
                                    }
                                }
                            }

                            if (k == nb_cell_3rdlist) k--;
                            goto MATCH;
                        }
                        else hamming_weight = sum;
                    }
                    //If the substring must be inserted in the last slot of the list
                    if (k == nb_cell_3rdlist) k--;

                    MATCH: if (pos_extra_filter3 == INT_MAX) pos_extra_filter3 = cc->nb_elem;

                    int a=0;
                    for (a=k; a<nb_cell_3rdlist; a++){
                        uint8_t word_tmp = cc->extra_filter3[a];
                        while (j<8){
                            if (((word_tmp >> j)&1) == 0){
                                if (a*8+j < cc->nb_elem) hamming_weight_0++;
                                else goto OUT_LOOP;
                            }
                            else goto OUT_LOOP;
                            j++;
                        }
                        j=0;
                    }

                    OUT_LOOP:res->pos_extra_filter3 = pos_extra_filter3;
                }
                else{
                    int size_line_children;
                    int pos_children;
                    int cpt_pv;
                    int cpt_node;
                    int nb_elem_in_pv;
                    int end, end_tmp;
                    int it, it_tmp;

                    int tmp_hamming_weight = hamming_weight;

                    uint8_t* children;

                    k = m*0xf8;
                    res->pos_children = k/NB_CHILDREN_PER_SKP;

                    (*func_on_types->count_Nodes_Children)(res->container, res->pos_children * NB_CHILDREN_PER_SKP, k, &(res->count_children), &(res->count_nodes));
                    cpt_node_tmp = (res->count_nodes += (*func_on_types->count_nodes)(res->container, 0, res->pos_children*NB_CHILDREN_PER_SKP));

                    nb_elem_in_pv = 0;
                    it = k - res->pos_children * NB_CHILDREN_PER_SKP;
                    it_tmp = it;

                    cpt_pv = res->count_children;
                    cpt_node = res->count_nodes;
                    pos_children = res->pos_children;

                    while (pos_children < nb_skp){

                        if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * NB_CHILDREN_PER_SKP;
                        else end = NB_CHILDREN_PER_SKP;

                        uc = &(((UC*)cc->children)[pos_children]);
                        size_line_children = func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot;
                        children = &(uc->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                        while (it < end){

                            it_tmp = it;

                            for (m=k; m<MIN(k+8, cc->nb_elem); m++){

                                if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, m)) == 0){
                                    tmp_hamming_weight += cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1;
                                    cpt_node++;
                                }
                                else{
                                    tmp_hamming_weight += children[cpt_pv*size_line_children] >> 7;
                                    cpt_pv += nb_elem_in_pv;
                                }

                                it++;
                            }

                            if (tmp_hamming_weight >= res->posFilter2){

                                uc_tmp = &(((UC*)cc->children)[res->pos_children]);
                                size_line_children = func_on_types->size_kmer_in_bytes_minus_1 + uc_tmp->size_annot;
                                children = &(uc_tmp->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                                if (res->pos_children == nb_skp - 1) end_tmp = cc->nb_elem - res->pos_children * NB_CHILDREN_PER_SKP;
                                else end_tmp = NB_CHILDREN_PER_SKP;

                                while (it_tmp < end_tmp){

                                    if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, k)) == 0){
                                        hamming_weight += cc->children_Node_container[res->count_nodes].UC_array.nb_children & 0x1;
                                        res->count_nodes++;
                                    }
                                    else{
                                        hamming_weight += children[res->count_children*size_line_children] >> 7;
                                        res->count_children += nb_elem_in_pv;
                                    }

                                    if (hamming_weight == res->posFilter2){
                                        pos_extra_filter3 = k;
                                        goto MATCH2;
                                    }

                                    k++;
                                    it_tmp++;
                                }
                            }
                            else{
                                hamming_weight = tmp_hamming_weight;
                                res->count_children = cpt_pv;
                                res->count_nodes = cpt_node;
                                res->pos_children = pos_children;
                                k = m;
                                //it_tmp = it;
                            }
                        }

                        it = 0;
                        it_tmp = 0;
                        cpt_pv = 0;
                        res->count_children = 0;

                        pos_children++;
                        res->pos_children = pos_children;
                    }

                    it_tmp = it;

                    MATCH2: if (pos_extra_filter3 == INT_MAX) pos_extra_filter3 = cc->nb_elem;

                    if (res->pos_children < nb_skp){

                        if (res->pos_children == nb_skp - 1) end_tmp = cc->nb_elem - res->pos_children * NB_CHILDREN_PER_SKP;
                        else end_tmp = NB_CHILDREN_PER_SKP;

                        it_tmp++;
                        k++;

                        //if (res->count_children >= ((UC*)cc->children)[res->pos_children].nb_children){
                        if (it_tmp >= end_tmp){
                            it = 0;
                            cpt_pv = 0;
                            pos_children = res->pos_children+1;
                        }
                        else{
                            it = it_tmp;
                            cpt_pv = res->count_children;
                            pos_children = res->pos_children;
                        }

                        cpt_node = res->count_nodes;

                        while (pos_children < nb_skp){

                            uc = &(((UC*)cc->children)[pos_children]);
                            size_line_children = func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot;
                            children = &(uc->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                            if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * NB_CHILDREN_PER_SKP;
                            else end = NB_CHILDREN_PER_SKP;

                            while (it < end){

                                if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, k)) == 0){
                                    if ((cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1) == 0) hamming_weight_0++;
                                    else goto OUT_LOOP2;
                                    cpt_node++;
                                }
                                else{
                                    if (children[cpt_pv*size_line_children] < 0x80) hamming_weight_0++;
                                    else goto OUT_LOOP2;
                                    cpt_pv += nb_elem_in_pv;
                                }

                                k++;
                                it++;
                            }

                            it = 0;
                            cpt_pv = 0;
                            pos_children++;
                        }
                    }

                    OUT_LOOP2:res->pos_extra_filter3 = pos_extra_filter3;
                }

                //We now compare the p_vs between start position and end position in filter3 to see if
                //the prefix we look for is present
                if (pos_extra_filter3 < cc->nb_elem){
                    if (s==8){ //if the length of p_vs is 8bits
                        uint8_t suffix = (res->substring[1] << 2) | (res->substring[2] >> 6); //compute p_v
                        for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){ //iterate between start and end position
                            if(cc->filter3[k] == suffix){ //if there is a match
                                res->presFilter3 = 1;
                                res->posFilter3 = k;
                                if (size_kmer != SIZE_SEED){ //if the current CC is not a leaf
                                    if ((*func_on_types->is_child)((void*)cc, k) == 1){
                                        res->bucket = k/NB_CHILDREN_PER_SKP;
                                        res->pos_sub_bucket = (*func_on_types->count_children)((void*)cc, (res->bucket)*NB_CHILDREN_PER_SKP, k);

                                        uc = &(((UC*)cc->children)[res->bucket]);

                                        res->link_child = &(uc->suffixes[res->pos_sub_bucket*(func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot)]);
                                        res->children_type_leaf = 1;
                                    }
                                    else{
                                        if (cpt_node_tmp == -1) res->link_child = &(cc->children_Node_container[(*func_on_types->count_nodes)((void*)cc, 0, k)]);
                                        else{
                                            res->link_child = &(cc->children_Node_container[MAX(res->count_nodes-1, 0) + (*func_on_types->count_nodes)((void*)cc, pos_extra_filter3, k)]);
                                        }

                                        if (cc->nb_Node_children == 0) ERROR( "presenceKmer(s=8): the count of nodes in this CC is 0, but some nodes are detected" )
                                    }
                                }
                                else{
                                    res->bucket = k/NB_CHILDREN_PER_SKP;
                                    res->pos_sub_bucket = k%NB_CHILDREN_PER_SKP;

                                    uc = &(((UC*)cc->children)[res->bucket]);

                                    res->link_child = &(uc->suffixes[res->pos_sub_bucket * uc->size_annot]);
                                }

                                return;
                            }
                            else if(cc->filter3[k] > suffix){
                                res->posFilter3 = k;
                                break;
                            }
                        }
                        if (k>pos_extra_filter3+hamming_weight_0) k--;
                        if (cc->filter3[k] < suffix) res->posFilter3 = k+1;
                    }
                    else {
                        uint8_t suffix = ((res->substring[1] & 0x3) << 2) | (res->substring[2] >> 6);
                        uint8_t tmp = 0;

                        for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){

                            if (IS_ODD(k)) tmp = cc->filter3[k/2] >> 4;
                            else tmp = cc->filter3[k/2] & 15;

                            if(tmp == suffix){

                                res->presFilter3 = 1;
                                res->posFilter3 = k;
                                //Need to determine which child
                                if (size_kmer != SIZE_SEED){
                                    if ((*func_on_types->is_child)((void*)cc, k)){

                                        res->bucket = k/NB_CHILDREN_PER_SKP;
                                        uc = &(((UC*)cc->children)[res->bucket]);

                                        res->pos_sub_bucket = (*func_on_types->count_children)((void*)cc, (res->bucket)*NB_CHILDREN_PER_SKP, k);
                                        res->link_child = &(uc->suffixes[res->pos_sub_bucket*(func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot)]);
                                        res->children_type_leaf = 1;
                                    }
                                    else{
                                        if (cpt_node_tmp == -1) res->link_child = &(cc->children_Node_container[(*func_on_types->count_nodes)((void*)cc, 0, k)]);
                                        else res->link_child = &(cc->children_Node_container[MAX(res->count_nodes-1, 0) + (*func_on_types->count_nodes)((void*)cc, pos_extra_filter3, k)]);

                                        if (cc->nb_Node_children == 0) ERROR ("presenceKmer(s=4): the count of nodes in this CC is 0, but some nodes are detected")
                                    }
                                }
                                else{
                                    res->bucket = k/NB_CHILDREN_PER_SKP;
                                    uc = &(((UC*)cc->children)[res->bucket]);

                                    res->pos_sub_bucket = k%NB_CHILDREN_PER_SKP;
                                    res->link_child = &(uc->suffixes[res->pos_sub_bucket * uc->size_annot]);
                                }

                                return;
                            }
                            else if(tmp > suffix){
                                res->posFilter3 = k;
                                break;
                            }
                        }

                        if (k>pos_extra_filter3+hamming_weight_0) k--;

                        if (IS_ODD(k)){
                            if((cc->filter3[k/2] >> 4) < suffix) res->posFilter3 = k+1;
                        }
                        else if((cc->filter3[k/2] & 15) < suffix) res->posFilter3 = k+1;
                    }
                }

                return;
            }
            else return;
        }
        while ((((CC*)node->CC_array)[i].type & 0x1) == 0);
    }

    //If the prefix was not found in any CC Bloom filters, we search in the UC of the node
    if (node->UC_array.suffixes != NULL){
        int nb_cell = func_on_types->size_kmer_in_bytes;
        int size_line = nb_cell+node->UC_array.size_annot;

        int k=0;
        for (k=0; k<(node->UC_array.nb_children >> 1)*size_line; k+=size_line){ //Iterate over every suffix stored
            if (memcmp(&(node->UC_array.suffixes[k]), kmer, nb_cell) == 0){ //if there is a match
                //Fill in the result_Presence structure
                res->container = &(node->UC_array);
                res->link_child = &(node->UC_array.suffixes[k]);
                res->pos_sub_bucket = k/size_line;
                res->container_is_UC = 1;
                res->posFilter2 = nb_cell;
                return;
            }
        }
    }

    return;
}

resultPresence* isKmerPresent(Node* restrict node, uint8_t* restrict kmer, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(node,"isKmerPresent()")
    ASSERT_NULL_PTR(kmer,"isKmerPresent()")
    ASSERT_NULL_PTR(func_on_types,"isKmerPresent()")

    int level = (size_kmer/SIZE_SEED)-1;
    int j=0;
    int nb_cell;
    int size_line;

    CC* cc;
    UC* uc;

    uint16_t nb_elt;

    uint8_t kmer_tmp[func_on_types[level].size_kmer_in_bytes];
    memcpy(kmer_tmp, kmer, func_on_types[level].size_kmer_in_bytes);

    __builtin_prefetch (&(func_on_types[level]), 0, 0);

    resultPresence* res = malloc(sizeof(resultPresence));
    ASSERT_NULL_PTR(res,"isKmerPresent()")
    presenceKmer(node, kmer_tmp, size_kmer, &(func_on_types[level]), res);

    if (size_kmer == SIZE_SEED) return res;
    else{
        int nb_cell_to_delete = 2;
        nb_cell = func_on_types[level].size_kmer_in_bytes;

        if (size_kmer == 45) nb_cell_to_delete++;

        for (j=0; j < nb_cell - nb_cell_to_delete; j++){
            kmer_tmp[j] = kmer_tmp[j+2] >> 2;
            if (j+3 < nb_cell) kmer_tmp[j] |= kmer_tmp[j+3] << 6;
        }

        kmer_tmp[j-1] &= func_on_types[level].mask_shift_kmer;

        if (res->link_child != NULL){
            if (res->children_type_leaf == 0){
                if (res->container_is_UC == 0){
                    resultPresence* res_tmp = isKmerPresent((Node*)res->link_child, kmer_tmp, size_kmer-SIZE_SEED, func_on_types);
                    free(res);
                    res = res_tmp;
                    return res;
                }
            }
            else{
                cc = (CC*)res->container;
                uc = &(((UC*)cc->children)[res->bucket]);

                nb_elt = (*func_on_types[level].getNbElts)(res->container, res->posFilter3);
                nb_cell = func_on_types[level].size_kmer_in_bytes_minus_1;
                size_line = nb_cell + uc->size_annot;

                res->link_child = NULL;
                res->container = NULL;

                if (size_kmer == 45){
                    for (j=res->pos_sub_bucket*size_line; j<(res->pos_sub_bucket+nb_elt)*size_line; j+=size_line){ //Iterate over every suffix stored
                        if (memcmp(&(uc->suffixes[j]), kmer_tmp, nb_cell*sizeof(uint8_t)) == 0){
                            res->link_child = &(uc->suffixes[j]);
                            res->container = uc;
                            res->pos_sub_bucket = j/size_line;
                            res->posFilter2 = nb_cell;
                            return res;
                        }
                    }
                }
                else{
                    for (j=res->pos_sub_bucket*size_line; j<(res->pos_sub_bucket+nb_elt)*size_line; j+=size_line){ //Iterate over every suffix stored
                        if (memcmp(&(uc->suffixes[j]), kmer_tmp, (nb_cell-1)*sizeof(uint8_t)) == 0){
                            if ((uc->suffixes[j+nb_cell-1] & 0x7f) == kmer_tmp[nb_cell-1]){
                                res->link_child = &(uc->suffixes[j]);
                                res->container = uc;
                                res->pos_sub_bucket = j/size_line;
                                res->posFilter2 = nb_cell;
                                return res;
                            }
                        }
                    }
                }
            }
        }
    }

    return res;
}
