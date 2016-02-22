#include "./../lib/insertNode.h"

/* ---------------------------------------------------------------------------------------------------------------
*  insertKmers(root, array_kmers, size_kmers, nb_kmers, id_genome, info_per_lvl, ann_inf, res, annot_sorted)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert kmers into the BFT
*  ---------------------------------------------------------------------------------------------------------------
*  root: ptr to the root of a BFT
*  array_kmers: array of kmers
*  size_kmers: length k of kmers in array_kmers
*  nb_kmers: number of kmers in array_kmers
*  id_genome: genome id to which belongs the kmers in array_kmers
*  info_per_lvl: ptr to info_per_level structure, used to manipulate compressed container field children_type
*  ann_inf: ptr to annotation_inform structure, used to make the transition between reading and modifying an annot
*  res: ptr to resultPresence structure, used to store information about the presence of a suffix prefix in a CC
*  ---------------------------------------------------------------------------------------------------------------
*/
void insertKmers(BFT_Root*  root, uint8_t*  array_kmers, int nb_kmers, uint32_t id_genome,
                 int size_id_genome){

    ASSERT_NULL_PTR(root,"insertKmers()")

    int i = 0;
    int curr_lvl = (root->k/NB_CHAR_SUF_PREF)-1;
    int nb_bytes = CEIL(root->k*2, SIZE_BITS_UINT_8T); //Nb bytes used to store a k-mer

    uint8_t kmer[nb_bytes];

    for (i=0; i<nb_kmers; i++){

        memcpy(kmer, &(array_kmers[i*nb_bytes]), nb_bytes * sizeof(uint8_t));

        insertKmer_Node(&(root->node), root, curr_lvl, &(array_kmers[i*nb_bytes]), root->k, kmer,
                        id_genome, size_id_genome, 0);
    }
}

void insertKmer_Node(Node*  node, BFT_Root*  root, int lvl_node, uint8_t*  suffix, int size_suffix,
                     uint8_t* kmer, uint32_t id_genome, int size_id_genome, int pos_CC_start_search){

    ASSERT_NULL_PTR(node,"insertKmer_Node()")
    ASSERT_NULL_PTR(suffix,"insertKmer_Node()")

    if (id_genome > NB_MAX_ID_GENOMES) ERROR ("insertKmer_Node(): unauthorized genome ID\n");

    if (node->UC_array.suffixes == NULL) initiateUC(&(node->UC_array));

    int node_nb_elem = -1;

    UC* uc;

    resultPresence* res = root->res;
    annotation_inform* ann_inf = root->ann_inf;
    annotation_array_elem* annot_sorted = root->comp_set_colors;

    //We search for the first prefix of the suffix contained in the parameter kmer
    //We want to start the search at the beginning of the node so 4th argument is 0
    presenceKmer(node, root, suffix, size_suffix, pos_CC_start_search, 1, res);

    //If the prefix is present, we insert the suffix into the corresponding child
    if (res->link_child != NULL){

        int nb_cell = root->info_per_lvl[lvl_node].size_kmer_in_bytes;
        int nb_elem;

        if (size_suffix == NB_CHAR_SUF_PREF){ //If the node is situated on the last level of the tree

            int size_substring = 0;

            if (res->container_is_UC == 0){ //If the container is not a UC

                uc = &(((UC*)((CC*)res->container)->children)[res->bucket]); //Define the array where the suffix is into CC->children

                //Determine the number of annotations in the array where the suffix is in CC->children
                if (res->bucket != CEIL(((CC*)res->container)->nb_elem, root->info_per_lvl[lvl_node].nb_ucs_skp)-1){
                    nb_elem = root->info_per_lvl[lvl_node].nb_ucs_skp;
                }
                else nb_elem = ((CC*)res->container)->nb_elem - res->bucket * root->info_per_lvl[lvl_node].nb_ucs_skp;
            }
            else { //if the container is a UC
                uc = (UC*)res->container;
                nb_elem = uc->nb_children >> 1;
                size_substring = nb_cell;
            }

            modify_annotations(root, uc, size_substring, nb_elem, res->pos_sub_bucket, id_genome,
                               size_id_genome, kmer, 1);
        }
        else { //If the node is not situated on the last level of the tree
            //Shift the suffix to delete the prefix of length NB_CHAR_SUF_PREF (this one is already resent in a container)
            int j=0;

            int nb_cell_to_delete = 2;
            if (size_suffix == 45) nb_cell_to_delete++;

            for (j=0; j < nb_cell - nb_cell_to_delete; j++){
                suffix[j] = suffix[j+2] >> 2;
                if (j+3 < nb_cell) suffix[j] |= suffix[j+3] << 6;
            }

            suffix[j-1] &= root->info_per_lvl[lvl_node].mask_shift_kmer;

            if (res->children_type_leaf == 0){ //If the container of result_presence is not a leaf (basically, a CC->children array)
                // If it is not a UC also, we continue insertion of the suffix into the child container, which can only be a node
                if (res->container_is_UC == 0){
                    insertKmer_Node((Node*)res->link_child, root, lvl_node-1, suffix, size_suffix-NB_CHAR_SUF_PREF, kmer,
                                    id_genome, size_id_genome, 0);
                }
                else{ // If it is a UC, we modify the annotation like usual (computation new size, realloc old annotation, modification of it)
                    uc = (UC*)res->container;
                    nb_elem = uc->nb_children >> 1;

                    modify_annotations(root, uc, nb_cell, nb_elem, res->pos_sub_bucket, id_genome,
                                       size_id_genome, kmer, 1);
                }
            }
            else{ //If the container of result_resence is a leaf (basically, a CC->children array), we modify it

                Node* new_node = insertKmer_Node_special(root, lvl_node, suffix, size_suffix-NB_CHAR_SUF_PREF, kmer,
                                                        id_genome, size_id_genome);

                if (new_node != NULL){
                    insertKmer_Node(new_node, root, lvl_node-1, suffix, size_suffix-NB_CHAR_SUF_PREF, kmer,
                                    id_genome, size_id_genome, 0);
                }
            }
        }
    }
    //Else if the prefix is present only as a false positive into a CC
    else if (res->presBF == 1){

        //We insert the prefix into the CC
        insertSP_CC(res, root, lvl_node, suffix, size_suffix, id_genome, size_id_genome);

        //If the CC, after insertion, stores NB_SUBSTRINGS_TRANSFORM prefixes, it is ready for recomputation
        if (((CC*)res->container)->nb_elem == root->info_per_lvl[lvl_node].tresh_suf_pref)
            transform_Filter2n3((CC*)res->container, 14, 4, &(root->info_per_lvl[lvl_node]));
    }
    //Else the prefix is not resent at all
    else{

        if (node->CC_array != NULL){

            node_nb_elem = 0;
            while ((((CC*)node->CC_array)[node_nb_elem].type & 0x1) == 0) node_nb_elem++;
            node_nb_elem++;

            if (((CC*)node->CC_array)[node_nb_elem-1].nb_elem < root->info_per_lvl[lvl_node].nb_kmers_uc){

                CC* cc_tmp = &(((CC*)node->CC_array)[node_nb_elem-1]);

                uint16_t hash1_v;
                uint16_t hash2_v;

                uint32_t substring_prefix;

                if (root->compressed > 0){
                    substring_prefix = (reverse_word_8(suffix[0]) << 10)
                                        | (reverse_word_8(suffix[1]) << 2)
                                        | (reverse_word_8(suffix[2]) >> 6);
                }
                else{
                    substring_prefix = ((reverse_word_8(suffix[0]) & 0x3f) << SIZE_BITS_UINT_8T) | reverse_word_8(suffix[1]);
                }

                hash1_v = root->hash_v[substring_prefix * 2] % root->info_per_lvl[lvl_node].modulo_hash;
                hash2_v = root->hash_v[substring_prefix * 2 + 1] % root->info_per_lvl[lvl_node].modulo_hash;

                cc_tmp->BF_filter2[hash1_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash1_v%SIZE_BITS_UINT_8T];
                cc_tmp->BF_filter2[hash2_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash2_v%SIZE_BITS_UINT_8T];

                res->pos_sub_bucket = 0;

                presenceKmer(node, root, suffix, size_suffix, node_nb_elem-1, 1, res);

                //We insert the prefix into the CC
                insertSP_CC(res, root, lvl_node, suffix, size_suffix, id_genome, size_id_genome);

                //If the CC, after insertion, stores NB_SUBSTRINGS_TRANSFORM prefixes, it is ready for recomputation
                if (((CC*)res->container)->nb_elem == root->info_per_lvl[lvl_node].tresh_suf_pref)
                    transform_Filter2n3((CC*)res->container, 14, 4, &(root->info_per_lvl[lvl_node]));
            }
            else{
                insertKmer_UC(&(node->UC_array), suffix, id_genome, size_id_genome,
                              size_suffix, res->pos_sub_bucket, ann_inf, annot_sorted);
            }
        }
        else {
            //We just insert the suffix in the UC
            //if (node->UC_array.suffixes == NULL) initiateUC(&(node->UC_array));

            insertKmer_UC(&(node->UC_array), suffix, id_genome, size_id_genome,
                          size_suffix, res->pos_sub_bucket, ann_inf, annot_sorted);
        }
    }

    //If the last container is UC_ptr and is full, it is transformed into a CC
    if ((node->UC_array.suffixes != NULL)
        && ((node->UC_array.nb_children >> 1) == root->info_per_lvl[lvl_node].nb_kmers_uc)){

        if (node->CC_array != NULL){

            if (node_nb_elem == -1){
                node_nb_elem = 0;
                while (IS_EVEN(((CC*)node->CC_array)[node_nb_elem].type)) node_nb_elem++;
                node_nb_elem++;
            }
        }
        else node_nb_elem = 0;

        node->CC_array = realloc(node->CC_array, (node_nb_elem+1)*sizeof(CC));
        ASSERT_NULL_PTR(node->CC_array,"insertKmer_Node()")

        CC* cc = &(((CC*)node->CC_array)[node_nb_elem]);

        initiateCC(cc, root->info_per_lvl[lvl_node].modulo_hash);

        if (node_nb_elem-1 >= 0) (&(((CC*)node->CC_array)[node_nb_elem-1]))->type &= 0xfffe;

        transform2CC(&(node->UC_array), cc, root, lvl_node, size_suffix);

        free(node->UC_array.suffixes);
        resetUC(&(node->UC_array));
    }

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  insertKmer_Node_special(pres, kmer, size_kmer, id_genome, unc_on_types, ann_inf)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert a suffix or modify a suffix annotation into a container (not a Node) which is a leaf in the tree.
*  ---------------------------------------------------------------------------------------------------------------
*  pres: ptr on resultPresence structure, contain information about the prefix presence/position in the container
*  kmer: ptr on the suffix to insert
*  size_kmer: size of the suffix to insert, in char.
*  id_genome: genome identity of the suffix to insert
*  info_per_lvl: ptr on info_per_level structure, contains information to manipulate CCs field CC->children_type
*  ann_inf: ptr on annotation_inform structure, used to make the transition between reading and modifying an annot
*  ---------------------------------------------------------------------------------------------------------------
*/
Node* insertKmer_Node_special(BFT_Root*  root, int lvl_cont, uint8_t*  suffix, int size_suffix,
                              uint8_t*  kmer, uint32_t id_genome, int size_id_genome){

    ASSERT_NULL_PTR(root->res->container,"insertKmer_Node_special()")
    ASSERT_NULL_PTR(root->info_per_lvl,"insertKmer_Node_special()")
    ASSERT_NULL_PTR(suffix,"insertKmer_Node_special()")

    resultPresence* pres = root->res;
    annotation_inform* ann_inf = root->ann_inf;
    annotation_array_elem* annot_sorted = root->comp_set_colors;

    CC* cc = (CC*)pres->container;
    UC* uc;

    uint8_t type = (cc->type >> 6) & 0x1;

    uint16_t nb_elt = getNbElts(cc, pres->posFilter3, type);

    int curr_lvl = lvl_cont - 1;
    int nb_cell = CEIL(size_suffix*2, SIZE_BITS_UINT_8T);
    int size_annot = ((UC*)cc->children)[pres->bucket].size_annot;
    int size_line = nb_cell+size_annot;

    int z=0;
    int pos = 0;
    int inside = 1;

    if (size_suffix == 36){
        pos = binary_search_UC_array((uint8_t*)pres->link_child, size_annot, 0, nb_elt - 1, suffix, nb_cell, 0xff);
        inside = memcmp(&(((uint8_t*)pres->link_child)[pos * size_line]), suffix, nb_cell * sizeof(uint8_t));
        if (inside < 0) pos++;
    }
    else{
        pos = binary_search_UC_array((uint8_t*)pres->link_child, size_annot, 0, nb_elt - 1, suffix, nb_cell, 0x7f);
        inside = memcmp(&(((uint8_t*)pres->link_child)[pos * size_line]), suffix, (nb_cell-1) * sizeof(uint8_t));

        if (inside < 0) pos++;
        else if (inside == 0){
            if ((((uint8_t*)pres->link_child)[pos * size_line + nb_cell - 1] & 0x7f) != suffix[nb_cell - 1]){
                inside = 1;
                if ((((uint8_t*)pres->link_child)[pos * size_line + nb_cell - 1] & 0x7f) < suffix[nb_cell - 1]) pos++;
            }
        }
    }

    if (inside != 0){
    //if (inside == 0){ //If there is no match in the suffixes comparison

        int pos_skp = pres->posFilter3/root->info_per_lvl[lvl_cont].nb_ucs_skp;

        if (nb_elt == root->info_per_lvl[curr_lvl].nb_kmers_uc){ //If the container (part of an array of CC->children) is at its max. capacity

            //Count the number of nodes after the position of insertion
            int count_Node_after_posFilter3;

            if (cc->nb_elem - pres->posFilter3+1 < pres->posFilter3+1){
                count_Node_after_posFilter3 = count_nodes(cc, pres->posFilter3+1, cc->nb_elem, type);
            }
            else count_Node_after_posFilter3 = cc->nb_Node_children - count_nodes(cc, 0, pres->posFilter3+1, type);

            //Allocate memory for a new node in CC->children_Node_container
            if (cc->nb_Node_children == 0) cc->children_Node_container = malloc(sizeof(Node));
            else cc->children_Node_container = realloc(cc->children_Node_container, (cc->nb_Node_children+1)*sizeof(Node));

            ASSERT_NULL_PTR(cc->children_Node_container,"insertKmer_Node_special()")

            int pos_Node = cc->nb_Node_children-count_Node_after_posFilter3;

            //Shift the nodes in CC->children_Node_container and initialize the new created node
            memmove(&(cc->children_Node_container[pos_Node+1]), &(cc->children_Node_container[pos_Node]), (count_Node_after_posFilter3)*sizeof(Node));
            Node* node = (Node*)&(cc->children_Node_container[pos_Node]);
            initiateNode(node);

            if (root->info_per_lvl[lvl_cont].level_min == 0){
                if ((((uint8_t*)(pres->link_child))[nb_cell-1] >> 7) == 1) node->UC_array.nb_children = 1;
                else node->UC_array.nb_children = 0;
            }

            //Create and initialize a new CC in this node
            node->CC_array = createCC(root->info_per_lvl[curr_lvl].modulo_hash);

            uc = &(((UC*)cc->children)[pos_skp]);

            //Transform the (suffixes+annotations) container into a CC (into the CC previously created)
            uint8_t** annot_extend = get_extend_annots(uc, nb_cell, uc->nb_children, pres->pos_sub_bucket, pres->pos_sub_bucket + nb_elt - 1);
            uint8_t** annot_complex = get_annots_cplx_nodes(uc, nb_cell, uc->nb_children, pres->pos_sub_bucket, pres->pos_sub_bucket + nb_elt - 1);

            transform2CC_from_arraySuffix((uint8_t*)pres->link_child, node->CC_array, root, curr_lvl, size_suffix, size_annot,
                                          annot_extend, annot_complex, uc->size_annot_cplx_nodes);

            if ((uc->nb_children - root->info_per_lvl[curr_lvl].nb_kmers_uc) > 0){

                delete_extend_annots(uc, nb_cell, uc->nb_children, pres->pos_sub_bucket, pres->pos_sub_bucket+nb_elt-1, 0, 1, 0);
                delete_annot_cplx_nodes(uc, nb_cell, uc->nb_children, pres->pos_sub_bucket, pres->pos_sub_bucket+nb_elt-1, 1, 1, 1);

                uc->nb_children -= root->info_per_lvl[curr_lvl].nb_kmers_uc;
            }
            else{
                free(uc->suffixes);
                initializeUC(uc);
                uc->size_annot = 1;
            }

            resetChildrenType(cc, pres->posFilter3, type);
            cc->nb_Node_children += 1;

            if (annot_extend != NULL) free(annot_extend);
            if (annot_complex != NULL) free(annot_complex);

            return node;
        }
        else{

            //If the container (part of an array of CC->children) is not at its max. capacity
            //We insert it into the container (part of an array of CC->children)
            uc = &(((UC*)cc->children)[pos_skp]);

            //Realloc the CC->children array where the suffix+annotation will be inserted
            //int pos_sub = pres->pos_sub_bucket+nb_elt;
            int pos_sub = pres->pos_sub_bucket + pos;
            int count_substrings_after_posFilter3 = uc->nb_children - pos_sub;
            int pos_sub_size_line;
            uint8_t* extend_annot = NULL;

            compute_best_mode(ann_inf, annot_sorted, NULL, 0, NULL, 0, id_genome, size_id_genome);

            pos_sub_size_line = pos_sub * size_line;

            if (ann_inf->min_size > size_annot){
                extend_annot = realloc_annotation(uc, nb_cell, uc->nb_children, ann_inf->min_size, 1, pos_sub);

                size_line = nb_cell + uc->size_annot;
                pos_sub_size_line = pos_sub*size_line;
            }
            else{
                z = uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                    + uc->nb_cplx_nodes * (uc->size_annot_cplx_nodes + SIZE_BYTE_CPLX_N);

                uc->suffixes = realloc(uc->suffixes, ((uc->nb_children+1) * size_line + z) * sizeof(uint8_t));
                ASSERT_NULL_PTR(uc->suffixes, "insertKmer_Node_special()")

                memmove(&(uc->suffixes[pos_sub_size_line+size_line]),
                        &(uc->suffixes[pos_sub_size_line]),
                        (count_substrings_after_posFilter3 * size_line + z)*sizeof(uint8_t));

                shift_extended_annot(uc, nb_cell, uc->nb_children+1, pos_sub);
                shift_annot_cplx_nodes(uc, nb_cell, uc->nb_children+1, pos_sub);
            }

            //Shift suffixes+annotations, insert the suffix, create the annotation
            memcpy(&(uc->suffixes[pos_sub_size_line]), suffix, nb_cell*sizeof(uint8_t));
            memset(&(uc->suffixes[pos_sub_size_line+nb_cell]), 0, uc->size_annot*sizeof(uint8_t));

            if ((pos == 0) && (root->info_per_lvl[lvl_cont].level_min == 0)
                && ((uc->suffixes[pos_sub_size_line + size_line + nb_cell - 1] >> 7) == 1)){

                uc->suffixes[pos_sub_size_line + nb_cell - 1] |= 0x80;
                uc->suffixes[pos_sub_size_line + size_line + nb_cell - 1] &= 0x7f;
            }

            modify_mode_annotation(ann_inf, &(uc->suffixes[pos_sub_size_line+nb_cell]),
                                   uc->size_annot, extend_annot, 1, id_genome, size_id_genome);

            if ((extend_annot != NULL) && (extend_annot[0] == 0))
                delete_extend_annots(uc, nb_cell, uc->nb_children, pos_sub, pos_sub, 0, 0, 1);

            type = addNewElt(cc, pres->posFilter3, cc->nb_elem, type);

            uc->nb_children++;
            cc->type = (cc->type & 0xffbf) | (type << 6);

            reinit_annotation_inform(ann_inf);
        }
    }
    else{
        uc = &(((UC*)cc->children)[pres->bucket]);

        modify_annotations(root, uc, nb_cell, uc->nb_children, pres->pos_sub_bucket + pos, id_genome,
                           size_id_genome, kmer, 1);
    }

    return NULL;
}
