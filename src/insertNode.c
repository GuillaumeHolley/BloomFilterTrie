#include "./../lib/insertNode.h"

/* ---------------------------------------------------------------------------------------------------------------
*  insertKmer_Node(node, kmer, size_kmer, id_genome, func_on_types, ann_inf)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert a suffix or modify a suffix annotation into a node of the tree.
*  ---------------------------------------------------------------------------------------------------------------
*  node: pointer on a Node structure to insert the suffix
*  kmer: pointer on the suffix to insert
*  size_kmer: size of the suffix to insert, in char.
*  id_genome: genome identity of the suffix to insert
*  func_on_types: ptr on ptrs_on_func structure, contains information to manipulate CCs field CC->children_type
*  ann_inf: ptr on annotation_inform structure, used to make the transition between reading and modifying an annot
*  ---------------------------------------------------------------------------------------------------------------
*/
void insertKmer_Node(Node* restrict node, Node* restrict root, uint8_t* restrict suffix, int size_suffix, uint8_t* kmer, int size_kmer, int id_genome,
                     ptrs_on_func* restrict func_on_types, annotation_inform* ann_inf, resultPresence* res, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(node,"insertKmer_Node()")
    ASSERT_NULL_PTR(suffix,"insertKmer_Node()")

    if (id_genome > NB_MAX_ID_GENOMES || id_genome < 0) ERROR ("insertKmer_Node(): unauthorized genome ID\n");

    if (node->CC_array == NULL) if (node->UC_array.suffixes == NULL) initiateUC(&(node->UC_array));

    int level = (size_suffix/SIZE_SEED)-1;

    UC* uc;

    __builtin_prefetch (&(func_on_types[level]), 0, 0);

    //We search for the first prefix of the suffix contained in the parameter kmer
    presenceKmer(node, suffix, size_suffix, &(func_on_types[level]), res);

    //If the prefix is present, we insert the suffix into the corresponding child
    if (res->link_child != NULL){
        int nb_cell = func_on_types[level].size_kmer_in_bytes;
        int nb_elem;

        if (size_suffix == SIZE_SEED){ //If the node is situated on the last level of the tree

            int size_substring = 0;

            if (res->container_is_UC == 0){ //If the container is not a UC

                uc = &(((UC*)((CC*)res->container)->children)[res->bucket]); //Define the array where the suffix is into CC->children

                //Determine the number of annotations in the array where the suffix is in CC->children
                if (res->bucket != CEIL(((CC*)res->container)->nb_elem,NB_CHILDREN_PER_SKP)-1) nb_elem = NB_CHILDREN_PER_SKP;
                else nb_elem = ((CC*)res->container)->nb_elem - res->bucket * NB_CHILDREN_PER_SKP;
            }
            else { //if the container is a UC
                uc = (UC*)res->container;
                nb_elem = uc->nb_children >> 1;
                size_substring = nb_cell;
            }

            modify_annotations(root, uc, size_substring, nb_elem, res->pos_sub_bucket, id_genome, kmer, size_kmer, func_on_types, ann_inf, annot_sorted, 1);
        }
        else { //If the node is not situated on the last level of the tree
            //Shift the suffix to delete the prefix of length SIZE_SEED (this one is already resent in a container)
            int j=0;

            int nb_cell_to_delete = 2;
            if (size_suffix == 45) nb_cell_to_delete++;

            for (j=0; j < nb_cell - nb_cell_to_delete; j++){
                suffix[j] = suffix[j+2] >> 2;
                if (j+3 < nb_cell) suffix[j] |= suffix[j+3] << 6;
            }

            suffix[j-1] &= func_on_types[level].mask_shift_kmer;

            if (res->children_type_leaf == 0){ //If the container of result_presence is not a leaf (basically, a CC->children array)
                // If it is not a UC also, we continue insertion of the suffix into the child container, which can only be a node
                if (res->container_is_UC == 0){
                    insertKmer_Node((Node*)res->link_child, root, suffix, size_suffix-SIZE_SEED, kmer, size_kmer, id_genome, func_on_types, ann_inf, res, annot_sorted);
                }
                else{ // If it is a UC, we modify the annotation like usual (computation new size, realloc old annotation, modification of it)
                    uc = (UC*)res->container;
                    nb_elem = uc->nb_children >> 1;

                    modify_annotations(root, uc, nb_cell, nb_elem, res->pos_sub_bucket, id_genome, kmer, size_kmer, func_on_types, ann_inf, annot_sorted, 1);
                }
            }
            else{ //If the container of result_resence is a leaf (basically, a CC->children array), we modify it

                Node* new_node = insertKmer_Node_special(root, res, suffix, size_suffix-SIZE_SEED, kmer, size_kmer, id_genome, func_on_types, ann_inf, annot_sorted);

                if (new_node != NULL)
                    insertKmer_Node(new_node, root, suffix, size_suffix-SIZE_SEED, kmer, size_kmer, id_genome, func_on_types, ann_inf, res, annot_sorted);
            }
        }
    }
    //Else if the prefix is resent only as a false positive into a CC
    else if (res->presBF == 1){

        //We insert the prefix into the CC
        insertSP_CC(res, suffix, size_suffix, id_genome, &(func_on_types[level]), ann_inf, annot_sorted);

        //If the CC, after insertion, stores NB_SUBSTRINGS_TRANSFORM prefixes, it is ready for recomputation
        if (((CC*)res->container)->nb_elem == NB_SUBSTRINGS_TRANSFORM) transform_Filter2n3((CC*)res->container, 14, 4, &(func_on_types[level]));
    }
    //Else the prefix is not resent at all
    else{

        //We just insert the suffix in the UC
        if (node->UC_array.suffixes == NULL) initiateUC(&(node->UC_array));
        insertKmer_UC(&(node->UC_array), suffix, id_genome, size_suffix, ann_inf, annot_sorted);
    }

    //If the last container is UC_ptr and is full, it is transformed into a CC
    if ((node->UC_array.suffixes != NULL) && ((node->UC_array.nb_children >> 1) == func_on_types[level].nb_kmers_per_cc)){
        int node_nb_elem = 0;

        if (node->CC_array != NULL){
            while ((((CC*)node->CC_array)[node_nb_elem].type & 0x1) == 0) node_nb_elem++;
            node_nb_elem++;
        }

        node->CC_array = realloc(node->CC_array, (node_nb_elem+1)*sizeof(CC));
        ASSERT_NULL_PTR(node->CC_array,"insertKmer_Node()")

        CC* cc = &(((CC*)node->CC_array)[node_nb_elem]);

        initiateCC(cc, MODULO_HASH);

        if (node_nb_elem-1 >= 0) (&(((CC*)node->CC_array)[node_nb_elem-1]))->type &= 0xfffe;

        transform2CC(&(node->UC_array), cc, size_suffix, &(func_on_types[level]));

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
*  func_on_types: ptr on ptrs_on_func structure, contains information to manipulate CCs field CC->children_type
*  ann_inf: ptr on annotation_inform structure, used to make the transition between reading and modifying an annot
*  ---------------------------------------------------------------------------------------------------------------
*/
Node* insertKmer_Node_special(Node* root, resultPresence* restrict pres, uint8_t* restrict suffix, int size_suffix, uint8_t* restrict kmer, int size_kmer, int id_genome,
                              ptrs_on_func* restrict func_on_types, annotation_inform* ann_inf, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(pres,"insertKmer_Node_special()")
    ASSERT_NULL_PTR(pres->container,"insertKmer_Node_special()")
    ASSERT_NULL_PTR(suffix,"insertKmer_Node_special()")
    ASSERT_NULL_PTR(func_on_types,"insertKmer_Node_special()")
    ASSERT_NULL_PTR(ann_inf,"insertKmer_Node_special()")

    CC* cc = (CC*)pres->container;
    UC* uc;

    int level = size_suffix/SIZE_SEED;
    uint16_t nb_elt = (*func_on_types[level].getNbElts)(pres->container, pres->posFilter3);

    int nb_cell = CEIL(size_suffix*2, SIZE_CELL);
    uint8_t size_annot = ((UC*)cc->children)[pres->bucket].size_annot;
    int size_line = nb_cell+size_annot;

    int z=0;
    int inside = 0;

    //Iterate over the suffixes + annotations to compare stored suffixes with the suffix we want to insert
    if (size_suffix == 36){
        for (z=0; z<nb_elt; z++){
            if (memcmp(&(((uint8_t*)pres->link_child)[z*size_line]), suffix, nb_cell*sizeof(uint8_t)) == 0){
                inside = 1;
                break;
            }
        }
    }
    else {
        for (z=0; z<nb_elt; z++){
            if ((memcmp(&(((uint8_t*)pres->link_child)[z*size_line]), suffix, (nb_cell-1)*sizeof(uint8_t)) == 0) &&
                ((((uint8_t*)pres->link_child)[z*size_line+nb_cell-1] & 0x7f) == suffix[nb_cell-1])){
                inside = 1;
                break;
            }
        }
    }

    if (inside == 0){ //If there is no match in the suffixes comparison

        if (nb_elt == func_on_types[level].nb_kmers_per_cc){ //If the container (part of an array of CC->children) is at its max. capacity

            int pos_skp = pres->posFilter3/NB_CHILDREN_PER_SKP;

            //Count the number of nodes after the position of insertion
            int count_Node_after_posFilter3 = (*func_on_types[level].count_nodes)(pres->container, pres->posFilter3+1, cc->nb_elem);

            //Allocate memory for a new node in CC->children_Node_container
            if (cc->nb_Node_children == 0) cc->children_Node_container = malloc(sizeof(Node));
            else cc->children_Node_container = realloc(cc->children_Node_container, (cc->nb_Node_children+1)*sizeof(Node));

            ASSERT_NULL_PTR(cc->children_Node_container,"insertKmer_Node_special()")

            int pos_Node = cc->nb_Node_children-count_Node_after_posFilter3;

            //Shift the nodes in CC->children_Node_container and initialize the new created node
            memmove(&(cc->children_Node_container[pos_Node+1]), &(cc->children_Node_container[pos_Node]), (count_Node_after_posFilter3)*sizeof(Node));
            Node* node = (Node*)&(cc->children_Node_container[pos_Node]);
            initiateNode(node);

            if (func_on_types[level].level_min == 0){
                if ((((uint8_t*)(pres->link_child))[nb_cell-1] >> 7) == 1) cc->children_Node_container[pos_Node].UC_array.nb_children = 1;
                else cc->children_Node_container[pos_Node].UC_array.nb_children = 0;
            }

            //Create and initialize a new CC in this node
            ((Node*)&(cc->children_Node_container[pos_Node]))->CC_array = createCC(MODULO_HASH);

            uc = &(((UC*)cc->children)[pos_skp]);

            //Transform the (suffixes+annotations) container into a CC (into the CC previously created)
            uint8_t** annot_extend = get_extend_annots(uc, nb_cell, uc->nb_children, pres->pos_sub_bucket, pres->pos_sub_bucket + nb_elt - 1);
            uint8_t** annot_complex = get_annots_cplx_nodes(uc, nb_cell, uc->nb_children, pres->pos_sub_bucket, pres->pos_sub_bucket + nb_elt - 1);

            transform2CC_from_arraySuffix((uint8_t*)pres->link_child,
                                       ((Node*)&(cc->children_Node_container[pos_Node]))->CC_array,
                                       size_suffix,
                                       size_annot,
                                       annot_extend,
                                       annot_complex,
                                       uc->size_annot_cplx_nodes,
                                       &(func_on_types[level-1]));

            if ((uc->nb_children - func_on_types[level].nb_kmers_per_cc) > 0){

                delete_extend_annots(uc, nb_cell, uc->nb_children, pres->pos_sub_bucket, pres->pos_sub_bucket+nb_elt-1, 0, 1, 0);
                delete_annot_cplx_nodes(uc, nb_cell, uc->nb_children, pres->pos_sub_bucket, pres->pos_sub_bucket+nb_elt-1, 1, 1, 1);

                uc->nb_children -= func_on_types[level].nb_kmers_per_cc;
            }
            else{
                free(uc->suffixes);
                initializeUC(uc);
                uc->size_annot = 1;
            }

            (*func_on_types[level].resetChildrenType)(pres->container, pres->posFilter3);
            cc->nb_Node_children += 1;

            if (annot_extend != NULL) free(annot_extend);
            if (annot_complex != NULL) free(annot_complex);

            return node;
        }
        else{ //If the container (part of an array of CC->children) is not at its max. capacity
            //We insert it into the container (part of an array of CC->children)
            int pos_skp = pres->posFilter3/NB_CHILDREN_PER_SKP;
            uc = &(((UC*)cc->children)[pos_skp]);

            //Realloc the CC->children array where the suffix+annotation will be inserted
            int pos_sub = pres->pos_sub_bucket+nb_elt;
            int count_substrings_after_posFilter3 = uc->nb_children - pos_sub;
            int pos_sub_size_line;
            uint8_t* extend_annot = NULL;

            compute_best_mode(ann_inf, annot_sorted, NULL, 0, NULL, 0, id_genome);

            pos_sub_size_line = pos_sub * size_line;

            if (ann_inf->min_size > size_annot){
                extend_annot = realloc_annotation(uc, nb_cell, uc->nb_children, ann_inf->min_size, 1, pos_sub);

                size_line = nb_cell + uc->size_annot;
                pos_sub_size_line = pos_sub*size_line;
            }
            else{
                z = uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + uc->nb_cplx_nodes * (uc->size_annot_cplx_nodes + SIZE_BYTE_CPLX_N);
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
            modify_mode_annotation(ann_inf, &(uc->suffixes[pos_sub_size_line+nb_cell]), uc->size_annot, extend_annot, 1, id_genome);

            if ((extend_annot != NULL) && (extend_annot[0] == 0)) delete_extend_annots(uc, nb_cell, uc->nb_children, pos_sub, pos_sub, 0, 0, 1);

            (*func_on_types[level].addNewElt)(pres->container, pres->posFilter3, cc->nb_elem);
            uc->nb_children++;

            reinit_annotation_inform(ann_inf);
        }
    }
    else{
        uc = &(((UC*)cc->children)[pres->bucket]);

        modify_annotations(root, uc, nb_cell, uc->nb_children, pres->pos_sub_bucket+z, id_genome, kmer, size_kmer, func_on_types, ann_inf, annot_sorted, 1);
    }

    return NULL;
}
