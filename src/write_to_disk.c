#include "./../lib/write_to_disk.h"

void write_Root(Root* restrict root, char* filename, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(root,"write_Root()")
    ASSERT_NULL_PTR(filename,"write_Root()")
    ASSERT_NULL_PTR(func_on_types,"write_Root()")

    FILE* file = fopen(filename, "w");
    ASSERT_NULL_PTR(file, "write_Root")

    int i = 0;

    uint16_t str_len = 0;

    fwrite(&(root->treshold_compression), sizeof(int), 1, file);
    fwrite(&(root->nb_genomes), sizeof(int), 1, file);
    fwrite(&(root->length_comp_set_colors), sizeof(int), 1, file);
    fwrite(&(root->k), sizeof(int), 1, file);

    if (root->comp_set_colors != NULL){

        for (i=0; i<root->length_comp_set_colors; i++){

            fwrite(&(root->comp_set_colors[i].last_index), sizeof(int64_t), 1, file);
            fwrite(&(root->comp_set_colors[i].size_annot), sizeof(uint8_t), 1, file);

            if (root->comp_set_colors[i].annot_array != NULL){

                if (i == 0){

                    fwrite(root->comp_set_colors[i].annot_array,
                           sizeof(uint8_t),
                           (root->comp_set_colors[i].last_index + 1) * root->comp_set_colors[i].size_annot,
                           file);
                }
                else{

                    fwrite(root->comp_set_colors[i].annot_array,
                           sizeof(uint8_t),
                           (root->comp_set_colors[i].last_index - root->comp_set_colors[i-1].last_index) * root->comp_set_colors[i].size_annot,
                           file);
                }
            }
        }
    }

    for (i=0; i<root->nb_genomes; i++){

        str_len = strlen(root->filenames[i])+1;

        fwrite(&str_len, sizeof(uint16_t), 1, file);
        fwrite(root->filenames[i], sizeof(uint8_t), str_len, file);
    }

    write_Node(&(root->node), file, root->k, func_on_types);

    fclose(file);

    return;
}

void write_Node(Node* restrict node, FILE* file, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(node,"write_Node()")
    ASSERT_NULL_PTR(func_on_types,"write_Node()")

    uint8_t flags = 0;

    int i = -1;

    if ((CC*)node->CC_array != NULL){
        do {
            i++;

            flags = CC_VERTEX | UC_PRESENT;

            fwrite(&flags, sizeof(uint8_t), 1, file);
            flags = 0;

            write_CC(&(((CC*)node->CC_array)[i]), file, size_kmer, func_on_types);
        }
        while ((((CC*)node->CC_array)[i].type & 0x1) == 0);
    }

    flags = UC_VERTEX | UC_PRESENT;
    fwrite(&flags, sizeof(uint8_t), 1, file);
    write_UC(&(node->UC_array), file, func_on_types[(size_kmer/SIZE_SEED)-1].size_kmer_in_bytes, node->UC_array.nb_children >> 1);

    return;
}

void write_UC(UC* uc, FILE* file, int size_substring, uint16_t nb_children){

    ASSERT_NULL_PTR(uc, "write_UC()")

    fwrite(&(uc->nb_children), sizeof(uint16_t), 1, file);
    fwrite(&(uc->size_annot), sizeof(uint8_t), 1, file);
    fwrite(&(uc->size_annot_cplx_nodes), sizeof(uint8_t), 1, file);
    fwrite(&(uc->nb_extended_annot), sizeof(uint16_t), 1, file);
    fwrite(&(uc->nb_cplx_nodes), sizeof(uint16_t), 1, file);

    if (uc->suffixes != NULL){
        fwrite(uc->suffixes, sizeof(uint8_t), nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot
               * SIZE_BYTE_EXT_ANNOT + uc->nb_cplx_nodes * (uc->size_annot_cplx_nodes + SIZE_BYTE_CPLX_N), file);
    }

    return;
}

void write_CC(CC* restrict cc, FILE* file, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(cc,"write_CC()")

    UC* uc;

    int i;

    int level = (size_kmer/SIZE_SEED)-1;
    int nb_skp = CEIL(cc->nb_elem, NB_CHILDREN_PER_SKP);

    uint16_t size_bf = cc->type >> 8;

    uint8_t s = (cc->type >> 2) & 0x3f;
    uint8_t p = SIZE_SEED*2-s;

    fwrite(&(cc->type), sizeof(uint16_t), 1, file);
    fwrite(&(cc->nb_elem), sizeof(uint16_t), 1, file);
    fwrite(&(cc->nb_Node_children), sizeof(uint16_t), 1, file);

    fwrite(cc->BF_filter2, sizeof(uint8_t), size_bf + (MASK_POWER_16[p]/SIZE_CELL), file);

    if (s == 8) fwrite(cc->filter3, sizeof(uint8_t), cc->nb_elem, file);
    else if (IS_ODD(cc->nb_elem)) fwrite(cc->filter3, sizeof(uint8_t), (cc->nb_elem/2)+1, file);
    else fwrite(cc->filter3, sizeof(uint8_t), cc->nb_elem/2, file);

    if (func_on_types[level].level_min == 1)
        fwrite(cc->extra_filter3, sizeof(uint8_t), CEIL(cc->nb_elem,SIZE_CELL), file);

    if (size_kmer != SIZE_SEED){

        if (((uint8_t*)cc->children_type)[0] == 8){
            fwrite(cc->children_type, sizeof(uint8_t), cc->nb_elem+1, file);
        }
        else fwrite(cc->children_type, sizeof(uint8_t), CEIL(cc->nb_elem,2)+1, file);

        for (i = 0; i < nb_skp; i++){

            uc = &(((UC*)cc->children)[i]);
            write_UC(uc, file, func_on_types[level].size_kmer_in_bytes_minus_1, uc->nb_children);
        }
    }
    else{
        for (i = 0; i < nb_skp; i++){

            if (i != nb_skp-1) write_UC(&(((UC*)cc->children)[i]), file, 0, NB_CHILDREN_PER_SKP);
            else write_UC(&(((UC*)cc->children)[i]), file, 0, cc->nb_elem - i * NB_CHILDREN_PER_SKP);
        }
    }

    for (i = 0; i < cc->nb_Node_children; i++)
        write_Node(&(cc->children_Node_container[i]), file, size_kmer-SIZE_SEED, func_on_types);

    return;
}

Root* read_Root(char* filename){

    ASSERT_NULL_PTR(filename,"read_Root()")

    int i = 0;

    int64_t tmp = 0;

    uint16_t str_len = 0;

    FILE* file = fopen(filename, "r");
    ASSERT_NULL_PTR(file, "read_Root() 1")

    Root* root = createRoot(NULL, 0, 0, 0);

    ptrs_on_func* func_on_types = NULL;

    if (fread(&(root->treshold_compression), sizeof(int), 1, file) != 1) ERROR("read_Root()")
    if (fread(&(root->nb_genomes), sizeof(int), 1, file) != 1) ERROR("read_Root()")
    if (fread(&(root->length_comp_set_colors), sizeof(int), 1, file) != 1) ERROR("read_Root()")
    if (fread(&(root->k), sizeof(int), 1, file) != 1) ERROR("read_Root()")

    func_on_types = create_ptrs_on_func(SIZE_SEED, root->k);

    if (root->length_comp_set_colors != 0){

        root->comp_set_colors = malloc(root->length_comp_set_colors * sizeof(annotation_array_elem));
        ASSERT_NULL_PTR(root->comp_set_colors, "read_Root() 2")

        for (i=0; i<root->length_comp_set_colors; i++){

            if (fread(&(root->comp_set_colors[i].last_index), sizeof(int64_t), 1, file) != 1) ERROR("read_Root()")
            if (fread(&(root->comp_set_colors[i].size_annot), sizeof(uint8_t), 1, file) != 1) ERROR("read_Root()")

            if (i == 0) tmp = (root->comp_set_colors[i].last_index + 1) * root->comp_set_colors[i].size_annot;
            else tmp = (root->comp_set_colors[i].last_index - root->comp_set_colors[i-1].last_index) * root->comp_set_colors[i].size_annot;

            if (tmp != 0){

                root->comp_set_colors[i].annot_array = malloc(tmp * sizeof(uint8_t));
                ASSERT_NULL_PTR(root->comp_set_colors[i].annot_array, "read_Root() 3")

                if (fread(root->comp_set_colors[i].annot_array, sizeof(uint8_t), tmp, file) != tmp)
                    ERROR("read_Root()")
            }
            else root->comp_set_colors[i].annot_array = NULL;
        }
    }
    else root->comp_set_colors = NULL;

    root->filenames = malloc(root->nb_genomes * sizeof(char*));
    ASSERT_NULL_PTR(root->filenames, "read_Root() 4")

    for (i=0; i<root->nb_genomes; i++){

        if (fread(&str_len, sizeof(uint16_t), 1, file) != 1) ERROR("read_Root()")

        root->filenames[i] = malloc(str_len * sizeof(char));
        ASSERT_NULL_PTR(root->filenames[i], "read_Root()")
        if (fread(root->filenames[i], sizeof(char), str_len, file) != str_len) ERROR("read_Root()")
    }

    read_Node(&(root->node), file, root->k, func_on_types);

    fclose(file);

    free(func_on_types);

    return root;
}

void read_Node(Node* restrict node, FILE* file, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(node,"read_Node()")
    ASSERT_NULL_PTR(func_on_types,"read_Node()")

    int i = 0;

    uint8_t flags = 0;

    uint16_t nb_children;

    if (fread(&flags, sizeof(uint8_t), 1, file) != 1) ERROR("read_Node()")

    if ((flags & CC_VERTEX) != 0){
        while ((flags & CC_VERTEX) != 0){

            if (i == 0) node->CC_array = malloc(sizeof(CC));
            else node->CC_array = realloc(node->CC_array, (i+1) * sizeof(CC));

            ASSERT_NULL_PTR(node->CC_array,"read_Node()")

            read_CC(&(((CC*)node->CC_array)[i]), file, size_kmer, func_on_types);

            if ((((CC*)node->CC_array)[i].type & 0x1) == 1){

                if ((flags & UC_PRESENT) != 0){

                    if (fread(&flags, sizeof(uint8_t), 1, file) != 1) ERROR("read_Node()")
                    if (fread(&nb_children, sizeof(uint16_t), 1, file) != 1) ERROR("read_Node()")

                    read_UC(&(node->UC_array), file, func_on_types[(size_kmer/SIZE_SEED)-1].size_kmer_in_bytes, nb_children >> 1);
                    node->UC_array.nb_children = nb_children;
                }

                break;
            }

            if (fread(&flags, sizeof(uint8_t), 1, file) != 1) ERROR("read_Node()")

            i++;
        }
    }
    else {

        if (fread(&nb_children, sizeof(uint16_t), 1, file) != 1) ERROR("read_Node()")

        node->CC_array = NULL;
        read_UC(&(node->UC_array), file, func_on_types[(size_kmer/SIZE_SEED)-1].size_kmer_in_bytes, nb_children >> 1);
        node->UC_array.nb_children = nb_children;
    }

    return;
}

void read_UC(UC* uc, FILE* file, int size_substring, uint16_t nb_children){

    ASSERT_NULL_PTR(uc, "read_UC()")

    size_t size_uc_suffixes;

    if (fread(&(uc->size_annot), sizeof(uint8_t), 1, file) != 1) ERROR("read_UC()")
    if (fread(&(uc->size_annot_cplx_nodes), sizeof(uint8_t), 1, file) != 1) ERROR("read_UC()")
    if (fread(&(uc->nb_extended_annot), sizeof(uint16_t), 1, file) != 1) ERROR("read_UC()")
    if (fread(&(uc->nb_cplx_nodes), sizeof(uint16_t), 1, file) != 1) ERROR("read_UC()")

    if (nb_children != 0){

        size_uc_suffixes = nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot
                            * SIZE_BYTE_EXT_ANNOT + uc->nb_cplx_nodes * (uc->size_annot_cplx_nodes + SIZE_BYTE_CPLX_N);

        uc->suffixes = malloc(size_uc_suffixes * sizeof(uint8_t));
        ASSERT_NULL_PTR(uc->suffixes, "read_UC()")
        if (fread(uc->suffixes, sizeof(uint8_t), size_uc_suffixes, file) != size_uc_suffixes) ERROR("read_UC()")
    }
    else uc->suffixes = NULL;

    return;
}

void read_CC(CC* restrict cc, FILE* file, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(cc,"read_CC() 0")

    UC* uc;

    int i, j, nb_skp, skipFilter2, skipFilter3;

    int level = (size_kmer/SIZE_SEED)-1;
    int inc = 0xf8/SIZE_CELL;
    int it_children = 0;
    int it_nodes = 0;
    int size_line = 0;
    int nb_elt = 0;
    int pos_uc = -1;

    uint16_t size_bf, nb_children;

    uint8_t s, p, tmp_uint8;
    uint8_t count_1 = 0;

    size_t tmp, bf_filter2;

    if (fread(&(cc->type), sizeof(uint16_t), 1, file) != 1) ERROR("read_CC() 1")
    if (fread(&(cc->nb_elem), sizeof(uint16_t), 1, file) != 1) ERROR("read_CC() 2")
    if (fread(&(cc->nb_Node_children), sizeof(uint16_t), 1, file) != 1) ERROR("read_CC() 3")

    nb_skp = CEIL(cc->nb_elem, NB_CHILDREN_PER_SKP);
    size_bf = cc->type >> 8;
    s = (cc->type >> 2) & 0x3f;
    p = SIZE_SEED*2-s;

    bf_filter2 = size_bf + (MASK_POWER_16[p]/SIZE_CELL); // BF + Filter2
    skipFilter3 = cc->nb_elem/0xf8;
    if (cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM) skipFilter2 = MASK_POWER_16[p]/0xf8; //SkipFilter2
    else skipFilter2 = 0;

    cc->BF_filter2 = malloc((bf_filter2 + skipFilter3 + skipFilter2) * sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->BF_filter2,"read_CC() 4")
    if (fread(cc->BF_filter2, sizeof(uint8_t), bf_filter2, file) != bf_filter2) ERROR("read_CC() 5")

    if (cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM){
        for (i = size_bf, j = bf_filter2; i < bf_filter2; i += inc, j++)
            cc->BF_filter2[j] = popcnt_8_par(cc->BF_filter2, i, i+inc);
    }

    if (s == 8) tmp = cc->nb_elem;
    else if (IS_ODD(cc->nb_elem)) tmp = (cc->nb_elem/2)+1;
    else tmp = cc->nb_elem/2;

    cc->filter3 = malloc(tmp * sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->filter3,"read_CC()6")
    if (fread(cc->filter3, sizeof(uint8_t), tmp, file) != tmp) ERROR("read_CC() 7")

    if (func_on_types[level].level_min == 1){
        tmp = CEIL(cc->nb_elem,SIZE_CELL);
        cc->extra_filter3 = malloc(tmp * sizeof(uint8_t));
        ASSERT_NULL_PTR(cc->extra_filter3,"read_CC() 8")
        if (fread(cc->extra_filter3, sizeof(uint8_t), tmp, file) != tmp) ERROR("read_CC() 9")
    }
    else cc->extra_filter3 = NULL;

    cc->children = malloc(nb_skp * sizeof(UC));
    ASSERT_NULL_PTR(cc->children,"read_CC() 10")

    if (cc->nb_Node_children != 0){
        cc->children_Node_container = malloc(cc->nb_Node_children * sizeof(Node));
        ASSERT_NULL_PTR(cc->children_Node_container,"read_CC() 11")
    }
    else cc->children_Node_container = NULL;

    if (size_kmer != SIZE_SEED){

        if (fread(&tmp_uint8, sizeof(uint8_t), 1, file) != 1) ERROR("read_CC() 12")

        if (tmp_uint8 == 8){
            cc->children_type = malloc((cc->nb_elem+1) * sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->children_type,"read_CC() 13")
            if (fread(&(((uint8_t*)cc->children_type)[1]), sizeof(uint8_t), (size_t)cc->nb_elem, file) != (size_t)cc->nb_elem)
                ERROR("read_CC() 14");
        }
        else{
            cc->children_type = malloc((CEIL(cc->nb_elem,2)+1) * sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->children_type,"read_CC() 15")
            if (fread(&(((uint8_t*)cc->children_type)[1]), sizeof(uint8_t), (size_t)CEIL(cc->nb_elem,2), file) != (size_t)CEIL(cc->nb_elem,2))
                ERROR("read_CC() 16")
        }

        ((uint8_t*)cc->children_type)[0] = tmp_uint8;

        for (i = 0; i < nb_skp; i++){

            if (fread(&nb_children, sizeof(uint16_t), 1, file) != 1) ERROR("read_CC() 17")

            uc = &(((UC*)cc->children)[i]);
            read_UC(uc, file, func_on_types[level].size_kmer_in_bytes_minus_1, nb_children);
            uc->nb_children = nb_children;
        }
    }
    else{
        cc->children_type = NULL;

        for (i = 0; i < nb_skp; i++){

            if (fread(&nb_children, sizeof(uint16_t), 1, file) != 1) ERROR("read_CC() 18")

            if (i != nb_skp-1) read_UC(&(((UC*)cc->children)[i]), file, 0, NB_CHILDREN_PER_SKP);
            else read_UC(&(((UC*)cc->children)[i]), file, 0, cc->nb_elem - i * NB_CHILDREN_PER_SKP);

            ((UC*)cc->children)[i].nb_children = nb_children;
        }
    }

    for (i = 0; i < cc->nb_Node_children; i++){
        initiateNode(&(cc->children_Node_container[i]));
        read_Node(&(cc->children_Node_container[i]), file, size_kmer-SIZE_SEED, func_on_types);
    }

    j = bf_filter2 + skipFilter2;

    if (func_on_types[level].level_min == 1){

        for (i = 0; i < skipFilter3*inc; i += inc){
            cc->BF_filter2[j] = popcnt_8_par(cc->extra_filter3, i, i+inc);
            j++;
        }
    }
    else{

        for (i=0; i<skipFilter3*0xf8; i++){

            if (i%NB_CHILDREN_PER_SKP == 0){
                pos_uc++;
                it_children = 0;
                uc = &(((UC*)cc->children)[pos_uc]);
                size_line = func_on_types[level].size_kmer_in_bytes_minus_1 + uc->size_annot;
            }

            if ((nb_elt = (*func_on_types[level].getNbElts)(cc, i)) == 0){
                count_1 += cc->children_Node_container[it_nodes].UC_array.nb_children & 0x1;
                it_nodes++;
            }
            else {
                count_1 += uc->suffixes[it_children * size_line + func_on_types[level].size_kmer_in_bytes_minus_1-1] >> 7;
                it_children += nb_elt;
            }

            if (i%0xf8 == 0xf8-1){
                cc->BF_filter2[j] = count_1;
                count_1 = 0;
                j++;
            }
        }
    }

    return;
}
