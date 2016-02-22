#include "./../lib/write_to_disk.h"

void write_kmers_2disk(BFT_Root* root, char* filename, bool compressed_output){

    ASSERT_NULL_PTR(root,"write_kmers_2disk()")
    ASSERT_NULL_PTR(filename,"write_kmers_2disk()")

    struct timeval tval_before, tval_after, tval_result;

    printf("\nExtraction of k-mers from the BFT to file %s\n\n", filename);

    gettimeofday(&tval_before, NULL);

    extract_kmers_to_disk(root, filename, compressed_output);

    gettimeofday(&tval_after, NULL);
    time_spent(&tval_before, &tval_after, &tval_result);
    printf("\nElapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
}

void write_BFT_Root_sparse(BFT_Root*  root, char* filename){

    ASSERT_NULL_PTR(root,"write_BFT_Root()")
    ASSERT_NULL_PTR(filename,"write_BFT_Root()")
    ASSERT_NULL_PTR(root->info_per_lvl,"write_BFT_Root()")

    FILE* file = fopen(filename, "w");
    ASSERT_NULL_PTR(file, "write_BFT_Root")

    int i = 0;

    uint16_t str_len = 0;

    fwrite(&(root->r1), sizeof(int), 1, file);
    fwrite(&(root->r2), sizeof(int), 1, file);
    fwrite(&(root->treshold_compression), sizeof(int), 1, file);
    fwrite(&(root->nb_genomes), sizeof(int), 1, file);
    fwrite(&(root->length_comp_set_colors), sizeof(int), 1, file);
    fwrite(&(root->k), sizeof(int), 1, file);
    fwrite(&(root->compressed), sizeof(uint8_t), 1, file);

    if (root->comp_set_colors != NULL){

        for (i=0; i<root->length_comp_set_colors; i++){

            fwrite(&(root->comp_set_colors[i].last_index), sizeof(int64_t), 1, file);
            fwrite(&(root->comp_set_colors[i].size_annot), sizeof(int), 1, file);

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

    for (i = 0; i < root->k / NB_CHAR_SUF_PREF; i++){
        fwrite(&(root->info_per_lvl[i].nb_bits_per_cell_skip_filter2), sizeof(int), 1, file);
        fwrite(&(root->info_per_lvl[i].nb_bits_per_cell_skip_filter3), sizeof(int), 1, file);
        fwrite(&(root->info_per_lvl[i].nb_ucs_skp), sizeof(int), 1, file);
        fwrite(&(root->info_per_lvl[i].nb_kmers_uc), sizeof(int), 1, file);
        fwrite(&(root->info_per_lvl[i].level_min), sizeof(int), 1, file);
        fwrite(&(root->info_per_lvl[i].modulo_hash), sizeof(int), 1, file);
        fwrite(&(root->info_per_lvl[i].tresh_suf_pref), sizeof(int), 1, file);
    }

    write_Node_sparse(&(root->node), root, (root->k / NB_CHAR_SUF_PREF) - 1, file, root->k);

    fclose(file);

    return;
}

void write_Node_sparse(Node*  node, BFT_Root* root, int lvl_node, FILE* file, int size_kmer){

    ASSERT_NULL_PTR(node,"write_Node()")
    ASSERT_NULL_PTR(root,"write_Node()")

    uint32_t count_cc = 0;

    write_UC_sparse(&(node->UC_array), root, file, root->info_per_lvl[(size_kmer/NB_CHAR_SUF_PREF)-1].size_kmer_in_bytes,
                    node->UC_array.nb_children >> 1);

    if ((CC*)node->CC_array != NULL){

        do {count_cc++;}
        while (IS_EVEN(((CC*)node->CC_array)[count_cc-1].type));

        fwrite(&count_cc, sizeof(uint32_t), 1, file);

        for (uint32_t i = 0; i < count_cc; i++) write_CC_sparse(&(((CC*)node->CC_array)[i]), root, lvl_node, file, size_kmer);
    }
    else fwrite(&count_cc, sizeof(uint32_t), 1, file);

    return;
}

void write_UC_sparse(UC* uc, BFT_Root* root, FILE* file, int size_substring, uint16_t nb_children){

    ASSERT_NULL_PTR(uc, "write_UC()")
    ASSERT_NULL_PTR(root, "write_UC()")

    //if (root->compressed > 0) create_annot_cplx_nodes_marked(uc, size_substring, nb_children);

    fwrite(&(uc->nb_children), sizeof(uint16_t), 1, file);
    fwrite(&(uc->size_annot), sizeof( UC_SIZE_ANNOT_T ), 1, file);
    fwrite(&(uc->size_annot_cplx_nodes), sizeof( UC_SIZE_ANNOT_CPLX_T ), 1, file);
    fwrite(&(uc->nb_extended_annot), sizeof(uint16_t), 1, file);
    fwrite(&(uc->nb_cplx_nodes), sizeof(uint16_t), 1, file);

    if (uc->suffixes != NULL){
        fwrite(uc->suffixes, sizeof(uint8_t), nb_children * (size_substring + uc->size_annot) + uc->nb_extended_annot
               * SIZE_BYTE_EXT_ANNOT + uc->nb_cplx_nodes * (uc->size_annot_cplx_nodes + SIZE_BYTE_CPLX_N), file);
    }

    return;
}

void write_CC_sparse(CC*  cc, BFT_Root* root, int lvl_cc, FILE* file, int size_kmer){

    ASSERT_NULL_PTR(cc,"write_CC()")

    UC* uc;

    int i;

    int level = (size_kmer/NB_CHAR_SUF_PREF)-1;
    int nb_skp = CEIL(cc->nb_elem, root->info_per_lvl[lvl_cc].nb_ucs_skp);

    uint16_t size_bf = cc->type >> 7;

    uint8_t s = (cc->type >> 1) & 0x1f;
    uint8_t p = NB_CHAR_SUF_PREF*2-s;

    fwrite(&(cc->type), sizeof(uint16_t), 1, file);
    fwrite(&(cc->nb_elem), sizeof(uint16_t), 1, file);
    fwrite(&(cc->nb_Node_children), sizeof(uint16_t), 1, file);

    fwrite(&(cc->BF_filter2[size_bf]), sizeof(uint8_t), MASK_POWER_16[p]/SIZE_BITS_UINT_8T, file);

    if (s == 8) fwrite(cc->filter3, sizeof(uint8_t), cc->nb_elem, file);
    else if (IS_ODD(cc->nb_elem)) fwrite(cc->filter3, sizeof(uint8_t), (cc->nb_elem/2)+1, file);
    else fwrite(cc->filter3, sizeof(uint8_t), cc->nb_elem/2, file);

    if (root->info_per_lvl[level].level_min == 1)
        fwrite(cc->extra_filter3, sizeof(uint8_t), CEIL(cc->nb_elem,SIZE_BITS_UINT_8T), file);

    if (size_kmer != NB_CHAR_SUF_PREF){

        if ((cc->type >> 6) & 0x1) fwrite(cc->children_type, sizeof(uint8_t), cc->nb_elem, file);
        else fwrite(cc->children_type, sizeof(uint8_t), CEIL(cc->nb_elem,2), file);

        for (i = 0; i < nb_skp; i++){
            uc = &(((UC*)cc->children)[i]);
            write_UC_sparse(uc, root, file, root->info_per_lvl[level].size_kmer_in_bytes_minus_1, uc->nb_children);
        }
    }
    else{
        for (i = 0; i < nb_skp; i++){
            if (i != nb_skp-1) write_UC_sparse(&(((UC*)cc->children)[i]), root, file, 0, root->info_per_lvl[lvl_cc].nb_ucs_skp);
            else write_UC_sparse(&(((UC*)cc->children)[i]), root, file, 0, cc->nb_elem - i * root->info_per_lvl[lvl_cc].nb_ucs_skp);
        }
    }

    for (i = 0; i < cc->nb_Node_children; i++)
        write_Node_sparse(&(cc->children_Node_container[i]), root, lvl_cc-1, file, size_kmer-NB_CHAR_SUF_PREF);

    return;
}

BFT_Root* read_BFT_Root_sparse(char* filename){

    ASSERT_NULL_PTR(filename,"read_BFT_Root()")

    int i = 0;

    uint16_t str_len = 0;

    int64_t tmp = 0;

    FILE* file = fopen(filename, "r");
    ASSERT_NULL_PTR(file, "read_BFT_Root() 1")

    BFT_Root* root = createBFT_Root(0, 0, 0);

    if (fread(&(root->r1), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
    if (fread(&(root->r2), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
    if (fread(&(root->treshold_compression), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
    if (fread(&(root->nb_genomes), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
    if (fread(&(root->length_comp_set_colors), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
    if (fread(&(root->k), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
    if (fread(&(root->compressed), sizeof(uint8_t), 1, file) != 1) ERROR("read_BFT_Root()")

    root->hash_v = create_hash_v_array(root->r1, root->r2);
    root->info_per_lvl = create_info_per_level(root->k);
    root->ann_inf = create_annotation_inform(root->nb_genomes);
    root->res = create_resultPresence();

    if (root->length_comp_set_colors != 0){

        root->comp_set_colors = malloc(root->length_comp_set_colors * sizeof(annotation_array_elem));
        ASSERT_NULL_PTR(root->comp_set_colors, "read_BFT_Root() 2")

        for (i=0; i<root->length_comp_set_colors; i++){

            if (fread(&(root->comp_set_colors[i].last_index), sizeof(int64_t), 1, file) != 1) ERROR("read_BFT_Root()")
            if (fread(&(root->comp_set_colors[i].size_annot), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")

            if (i == 0) tmp = (root->comp_set_colors[i].last_index + 1) * root->comp_set_colors[i].size_annot;
            else tmp = (root->comp_set_colors[i].last_index - root->comp_set_colors[i-1].last_index) * root->comp_set_colors[i].size_annot;

            if (tmp != 0){

                root->comp_set_colors[i].annot_array = malloc(tmp * sizeof(uint8_t));
                ASSERT_NULL_PTR(root->comp_set_colors[i].annot_array, "read_BFT_Root() 3")

                if (fread(root->comp_set_colors[i].annot_array, sizeof(uint8_t), tmp, file) != tmp)
                    ERROR("read_BFT_Root()")
            }
            else root->comp_set_colors[i].annot_array = NULL;
        }
    }
    else root->comp_set_colors = NULL;

    root->filenames = malloc(root->nb_genomes * sizeof(char*));
    ASSERT_NULL_PTR(root->filenames, "read_BFT_Root() 4")

    for (i=0; i<root->nb_genomes; i++){

        if (fread(&str_len, sizeof(uint16_t), 1, file) != 1) ERROR("read_BFT_Root()")

        root->filenames[i] = malloc(str_len * sizeof(char));
        ASSERT_NULL_PTR(root->filenames[i], "read_BFT_Root()")
        if (fread(root->filenames[i], sizeof(char), str_len, file) != str_len) ERROR("read_BFT_Root()")
    }

    for (i = 0; i < root->k / NB_CHAR_SUF_PREF; i++){
        if (fread(&(root->info_per_lvl[i].nb_bits_per_cell_skip_filter2), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
        if (fread(&(root->info_per_lvl[i].nb_bits_per_cell_skip_filter3), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
        if (fread(&(root->info_per_lvl[i].nb_ucs_skp), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
        if (fread(&(root->info_per_lvl[i].nb_kmers_uc), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
        if (fread(&(root->info_per_lvl[i].level_min), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
        if (fread(&(root->info_per_lvl[i].modulo_hash), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")
        if (fread(&(root->info_per_lvl[i].tresh_suf_pref), sizeof(int), 1, file) != 1) ERROR("read_BFT_Root()")

        root->info_per_lvl[i].nb_bytes_per_cell_skip_filter2 = CEIL(root->info_per_lvl[i].nb_bits_per_cell_skip_filter2, SIZE_BITS_UINT_8T);
        root->info_per_lvl[i].nb_bytes_per_cell_skip_filter3 = CEIL(root->info_per_lvl[i].nb_bits_per_cell_skip_filter3, SIZE_BITS_UINT_8T);
    }

    read_Node_sparse(&(root->node), root, (root->k / NB_CHAR_SUF_PREF) - 1, file, root->k);

    fclose(file);

    return root;
}

void read_Node_sparse(Node*  node, BFT_Root* root, int lvl_node, FILE* file, int size_kmer){

    ASSERT_NULL_PTR(node,"read_Node()\n")
    ASSERT_NULL_PTR(root,"read_Node()\n")

    uint32_t count_ccs = 0;

    if (fread(&(node->UC_array.nb_children), sizeof(uint16_t), 1, file) != 1) ERROR("read_Node_sparse()\n")

    read_UC_sparse(&(node->UC_array), file, root->info_per_lvl[(size_kmer/NB_CHAR_SUF_PREF)-1].size_kmer_in_bytes,
                   node->UC_array.nb_children >> 1);

    if (fread(&count_ccs, sizeof(uint32_t), 1, file) != 1) ERROR("read_Node_sparse()\n")

    if (count_ccs){

        node->CC_array = malloc(count_ccs * sizeof(CC));
        ASSERT_NULL_PTR(node->CC_array,"read_Node_sparse()\n")

        for (uint32_t i = 0; i < count_ccs; i++) read_CC_sparse(&(((CC*)node->CC_array)[i]), root, lvl_node, file, size_kmer);
    }

    return;
}

void read_UC_sparse(UC* uc, FILE* file, int size_substring, uint16_t nb_children){

    ASSERT_NULL_PTR(uc, "read_UC()")

    size_t size_uc_suffixes;

    if (fread(&(uc->size_annot), sizeof( UC_SIZE_ANNOT_T ), 1, file) != 1) ERROR("read_UC()")
    if (fread(&(uc->size_annot_cplx_nodes), sizeof( UC_SIZE_ANNOT_CPLX_T ), 1, file) != 1) ERROR("read_UC()")
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

void read_CC_sparse(CC*  cc, BFT_Root* root, int lvl_cc, FILE* file, int size_kmer){

    ASSERT_NULL_PTR(cc,"read_CC() 0")

    UC* uc;

    int i, j, nb_skp, first_bit, skipFilter2, skipFilter3;

    int size_line = 0;
    int nb_elt = 0;
    int it_filter3 = 0;

    uint32_t current_sp, current_sp_tmp;

    uint16_t size_bf, it_filter2, nb_children, hash1_v, hash2_v;

    uint8_t s, p;
    uint8_t count_1 = 0;
    uint8_t type;

    size_t tmp, bf_filter2;

    if (fread(&(cc->type), sizeof(uint16_t), 1, file) != 1) ERROR("read_CC() 1")
    if (fread(&(cc->nb_elem), sizeof(uint16_t), 1, file) != 1) ERROR("read_CC() 2")
    if (fread(&(cc->nb_Node_children), sizeof(uint16_t), 1, file) != 1) ERROR("read_CC() 3")

    nb_skp = CEIL(cc->nb_elem, root->info_per_lvl[lvl_cc].nb_ucs_skp);
    type = (cc->type >> 6) & 0x1;
    size_bf = cc->type >> 7;
    s = (cc->type >> 1) & 0x1f;
    p = NB_CHAR_SUF_PREF*2-s;

    bf_filter2 = size_bf + (MASK_POWER_16[p]/SIZE_BITS_UINT_8T); // BF + Filter2
    skipFilter3 = cc->nb_elem/root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter3;
    if (cc->nb_elem >= root->info_per_lvl[lvl_cc].tresh_suf_pref) skipFilter2 = MASK_POWER_16[p]/root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter2; //SkipFilter2
    else skipFilter2 = 0;

    cc->BF_filter2 = calloc(bf_filter2 + skipFilter2 + skipFilter3, sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->BF_filter2,"read_CC() 4")
    if (fread(&(cc->BF_filter2[size_bf]), sizeof(uint8_t), bf_filter2-size_bf, file) != bf_filter2 - size_bf)
        ERROR("read_CC() 5")

    if (cc->nb_elem >= root->info_per_lvl[lvl_cc].tresh_suf_pref){
        for (i = size_bf, j = bf_filter2; i < bf_filter2; i += root->info_per_lvl[lvl_cc].nb_bytes_per_cell_skip_filter2, j++)
            cc->BF_filter2[j] = popcnt_8_par(cc->BF_filter2, i, i + root->info_per_lvl[lvl_cc].nb_bytes_per_cell_skip_filter2);
    }

    if (s == 8) tmp = cc->nb_elem;
    else if (IS_ODD(cc->nb_elem)) tmp = (cc->nb_elem/2)+1;
    else tmp = cc->nb_elem/2;

    cc->filter3 = malloc(tmp * sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->filter3,"read_CC()6")
    if (fread(cc->filter3, sizeof(uint8_t), tmp, file) != tmp) ERROR("read_CC() 7")

    if (root->info_per_lvl[lvl_cc].level_min == 1){
        tmp = CEIL(cc->nb_elem,SIZE_BITS_UINT_8T);
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

    if (size_kmer != NB_CHAR_SUF_PREF){

        if (type){
            cc->children_type = malloc(cc->nb_elem * sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->children_type,"read_CC() 13")
            if (fread(cc->children_type, sizeof(uint8_t), (size_t)cc->nb_elem, file) != (size_t)cc->nb_elem)
                ERROR("read_CC() 14");
        }
        else{
            cc->children_type = malloc(CEIL(cc->nb_elem,2) * sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->children_type,"read_CC() 15")
            if (fread(cc->children_type, sizeof(uint8_t), (size_t)CEIL(cc->nb_elem,2), file) != (size_t)CEIL(cc->nb_elem,2))
                ERROR("read_CC() 16")
        }

        for (i = 0; i < nb_skp; i++){

            if (fread(&nb_children, sizeof(uint16_t), 1, file) != 1) ERROR("read_CC() 17")

            uc = &(((UC*)cc->children)[i]);
            read_UC_sparse(uc, file, root->info_per_lvl[lvl_cc].size_kmer_in_bytes_minus_1, nb_children);
            uc->nb_children = nb_children;
        }
    }
    else{
        cc->children_type = NULL;

        for (i = 0; i < nb_skp; i++){

            if (fread(&nb_children, sizeof(uint16_t), 1, file) != 1) ERROR("read_CC() 18")

            if (i != nb_skp-1) read_UC_sparse(&(((UC*)cc->children)[i]), file, 0, root->info_per_lvl[lvl_cc].nb_ucs_skp);
            else read_UC_sparse(&(((UC*)cc->children)[i]), file, 0, cc->nb_elem - i * root->info_per_lvl[lvl_cc].nb_ucs_skp);

            ((UC*)cc->children)[i].nb_children = nb_children;
        }
    }

    for (i = 0; i < cc->nb_Node_children; i++){
        initiateNode(&(cc->children_Node_container[i]));
        read_Node_sparse(&(cc->children_Node_container[i]), root, lvl_cc-1, file, size_kmer-NB_CHAR_SUF_PREF);
    }

    j = bf_filter2 + skipFilter2;

    if (root->info_per_lvl[lvl_cc].level_min == 1){

        for (i = 0; i < skipFilter3 * root->info_per_lvl[lvl_cc].nb_bytes_per_cell_skip_filter3; i += root->info_per_lvl[lvl_cc].nb_bytes_per_cell_skip_filter3, j++)
            cc->BF_filter2[j] = popcnt_8_par(cc->extra_filter3, i, i + root->info_per_lvl[lvl_cc].nb_bytes_per_cell_skip_filter3);

        for (it_filter2=0; it_filter2<MASK_POWER_16[p]; it_filter2++){

            if ((cc->BF_filter2[size_bf+it_filter2/SIZE_BITS_UINT_8T] & (MASK_POWER_8[it_filter2%SIZE_BITS_UINT_8T])) != 0){

                first_bit = 1;
                current_sp = ((uint32_t)it_filter2) << s;

                while((it_filter3 < cc->nb_elem) &&
                      (((cc->extra_filter3[it_filter3/SIZE_BITS_UINT_8T] & MASK_POWER_8[it_filter3%SIZE_BITS_UINT_8T]) == 0) || (first_bit == 1))){

                    if (p == 10) current_sp_tmp = current_sp | cc->filter3[it_filter3];
                    else if (IS_ODD(it_filter3)) current_sp_tmp = current_sp | (cc->filter3[it_filter3/2] >> 4);
                    else current_sp_tmp = current_sp | (cc->filter3[it_filter3/2] & 0xf);

                    if (root->compressed <= 0) current_sp_tmp >>= 4;

                    hash1_v = root->hash_v[current_sp_tmp * 2] % root->info_per_lvl[lvl_cc].modulo_hash;
                    hash2_v = root->hash_v[current_sp_tmp * 2 + 1] % root->info_per_lvl[lvl_cc].modulo_hash;

                    cc->BF_filter2[hash1_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash1_v%SIZE_BITS_UINT_8T];
                    cc->BF_filter2[hash2_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash2_v%SIZE_BITS_UINT_8T];

                    it_filter3++;
                    first_bit=0;
                }
            }
        }
    }
    else{

        int lim_skipFilter3 = skipFilter3 * root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter3;

        int it = 0;
        int end = 0;
        int cpt_pv = 0;
        int it_node = 0;
        int it_children_bucket = 0;
        int nb_skp = CEIL(cc->nb_elem, root->info_per_lvl[lvl_cc].nb_ucs_skp);
        int nb_cell_children = root->info_per_lvl[lvl_cc].size_kmer_in_bytes_minus_1-1;

        for (it_filter2=0; it_filter2 < MASK_POWER_16[p]; it_filter2++){

            if ((cc->BF_filter2[size_bf+it_filter2/SIZE_BITS_UINT_8T] & (MASK_POWER_8[it_filter2%SIZE_BITS_UINT_8T])) != 0){

                first_bit = 1;

                current_sp = ((uint32_t)it_filter2) << s;

                while (it_children_bucket < nb_skp){

                    if (it_children_bucket == nb_skp - 1) end = cc->nb_elem - it_children_bucket * root->info_per_lvl[lvl_cc].nb_ucs_skp;
                    else end = root->info_per_lvl[lvl_cc].nb_ucs_skp;

                    uc = &(((UC*)cc->children)[it_children_bucket]);
                    size_line = root->info_per_lvl[lvl_cc].size_kmer_in_bytes_minus_1 + uc->size_annot;

                    while (it < end){

                        if ((nb_elt = getNbElts(cc, it_filter3, type)) == 0){

                            if (((cc->children_Node_container[it_node].UC_array.nb_children & 0x1) == 0) || (first_bit == 1)){

                                if ((first_bit == 1) && (it_filter3 < lim_skipFilter3)) count_1++;
                                first_bit=0;

                                if (p == 10) current_sp_tmp = current_sp | cc->filter3[it_filter3];
                                else if (IS_ODD(it_filter3)) current_sp_tmp = current_sp | (cc->filter3[it_filter3/2] >> 4);
                                else current_sp_tmp = current_sp | (cc->filter3[it_filter3/2] & 0xf);

                                if (root->compressed <= 0) current_sp_tmp >>= 4;

                                it_node++;
                            }
                            else goto OUT_LOOP;
                        }
                       else{
                            if (((uc->suffixes[cpt_pv*size_line+nb_cell_children] >> 7) == 0)  || (first_bit == 1)){

                                if ((first_bit == 1) && (it_filter3 < lim_skipFilter3)) count_1++;
                                first_bit=0;

                                if (p == 10) current_sp_tmp = current_sp | cc->filter3[it_filter3];
                                else if (IS_ODD(it_filter3)) current_sp_tmp = current_sp | (cc->filter3[it_filter3/2] >> 4);
                                else current_sp_tmp = current_sp | (cc->filter3[it_filter3/2] & 0xf);

                                if (root->compressed <= 0) current_sp_tmp >>= 4;

                                cpt_pv += nb_elt;
                            }
                            else goto OUT_LOOP;
                        }

                        hash1_v = root->hash_v[current_sp_tmp * 2] % root->info_per_lvl[lvl_cc].modulo_hash;
                        hash2_v = root->hash_v[current_sp_tmp * 2 + 1] % root->info_per_lvl[lvl_cc].modulo_hash;

                        cc->BF_filter2[hash1_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash1_v%SIZE_BITS_UINT_8T];
                        cc->BF_filter2[hash2_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash2_v%SIZE_BITS_UINT_8T];

                        if ((it_filter3 % root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter3 == root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter3 - 1)
                            && (it_filter3 < lim_skipFilter3)){
                            cc->BF_filter2[j] = count_1;
                            count_1 = 0;
                            j++;
                        }

                        it++;
                        it_filter3++;
                    }

                    it = 0;
                    cpt_pv = 0;
                    it_children_bucket++;
                }
            }

            OUT_LOOP: continue;
        }
    }

    return;
}
