#include "merge.h"

/*void merging_BFT(char* prefix_bft1, char* prefix_bft2, char* output_prefix, int cut_lvl, bool packed_in_subtries){

    ASSERT_NULL_PTR(prefix_bft1, "merging_BFT()\n")
    ASSERT_NULL_PTR(prefix_bft2, "merging_BFT()\n")

    UC* uc = NULL;
    CC* cc = NULL;
    Node* node = NULL;

    BFT** bfts_insert;

    BFT* bft1 = l_read_BFT_Root(prefix_bft1, cut_lvl);
    BFT* bft2 = l_read_BFT_Root(prefix_bft2, cut_lvl);

    bool overlap_genome_ids = are_genomes_ids_overlapping(bft1, bft2);

    int nb_ccs_bft1 = 0;
    int lvl_root = (bft1->k / NB_CHAR_SUF_PREF) - 1;
    int nb_bytes = CEIL(bft1->k * 2, SIZE_BITS_UINT_8T);
    int len_filename_cc = strlen(prefix_bft1);
    int len_output_filename = strlen(output_prefix);
    int old_nb_genome_ids_bft = bft1->nb_genomes;

    int i, j, nb_threads, thread_id;
    int len_filename_container;
    int len_output_filename_cont;
    int nb_skp;

    int* presence_bf = get_bf_presence_per_cc(bft1);

    BFT_annotation** bft_annot;

    FILE** files_cc;

    FILE* file_tmp;

    uint8_t** kmer_comp;
    uint8_t** kmer_comp_cpy;

    char** new_filename;

    char** output_filename = malloc(sizeof(char*));
    ASSERT_NULL_PTR(output_filename, "merging_BFT() 4\n");

    output_filename[0] = malloc((len_output_filename + 30) * sizeof(char));
    ASSERT_NULL_PTR(output_filename[0], "merging_BFT() 4\n");

    strcpy(output_filename[0], output_prefix);
    strcpy(&(output_filename[0][len_output_filename]), "_c");

    if (bft1->node.CC_array != NULL){
        do {nb_ccs_bft1++;}
        while (IS_EVEN(((CC*)bft1->node.CC_array)[nb_ccs_bft1-1].type));
    }

    files_cc = malloc((nb_ccs_bft1 + 1) * sizeof(FILE*));
    ASSERT_NULL_PTR(files_cc, "merging_BFT() 6\n");

    for (j = 0; j <= nb_ccs_bft1; j++){
        sprintf(&(output_filename[0][len_output_filename+2]), "%d", j);

        files_cc[j] = fopen(output_filename[0], "w");
        ASSERT_NULL_PTR(files_cc[j], "merging_BFT() 7\n");
    }

    free(output_filename[0]);
    free(output_filename);

    l_iterate_over_kmers(bft2, cut_lvl, packed_in_subtries, prefix_bft2, l_insert_kmer, nb_bytes, files_cc, presence_bf, nb_ccs_bft1);

    for (j = 0; j <= nb_ccs_bft1; j++) fclose(files_cc[j]);

    if (overlap_genome_ids) add_genomes_BFT_Root(bft2->nb_genomes - 1, &(bft2->filenames[1]), bft1);
    else add_genomes_BFT_Root(bft2->nb_genomes, bft2->filenames, bft1);

    if (bft1->node.CC_array != NULL){

        #pragma omp parallel                                                                                                        \
        shared(bft1, nb_ccs_bft1, new_filename, output_filename, len_filename_cc, nb_bytes, kmer_comp, kmer_comp_cpy, files_cc,     \
               bft_annot, old_nb_genome_ids_bft, overlap_genome_ids, lvl_root, nb_threads, bfts_insert)                             \
        private(i, j, uc, cc, node, nb_skp, file_tmp, len_filename_container, len_output_filename_cont, thread_id)
        {

            #pragma omp single
            {
                nb_threads = omp_get_num_threads();

                bfts_insert = malloc(nb_threads * sizeof(BFT*));
                ASSERT_NULL_PTR(bfts_insert, "merging_BFT() 3\n");

                new_filename = malloc(nb_threads * sizeof(char*));
                ASSERT_NULL_PTR(new_filename, "merging_BFT() 3\n");

                output_filename = malloc(nb_threads * sizeof(char*));
                ASSERT_NULL_PTR(output_filename, "merging_BFT() 3\n");

                kmer_comp = malloc(nb_threads * sizeof(uint8_t*));
                ASSERT_NULL_PTR(kmer_comp, "merging_BFT() 3\n");

                kmer_comp_cpy = malloc(nb_threads * sizeof(uint8_t*));
                ASSERT_NULL_PTR(kmer_comp_cpy, "merging_BFT() 3\n");

                bft_annot = malloc(nb_threads * sizeof(BFT_annotation*));
                ASSERT_NULL_PTR(bft_annot, "merging_BFT() 3\n");

                for (i = 0; i < nb_threads; i++){

                    bfts_insert[i] = copy_BFT_Root(bft1);

                    kmer_comp[i] = malloc(nb_bytes * sizeof(uint8_t));
                    ASSERT_NULL_PTR(kmer_comp[i], "merging_BFT() 1\n");

                    kmer_comp_cpy[i] = malloc(nb_bytes * sizeof(uint8_t));
                    ASSERT_NULL_PTR(kmer_comp_cpy[i], "merging_BFT() 1\n");

                    new_filename[i] = malloc((len_filename_cc + 30) * sizeof(char));
                    ASSERT_NULL_PTR(new_filename[i], "merging_BFT() 3\n");

                    strcpy(new_filename[i], prefix_bft1);
                    new_filename[i][len_filename_cc] = '_';

                    output_filename[i] = malloc((len_output_filename + 30) * sizeof(char));
                    ASSERT_NULL_PTR(output_filename[i], "merging_BFT() 3\n");

                    strcpy(output_filename[i], output_prefix);
                    output_filename[i][len_output_filename] = '_';

                    bft_annot[i] = create_BFT_annotation();

                    bft_annot[i]->annot = calloc(SIZE_MAX_BYTE_ANNOT, sizeof(uint8_t));
                    ASSERT_NULL_PTR(bft_annot[i]->annot, "merging_BFT() 5\n");
                }
            }

            #pragma omp for
            for (i = 0; i < nb_ccs_bft1; i++){

                thread_id = omp_get_thread_num();

                cc = &(((CC*)bfts_insert[thread_id]->node.CC_array)[i]);
                nb_skp = CEIL(cc->nb_elem, bfts_insert[thread_id]->info_per_lvl[lvl_root].nb_ucs_skp);

                sprintf(&(new_filename[thread_id][len_filename_cc+1]), "%d", i);
                len_filename_container = strlen(new_filename[thread_id]);

                if (!packed_in_subtries) new_filename[thread_id][len_filename_container] = '_';

                output_filename[thread_id][len_output_filename+1] = 'c';
                sprintf(&(output_filename[thread_id][len_output_filename+2]), "%d", i);

                files_cc[i] = fopen(output_filename[thread_id], "r");
                ASSERT_NULL_PTR(files_cc[i], "merging_BFT() 8\n");

                sprintf(&(output_filename[thread_id][len_output_filename+1]), "%d", i);
                len_output_filename_cont = strlen(output_filename[thread_id]);

                if (!packed_in_subtries) output_filename[thread_id][len_output_filename_cont] = '_';
                else{
                    file_tmp = fopen(new_filename[thread_id], "r");
                    ASSERT_NULL_PTR(file_tmp, "merging_BFT() 9\n");
                }

                for (j = 0; j < nb_skp; j++){

                    uc = &(((UC*)cc->children)[j]);

                    if (!packed_in_subtries){
                        sprintf(&(new_filename[thread_id][len_filename_container+1]), "%d", j);

                        file_tmp = fopen(new_filename[thread_id], "r");
                        ASSERT_NULL_PTR(file_tmp, "merging_BFT() 9.1\n");
                    }

                    if (fread(&(uc->nb_children), sizeof(uint16_t), 1, file_tmp) != 1) ERROR("merging_BFT() 10\n")

                    if (lvl_root){
                        read_UC_sparse(uc, bfts_insert[thread_id]->ann_inf, file_tmp,
                                       bfts_insert[thread_id]->info_per_lvl[lvl_root].size_kmer_in_bytes_minus_1, uc->nb_children);
                    }
                    else if (j != nb_skp-1){
                        read_UC_sparse(uc, bfts_insert[thread_id]->ann_inf, file_tmp, 0,
                                       bfts_insert[thread_id]->info_per_lvl[lvl_root].nb_ucs_skp);
                    }
                    else {
                        read_UC_sparse(uc, bfts_insert[thread_id]->ann_inf, file_tmp, 0,
                                       cc->nb_elem - j * bfts_insert[thread_id]->info_per_lvl[lvl_root].nb_ucs_skp);
                    }

                    if (!packed_in_subtries) fclose(file_tmp);
                }

                for (j = 0; j < cc->nb_Node_children; j++){

                    node = &(cc->children_Node_container[j]);

                    if (!packed_in_subtries){
                        sprintf(&(new_filename[thread_id][len_filename_container+1]), "%d", j + nb_skp);

                        file_tmp = fopen(new_filename[thread_id], "r");
                        ASSERT_NULL_PTR(file_tmp, "merging_BFT() 11\n");
                    }

                    l_read_Node(node, bfts_insert[thread_id], lvl_root-1, cut_lvl, file_tmp, new_filename[thread_id],
                                bfts_insert[thread_id]->k - NB_CHAR_SUF_PREF);

                    if (!packed_in_subtries) fclose(file_tmp);
                }

                while (fread(kmer_comp[thread_id], sizeof(uint8_t), nb_bytes, files_cc[i]) == nb_bytes){
                    if (fread(&(bft_annot[thread_id]->size_annot), sizeof(int), 1, files_cc[i]) != 1) ERROR("merging_BFT() 13\n")
                    if (fread(bft_annot[thread_id]->annot, sizeof(uint8_t), bft_annot[thread_id]->size_annot, files_cc[i])
                        != bft_annot[thread_id]->size_annot) ERROR("merging_BFT() 14\n")

                    l_insert_kmer_bis(bfts_insert[thread_id], lvl_root, kmer_comp[thread_id], kmer_comp_cpy[thread_id],
                                       nb_bytes, old_nb_genome_ids_bft, overlap_genome_ids, bft_annot[thread_id], i);

                    memset(bft_annot[thread_id]->annot, 0, bft_annot[thread_id]->size_annot * sizeof(uint8_t));
                }

                fclose(files_cc[i]);

                if (packed_in_subtries){
                    fclose(file_tmp);

                    file_tmp = fopen(output_filename[thread_id], "w");
                    ASSERT_NULL_PTR(file_tmp, "merging_BFT() 15\n");
                }

                cc = &(((CC*)bfts_insert[thread_id]->node.CC_array)[i]);
                nb_skp = CEIL(cc->nb_elem, bfts_insert[thread_id]->info_per_lvl[lvl_root].nb_ucs_skp);

                if (i < nb_ccs_bft1 - 1){

                    for (j = 0; j < nb_skp; j++){

                        uc = &(((UC*)cc->children)[j]);

                        if (!packed_in_subtries){
                            sprintf(&(output_filename[thread_id][len_output_filename_cont+1]), "%d", j);

                            file_tmp = fopen(output_filename[thread_id], "w");
                            ASSERT_NULL_PTR(file_tmp, "merging_BFT() 15\n");
                        }

                        if (lvl_root){
                            write_UC_sparse(uc, bfts_insert[thread_id], file_tmp,
                                            bfts_insert[thread_id]->info_per_lvl[lvl_root].size_kmer_in_bytes_minus_1,
                                            uc->nb_children, true);
                        }
                        else if (j != nb_skp-1){
                            write_UC_sparse(uc, bfts_insert[thread_id], file_tmp, 0,
                                            bfts_insert[thread_id]->info_per_lvl[lvl_root].nb_ucs_skp, false);
                        }
                        else{
                            write_UC_sparse(uc, bfts_insert[thread_id], file_tmp, 0,
                                            cc->nb_elem - j * bfts_insert[thread_id]->info_per_lvl[lvl_root].nb_ucs_skp, false);
                        }

                        if (!packed_in_subtries) fclose(file_tmp);
                        if (uc->suffixes != NULL) free(uc->suffixes);
                    }

                    for (j = 0; j < cc->nb_Node_children; j++){

                        node = &(cc->children_Node_container[j]);

                        if (!packed_in_subtries){
                            sprintf(&(output_filename[thread_id][len_output_filename_cont+1]), "%d", j + nb_skp);

                            file_tmp = fopen(output_filename[thread_id], "w");
                            ASSERT_NULL_PTR(file_tmp, "merging_BFT() 16\n");
                        }

                        l_write_Node(node, bfts_insert[thread_id], lvl_root-1, bfts_insert[thread_id]->k - NB_CHAR_SUF_PREF,
                                     cut_lvl, true, file_tmp, output_filename[thread_id]);

                        if (!packed_in_subtries) fclose(file_tmp);

                        freeNode(node, lvl_root-1, bfts_insert[thread_id]->info_per_lvl);
                    }

                    if (packed_in_subtries) fclose(file_tmp);

                    //printf("Merging CC %d finished\n", i);
                }
            }

            #pragma omp single
            {
                for (i = 1; i < nb_threads; i++){
                    free(kmer_comp[i]);
                    free(kmer_comp_cpy[i]);
                    free(new_filename[i]);
                    free(output_filename[i]);
                    free(bft_annot[i]->annot);
                    free(bft_annot[i]);
                }
            }
        }

        output_filename[0][len_output_filename+1] = 'c';
        sprintf(&(output_filename[0][len_output_filename+2]), "%d", nb_ccs_bft1);

        files_cc[nb_ccs_bft1] = fopen(output_filename[0], "r");
        ASSERT_NULL_PTR(files_cc[nb_ccs_bft1], "merging_BFT() 17\n");

        i = nb_ccs_bft1 - 1;

        while (fread(kmer_comp[0], sizeof(uint8_t), nb_bytes, files_cc[nb_ccs_bft1]) == nb_bytes){
            if (fread(&(bft_annot[0]->size_annot), sizeof(int), 1, files_cc[nb_ccs_bft1]) != 1) ERROR("merging_BFT() 19\n")
            if (fread(bft_annot[0]->annot, sizeof(uint8_t), bft_annot[0]->size_annot, files_cc[nb_ccs_bft1])
                != bft_annot[0]->size_annot) ERROR("merging_BFT() 20\n")

            l_insert_kmer_bis(bft1, lvl_root, kmer_comp[0], kmer_comp_cpy[0], nb_bytes, old_nb_genome_ids_bft,
                               overlap_genome_ids, bft_annot[0], i);

            memset(bft_annot[0]->annot, 0, bft_annot[0]->size_annot * sizeof(uint8_t));
        }

        fclose(files_cc[nb_ccs_bft1]);

        do {
            cc = &(((CC*)bft1->node.CC_array)[i]);
            nb_skp = CEIL(cc->nb_elem, bft1->info_per_lvl[lvl_root].nb_ucs_skp);

            sprintf(&(output_filename[0][len_output_filename+1]), "%d", i);
            len_output_filename_cont = strlen(output_filename[0]);

            if (!packed_in_subtries) output_filename[0][len_output_filename_cont] = '_';
            else{
                file_tmp = fopen(output_filename[0], "w");
                ASSERT_NULL_PTR(file_tmp, "merging_BFT() 9\n");
            }

            for (j = 0; j < nb_skp; j++){

                uc = &(((UC*)cc->children)[j]);

                if (!packed_in_subtries){
                    sprintf(&(output_filename[0][len_output_filename_cont+1]), "%d", j);

                    file_tmp = fopen(output_filename[0], "w");
                    ASSERT_NULL_PTR(file_tmp, "merging_BFT() 21\n");
                }

                if (lvl_root){
                    write_UC_sparse(uc, bft1, file_tmp, bft1->info_per_lvl[lvl_root].size_kmer_in_bytes_minus_1,
                                    uc->nb_children, true);
                }
                else if (j != nb_skp-1) write_UC_sparse(uc, bft1, file_tmp, 0, bft1->info_per_lvl[lvl_root].nb_ucs_skp, false);
                else write_UC_sparse(uc, bft1, file_tmp, 0, cc->nb_elem - j * bft1->info_per_lvl[lvl_root].nb_ucs_skp, false);

                if (!packed_in_subtries) fclose(file_tmp);
                if (uc->suffixes != NULL) free(uc->suffixes);
            }

            for (j = 0; j < cc->nb_Node_children; j++){

                node = &(cc->children_Node_container[j]);

                if (!packed_in_subtries){
                    sprintf(&(output_filename[0][len_output_filename_cont+1]), "%d", j + nb_skp);

                    file_tmp = fopen(output_filename[0], "w");
                    ASSERT_NULL_PTR(file_tmp, "merging_BFT() 22\n");
                }

                l_write_Node(node, bft1, lvl_root-1, bft1->k - NB_CHAR_SUF_PREF, cut_lvl, true, file_tmp,
                             output_filename[0]);

                if (!packed_in_subtries) fclose(file_tmp);

                freeNode(node, lvl_root-1, bft1->info_per_lvl);
            }

            if (packed_in_subtries) fclose(file_tmp);

            i++;
        }
        while (IS_EVEN(((CC*)bft1->node.CC_array)[i-1].type));

        l_write_BFT_Root(bft1, output_prefix, cut_lvl, false);
    }

    free(kmer_comp[0]);
    free(kmer_comp_cpy[0]);
    free(bft_annot[0]->annot);
    free(bft_annot[0]);
    free(new_filename[0]);
    free(output_filename[0]);
    free(kmer_comp);
    free(kmer_comp_cpy);
    free(bft_annot);
    free(new_filename);
    free(output_filename);
    free(presence_bf);
    free(files_cc);
}

size_t l_insert_kmer(BFT_kmer* kmer, BFT* graph, va_list args){

    uint32_t substring_prefix;

    uint8_t kmer_cpy[SIZE_BYTES_SUF_PREF];

    int id_file;

    int nb_bytes = va_arg(args, int);

    FILE** files_cc = va_arg(args, FILE**);

    int* presence_bf = va_arg(args, int*);
    int nb_ccs = va_arg(args, int);

    BFT_annotation* bft_annot = get_annotation(kmer);

    if ((bft_annot->annot[0] & 0x3) == 3){

        uint32_t position = bft_annot->annot[0] >> 2;
        int i = 0;

        while ((i < bft_annot->size_annot) && (IS_ODD(bft_annot->annot[i]))){
            position |= ((uint32_t)(bft_annot->annot[i] >> 1)) << (6+(i-1)*7);
            i++;
        }

        if ((i >= bft_annot->size_annot) && (bft_annot->annot_ext != NULL)){
            if (IS_ODD(bft_annot->annot_ext[0])){
                position |= ((uint32_t)(bft_annot->annot_ext[0] >> 1)) << (6 + (i + bft_annot->size_annot - 1) * 7);
                i++;
            }
        }

        bft_annot->annot = extract_from_annotation_array_elem(graph->comp_set_colors, position, &bft_annot->size_annot);

        i = decomp_annotation(graph->ann_inf, bft_annot->annot, bft_annot->size_annot, NULL, 0, false);

        if (i != 0){
            bft_annot->annot = graph->ann_inf->annotation;
            bft_annot->size_annot = i;
        }

        bft_annot->annot_ext = NULL;
    }

    kmer_cpy[0] = reverse_word_8(kmer->kmer_comp[0]);
    kmer_cpy[1] = reverse_word_8(kmer->kmer_comp[1]);
    kmer_cpy[2] = reverse_word_8(kmer->kmer_comp[2]) & 0xc0;

    substring_prefix = (kmer_cpy[0] << 10) | (kmer_cpy[1] << 2) | (kmer_cpy[2] >> 6);

    if (presence_bf[substring_prefix] == -1) id_file = nb_ccs;
    else id_file = presence_bf[substring_prefix];

    if (fwrite(kmer->kmer_comp, sizeof(uint8_t), nb_bytes, files_cc[id_file]) != nb_bytes) ERROR("l_insert_kmer2()\n")

    if (bft_annot->annot_ext == NULL){
        if (fwrite(&(bft_annot->size_annot), sizeof(int), 1, files_cc[id_file]) != 1) ERROR("l_insert_kmer2()\n")
        if (fwrite(bft_annot->annot, sizeof(uint8_t), bft_annot->size_annot, files_cc[id_file]) != bft_annot->size_annot) ERROR("l_insert_kmer2()\n")
    }
    else{
        bft_annot->size_annot++;
        if (fwrite(&(bft_annot->size_annot), sizeof(int), 1, files_cc[id_file]) != 1) ERROR("l_insert_kmer2()\n")
        bft_annot->size_annot--;
        if (fwrite(bft_annot->annot, sizeof(uint8_t), bft_annot->size_annot, files_cc[id_file]) != bft_annot->size_annot) ERROR("l_insert_kmer2()\n")
        if (fwrite(bft_annot->annot_ext, sizeof(uint8_t), 1, files_cc[id_file]) != 1) ERROR("l_insert_kmer2()\n")
    }

    memset(graph->ann_inf->annotation, 0, graph->ann_inf->size_annot * sizeof(uint8_t));
    graph->ann_inf->size_annot = 0;

    free_BFT_annotation(bft_annot);
    return 1;
}

void l_insert_kmer_bis(BFT* bft, int lvl_root, uint8_t* kmer_comp, uint8_t* kmer_comp_cpy, int nb_bytes, int old_nb_genome_ids_bft,
                          bool overlap_genome_ids, BFT_annotation* bft_annot, int pos_start_search){

    uint32_t id;
    uint32_t size_id;

    uint32_t pow2_imin = 0;

    uint8_t* annot;
    uint8_t* annot_ext;
    uint8_t* annot_cplx;

    int size_annot;
    int size_annot_cplx;

    resultPresence* res;

    UC* uc;

    uint32_t i = 1;

    uint32_t* genome_ids = get_list_id_genomes(bft_annot, bft);

    memcpy(kmer_comp_cpy, kmer_comp, nb_bytes * sizeof(uint8_t));
    res = isKmerPresent(&(bft->node), bft, lvl_root, kmer_comp_cpy, bft->k);

    if (res->link_child == NULL){

        id = genome_ids[i] + old_nb_genome_ids_bft;
        if (overlap_genome_ids) id--;

        if (id >= pow2_imin){
            pow2_imin = round_up_next_highest_power2(id);
            size_id = get_nb_bytes_power2_annot_bis(id, pow2_imin);
        }

        memcpy(kmer_comp_cpy, kmer_comp, nb_bytes * sizeof(uint8_t));

        if (pos_start_search == -1) insertKmer_Node(&(bft->node), bft, lvl_root, kmer_comp_cpy, bft->k, kmer_comp, id, size_id, 0);
        else insertKmer_Node(&(bft->node), bft, lvl_root, kmer_comp_cpy, bft->k, kmer_comp, id, size_id, pos_start_search);

        i++;

        if (i <= genome_ids[0]){
            free(res);
            memcpy(kmer_comp_cpy, kmer_comp, nb_bytes * sizeof(uint8_t));
            res = isKmerPresent(&(bft->node), bft, lvl_root, kmer_comp_cpy, bft->k);
        }
    }

    for (; i <= genome_ids[0]; i++){

        id = genome_ids[i] + old_nb_genome_ids_bft;
        if (overlap_genome_ids) id--;

        if (id >= pow2_imin){
            pow2_imin = round_up_next_highest_power2(id);
            size_id = get_nb_bytes_power2_annot_bis(id, pow2_imin);
        }

        if (res->posFilter2 != 0){
            uc = (UC*)res->container;
            get_annot(uc, &annot, &annot_ext, &annot_cplx, &size_annot, &size_annot_cplx,
                           res->posFilter2, res->posFilter3, res->pos_sub_bucket);
        }
        else{
            uc = &(((UC*)((CC*)res->container)->children)[res->bucket]);
            get_annot(uc, &annot, &annot_ext, &annot_cplx, &size_annot, &size_annot_cplx,
                           res->posFilter2, res->posFilter3, res->pos_sub_bucket);
        }

        compute_best_mode(bft->ann_inf, bft->comp_set_colors, annot, size_annot, annot_ext, 1, id, size_id);

        if (bft->ann_inf->last_added != id){

            if (bft->ann_inf->min_size > uc->size_annot){
                if ((annot_ext == NULL) || (bft->ann_inf->min_size > uc->size_annot+1))
                    annot_ext = realloc_annotation(uc, res->posFilter2, res->posFilter3, bft->ann_inf->min_size, 0, res->pos_sub_bucket);
            }

            modify_mode_annotation(bft->ann_inf, &(uc->suffixes[res->pos_sub_bucket * (res->posFilter2 + uc->size_annot) + res->posFilter2]),
                                   uc->size_annot, annot_ext, 1, id, size_id);

            if ((annot_ext != NULL) && (annot_ext[0] == 0))
                delete_extend_annots(uc, res->posFilter2, res->posFilter3, res->pos_sub_bucket, res->pos_sub_bucket, 0, 0, 1);
        }

        reinit_annotation_inform(bft->ann_inf);
    }

    free(res);
    free(genome_ids);

    return;
}
*/
