#include "file_io.h"

void compress_annotations_disk(BFT_Root* bft, char* filename_bft){

    ASSERT_NULL_PTR(bft, "compress_annotations_disk()\n")
    ASSERT_NULL_PTR(filename_bft, "compress_annotations_disk()\n")

    Pvoid_t comp_annots = (PWord_t)NULL;
    Word_t Rc_word;

    memory_Used* mem;

    bool is_compressed = (bft->compressed == 1 ? true : false);

    int len_longest_annot;
    int lvl_bft = (bft->k / NB_CHAR_SUF_PREF) - 1;

    char* filename_bft_tmp = malloc((strlen(filename_bft) + 50) * sizeof(char));
    ASSERT_NULL_PTR(filename_bft_tmp, "compress_annotations_disk()\n")

    strcpy(filename_bft_tmp, filename_bft);
    strcpy(&filename_bft_tmp[strlen(filename_bft)], "_annots");

    mem = printMemoryUsedFromNode(&(bft->node), lvl_bft, bft->k, bft->info_per_lvl);
    len_longest_annot = (int) MAX(mem->size_biggest_annot+1, getMaxSize_annotation_array_elem(bft->comp_set_colors));
    free(mem);

    if (is_compressed){

        bft->compressed = 0;
        write_BFT_Root(bft, filename_bft, false);
        bft->compressed = 1;
    }
    else write_BFT_Root(bft, filename_bft, false);

    load_annotation_from_Node(&(bft->node), lvl_bft, bft->k, len_longest_annot,
                              bft->info_per_lvl, &comp_annots, bft->comp_set_colors, bft->ann_inf); //Load annot in a compressed way

    freeNode(&bft->node, lvl_bft, bft->info_per_lvl);
    free_annotation_array_elem(&(bft->comp_set_colors), &(bft->length_comp_set_colors));

    sort_annotations3(&comp_annots, len_longest_annot);
    write_partial_comp_set_colors(filename_bft_tmp, &comp_annots, len_longest_annot); //Write new annots compressed

    read_BFT_replace_comp_annots_bis(bft, filename_bft, filename_bft_tmp, &comp_annots, len_longest_annot, is_compressed);
    if (remove(filename_bft)) printf("Warning: Could not remove temporary file.\n");

    if (is_compressed) bft->compressed = 1;

    #if defined (_WORDx86)
        Word_t * PValue;

        uint8_t* it_index = calloc((len_longest_annot + CEIL(len_longest_annot, SIZE_BITS_UINT_8T - 1) + 4), sizeof(uint8_t));
        ASSERT_NULL_PTR(it_index, "compressKmers_from_KmerFiles()\n");

        JSLF(PValue, comp_annots, it_index);

        while (PValue != NULL){
            free(*PValue);
            JSLN(PValue, comp_annots, it_index);
        }

        free(it_index);
    #endif

    JSLFA(Rc_word, comp_annots);

    free_annotation_array_elem(&(bft->comp_set_colors), &(bft->length_comp_set_colors));
    read_annotation_array_elem(filename_bft_tmp, &(bft->comp_set_colors), &(bft->length_comp_set_colors));

    if (remove(filename_bft_tmp)) printf("Warning: Could not remove temporary file.\n");

    free(filename_bft_tmp);

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  insert_Genomes_from_KmerFiles(root, filenames, binary_files, size_kmer, ptr_ptr_annot_sorted)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert k-mers from k-mer files into a BFT
*  ---------------------------------------------------------------------------------------------------------------
*  root: ptr to the root of a BFT
*  filenames: array of filenames. The files contains the k-mers to insert.
*  binary_files: Indicate if the files contains k-mers (ASCII) or compressed k-mers (2 bits per nuc.)
*  size_kmer: length k of k-mers in files
*  ---------------------------------------------------------------------------------------------------------------
*/
void insert_Genomes_from_KmerFiles(BFT_Root* root, int nb_files, char** filenames, int binary_files, char* filename_bft){

    ASSERT_NULL_PTR(root,"insert_Genomes_from_KmerFiles()")
    ASSERT_NULL_PTR(filenames,"insert_Genomes_from_KmerFiles()")

    struct timeval tval_before, tval_after, tval_last, tval_result;
    gettimeofday(&tval_before, NULL);
    tval_last = tval_before;

    FILE* file;

    uint8_t* array_kmers = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(array_kmers,"insert_Genomes_from_KmerFiles()")

    char* line = calloc(100, sizeof(char));
    ASSERT_NULL_PTR(line,"insert_Genomes_from_KmerFiles()")

    char* str_tmp;

    int i = 0, j = 0, k = 0;
    int size_id_genome = 0;
    int lvl_root = (root->k / NB_CHAR_SUF_PREF) - 1;
    int nb_bytes_kmer = CEIL(root->k*2, SIZE_BITS_UINT_8T);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer;

    size_t return_fread;

    uint64_t kmers_read;

    for (i = 0; i < nb_files; i++){ //For each file in input

        k = 0;
        j = 0;
        kmers_read = 0;

        str_tmp = basename(filenames[i]);
        add_genomes_BFT_Root(1, &str_tmp, root);

        size_id_genome = get_nb_bytes_power2_annot(root->nb_genomes-1);

        file = fopen(filenames[i], "r");
        ASSERT_NULL_PTR(file,"insert_Genomes_from_KmerFiles()")

        printf("\nFile %d: %s\n\n", root->nb_genomes-1, filenames[i]);

        if (binary_files){

            if (fgets(line, 100, file) != NULL) k = atoi(line);
            else ERROR("Cannot read header of the file")

            if (fgets(line, 100, file) != NULL) printf("%d %d-mers in the file\n\n", atoi(line), k);
            else ERROR("Cannot read header of the file")

            while ((!ferror(file)) && (!feof(file))){

                return_fread = fread(array_kmers, (size_t)nb_bytes_kmer, (size_t)nb_kmer_in_buf, file);

                insertKmers(root, array_kmers, return_fread, root->nb_genomes-1, size_id_genome);

                memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));

                if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+return_fread)%PRINT_EVERY_X_KMERS)){
                    printf("%" PRIu64 " kmers read\n", kmers_read+return_fread);
                }

                kmers_read += return_fread;
            }
        }
        else {
            while (fgets(line, 100, file) != NULL){

                if (parseKmerCount(line, root->k, array_kmers, k) == 1){

                    k += nb_bytes_kmer;
                    j++;

                    if (j == nb_kmer_in_buf){
                        insertKmers(root, array_kmers, nb_kmer_in_buf, root->nb_genomes-1, size_id_genome);

                        j = 0;
                        k = 0;
                        memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));

                        if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+nb_kmer_in_buf)%PRINT_EVERY_X_KMERS)){
                            printf("%" PRIu64 " kmers read\n", kmers_read+nb_kmer_in_buf);
                        }

                        kmers_read += nb_kmer_in_buf;
                    }
                }
            }

            insertKmers(root, array_kmers, j, root->nb_genomes-1, size_id_genome);
            kmers_read += j;

            memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));
        }

        fclose(file);

        if (root->treshold_compression && (root->nb_genomes - 1 > 5) && ((root->nb_genomes - 1) % root->treshold_compression == 0))
            compress_annotations_disk(root, filename_bft);

        gettimeofday(&tval_after, NULL);

        time_spent(&tval_last, &tval_after, &tval_result);
        printf("\nElapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

        time_spent(&tval_before, &tval_after, &tval_result);
        printf("Total elapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

        printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);
        printf("Current memory: %llu mb\n", ((unsigned long long int)getCurrentRSS())/1024);

        tval_last = tval_after;
    }

    memory_Used* mem = printMemoryUsedFromNode(&(root->node), lvl_root, root->k, root->info_per_lvl);
    printMemory(mem);

    free(mem);
    free(line);
    free(array_kmers);

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  insert_Genomes_from_FASTxFiles(root, filenames, size_kmer, ptr_ptr_annot_sorted)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert k-mers from FASTx files into a BFT
*  ---------------------------------------------------------------------------------------------------------------
*  root: ptr to the root of a BFT
*  filenames: array of FASTx filenames
*  size_kmer: length k of k-mers to extract from the FASTx files
*  ---------------------------------------------------------------------------------------------------------------
*/
void insert_Genomes_from_FASTxFiles(BFT_Root* root, int nb_files, char** filenames){

    ASSERT_NULL_PTR(root,"insert_Genomes_from_FASTxFiles()")
    ASSERT_NULL_PTR(filenames,"insert_Genomes_from_FASTxFiles()")

    /*struct timeval tval_before, tval_after, tval_last, tval_result;
    gettimeofday(&tval_before, NULL);
    tval_last = tval_before;

    int i = 0;
    int size_buf_tmp = 0; //How many characters are stored in buf_tmp
    int nb_kmers_buf = 0;
    int size_id_genome = 0;
    int length_comp_set_colors_tmp = 0;
    int nb_cell_kmer = CEIL(size_kmer*2, SIZE_BITS_UINT_8T); //Size of kmers in bytes

    annotation_array_elem* comp_set_colors_tmp = NULL;

    Pvoid_t PJArray = (PWord_t)NULL;
    Word_t Rc_word;

    char* str_tmp;

    char* buf_tmp = calloc((size_kmer-1)*2, sizeof(char)); //Allocate temporary buffer
    ASSERT_NULL_PTR(buf_tmp,"insert_Genomes_from_FASTxFiles()")

    uint8_t* tab_kmers = calloc(SIZE_BUFFER*nb_cell_kmer, sizeof(uint8_t)); //Allocate buffer for kmers
    ASSERT_NULL_PTR(tab_kmers,"insert_Genomes_from_FASTxFiles()")

    uint64_t kmers_read = 0;
    uint64_t tmp_kmers_read = 0;

    for (i = 0; i < nb_files; i++){ //For each file in input

        size_buf_tmp = 0;
        kmers_read = 0;
        tmp_kmers_read = 0;
        nb_kmers_buf = 0;

        str_tmp = basename(filenames[i]);
        add_genomes_BFT_Root(1, &str_tmp, root);

        size_id_genome = get_nb_bytes_power2_annot(root->nb_genomes-1);

        int fp = open(filenames[i], O_RDONLY); //Open it
        kseq_t *seq = kseq_init(fp); //initialize the parser for this file
        int size_seq = kseq_read(seq, -1); //Start reading file, seq contains a buffer with a part of a sequence from the file

        printf("\nFile : %s\n\n", filenames[i]);

        while (size_seq > -1) { //While the end of the file is not reached

            if (size_seq > 0) size_buf_tmp = 0; //New sequence

            int current_buf_length = seq->seq.l - seq->seq.z; //Number of characters put into the seq buffer

            if (current_buf_length > 0){ //If the seq buffer is not empty

                nb_kmers_buf = MAX(current_buf_length-size_kmer+1, 0); //Number of kmers it is possible to read in seq buffer

                if (size_buf_tmp == 0){ //If the number of characters in the temporary buffer is 0
                    if (nb_kmers_buf != 0){ //If there is at least one kmer in the seq buffer
                        memcpy(buf_tmp, &(seq->seq.s[nb_kmers_buf]), size_kmer-1); //Copy the last size_kmer-1 characters of seq-buffer into buf_tmp
                        size_buf_tmp = (size_kmer-1);
                    }
                    else{
                        memcpy(buf_tmp, &(seq->seq.s[0]), current_buf_length); //Copy the content of seq buffer into buf_tmp
                        size_buf_tmp = current_buf_length;
                    }
                }
                else { //If the number of characters in the temporary buffer is not 0
                    //Insertion of kmers overlapping the last buffer and the current one (they are in buf_tmp)
                    int size_to_copy = MIN(size_kmer-1, current_buf_length);
                    memcpy(&(buf_tmp[size_buf_tmp]), seq->seq.s, size_to_copy);
                    size_buf_tmp += size_to_copy;

                    int nb_kmers = size_buf_tmp - size_kmer + 1;
                    if (nb_kmers > 0){
                        parseSequenceBuffer(buf_tmp, tab_kmers, &nb_kmers, size_kmer, nb_cell_kmer); //Read buf_tmp, extract the kmers in tab_kmers
                        insertKmers(root, tab_kmers, size_kmer, nb_kmers, root->nb_genomes-1, size_id_genome, 0); //Insert the kmers into the tree
                        memset(tab_kmers, 0, nb_kmers*nb_cell_kmer*sizeof(uint8_t)); //Reinit tab_kmers
                        tmp_kmers_read = nb_kmers;
                    }
                    else tmp_kmers_read = 0;

                    if (nb_kmers_buf != 0){
                        memcpy(buf_tmp, &(seq->seq.s[nb_kmers_buf]), size_kmer-1);
                        size_buf_tmp = size_kmer-1;
                    }
                    else{
                        memcpy(buf_tmp, &(seq->seq.s[0]), current_buf_length);
                        size_buf_tmp = current_buf_length;
                    }
                }

                //Extraction of buffer's kmers. Insertion in the tree.
                if (nb_kmers_buf > 0){
                    parseSequenceBuffer(seq->seq.s, tab_kmers, &nb_kmers_buf, size_kmer, nb_cell_kmer);
                    insertKmers(root, tab_kmers, size_kmer, nb_kmers_buf, root->nb_genomes-1, size_id_genome, 0);
                    memset(tab_kmers, 0, nb_kmers_buf*nb_cell_kmer*sizeof(uint8_t));
                    tmp_kmers_read += nb_kmers_buf;
                }

                //Display how many kmers were read
                if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+tmp_kmers_read)%PRINT_EVERY_X_KMERS))
                    printf("%" PRIu64 " kmers read\n", kmers_read+tmp_kmers_read);

                kmers_read += tmp_kmers_read;
            }

            size_seq = kseq_read(seq, size_seq);
        }

        if (root->treshold_compression != 0){
            if ((root->nb_genomes-1 > 5) && ((root->nb_genomes-1)%root->treshold_compression == 0)){

                load_annotation_from_Node(&(root->node), size_kmer, info_per_lvl, &PJArray, root->comp_set_colors);

                comp_set_colors_tmp = root->comp_set_colors;

                length_comp_set_colors_tmp = root->length_comp_set_colors;

                memory_Used* mem = printMemoryUsedFromNode(&(root->node), lvl_root, root->k, root->info_per_lvl);
                root->comp_set_colors = sort_annotations(&PJArray, &(root->length_comp_set_colors), mem->size_biggest_annot);
                free(mem);

                compress_annotation_from_Node(&(root->node), size_kmer, info_per_lvl, &PJArray, comp_set_colors_tmp);

                free_annotation_array_elem(&comp_set_colors_tmp, length_comp_set_colors_tmp);

                #if defined (_WORDx86)
                    Word_t * PValue;

                    uint8_t* it_index = calloc((len_longest_annot + CEIL(len_longest_annot, SIZE_BITS_UINT_8T) + 4), sizeof(uint8_t));
                    ASSERT_NULL_PTR(it_index, "sort_annotations()");

                    JSLF(PValue, PJArray, it_index);

                    while (PValue != NULL){
                        free(*PValue);
                        JSLN(PValue, PJArray, it_index);
                    }

                    free(it_index);
                #endif

                JSLFA(Rc_word, PJArray);
            }
        }

        kseq_destroy(seq);
        close(fp);

        gettimeofday(&tval_after, NULL);

        time_spent(&tval_last, &tval_after, &tval_result);
        printf("\nElapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

        time_spent(&tval_before, &tval_after, &tval_result);
        printf("Total elapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

        printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);
        printf("Current memory: %llu mb\n", ((unsigned long long int)getCurrentRSS())/1024);

        tval_last = tval_after;
    }

    memory_Used* mem = printMemoryUsedFromNode(&(root->node), size_kmer, info_per_lvl);
    printMemory(mem);

    free(mem);
    free(buf_tmp);
    free(tab_kmers);*/

    return;
}

int queryBFT_kmerPresences_from_KmerFiles(BFT_Root* root, char* query_filename, int binary_file, char* output_filename){

    ASSERT_NULL_PTR(root,"queryBFT_kmerPresences_from_KmerFiles()")
    ASSERT_NULL_PTR(query_filename,"queryBFT_kmerPresences_from_KmerFiles()")

    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    const char comma = ',';

    //int annot_present;
    int size_annot;
    int size_annot_cplx;
    int size_annot_res;

    int i = 0;
    int j = 0;
    int k = 0;
    int nb_kmers_present = 0;
    int nb_bytes_kmer = CEIL(root->k*2, SIZE_BITS_UINT_8T);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer;
    int lvl_root = (root->k / NB_CHAR_SUF_PREF) - 1;

    uint64_t kmers_read = 0;

    FILE* file_query;
    FILE* file_output;

    resultPresence* res;

    size_t return_fread;

    uint8_t* annot;
    uint8_t* annot_ext;
    uint8_t* annot_cplx;

    uint8_t* annot_res = calloc(CEIL(root->nb_genomes+2, SIZE_BITS_UINT_8T), sizeof(uint8_t));
    ASSERT_NULL_PTR(annot_res,"queryBFT_kmerPresences_from_KmerFiles()")

    uint8_t* array_kmers = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(array_kmers,"queryBFT_kmerPresences_from_KmerFiles()")

    char* line = calloc(100, sizeof(char));
    ASSERT_NULL_PTR(line,"queryBFT_kmerPresences_from_KmerFiles()")

    file_query = fopen(query_filename, "r");
    ASSERT_NULL_PTR(file_query,"queryBFT_kmerPresences_from_KmerFiles()")

    file_output = fopen(output_filename, "w");
    ASSERT_NULL_PTR(file_output,"queryBFT_kmerPresences_from_KmerFiles()")

    printf("\nQuerying BFT for k-mers in %s\n\n", query_filename);

    for (i=0; i<root->nb_genomes-1; i++){
        fwrite(root->filenames[i], sizeof(char), strlen(root->filenames[i])-1, file_output);
        fwrite(&comma, sizeof(char), 1, file_output);
    }

    fwrite(root->filenames[i], sizeof(char), strlen(root->filenames[i]), file_output);

    if (binary_file){

        if (fgets(line, 100, file_query) == NULL) ERROR("Cannot read header of the queries file")
        if (fgets(line, 100, file_query) == NULL) ERROR("Cannot read header of the queries file")

        while ((!ferror(file_query)) && (!feof(file_query))){

            return_fread = fread(array_kmers, (size_t)nb_bytes_kmer, (size_t)nb_kmer_in_buf, file_query);

            for (k=0; k<(int)return_fread; k++){

                res = isKmerPresent(&(root->node), root, lvl_root, &(array_kmers[k*nb_bytes_kmer]), root->k);

                if (res->link_child != NULL){

                    if (res->posFilter2 != 0){
                        get_annot((UC*)res->container, &annot, &annot_ext, &annot_cplx, &size_annot,
                                        &size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);
                    }
                    else{
                        get_annot(&(((UC*)((CC*)res->container)->children)[res->bucket]), &annot, &annot_ext,
                                       &annot_cplx, &size_annot, &size_annot_cplx, res->posFilter2, res->posFilter3,
                                       res->pos_sub_bucket);
                    }

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

                    printAnnotation_CSV(file_output, annot_res, size_annot_res, NULL, 0, root->nb_genomes-1, root->comp_set_colors);

                    nb_kmers_present++;
                }
                else {
                    annot_res[0] = 0;
                    printAnnotation_CSV(file_output, annot_res, 1, NULL, 0, root->nb_genomes-1, root->comp_set_colors);
                }

                free(res);
            }

            //if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+return_fread)%PRINT_EVERY_X_KMERS))
            //    printf("%" PRIu64 " kmers read\n", kmers_read+return_fread);

            kmers_read += return_fread;

            memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));
        }
    }
    else{

        while (fgets(line, 100, file_query) != NULL){

            if (parseKmerCount(line, root->k, array_kmers, k) == 1){
                k += nb_bytes_kmer;
                j++;

                if (j == nb_kmer_in_buf){

                    for (i=0; i<nb_kmer_in_buf; i++){

                        res = isKmerPresent(&(root->node), root, lvl_root, &(array_kmers[i*nb_bytes_kmer]), root->k);

                        if (res->link_child != NULL){

                            if (res->posFilter2 != 0){
                                get_annot((UC*)res->container, &annot, &annot_ext, &annot_cplx, &size_annot,
                                                &size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);
                            }
                            else{
                                get_annot(&(((UC*)((CC*)res->container)->children)[res->bucket]), &annot, &annot_ext,
                                               &annot_cplx, &size_annot, &size_annot_cplx, res->posFilter2, res->posFilter3,
                                               res->pos_sub_bucket);
                            }

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

                            printAnnotation_CSV(file_output, annot_res, size_annot_res, NULL, 0, root->nb_genomes-1, root->comp_set_colors);

                            nb_kmers_present++;
                        }
                        else {
                            annot_res[0] = 0;
                            printAnnotation_CSV(file_output, annot_res, 1, NULL, 0, root->nb_genomes-1, root->comp_set_colors);
                        }

                        free(res);
                    }

                    j = 0;
                    k = 0;
                    memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));

                    //if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+nb_kmer_in_buf)%PRINT_EVERY_X_KMERS))
                    //    printf("%" PRIu64 " kmers read\n", kmers_read+nb_kmer_in_buf);

                    kmers_read += nb_kmer_in_buf;
                }
            }
        }

        for (i=0; i<j; i++){

            res = isKmerPresent(&(root->node), root, lvl_root, &(array_kmers[i*nb_bytes_kmer]), root->k);

            if (res->link_child != NULL){

                if (res->posFilter2 != 0){
                    get_annot((UC*)res->container, &annot, &annot_ext, &annot_cplx, &size_annot,
                                    &size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);
                }
                else{
                    get_annot(&(((UC*)((CC*)res->container)->children)[res->bucket]), &annot, &annot_ext,
                                   &annot_cplx, &size_annot, &size_annot_cplx, res->posFilter2, res->posFilter3,
                                   res->pos_sub_bucket);
                }

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

                printAnnotation_CSV(file_output, annot_res, size_annot_res, NULL, 0, root->nb_genomes-1, root->comp_set_colors);

                nb_kmers_present++;
            }
            else {
                annot_res[0] = 0;
                printAnnotation_CSV(file_output, annot_res, 1, NULL, 0, root->nb_genomes-1, root->comp_set_colors);
            }

            free(res);
        }

        memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));
    }

    fclose(file_query);
    fclose(file_output);

    free(array_kmers);
    free(annot_res);
    free(line);

    gettimeofday(&tval_after, NULL);
    time_spent(&tval_before, &tval_after, &tval_result);

    printf("\nElapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

    printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);
    printf("Current memory: %llu mb\n", ((unsigned long long int)getCurrentRSS())/1024);

    return nb_kmers_present;
}

int queryBFT_kmerBranching_from_KmerFiles(BFT_Root* root, char* query_filename, int binary_file){

    ASSERT_NULL_PTR(root,"queryBFT_kmerBranching_from_KmerFiles()")
    ASSERT_NULL_PTR(query_filename,"queryBFT_kmerBranching_from_KmerFiles()")

    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    int i = 0;
    int j = 0;
    int k = 0;
    int count_branching_node = 0;
    int lvl_root = (root->k / NB_CHAR_SUF_PREF) - 1;
    int nb_bytes_kmer = CEIL(root->k*2, SIZE_BITS_UINT_8T);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer;

    uint64_t kmers_read = 0;

    FILE* file;

    size_t return_fread;

    uint8_t* array_kmers = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(array_kmers,"queryBFT_kmerBranching_from_KmerFiles()")

    char* line = calloc(100, sizeof(char));
    ASSERT_NULL_PTR(line,"queryBFT_kmerBranching_from_KmerFiles()")

    root->skip_sp = build_skip_nodes(&(root->node));

    file = fopen(query_filename, "r");
    ASSERT_NULL_PTR(file,"queryBFT_kmerBranching_from_KmerFiles()")

    printf("\nQuerying BFT for branching k-mers in %s\n\n", query_filename);

    if (binary_file){

        if (fgets(line, 100, file) == NULL) ERROR("Cannot read header of the file")
        if (fgets(line, 100, file) == NULL) ERROR("Cannot read header of the file")

        while ((!ferror(file)) && (!feof(file))){

            return_fread = fread(array_kmers, (size_t)nb_bytes_kmer, (size_t)nb_kmer_in_buf, file);

            for (k=0; k<(int)return_fread; k++){

                if (isBranchingRight(&(root->node), root, lvl_root, &(array_kmers[k*nb_bytes_kmer]), root->k) > 1){
                    count_branching_node++;
                }
                else if (isBranchingLeft(&(root->node), root, lvl_root, &(array_kmers[k*nb_bytes_kmer]), root->k) > 1){
                    count_branching_node++;
                }
            }

            if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+return_fread)%PRINT_EVERY_X_KMERS))
                printf("%" PRIu64 " kmers read\n", kmers_read+return_fread);

            kmers_read += return_fread;

            memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));
        }
    }
    else{

        while (fgets(line, 100, file) != NULL){

            if (parseKmerCount(line, root->k, array_kmers, k) == 1){
                k += nb_bytes_kmer;
                j++;

                if (j == nb_kmer_in_buf){

                    for (i=0; i<nb_kmer_in_buf; i++){

                        if (isBranchingRight(&(root->node), root, lvl_root, &(array_kmers[i*nb_bytes_kmer]), root->k) > 1){
                            count_branching_node++;
                        }
                        else if (isBranchingLeft(&(root->node), root, lvl_root, &(array_kmers[i*nb_bytes_kmer]), root->k) > 1){
                            count_branching_node++;
                        }
                    }

                    j = 0;
                    k = 0;
                    memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));

                    if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+nb_kmer_in_buf)%PRINT_EVERY_X_KMERS))
                        printf("%" PRIu64 " kmers read\n", kmers_read+nb_kmer_in_buf);

                    kmers_read += nb_kmer_in_buf;
                }
            }
        }

        for (i=0; i<j; i++){

            if (isBranchingRight(&(root->node), root, lvl_root, &(array_kmers[i*nb_bytes_kmer]), root->k) > 1){
                count_branching_node++;
            }
            else if (isBranchingLeft(&(root->node), root, lvl_root, &(array_kmers[i*nb_bytes_kmer]), root->k) > 1){
                count_branching_node++;
            }
        }

        memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));
    }

    fclose(file);

    free(array_kmers);
    free(line);

    if (root->skip_sp != NULL) free_skip_nodes(&(root->node), root->skip_sp);

    gettimeofday(&tval_after, NULL);
    time_spent(&tval_before, &tval_after, &tval_result);

    printf("\nElapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

    printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);
    printf("Current memory: %llu mb\n", ((unsigned long long int)getCurrentRSS())/1024);

    return count_branching_node;
}

/*void par_insert_Genomes_from_KmerFiles(int nb_files, char** filenames, int binary_files, int size_kmer,
                                       int treshold_compression, char* prefix_output, int cut_lvl, int memory_limit){

    ASSERT_NULL_PTR(filenames,"insert_Genomes_from_KmerFiles()")

    struct timeval tval_before, tval_after, tval_result;

    Pvoid_t* PJArray;

    Word_t Rc_word;

    annotation_array_elem* comp_set_colors_tmp;

    BFT_Root** root;

    FILE** file;

    uint8_t** array_kmers;

    char** line;
    char** output_filename;
    char** output_filename2;

    char* str_tmp;
    char* output_filename3;

    int* nb_bfts_on_disk;

    int steps = 2;

    int i, j, k;
    int lvl_root;
    int nb_bytes_kmer;
    int nb_kmer_in_buf;
    int it_thread, thread_id;
    int len_longest_annot;
    int len_output_filename;
    int len_output_filename2;
    int nb_merging;
    int nb_threads;
    int size_id_genome;
    int it_bft_thread;
    int length_comp_set_colors_tmp;

    size_t return_fread;

    uint64_t kmers_read;

    //if (memory_limit > 0) omp_set_num_threads(1);
    //else memory_limit = INT_MAX;

    omp_set_num_threads(1);

    #pragma omp parallel                                                                                                                    \
    shared(line, array_kmers, file, tval_before, tval_after, tval_result, PJArray, nb_threads, nb_bfts_on_disk,                             \
           lvl_root, nb_bytes_kmer, nb_kmer_in_buf, root, output_filename, output_filename2, len_output_filename,)                          \
    private(i, j, k, thread_id, Rc_word, comp_set_colors_tmp, str_tmp, size_id_genome, return_fread, it_thread,                             \
            length_comp_set_colors_tmp, len_longest_annot, kmers_read, len_output_filename2)
    {

        #pragma omp single
        {
            gettimeofday(&tval_before, NULL);

            nb_threads = omp_get_num_threads();

            length_comp_set_colors_tmp = 0;
            size_id_genome = 0;

            comp_set_colors_tmp = NULL;

            len_output_filename = strlen(prefix_output);

            nb_bfts_on_disk = calloc(nb_threads, sizeof(int));
            ASSERT_NULL_PTR(nb_bfts_on_disk, "par_insert_Genomes_from_KmerFiles()\n")

            output_filename = malloc(nb_threads * sizeof(char*));
            ASSERT_NULL_PTR(output_filename, "par_insert_Genomes_from_KmerFiles()\n")

            output_filename2 = malloc(nb_threads * sizeof(char*));
            ASSERT_NULL_PTR(output_filename2, "par_insert_Genomes_from_KmerFiles()\n")

            output_filename3 = malloc((len_output_filename + 30) * sizeof(char));
            ASSERT_NULL_PTR(output_filename3, "merging_BFT()\n");

            strcpy(output_filename3, prefix_output);
            strcpy(&output_filename3[strlen(output_filename3)], "_tmp");

            line = malloc(nb_threads * sizeof(char*));
            ASSERT_NULL_PTR(line, "par_insert_Genomes_from_KmerFiles()\n")

            array_kmers = malloc(nb_threads * sizeof(uint8_t*));
            ASSERT_NULL_PTR(array_kmers, "par_insert_Genomes_from_KmerFiles()\n")

            file = malloc(nb_threads * sizeof(FILE*));
            ASSERT_NULL_PTR(file, "par_insert_Genomes_from_KmerFiles()\n")

            PJArray = malloc(nb_threads * sizeof(PWord_t));
            ASSERT_NULL_PTR(PJArray, "par_insert_Genomes_from_KmerFiles()\n")

            root = malloc(nb_threads * sizeof(BFT_Root*));
            ASSERT_NULL_PTR(root, "par_insert_Genomes_from_KmerFiles()\n")

            for (it_thread = 0; it_thread < nb_threads; it_thread++){

                line[it_thread] = calloc(100, sizeof(char));
                ASSERT_NULL_PTR(line[it_thread], "par_insert_Genomes_from_KmerFiles()\n")

                array_kmers[it_thread] = calloc(SIZE_BUFFER, sizeof(uint8_t));
                ASSERT_NULL_PTR(array_kmers[it_thread], "par_insert_Genomes_from_KmerFiles()\n")

                output_filename[it_thread] = malloc((len_output_filename + 30) * sizeof(char));
                ASSERT_NULL_PTR(output_filename[it_thread], "merging_BFT()\n");

                output_filename2[it_thread] = malloc((len_output_filename + 30) * sizeof(char));
                ASSERT_NULL_PTR(output_filename2[it_thread], "merging_BFT()\n");

                strcpy(output_filename[it_thread], prefix_output);
                strcpy(output_filename2[it_thread], prefix_output);

                PJArray[it_thread] = (PWord_t)NULL;

                root[it_thread] = createBFT_Root(size_kmer, treshold_compression, 0);
            }

            lvl_root = (root[0]->k / NB_CHAR_SUF_PREF) - 1;
            nb_bytes_kmer = CEIL(root[0]->k * 2, SIZE_BITS_UINT_8T);
            nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer;
        }

        #pragma omp for
        for (i = 0; i < nb_files; i++){ //For each file in input

            thread_id = omp_get_thread_num();

            kmers_read = 0;
            k = 0;
            j = 0;

            str_tmp = basename(filenames[i]);
            add_genomes_BFT_Root(1, &str_tmp, root[thread_id]);

            size_id_genome = get_nb_bytes_power2_annot(root[thread_id]->nb_genomes-1);

            file[thread_id] = fopen(filenames[i], "r");
            ASSERT_NULL_PTR(file[thread_id],"insert_Genomes_from_KmerFiles()")

            printf("Processing file %s\n", filenames[i]);

            if (binary_files){

                if (fgets(line[thread_id], 100, file[thread_id]) != NULL) k = atoi(line[thread_id]);
                else ERROR("Cannot read header of the file")

                if (fgets(line[thread_id], 100, file[thread_id]) == NULL) ERROR("Cannot read header of the file")

                while ((return_fread = fread(array_kmers[thread_id], nb_bytes_kmer, nb_kmer_in_buf, file[thread_id])) == nb_kmer_in_buf)
                {
                    insertKmers(root[thread_id], array_kmers[thread_id], return_fread, root[thread_id]->nb_genomes-1, size_id_genome);

                    memset(array_kmers[thread_id], 0, SIZE_BUFFER * sizeof(uint8_t));

                    if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+return_fread)%PRINT_EVERY_X_KMERS)){
                    //    printf("%" PRIu64 " kmers read\n", kmers_read+return_fread);

                        if (((unsigned long long int)getCurrentRSS())/1024 >= memory_limit){

                            sprintf(&(output_filename[thread_id][len_output_filename]), "%d", thread_id);
                            len_output_filename2 = strlen(output_filename[thread_id]);
                            output_filename[thread_id][len_output_filename2] = '_';

                            sprintf(&(output_filename[thread_id][len_output_filename2+1]), "%d", nb_bfts_on_disk[thread_id]);
                            nb_bfts_on_disk[thread_id]++;

                            write_BFT_Root_sparse(root[thread_id], output_filename[thread_id], false);
                            freeBFT_Root(root[thread_id]);

                            root[thread_id] = createBFT_Root(size_kmer, treshold_compression, 0);
                            str_tmp = basename(filenames[i]);
                            add_genomes_BFT_Root(1, &str_tmp, root[thread_id]);
                        }
                    }

                    kmers_read += return_fread;
                    return_fread = 0;
                }

                insertKmers(root[thread_id], array_kmers[thread_id], return_fread, root[thread_id]->nb_genomes-1, size_id_genome);
            }
            else {
                while (fgets(line[thread_id], 100, file[thread_id]) != NULL){

                    if (parseKmerCount(line[thread_id], root[thread_id]->k, array_kmers[thread_id], k) == 1){
                        k += nb_bytes_kmer;
                        j++;

                        if (j == nb_kmer_in_buf){
                            insertKmers(root[thread_id], array_kmers[thread_id], nb_kmer_in_buf, root[thread_id]->nb_genomes-1, size_id_genome);

                            j = 0;
                            k = 0;
                            memset(array_kmers[thread_id], 0, SIZE_BUFFER * sizeof(uint8_t));

                            if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+nb_kmer_in_buf)%PRINT_EVERY_X_KMERS)){
                                //printf("%" PRIu64 " kmers read\n", kmers_read+nb_kmer_in_buf);

                                if (((unsigned long long int)getCurrentRSS())/1024 >= memory_limit){

                                    sprintf(&(output_filename[thread_id][len_output_filename]), "%d", thread_id);
                                    len_output_filename2 = strlen(output_filename[thread_id]);
                                    output_filename[thread_id][len_output_filename2] = '_';

                                    sprintf(&(output_filename[thread_id][len_output_filename2+1]), "%d", nb_bfts_on_disk[thread_id]);
                                    nb_bfts_on_disk[thread_id]++;

                                    write_BFT_Root_sparse(root[thread_id], output_filename[thread_id], false);
                                    freeBFT_Root(root[thread_id]);

                                    root[thread_id] = createBFT_Root(size_kmer, treshold_compression, 0);
                                    str_tmp = basename(filenames[i]);
                                    add_genomes_BFT_Root(1, &str_tmp, root[thread_id]);
                                }
                            }

                            kmers_read += nb_kmer_in_buf;
                        }
                    }
                }

                insertKmers(root[thread_id], array_kmers[thread_id], j, root[thread_id]->nb_genomes-1, size_id_genome);
                kmers_read += j;

                memset(array_kmers[thread_id], 0, SIZE_BUFFER * sizeof(uint8_t));
            }

            fclose(file[thread_id]);

            if (root[thread_id]->treshold_compression != 0){
                if ((root[thread_id]->nb_genomes-1 > 5) && ((root[thread_id]->nb_genomes-1)%root[thread_id]->treshold_compression == 0)){

                    memory_Used* mem = printMemoryUsedFromNode(&(root[thread_id]->node), lvl_root, root[thread_id]->k, root[thread_id]->info_per_lvl);
                    len_longest_annot = MAX(mem->size_biggest_annot+1, getMaxSize_annotation_array_elem(root[thread_id]->comp_set_colors));
                    free(mem);

                    load_annotation_from_Node(&(root[thread_id]->node), lvl_root, root[thread_id]->k, len_longest_annot, root[thread_id]->info_per_lvl,
                                              &(PJArray[thread_id]), root[thread_id]->comp_set_colors, root[thread_id]->ann_inf);

                    comp_set_colors_tmp = root[thread_id]->comp_set_colors;
                    length_comp_set_colors_tmp = root[thread_id]->length_comp_set_colors;

                    root[thread_id]->comp_set_colors = sort_annotations(&(PJArray[thread_id]), &(root[thread_id]->length_comp_set_colors), len_longest_annot);

                    if (root[thread_id]->comp_set_colors != NULL){
                        compress_annotation_from_Node(&(root[thread_id]->node), lvl_root, root[thread_id]->k, root[thread_id]->info_per_lvl, &(PJArray[thread_id]),
                                                      comp_set_colors_tmp, root[thread_id]->ann_inf);

                        free_annotation_array_elem(&comp_set_colors_tmp, &length_comp_set_colors_tmp);
                    }

                    #if defined (_WORDx86)
                        Word_t * PValue;

                        uint8_t* it_index = calloc((len_longest_annot + CEIL(len_longest_annot, SIZE_BITS_UINT_8T) + 4), sizeof(uint8_t));
                        ASSERT_NULL_PTR(it_index, "sort_annotations()");

                        JSLF(PValue, PJArray[thread_id], it_index);

                        while (PValue != NULL){
                            free(*PValue);
                            JSLN(PValue, PJArray[thread_id], it_index);
                        }

                        free(it_index);
                    #endif

                    JSLFA(Rc_word, PJArray[thread_id]);
                }
            }

            gettimeofday(&tval_after, NULL);

            time_spent(&tval_before, &tval_after, &tval_result);
            printf("Total elapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

            printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);
            printf("Current memory: %llu mb\n", ((unsigned long long int)getCurrentRSS())/1024);
        }

        #pragma omp for
        for (it_thread = 0; it_thread < nb_threads; it_thread++){

            sprintf(&(output_filename[it_thread][len_output_filename]), "%d", it_thread);
            len_output_filename2 = strlen(output_filename[it_thread]);
            output_filename[it_thread][len_output_filename2] = '_';

            sprintf(&(output_filename[it_thread][len_output_filename2+1]), "%d", nb_bfts_on_disk[it_thread]);
            nb_bfts_on_disk[it_thread]++;

            write_BFT_Root_sparse(root[it_thread], output_filename[it_thread], false);

            free(line[it_thread]);
            free(array_kmers[it_thread]);
            freeBFT_Root(root[it_thread]);
        }
    }

    free(PJArray);
    free(file);
    free(array_kmers);
    free(line);
    free(root);

    omp_set_nested(false);
    omp_set_num_threads(1);

    if (omp_get_num_threads() == 1){

        strcpy(&output_filename[0][len_output_filename], "0_0");
        strcpy(&output_filename2[0][len_output_filename], "0_0_pkd");
        read_cut_BFT_Root(output_filename[0], output_filename2[0], cut_lvl, true);
        strcpy(&output_filename[0][len_output_filename], "0_0_pkd");

        for (it_thread = 0; it_thread < nb_threads; it_thread++){

            sprintf(&(output_filename2[it_thread][len_output_filename]), "%d", it_thread);
            len_output_filename2 = strlen(output_filename2[it_thread]);
            output_filename2[it_thread][len_output_filename2] = '_';

            for (it_bft_thread = 0; it_bft_thread < nb_bfts_on_disk[it_thread]; it_bft_thread++){

                if (it_thread || it_bft_thread){
                    sprintf(&(output_filename2[it_thread][len_output_filename2+1]), "%d", it_bft_thread);

                    read_cut_BFT_Root(output_filename2[it_thread], output_filename3, cut_lvl, true);

                    printf("%s - %s\n", output_filename[0], output_filename2[it_thread]);

                    merging_BFT(output_filename[0], output_filename3, output_filename[0], cut_lvl, true);
                }
            }

            gettimeofday(&tval_after, NULL);

            time_spent(&tval_before, &tval_after, &tval_result);
            printf("Total elapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

            printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);
            printf("Current memory: %llu mb\n", ((unsigned long long int)getCurrentRSS())/1024);
        }

        free(nb_bfts_on_disk);
    }
    else{
        #pragma omp parallel                                                                                                            \
        shared(output_filename, output_filename2, len_output_filename, nb_bfts_on_disk, tval_before, tval_after, tval_result)           \
        private(it_thread, it_bft_thread, nb_merging, steps, len_output_filename2)
        {
            #pragma omp for
            for (it_thread = 0; it_thread < nb_threads; it_thread++){

                steps = 2;
                nb_merging = nb_bfts_on_disk[it_thread] - 1;

                sprintf(&(output_filename[it_thread][len_output_filename]), "%d", it_thread);
                sprintf(&(output_filename2[it_thread][len_output_filename]), "%d", it_thread);

                len_output_filename2 = strlen(output_filename[it_thread]);

                output_filename[it_thread][len_output_filename2] = '_';
                output_filename2[it_thread][len_output_filename2] = '_';

                while (nb_merging > 0){

                    for (it_bft_thread = 0; it_bft_thread + steps/2 < nb_bfts_on_disk[it_thread]; it_bft_thread += steps){

                        sprintf(&(output_filename[it_thread][len_output_filename2+1]), "%d", it_bft_thread);
                        sprintf(&(output_filename2[it_thread][len_output_filename2+1]), "%d", it_bft_thread + steps/2);

                        printf("%s - %s\n", output_filename[it_thread], output_filename2[it_thread]);

                        merging_BFT(output_filename[it_thread], output_filename2[it_thread], output_filename[it_thread], cut_lvl, false);

                        nb_merging--;
                    }

                    steps *= 2;
                }

                gettimeofday(&tval_after, NULL);

                time_spent(&tval_before, &tval_after, &tval_result);
                printf("Total elapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

                printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);
                printf("Current memory: %llu mb\n", ((unsigned long long int)getCurrentRSS())/1024);
            }
        }

        free(nb_bfts_on_disk);

        steps = 2;
        nb_merging = nb_threads - 1;

        while (nb_merging > 0){

            #pragma omp parallel shared(nb_merging, steps, output_filename, output_filename2) private(it_thread, len_output_filename2)
            {
                #pragma omp for
                for (it_thread = 0; it_thread < nb_threads; it_thread += steps){

                    if (it_thread + steps/2 < nb_threads){

                        sprintf(&(output_filename[it_thread][len_output_filename]), "%d", it_thread);
                        sprintf(&(output_filename2[it_thread][len_output_filename]), "%d", it_thread + steps/2);

                        len_output_filename2 = strlen(output_filename[it_thread]);

                        strcpy(&output_filename[it_thread][len_output_filename2], "_0");
                        strcpy(&output_filename2[it_thread][len_output_filename2], "_0");

                        merging_BFT(output_filename[it_thread], output_filename2[it_thread], output_filename[it_thread], cut_lvl, false);

                        nb_merging--;
                    }
                }
            }

            steps *= 2;

            gettimeofday(&tval_after, NULL);

            time_spent(&tval_before, &tval_after, &tval_result);
            printf("Total elapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

            printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);
            printf("Current memory: %llu mb\n", ((unsigned long long int)getCurrentRSS())/1024);
        }
    }

    return;
}*/
