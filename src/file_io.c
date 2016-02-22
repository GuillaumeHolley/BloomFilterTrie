#include "./../lib/file_io.h"

extern void time_spent(struct timeval *start_time, struct timeval *end_time, struct timeval *resulting_time);

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
void insert_Genomes_from_KmerFiles(BFT_Root* root, int nb_files, char** filenames, int binary_files){

    ASSERT_NULL_PTR(root,"insert_Genomes_from_KmerFiles()")
    ASSERT_NULL_PTR(filenames,"insert_Genomes_from_KmerFiles()")

    struct timeval tval_before, tval_after, tval_last, tval_result;
    gettimeofday(&tval_before, NULL);
    tval_last = tval_before;

    FILE* file;

    Pvoid_t PJArray = (PWord_t)NULL;
    Word_t Rc_word;

    annotation_array_elem* comp_set_colors_tmp = NULL;

    uint8_t* array_kmers = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(array_kmers,"insert_Genomes_from_KmerFiles()")

    char* line = calloc(100, sizeof(char));
    ASSERT_NULL_PTR(line,"insert_Genomes_from_KmerFiles()")

    char* str_tmp;

    int i = 0;
    int j = 0;
    int k = 0;
    int size_id_genome = 0;
    int lvl_root = (root->k / NB_CHAR_SUF_PREF) - 1;
    int nb_bytes_kmer = CEIL(root->k*2, SIZE_BITS_UINT_8T);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer;

    int length_comp_set_colors_tmp = 0;

    int len_longest_annot;

    size_t return_fread;

    uint64_t kmers_read;

    for (i = 0; i < nb_files; i++){ //For each file in input

        kmers_read = 0;
        k = 0;
        j = 0;

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

        if (root->treshold_compression != 0){
            if ((root->nb_genomes-1 > 5) && ((root->nb_genomes-1)%root->treshold_compression == 0)){

                memory_Used* mem = printMemoryUsedFromNode(&(root->node), lvl_root, root->k, root->info_per_lvl);
                len_longest_annot = MAX(mem->size_biggest_annot+1, getMaxSize_annotation_array_elem(root->comp_set_colors));
                free(mem);

                load_annotation_from_Node(&(root->node), lvl_root, root->k, len_longest_annot, root->info_per_lvl,
                                          &PJArray, root->comp_set_colors, root->ann_inf);

                comp_set_colors_tmp = root->comp_set_colors;
                length_comp_set_colors_tmp = root->length_comp_set_colors;

                root->comp_set_colors = sort_annotations(&PJArray, &(root->length_comp_set_colors), len_longest_annot);

                if (root->comp_set_colors != NULL){
                    compress_annotation_from_Node(&(root->node), lvl_root, root->k, root->info_per_lvl, &PJArray, comp_set_colors_tmp,
                                                  root->ann_inf);

                    free_annotation_array_elem(comp_set_colors_tmp, length_comp_set_colors_tmp);
                }

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

                free_annotation_array_elem(comp_set_colors_tmp, length_comp_set_colors_tmp);

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
