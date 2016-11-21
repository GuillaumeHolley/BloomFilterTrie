#include "compression.h"

uint32_t compress_FASTx_files(char* filename_read_1, char* filename_read_2, bool pair_ended,
                             int size_seed, int size_minimizer, int size_kmer, uint32_t prev_nb_part,
                             char* output_prefix, bool compress_shift, bool recycle_paths,
                             BFT_Root* root_no_iupac, BFT_Root* root_iupac){

    uint32_t nb_parts;

    int len_output_prefix = strlen(output_prefix);

    char* output_prefix_buffer = malloc((len_output_prefix + 50) * sizeof(char));
    ASSERT_NULL_PTR(output_prefix_buffer, "min_simReads_max_simkmers() 2\n")

    strcpy(output_prefix_buffer, output_prefix);

    prepare_shuffling_dictionary();

    printf("\n[BFTC step 1/3: Maximizing k-mer similarity]\n\n");

    min_simReads_max_simkmers2(filename_read_1, filename_read_2, pair_ended,
                              output_prefix, size_kmer, size_seed, size_minimizer,
                              compress_shift, root_no_iupac, root_iupac);

    strcpy(&output_prefix_buffer[len_output_prefix], "_span_sup_reads");
    if (remove(output_prefix_buffer)) printf("Warning: could not delete temporary file.\n");

    printf("[BFTC step 2/3: Minimizing partition integers]\n\n");

    if (prev_nb_part) prev_nb_part++;

    nb_parts = compressKmers_from_KmerFiles_bis3(output_prefix, size_seed, size_kmer, prev_nb_part);

    printf("[BFTC step 3/3: G-DBG encoding]\n\n");

    insert_kmers_partitions2(output_prefix, size_seed, size_kmer, recycle_paths, root_no_iupac, root_iupac);

    strcpy(&output_prefix_buffer[len_output_prefix], "_kmers_tmp");
    if (remove(output_prefix_buffer)) printf("Warning: could not delete temporary file.\n");

    strcpy(&output_prefix_buffer[len_output_prefix], "_parts_count_tmp");
    if (remove(output_prefix_buffer))  printf("Warning: could not delete temporary file.\n");

    free(output_prefix_buffer);

    return nb_parts;
}

void decompress_FASTx_file(char* output_prefix, bool pair_ended, int size_kmer, int size_seed){

    prepare_shuffling_dictionary();

    int len_output_prefix = strlen(output_prefix);

    char* output_prefix_buffer = malloc((len_output_prefix + 50) * sizeof(char));
    ASSERT_NULL_PTR(output_prefix_buffer, "min_simReads_max_simkmers() 2\n")

    strcpy(output_prefix_buffer, output_prefix);

    decompress2(output_prefix_buffer, pair_ended, size_kmer, size_seed);
    extract_reads(output_prefix_buffer, pair_ended);

    strcpy(&output_prefix_buffer[len_output_prefix], ".super_reads");
    if (remove(output_prefix_buffer)) printf("Warning: could not delete temporary file.\n");

    //strcpy(&output_prefix_buffer[len_output_prefix], ".reads");
    //if (remove(output_prefix_buffer)) printf("Warning: could not delete temporary file.\n");

    if (pair_ended){
        output_prefix_buffer[len_output_prefix] = '\0';
        reorder_paired_reads(output_prefix_buffer);
    }

    free(output_prefix_buffer);

    return;
}

size_t insert_kmer_into_bf_from_graph(BFT_kmer* kmer, BFT* graph, va_list args){

    bool graph_is_iupac = (bool)va_arg(args, int);

    BloomFilter* bloom_filter = va_arg(args, BloomFilter*);

    uint64_t hash_tmp;

    if (graph_is_iupac) kmer_iupac_comp_to_ascii(kmer->kmer_comp, graph->k/2, kmer->kmer);

    const uint64_t hash_v = XXH64(kmer->kmer, strlen(kmer->kmer), bloom_filter->seed_hash1);
    const uint64_t hash_v2 = XXH64(kmer->kmer, strlen(kmer->kmer), bloom_filter->seed_hash2);

    uint8_t* curr_block = &(bloom_filter->bf[(hash_v % bloom_filter->nb_blocks) * bloom_filter->size_block_bytes]);
    __builtin_prefetch(curr_block, 0, 3);

    for (int i = 1; i <= bloom_filter->nb_hash; i++){
        hash_tmp = (hash_v + i * hash_v2) % bloom_filter->size_block_bits;
        curr_block[hash_tmp/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash_tmp%SIZE_BITS_UINT_8T];
    }

    return 1;
}

size_t insert_kmer_into_bf(BloomFilter* bloom_filter, char* kmer, int length_kmer){

    uint64_t hash_tmp;

    const uint64_t hash_v = XXH64(kmer, length_kmer, bloom_filter->seed_hash1);
    const uint64_t hash_v2 = XXH64(kmer, length_kmer, bloom_filter->seed_hash2);

    uint8_t* curr_block = &(bloom_filter->bf[(hash_v % bloom_filter->nb_blocks) * bloom_filter->size_block_bytes]);
    __builtin_prefetch(curr_block, 0, 3);

    for (int i = 1; i <= bloom_filter->nb_hash; i++){
        hash_tmp = (hash_v + i * hash_v2) % bloom_filter->size_block_bits;
        curr_block[hash_tmp/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash_tmp%SIZE_BITS_UINT_8T];
    }

    return 1;
}

BloomFilter* create_BloomFilter(uint64_t nb_elem, double fp_rate_max){

    BloomFilter* bloom_filter = malloc(sizeof(BloomFilter));
    ASSERT_NULL_PTR(bloom_filter, "create_BloomFilter()\n")

    bloom_filter->size_block_bytes = SIZE_BLOCK_BYTES;
    bloom_filter->size_block_bits = bloom_filter->size_block_bytes * SIZE_BITS_UINT_8T;

    bloom_filter->fp_rate_bf = fp_rate_max;

    bloom_filter->nb_bits_bf = ceil((nb_elem * log(bloom_filter->fp_rate_bf)) / log(1.0 / (pow(2.0, log(2.0)))));

    /*max_memory_in_mb *= 8388608;

    if (bloom_filter->nb_bits_bf * 1.3 > max_memory_in_mb){
        while ((bloom_filter->nb_bits_bf * 1.3 > max_memory_in_mb) && (bloom_filter->fp_rate_bf > 0.001)){
            bloom_filter->fp_rate_bf *= 1.5;
            bloom_filter->nb_bits_bf = ceil((nb_elem * log(bloom_filter->fp_rate_bf)) / log(1.0 / (pow(2.0, log(2.0)))));
        }
    }
    else{
        while ((bloom_filter->nb_bits_bf * 1.3 < max_memory_in_mb) && (bloom_filter->fp_rate_bf > 0.001)){
            bloom_filter->fp_rate_bf /= 1.5;
            bloom_filter->nb_bits_bf = ceil((nb_elem * log(bloom_filter->fp_rate_bf)) / log(1.0 / (pow(2.0, log(2.0)))));
        }
    }*/

    bloom_filter->nb_hash = round(log(2.0) * bloom_filter->nb_bits_bf / nb_elem);
    bloom_filter->nb_bits_bf *= 1.3;

    bloom_filter->nb_blocks = CEIL(bloom_filter->nb_bits_bf, bloom_filter->size_block_bits);
    bloom_filter->nb_bytes_bf = bloom_filter->nb_blocks * bloom_filter->size_block_bytes;

    bloom_filter->seed_hash1 = rand();
    while (bloom_filter->seed_hash1 == (bloom_filter->seed_hash2 = rand())){};

    bloom_filter->bf = calloc(bloom_filter->nb_bytes_bf, sizeof(uint8_t));
    ASSERT_NULL_PTR(bloom_filter->bf,"create_BloomFilter()\n")

    return bloom_filter;
}

void free_BloomFilter(BloomFilter* bloom_filter){

    ASSERT_NULL_PTR(bloom_filter,"free_BloomFilter()\n")

    free(bloom_filter->bf);
    free(bloom_filter);

    return;
}

int64_t binning_reads(char* filename_reads_mate_1, char* filename_reads_mate_2, bool pair_ended, char* prefix_bin_name, int size_minimizer){

    //printf("binning_reads()\n");

    kseq_t* seq_read;

    int i, j;
    int length_read;
    int fp_bin;
    int fp_read;
    int pos_curr_minimzer;

    int minimizer_shift = 0;
    int length_buffer_read_rev = SIZE_BUFFER;
    int length_prefix_bin_name = strlen(prefix_bin_name);
    int size_minimizer_bytes = CEIL(size_minimizer * 2, SIZE_BITS_UINT_8T);

    uint64_t minimizer;

    int64_t nb_reads = 0;
    int64_t nb_reads_return = 0;

    bool min_is_rev;
    bool pair_ended_tmp = pair_ended;

    char saved_char;
    char nl = '\n';

    char* buffer;

    char* curr_minimizer;

    char buffer_int[16] = {0};
    char* info_read_non_rev = " 0";
    char* info_read_rev = " 1";

    uint8_t* id_minimizer = malloc(size_minimizer_bytes * sizeof(uint8_t));
    ASSERT_NULL_PTR(id_minimizer, "binning_reads()\n")

    char* bin_name = malloc((length_prefix_bin_name + size_minimizer + 100) * sizeof(char));
    ASSERT_NULL_PTR(bin_name, "binning_reads()\n")

    char* buffer_read_rev = malloc(length_buffer_read_rev * sizeof(char));
    ASSERT_NULL_PTR(buffer_read_rev, "binning_reads()\n")

    strcpy(bin_name, prefix_bin_name);

    if (size_minimizer % 4) minimizer_shift = SIZE_BITS_UINT_8T - ((size_minimizer % 4) * 2);

    fp_read = open(filename_reads_mate_1, O_RDONLY);
    seq_read = kseq_init(fp_read);

    START_BINNING:

    while (kseq_read(seq_read) > -1){

        buffer = seq_read->seq.s;
        length_read = strlen(seq_read->seq.s);

        curr_minimizer = NULL;
        pos_curr_minimzer = -1;

        min_is_rev = false;

        //What if the read AND its rc contains only minimizers with IUPAC characters ?
        for (i = 0; i <= length_read - size_minimizer; i++){

            saved_char = seq_read->seq.s[i + size_minimizer];
            seq_read->seq.s[i + size_minimizer] = '\0';

            if (is_substring_IUPAC(&(seq_read->seq.s[i])) == 0){

                if ((strstr(&(seq_read->seq.s[i]), "AAA") == NULL) && (strstr(&(seq_read->seq.s[i]), "CCC") == NULL)
                    && (strstr(&(seq_read->seq.s[i]), "GGG") == NULL) && (strstr(&(seq_read->seq.s[i]), "TTT") == NULL)){

                        pos_curr_minimzer = i;
                        curr_minimizer = &(seq_read->seq.s[i]);
                        seq_read->seq.s[i + size_minimizer] = saved_char;
                        break;
                }
            }

            seq_read->seq.s[i + size_minimizer] = saved_char;
        }

        //What if the read AND its rc contains only minimizers with IUPAC characters ?
        for (i++; i <= length_read - size_minimizer; i++){

            saved_char = seq_read->seq.s[i + size_minimizer];
            seq_read->seq.s[i + size_minimizer] = '\0';

            if (is_substring_IUPAC(&(seq_read->seq.s[i])) == 0){

                if (memcmp(&(seq_read->seq.s[i]), curr_minimizer, size_minimizer * sizeof(char)) < 0){

                    if ((strstr(&(seq_read->seq.s[i]), "AAA") == NULL) && (strstr(&(seq_read->seq.s[i]), "CCC") == NULL)
                        && (strstr(&(seq_read->seq.s[i]), "GGG") == NULL) && (strstr(&(seq_read->seq.s[i]), "TTT") == NULL)){

                        curr_minimizer = &(seq_read->seq.s[i]);
                        pos_curr_minimzer = i;
                    }
                }
            }

            seq_read->seq.s[i + size_minimizer] = saved_char;
        }

        if (length_read + 1 > length_buffer_read_rev){

            length_buffer_read_rev = length_read + 1;

            buffer_read_rev = realloc(buffer_read_rev, length_buffer_read_rev * sizeof(char));
            ASSERT_NULL_PTR(buffer_read_rev, "binning_reads()\n")
        }

        reverse_complement(seq_read->seq.s, buffer_read_rev, length_read);

        if (pos_curr_minimzer == -1){

            for (i = 0; i <= length_read - size_minimizer; i++){

                saved_char = buffer_read_rev[i + size_minimizer];
                buffer_read_rev[i + size_minimizer] = '\0';

                if (is_substring_IUPAC(&buffer_read_rev[i]) == 0){

                    if ((strstr(&buffer_read_rev[i], "AAA") == NULL) && (strstr(&buffer_read_rev[i], "CCC") == NULL)
                        && (strstr(&buffer_read_rev[i], "GGG") == NULL) && (strstr(&buffer_read_rev[i], "TTT") == NULL)){

                        pos_curr_minimzer = i;
                        curr_minimizer = &buffer_read_rev[i];
                        buffer_read_rev[i + size_minimizer] = saved_char;
                        min_is_rev = true;
                        break;
                    }
                }

                buffer_read_rev[i + size_minimizer] = saved_char;
            }

            i++;
        }
        else i = 0;

        for (; i <= length_read - size_minimizer; i++){

            saved_char = buffer_read_rev[i + size_minimizer];
            buffer_read_rev[i + size_minimizer] = '\0';

            if (is_substring_IUPAC(&buffer_read_rev[i]) == 0){

                if (memcmp(&buffer_read_rev[i], curr_minimizer, size_minimizer * sizeof(char)) < 0){

                    if ((strstr(&buffer_read_rev[i], "AAA") == NULL) && (strstr(&buffer_read_rev[i], "CCC") == NULL)
                        && (strstr(&buffer_read_rev[i], "GGG") == NULL) && (strstr(&buffer_read_rev[i], "TTT") == NULL)){

                        curr_minimizer = &buffer_read_rev[i];
                        pos_curr_minimzer = i;
                        buffer = buffer_read_rev;
                        min_is_rev = true;
                    }
                }
            }

            buffer_read_rev[i + size_minimizer] = saved_char;
        }

        if (pos_curr_minimzer == -1) minimizer = 0;
        else {

            memset(id_minimizer, 0, size_minimizer_bytes * sizeof(uint8_t));
            parseKmerCount(curr_minimizer, size_minimizer, id_minimizer, 0);

            for (j = 0, minimizer = 0; j < size_minimizer_bytes; j++)
                minimizer = (minimizer << SIZE_BITS_UINT_8T) | reverse_word_8(id_minimizer[j]);

            minimizer >>= minimizer_shift;
        }

        sprintf(buffer_int, "_%" PRIu64, minimizer);
        strcpy(&bin_name[length_prefix_bin_name], buffer_int);

        fp_bin = open(bin_name, O_WRONLY | O_APPEND | O_CREAT, S_IRUSR | S_IWUSR);
        if (fp_bin == -1) ERROR("binning_reads(): Could not open/create the file corrsponding to a bin")

        sprintf(buffer_int, "%d ", pos_curr_minimzer);
        write(fp_bin, buffer_int, strlen(buffer_int));

        write(fp_bin, buffer, length_read);

        if (min_is_rev) write(fp_bin, info_read_rev, 2);
        else write(fp_bin, info_read_non_rev, 2);

        if (pair_ended){
            sprintf(buffer_int, " %" PRId64, nb_reads);
            write(fp_bin, buffer_int, strlen(buffer_int));
        }

        write(fp_bin, &nl, 1);

        close(fp_bin);

        nb_reads++;
        //if (nb_reads % 1000000 == 0) printf("%" PRIu64 " reads binned\n", nb_reads);
    }

    kseq_destroy(seq_read);
    close(fp_read);

    if (pair_ended_tmp){

        nb_reads_return = nb_reads;

        fp_read = open(filename_reads_mate_2, O_RDONLY);
        seq_read = kseq_init(fp_read);

        pair_ended_tmp = false;

        goto START_BINNING;
    }

    free(bin_name);
    free(id_minimizer);
    free(buffer_read_rev);

    if (pair_ended){

        if (nb_reads_return != nb_reads / 2)
            ERROR("binning_reads(): Number of reads in files containing pair-ended reads is not the same.\n")

        return nb_reads_return;
    }

    return nb_reads;
}

int memcmp_mismatch(const uint8_t * ptr1, const uint8_t * ptr2, size_t num, int nb_max_mismatch, Mismatch* mis, int id_read){

    ASSERT_NULL_PTR(ptr1, "memcmp_mismatch()\n")
    ASSERT_NULL_PTR(ptr2, "memcmp_mismatch()\n")

    int nb_mismatch = 0;

    for (int it = 0; it < num; it++){

        if (ptr1[it] != ptr2[it]){

            if ((nb_mismatch + 1) <= nb_max_mismatch){

                mis[nb_mismatch].position = it;
                mis[nb_mismatch].mismatch_char = ptr2[it];
                mis[nb_mismatch].id_read = id_read;

                nb_mismatch++;
            }
            else return -1;
        }
    }

    return nb_mismatch;
}

void create_super_reads(char* prefix_bin_name, bool pair_ended, int size_minimizer){

    //printf("create_super_reads()\n");

    Pvoid_t reads = (Pvoid_t) NULL;
    PWord_t reads_value;
    Word_t  reads_flag;

    Word_t pos_minimizer;
    Word_t pos_minimizer_tmp;

    FILE* file_bin;

    List* list_reads;
    List* list_pos_length_occ;

    ListNode* cur = NULL;

    Read_count* read_count;
    Mismatch* mis;
    Pos_length_occ* pos_length_occ;

    int rev;
    int it_list;
    int it_pos_length_occ;
    int list_count;
    int length_read;
    int length_delimiter;
    int length_super_read;
    int reads_delete;
    int length_after_minimizer;

    int it;
    int match;
    int len_array_read;
    int len_read_count;
    int nb_mis_max = 5;

    int size_array_reads = 1;
    int length_prefix_bin_name = strlen(prefix_bin_name);

    size_t size_buffer_read = SIZE_BUFFER;

    uint64_t it_buffer;

    char buffer_int[64] = {0};

    char nl = '\n';
    char space = ' ';

    char* delimiter;
    char* read;

    char* bin_name = malloc((length_prefix_bin_name + size_minimizer + 100) * sizeof(char));
    ASSERT_NULL_PTR(bin_name, "binning_reads()\n")

    char* buffer_read = malloc(size_buffer_read * sizeof(char));
    ASSERT_NULL_PTR(buffer_read, "binning_reads()\n")

    char** array_reads = malloc(size_array_reads * sizeof(char*));
    ASSERT_NULL_PTR(array_reads, "binning_reads()\n")

    mis = malloc(nb_mis_max * sizeof(Mismatch));
    ASSERT_NULL_PTR(mis, "binning_reads()\n")

    strcpy(bin_name, prefix_bin_name);

    for (it_buffer = 0; it_buffer < pow(4, size_minimizer); it_buffer++){

        sprintf(buffer_int, "_%" PRIu64, it_buffer);
        strcpy(&bin_name[length_prefix_bin_name], buffer_int);

        if ((file_bin = fopen(bin_name, "r")) != NULL){

            //Index reads in memory (position_minimizer, list of reads)
            while (getline(&buffer_read, &size_buffer_read, file_bin) != -1){

                pos_minimizer = custom_atoi(buffer_read, &delimiter, ' ');
                delimiter++;

                length_delimiter = strlen(delimiter);

                if (delimiter[length_delimiter-1] == '\n'){
                    length_delimiter--;
                    delimiter[length_delimiter] = '\0';
                }

                read = malloc((length_delimiter + 1) * sizeof(char));
                ASSERT_NULL_PTR(read, "binning_reads()\n")

                strcpy(read, delimiter);

                JLG(reads_value, reads, pos_minimizer);

                if (reads_value == NULL){

                    list_reads = List_create();

                    JLI(reads_value, reads, pos_minimizer);
                    *reads_value = (Word_t) list_reads;
                }
                else list_reads = (List*) *reads_value;

                List_push(list_reads, read);
            }

            fclose(file_bin);

            //Delete duplicates
            pos_minimizer_tmp = 0;
            JLF(reads_value, reads, pos_minimizer_tmp);

            while (reads_value != NULL){

                list_reads = (List*) *reads_value;
                list_count = list_reads->count;

                if (list_count > size_array_reads){

                    size_array_reads = list_count;

                    array_reads = realloc(array_reads, size_array_reads * sizeof(char*));
                    ASSERT_NULL_PTR(array_reads, "binning_reads()\n")
                }

                for (cur = list_reads->first, it_list = 0; cur != NULL; cur = cur->next, it_list++)
                    array_reads[it_list] = (char*) cur->value;

                qsort(array_reads, list_count, sizeof(char*), string_cmp);

                List_destroy(list_reads);
                list_reads = List_create();

                read_count = create_Read_count();

                read_count->read = array_reads[0];
                read_count->occ_count = 1;

                List_push(list_reads, read_count);

                if (!pair_ended){

                    for (it_list = 1; it_list < list_count; it_list++){

                        len_array_read = strlen(array_reads[it_list]);
                        len_read_count = strlen(read_count->read);

                        if (len_array_read == len_read_count){

                            match = memcmp_mismatch((uint8_t*) read_count->read, (uint8_t*) array_reads[it_list], len_read_count - 1,
                                                    nb_mis_max, mis, read_count->occ_count);
                        }
                        else if ((match = strcmp(read_count->read, array_reads[it_list])) != 0) match = -1;

                        if ((match >= 0) && (match <= nb_mis_max) && (read_count->read[len_read_count-1] == array_reads[it_list][len_array_read-1])){

                            free(array_reads[it_list]);

                            read_count->occ_count += 1;

                            if (match){

                                read_count->list_mismatches = realloc(read_count->list_mismatches, (read_count->nb_mismatches + match) * sizeof(Mismatch));
                                ASSERT_NULL_PTR(read_count->list_mismatches, "binning_reads()\n")

                                memcpy(&(read_count->list_mismatches[read_count->nb_mismatches]), mis, match * sizeof(Mismatch));

                                read_count->nb_mismatches += match;
                            }
                        }
                        else{

                            read_count = create_Read_count();

                            read_count->read = array_reads[it_list];
                            read_count->occ_count = 1;

                            List_push(list_reads, read_count);
                        }
                    }
                }
                else{

                    read_count->id_occ = malloc(sizeof(int64_t));
                    ASSERT_NULL_PTR(read_count->id_occ, "binning_reads()\n")

                    delimiter = strchr(strchr(read_count->read, ' ') + 1, ' ');
                    delimiter[0] = '\0';
                    delimiter++;

                    read_count->id_occ[0] = (int64_t) strtoll(delimiter, &delimiter, 10);

                    for (it_list = 1; it_list < list_count; it_list++){

                        delimiter = strchr(strchr(array_reads[it_list], ' ') + 1, ' ');
                        delimiter[0] = '\0';
                        delimiter++;

                        len_array_read = strlen(array_reads[it_list]);
                        len_read_count = strlen(read_count->read);

                        if (len_array_read == len_read_count){

                            match = memcmp_mismatch((uint8_t*) read_count->read, (uint8_t*) array_reads[it_list], len_read_count - 1,
                                                    nb_mis_max, mis, read_count->occ_count);
                        }
                        else if ((match = strcmp(read_count->read, array_reads[it_list])) != 0) match = -1;

                        if ((match >= 0) && (match <= nb_mis_max) && (read_count->read[len_read_count-1] == array_reads[it_list][len_array_read-1])){

                            read_count->occ_count += 1;

                            read_count->id_occ = realloc(read_count->id_occ, read_count->occ_count * sizeof(int64_t));
                            ASSERT_NULL_PTR(read_count->id_occ, "binning_reads()\n")

                            read_count->id_occ[read_count->occ_count - 1] = (int64_t) strtoll(delimiter, &delimiter, 10);

                            if (match){

                                read_count->list_mismatches = realloc(read_count->list_mismatches, (read_count->nb_mismatches + match) * sizeof(Mismatch));
                                ASSERT_NULL_PTR(read_count->list_mismatches, "binning_reads()\n")

                                memcpy(&(read_count->list_mismatches[read_count->nb_mismatches]), mis, match * sizeof(Mismatch));

                                read_count->nb_mismatches += match;
                            }

                            free(array_reads[it_list]);
                        }
                        else {

                            read_count = create_Read_count();

                            read_count->read = array_reads[it_list];
                            read_count->occ_count = 1;

                            read_count->id_occ = malloc(sizeof(int64_t));
                            ASSERT_NULL_PTR(read_count->id_occ, "binning_reads()\n")

                            read_count->id_occ[0] = (int64_t) strtoll(delimiter, &delimiter, 10);

                            List_push(list_reads, read_count);
                        }
                    }
                }

                //Replace the list of duplicated reads by the list of unique reads
                *reads_value = (Word_t) list_reads;

                JLN(reads_value, reads, pos_minimizer_tmp);
            }

            //Create super reads
            file_bin = fopen(bin_name, "w");
            ASSERT_NULL_PTR(file_bin, "binning_reads()\n")

            pos_minimizer = INT_MAX;
            JLL(reads_value, reads, pos_minimizer);

            while (reads_value != NULL){

                list_reads = (List*) *reads_value;
                read_count = (Read_count*) list_reads->first->value;

                delimiter = strchr(read_count->read, ' ');
                delimiter[0] = '\0';
                delimiter++;

                rev = custom_atoi(delimiter, &delimiter, '\0');

                length_read = strlen(read_count->read);
                length_super_read = length_read;

                strcpy(buffer_read, read_count->read);

                list_pos_length_occ = List_create();

                pos_length_occ = create_set_Pos_length_occ_from_Read_count(0, length_read, rev, read_count);
                List_push(list_pos_length_occ, pos_length_occ);

                free(read_count->read);
                free(read_count);

                List_remove(list_reads, list_reads->first);

                if (list_reads->count == 0){
                    List_destroy(list_reads);
                    JLD(reads_delete, reads, pos_minimizer);
                }

                pos_minimizer_tmp = pos_minimizer;

                sprintf(buffer_int, "%d ", (int) pos_minimizer);
                fwrite(buffer_int, sizeof(char), strlen(buffer_int), file_bin);
                fwrite(buffer_read, sizeof(char), length_read, file_bin);

                length_after_minimizer = length_read - pos_minimizer - size_minimizer;

                JLL(reads_value, reads, pos_minimizer_tmp);

                while (reads_value != NULL){ //Iterate over reads of index

                    list_reads = (List*) *reads_value;

                    for (cur = list_reads->first; cur != NULL; cur = cur->next){

                        read_count = (Read_count*) cur->value;

                        if (strcmp(buffer_read, read_count->read)){ // If reads are not the same

                            match = memcmp_mismatch((uint8_t*) &buffer_read[pos_minimizer - pos_minimizer_tmp], (uint8_t*) read_count->read,
                                                    (pos_minimizer_tmp + length_read - pos_minimizer) * sizeof(char), nb_mis_max, mis, read_count->occ_count);

                            if ((match >= 0) && (match <= nb_mis_max)){

                                delimiter = strchr(read_count->read, ' ');
                                delimiter[0] = '\0';
                                delimiter++;

                                rev = custom_atoi(delimiter, &delimiter, '\0');

                                memmove(buffer_read, &buffer_read[pos_minimizer - pos_minimizer_tmp], (pos_minimizer_tmp + length_read - pos_minimizer) * sizeof(char));
                                strcpy(&buffer_read[pos_minimizer_tmp + length_read - pos_minimizer], &(read_count->read[pos_minimizer_tmp + length_read - pos_minimizer]));

                                length_read = strlen(buffer_read);
                                pos_minimizer = pos_minimizer_tmp;

                                fwrite(&buffer_read[pos_minimizer + size_minimizer + length_after_minimizer], sizeof(char),
                                       length_read - pos_minimizer - size_minimizer - length_after_minimizer, file_bin);

                                length_super_read += length_read - pos_minimizer - size_minimizer - length_after_minimizer;
                                length_after_minimizer = length_read - pos_minimizer - size_minimizer;

                                pos_length_occ = create_set_Pos_length_occ_from_Read_count(length_super_read - length_read, length_read, rev, read_count);

                                if (match){

                                    pos_length_occ->list_mismatches = realloc(pos_length_occ->list_mismatches, (pos_length_occ->nb_mismatches + match) * sizeof(Mismatch));
                                    ASSERT_NULL_PTR(pos_length_occ->list_mismatches, "create_super_reads()\n")

                                    memcpy(&(pos_length_occ->list_mismatches[pos_length_occ->nb_mismatches]), mis, match * sizeof(Mismatch));

                                    pos_length_occ->nb_mismatches += match;
                                }

                                List_push(list_pos_length_occ, pos_length_occ);

                                free(read_count->read);
                                free(read_count);

                                List_remove(list_reads, cur);

                                if (list_reads->count == 0){
                                    List_destroy(list_reads);
                                    JLD(reads_delete, reads, pos_minimizer_tmp);
                                }

                                break;
                            }
                        }
                    }

                    JLP(reads_value, reads, pos_minimizer_tmp);
                }

                //Dump pos_length_occ

                fwrite(&space, sizeof(char), 1, file_bin);

                for (cur = list_pos_length_occ->first; cur != NULL; cur = cur->next){

                    pos_length_occ = (Pos_length_occ*) cur->value;

                    qsort(pos_length_occ->list_mismatches, pos_length_occ->nb_mismatches, sizeof(Mismatch), mismatch_cmp);

                    sprintf(buffer_int, "%d,", pos_length_occ->position);
                    fwrite(buffer_int, sizeof(char), strlen(buffer_int), file_bin);

                    sprintf(buffer_int, "%d,", pos_length_occ->length);
                    fwrite(buffer_int, sizeof(char), strlen(buffer_int), file_bin);

                    sprintf(buffer_int, "%d,", pos_length_occ->rev);
                    fwrite(buffer_int, sizeof(char), strlen(buffer_int), file_bin);

                    sprintf(buffer_int, "%d,", pos_length_occ->occ_count);
                    fwrite(buffer_int, sizeof(char), strlen(buffer_int), file_bin);

                    sprintf(buffer_int, "%d,", pos_length_occ->nb_mismatches);
                    fwrite(buffer_int, sizeof(char), strlen(buffer_int), file_bin);

                    if (pair_ended){

                        for (it_pos_length_occ = 0; it_pos_length_occ < pos_length_occ->occ_count; it_pos_length_occ++){

                            sprintf(buffer_int, "%" PRId64 ",", pos_length_occ->id_occ[it_pos_length_occ]);
                            fwrite(buffer_int, sizeof(char), strlen(buffer_int), file_bin);

                        }
                    }

                    for (it = 0; it < pos_length_occ->nb_mismatches; it++){

                        sprintf(buffer_int, "%d,%c,%d,",
                                pos_length_occ->list_mismatches[it].position,
                                pos_length_occ->list_mismatches[it].mismatch_char,
                                pos_length_occ->list_mismatches[it].id_read);

                        fwrite(buffer_int, sizeof(char), strlen(buffer_int), file_bin);
                    }

                    free(pos_length_occ->id_occ);
                    free(pos_length_occ->list_mismatches);
                    free(pos_length_occ);
                }

                fwrite(&nl, sizeof(char), 1, file_bin);

                List_destroy(list_pos_length_occ);

                pos_minimizer = INT_MAX;
                JLL(reads_value, reads, pos_minimizer);
            }

            fclose(file_bin);

            pos_minimizer_tmp = 0;
            JLF(reads_value, reads, pos_minimizer_tmp);

            while (reads_value != NULL){

                list_reads = (List*) *reads_value;

                for (cur = list_reads->first; cur != NULL; cur = cur->next){

                    read_count = (Read_count*) cur->value;

                    free(read_count->read);
                    free(read_count->list_mismatches);
                    free(read_count->id_occ);
                    free(read_count);
                }

                List_destroy(list_reads);

                JLN(reads_value, reads, pos_minimizer_tmp);
            }

            JLFA(reads_flag, reads);
        }
    }

    free(mis);
    free(array_reads);
    free(bin_name);
    free(buffer_read);

    return;
}

uint64_t create_spanning_super_reads(char* prefix_bin_name, bool pair_ended, int64_t nb_reads_per_file, int size_minimizer){

    //printf("create_spanning_super_reads()\n");

    FILE* file_bin;
    FILE* file_new_bin;
    FILE* file_span_super_reads;

    FILE* file_span_super_reads_mate;
    FILE* file_span_super_reads_mate_info;
    FILE* file_span_super_reads_pos_mate;
    FILE* file_span_super_reads_pos_mate_tmp;
    FILE* file_span_super_reads_pos;
    FILE* file_span_super_reads_length;
    FILE* file_span_super_reads_occ;
    FILE* file_span_super_reads_rev;
    FILE* file_span_super_reads_pos_mis;
    FILE* file_span_super_reads_char_mis;
    FILE* file_span_super_reads_id_mis;
    FILE* file_span_super_reads_nb_mis;

    List* list_pos_length_occ;

    ListNode* cur = NULL;

    Pos_length_occ* pos_length_occ;

    int match;
    int length_spanning;
    int list_count;

    int it_id_occ;
    int length_read;
    int length_read_tmp;
    int length_after_minimizer;
    int pos_minimizer;
    int pos_minimizer_tmp;
    int diff_pos_minimizers;
    int pos_delimiter_read_tmp;

    int length_read_info = 0;
    int minimizer_shift = 0;
    int length_prefix_bin_name = strlen(prefix_bin_name);
    int size_minimizer_bytes = CEIL(size_minimizer * 2, SIZE_BITS_UINT_8T);

    long int pos_start_line;
    long int pos_curr_line;

    long int pos_start_line_tmp;
    long int pos_curr_line_tmp;

    bool continue_loop;
    bool first_pos_length_occ;
    bool is_rev;

    size_t size_buffer_read = SIZE_BUFFER;
    size_t size_buffer_read_tmp = SIZE_BUFFER;
    size_t size_buffer_read_rev = SIZE_BUFFER;

    size_t size_array_pos_length_occ = 1;

    uint32_t i;
    uint32_t position, position_tmp;
    uint32_t length;
    uint32_t occ_count;
    uint32_t rev;

    int64_t nb_reads = 0;
    int64_t nb_reads_tmp = 0;

    int64_t id_occ;

    uint64_t opened_file_minimizer;
    uint64_t minimizer;
    uint64_t it_buffer;
    uint64_t tot_length_seq;

    uint8_t read_info = 0;
    uint8_t mate_info = 0;

    char saved_char;

    char nl = '\n';
    char eol = '\0';
    char read_used = 'X';
    char buffer_int[16] = {0};

    char* delimiter_read;
    char* delimiter_read_tmp;

    char* delimiter_meta;
    char* delimiter_meta_tmp;

    char* read;

    char* curr_minimizer;

    uint8_t* id_minimizer = malloc(size_minimizer_bytes * sizeof(uint8_t));
    ASSERT_NULL_PTR(id_minimizer, "create_spanning_super_reads()\n")

    char* bin_name = malloc((length_prefix_bin_name + size_minimizer + 100) * sizeof(char));
    ASSERT_NULL_PTR(bin_name, "create_spanning_super_reads()\n")

    char* buffer_read = malloc(size_buffer_read * sizeof(char));
    ASSERT_NULL_PTR(buffer_read, "create_spanning_super_reads()\n")

    char* buffer_read_tmp = malloc(size_buffer_read_tmp * sizeof(char));
    ASSERT_NULL_PTR(buffer_read_tmp, "create_spanning_super_reads()\n")

    char* buffer_read_rev = malloc(size_buffer_read_rev * sizeof(char));
    ASSERT_NULL_PTR(buffer_read_rev, "create_spanning_super_reads()\n")

    Pos_length_occ** array_pos_length_occ = malloc(size_array_pos_length_occ * sizeof(Pos_length_occ*));
    ASSERT_NULL_PTR(array_pos_length_occ, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads");

    file_span_super_reads = fopen(bin_name, "w");
    ASSERT_NULL_PTR(file_span_super_reads, "create_spanning_super_reads()\n")

    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_pos");

    file_span_super_reads_pos = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_pos, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_length");

    file_span_super_reads_length = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_length, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_occ");

    file_span_super_reads_occ = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_occ, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_rev");

    file_span_super_reads_rev = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_rev, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_pos_mis");

    file_span_super_reads_pos_mis = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_pos_mis, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_char_mis");

    file_span_super_reads_char_mis = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_char_mis, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_id_mis");

    file_span_super_reads_id_mis = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_id_mis, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_nb_mis");

    file_span_super_reads_nb_mis = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_nb_mis, "create_spanning_super_reads()\n")

    if (pair_ended){

        strcpy(bin_name, prefix_bin_name);
        strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_tmp");

        file_span_super_reads_pos_mate = fopen(bin_name, "wb");
        ASSERT_NULL_PTR(file_span_super_reads_pos_mate, "create_spanning_super_reads()\n")
    }

    if (size_minimizer % 4) minimizer_shift = SIZE_BITS_UINT_8T - ((size_minimizer % 4) * 2);

    for (it_buffer = 0; it_buffer < pow(4, size_minimizer); it_buffer++){

        sprintf(buffer_int, "_%" PRIu64, it_buffer);
        strcpy(&bin_name[length_prefix_bin_name], buffer_int);

        opened_file_minimizer = it_buffer;

        if ((file_bin = fopen(bin_name, "r+")) != NULL){

            pos_start_line = ftell(file_bin);

            while (getline(&buffer_read, &size_buffer_read, file_bin) != -1){

                if (buffer_read[0] != read_used){

                    length_spanning = 0;

                    pos_minimizer = custom_atoi(buffer_read, &delimiter_read, ' ') + size_minimizer;
                    delimiter_read[0] = eol;
                    delimiter_read++;

                    delimiter_meta = strchr(delimiter_read, ' ');
                    delimiter_meta[0] = eol;
                    delimiter_meta++;

                    if (delimiter_meta[strlen(delimiter_meta)-1] == nl) delimiter_meta[strlen(delimiter_meta)-1] = eol;

                    read = delimiter_read;
                    length_read = strlen(delimiter_read);

                    is_rev = false;

                    curr_minimizer = &read[pos_minimizer];

                    //What if the read AND its rc contains only minimizers with IUPAC characters ?
                    for (i = custom_atoi(buffer_read, &delimiter_read_tmp, eol) + 1; i <= length_read - size_minimizer; i++){

                        saved_char = read[i + size_minimizer];
                        read[i + size_minimizer] = eol;

                        if (is_substring_IUPAC(&read[i]) == 0){

                            if (memcmp(&read[i], curr_minimizer, size_minimizer * sizeof(char)) < 0){

                                if ((strstr(&read[i], "AAA") == NULL) && (strstr(&read[i], "CCC") == NULL)
                                    && (strstr(&read[i], "GGG") == NULL) && (strstr(&read[i], "TTT") == NULL)){

                                    curr_minimizer = &read[i];
                                    pos_minimizer = i;
                                }
                            }
                        }

                        read[i + size_minimizer] = saved_char;
                    }

                    if (size_buffer_read_rev < size_buffer_read){

                        size_buffer_read_rev = size_buffer_read;

                        buffer_read_rev = realloc(buffer_read_rev, size_buffer_read_rev * sizeof(char));
                        ASSERT_NULL_PTR(buffer_read_rev, "create_spanning_super_reads()\n")
                    }

                    //if (i > length_read - size_minimizer) continue;

                    reverse_complement(delimiter_read, buffer_read_rev, length_read);

                    for (i = 0; i < length_read - size_minimizer; i++){

                        saved_char = buffer_read_rev[i + size_minimizer];
                        buffer_read_rev[i + size_minimizer] = eol;

                        if (is_substring_IUPAC(&buffer_read_rev[i]) == 0){

                            if (memcmp(&buffer_read_rev[i], curr_minimizer, size_minimizer * sizeof(char)) < 0){

                                if ((strstr(&buffer_read_rev[i], "AAA") == NULL) && (strstr(&buffer_read_rev[i], "CCC") == NULL)
                                    && (strstr(&buffer_read_rev[i], "GGG") == NULL) && (strstr(&buffer_read_rev[i], "TTT") == NULL)){

                                    curr_minimizer = &buffer_read_rev[i];
                                    pos_minimizer = i;
                                    read = buffer_read_rev;
                                    is_rev = true;
                                }
                            }
                        }

                        buffer_read_rev[i + size_minimizer] = saved_char;
                    }

                    memset(id_minimizer, 0, size_minimizer_bytes * sizeof(uint8_t));
                    parseKmerCount(curr_minimizer, size_minimizer, id_minimizer, 0);

                    for (i = 0, minimizer = 0; i < size_minimizer_bytes; i++)
                        minimizer = (minimizer << SIZE_BITS_UINT_8T) | reverse_word_8(id_minimizer[i]);

                    minimizer >>= minimizer_shift;

                    if (minimizer != opened_file_minimizer){

                        sprintf(buffer_int, "_%" PRIu64, minimizer);
                        strcpy(&bin_name[length_prefix_bin_name], buffer_int);

                        while ((file_new_bin = fopen(bin_name, "r+")) != NULL){

                            continue_loop = false;

                            pos_start_line_tmp = ftell(file_new_bin);

                            while (getline(&buffer_read_tmp, &size_buffer_read_tmp, file_new_bin) != -1){

                                if (buffer_read_tmp[0] != read_used){

                                    pos_minimizer_tmp = custom_atoi(buffer_read_tmp, &delimiter_read_tmp, ' ');
                                    delimiter_read_tmp[0] = eol;
                                    delimiter_read_tmp++;

                                    delimiter_meta_tmp = strchr(delimiter_read_tmp, ' ');
                                    delimiter_meta_tmp[0] = eol;
                                    delimiter_meta_tmp++;

                                    if (delimiter_meta_tmp[strlen(delimiter_meta_tmp)-1] == nl) delimiter_meta_tmp[strlen(delimiter_meta_tmp)-1] = eol;

                                    length_read_tmp = strlen(delimiter_read_tmp);

                                    if (pos_minimizer_tmp <= pos_minimizer){

                                        //Compute matches
                                        match = memcmp(&read[pos_minimizer], &delimiter_read_tmp[pos_minimizer_tmp],
                                                          MIN(length_read - pos_minimizer, length_read_tmp - pos_minimizer_tmp) * sizeof(char));

                                        if (match == 0){

                                            match = memcmp(&read[pos_minimizer - pos_minimizer_tmp], delimiter_read_tmp, pos_minimizer_tmp * sizeof(char));

                                            if (match == 0){

                                                length_spanning++;

                                                if (length_spanning == 1){

                                                    if (!is_rev) memmove(buffer_read, read, (length_read + 1) * sizeof(char));
                                                    else strcpy(buffer_read, read);

                                                    length_after_minimizer = length_read - pos_minimizer - size_minimizer;

                                                    pos_curr_line = ftell(file_bin);

                                                    fseek(file_bin, pos_start_line, SEEK_SET);
                                                    fwrite(&read_used, sizeof(char), 1, file_bin);
                                                    fseek(file_bin, pos_curr_line, SEEK_SET);

                                                    list_pos_length_occ = List_create();

                                                    while (delimiter_meta[0] != eol){

                                                        pos_length_occ = parse_pos_length_occ(&delimiter_meta, pair_ended);
                                                        List_push(list_pos_length_occ, pos_length_occ);
                                                    }

                                                    for (cur = list_pos_length_occ->first; (cur != NULL) && is_rev; cur = cur->next){

                                                        pos_length_occ = cur->value;
                                                        pos_length_occ->rev = pos_length_occ->rev ? 0 : 1;
                                                        pos_length_occ->position = length_read - pos_length_occ->position - pos_length_occ->length;

                                                        for (it_id_occ = 0; it_id_occ < pos_length_occ->nb_mismatches; it_id_occ++){

                                                            pos_length_occ->list_mismatches[it_id_occ].position = pos_length_occ->length - pos_length_occ->list_mismatches[it_id_occ].position - 1;
                                                            pos_length_occ->list_mismatches[it_id_occ].mismatch_char = reverse_complement_char(pos_length_occ->list_mismatches[it_id_occ].mismatch_char);

                                                            if (pos_length_occ->list_mismatches[it_id_occ].id_read < pos_length_occ->occ_count){
                                                                pos_length_occ->list_mismatches[it_id_occ].id_read = (pos_length_occ->occ_count - 1)
                                                                                                                    - pos_length_occ->list_mismatches[it_id_occ].id_read;
                                                            }
                                                        }
                                                    }
                                                }

                                                pos_curr_line_tmp = ftell(file_new_bin);

                                                fseek(file_new_bin, pos_start_line_tmp, SEEK_SET);
                                                fwrite(&read_used, sizeof(char), 1, file_new_bin);
                                                fseek(file_new_bin, pos_curr_line_tmp, SEEK_SET);

                                                diff_pos_minimizers = pos_minimizer - pos_minimizer_tmp;

                                                while (delimiter_meta_tmp[0] != eol){

                                                    pos_length_occ = parse_pos_length_occ(&delimiter_meta_tmp, pair_ended);
                                                    pos_length_occ->position += diff_pos_minimizers;

                                                    List_push(list_pos_length_occ, pos_length_occ);
                                                }

                                                pos_delimiter_read_tmp = pos_minimizer_tmp + size_minimizer + length_after_minimizer;

                                                if (pos_delimiter_read_tmp < length_read_tmp){

                                                    if (size_buffer_read < length_read + length_read_tmp - pos_delimiter_read_tmp + 1){

                                                        size_buffer_read = length_read + length_read_tmp - pos_delimiter_read_tmp + 1;

                                                        buffer_read = realloc(buffer_read, size_buffer_read * sizeof(char));
                                                        ASSERT_NULL_PTR(buffer_read, "create_spanning_super_reads()\n")
                                                    }

                                                    strcpy(&buffer_read[length_read], &delimiter_read_tmp[pos_delimiter_read_tmp]);
                                                    length_read += length_read_tmp - pos_delimiter_read_tmp;
                                                }

                                                is_rev = false;
                                                read = buffer_read;
                                                pos_minimizer_tmp = pos_minimizer + 1;
                                                curr_minimizer = &buffer_read[pos_minimizer];

                                                //What if the read AND its rc contains only minimizers with IUPAC characters ?
                                                for (i = pos_minimizer_tmp; i <= length_read - size_minimizer; i++){

                                                    saved_char = buffer_read[i + size_minimizer];
                                                    buffer_read[i + size_minimizer] = eol;

                                                    if (is_substring_IUPAC(&buffer_read[i]) == 0){

                                                        if (memcmp(&buffer_read[i], curr_minimizer, size_minimizer * sizeof(char)) < 0){

                                                            if ((strstr(&buffer_read[i], "AAA") == NULL) && (strstr(&buffer_read[i], "CCC") == NULL)
                                                                && (strstr(&buffer_read[i], "GGG") == NULL) && (strstr(&buffer_read[i], "TTT") == NULL)){

                                                                curr_minimizer = &buffer_read[i];
                                                                pos_minimizer = i;
                                                            }
                                                        }
                                                    }

                                                    buffer_read[i + size_minimizer] = saved_char;
                                                }

                                                if (size_buffer_read_rev < size_buffer_read){

                                                    size_buffer_read_rev = size_buffer_read;

                                                    buffer_read_rev = realloc(buffer_read_rev, size_buffer_read_rev * sizeof(char));
                                                    ASSERT_NULL_PTR(buffer_read_rev, "create_spanning_super_reads()\n")
                                                }

                                                //if (i > length_read - size_minimizer) continue;

                                                reverse_complement(buffer_read, buffer_read_rev, length_read);

                                                for (i = 0; i <= length_read - size_minimizer; i++){

                                                    saved_char = buffer_read_rev[i + size_minimizer];
                                                    buffer_read_rev[i + size_minimizer] = eol;

                                                    if (is_substring_IUPAC(&buffer_read_rev[i]) == 0){

                                                        if (memcmp(&buffer_read_rev[i], curr_minimizer, size_minimizer * sizeof(char)) < 0){

                                                            if ((strstr(&buffer_read_rev[i], "AAA") == NULL) && (strstr(&buffer_read_rev[i], "CCC") == NULL)
                                                                && (strstr(&buffer_read_rev[i], "GGG") == NULL) && (strstr(&buffer_read_rev[i], "TTT") == NULL)){

                                                                curr_minimizer = &buffer_read_rev[i];
                                                                pos_minimizer = i;
                                                                read = buffer_read_rev;
                                                                is_rev = true;
                                                            }
                                                        }
                                                    }

                                                    buffer_read_rev[i + size_minimizer] = saved_char;
                                                }

                                                if (is_rev){

                                                    strcpy(buffer_read, read);

                                                    for (cur = list_pos_length_occ->first; cur != NULL; cur = cur->next){

                                                        pos_length_occ = cur->value;
                                                        pos_length_occ->rev = (pos_length_occ->rev ? 0 : 1);
                                                        pos_length_occ->position = length_read - pos_length_occ->position - pos_length_occ->length;

                                                        for (it_id_occ = 0; it_id_occ < pos_length_occ->nb_mismatches; it_id_occ++){

                                                            pos_length_occ->list_mismatches[it_id_occ].position = pos_length_occ->length - pos_length_occ->list_mismatches[it_id_occ].position - 1;
                                                            pos_length_occ->list_mismatches[it_id_occ].mismatch_char = reverse_complement_char(pos_length_occ->list_mismatches[it_id_occ].mismatch_char);

                                                            if (pos_length_occ->list_mismatches[it_id_occ].id_read < pos_length_occ->occ_count)
                                                                pos_length_occ->list_mismatches[it_id_occ].id_read = (pos_length_occ->occ_count - 1) - pos_length_occ->list_mismatches[it_id_occ].id_read;
                                                        }
                                                    }
                                                }

                                                length_after_minimizer = length_read - pos_minimizer - size_minimizer;

                                                memset(id_minimizer, 0, size_minimizer_bytes * sizeof(uint8_t));
                                                parseKmerCount(curr_minimizer, size_minimizer, id_minimizer, 0);

                                                for (i = 0, minimizer = 0; i < size_minimizer_bytes; i++)
                                                    minimizer = (minimizer << SIZE_BITS_UINT_8T) | reverse_word_8(id_minimizer[i]);

                                                minimizer >>= minimizer_shift;

                                                sprintf(buffer_int, "_%" PRIu64, minimizer);
                                                strcpy(&bin_name[length_prefix_bin_name], buffer_int);

                                                if (minimizer != opened_file_minimizer) continue_loop = true;

                                                break;
                                            }
                                        }
                                    }
                                }

                                pos_start_line_tmp = ftell(file_new_bin);
                            }

                            fclose(file_new_bin);

                            if (!continue_loop) break;
                        }

                        if (length_spanning){

                            fwrite(read, sizeof(char), length_read, file_span_super_reads);
                            fwrite(&nl, sizeof(char), 1, file_span_super_reads);

                            //Sort read position
                            list_count = list_pos_length_occ->count;

                            if (list_count > size_array_pos_length_occ){

                                size_array_pos_length_occ = list_count;

                                array_pos_length_occ = realloc(array_pos_length_occ, size_array_pos_length_occ * sizeof(Pos_length_occ*));
                                ASSERT_NULL_PTR(array_pos_length_occ, "binning_reads()\n")
                            }

                            for (cur = list_pos_length_occ->first, i = 0; cur != NULL; cur = cur->next, i++)
                                array_pos_length_occ[i] = (Pos_length_occ*) cur->value;

                            List_destroy(list_pos_length_occ);

                            qsort(array_pos_length_occ, list_count, sizeof(Pos_length_occ*), pos_length_occ_cmp);

                            first_pos_length_occ = true;
                            position = 0;

                            for (i = 0; i < list_count; i++){

                                pos_length_occ = array_pos_length_occ[i];

                                position = (uint32_t) pos_length_occ->position - position;
                                occ_count = (uint32_t) pos_length_occ->occ_count;
                                rev = (uint32_t) pos_length_occ->rev;
                                length = ((uint32_t) pos_length_occ->length) << 1;

                                if (first_pos_length_occ){
                                    length |= 0x1;
                                    first_pos_length_occ = false;
                                }

                                read_info = (read_info << 1) | ((uint8_t) rev);
                                length_read_info++;

                                if (length_read_info == SIZE_BITS_UINT_8T){
                                    fwrite(&read_info, sizeof(uint8_t), 1, file_span_super_reads_rev);
                                    read_info = 0;
                                    length_read_info = 0;
                                }

                                fwrite(&position, sizeof(uint32_t), 1, file_span_super_reads_pos);
                                fwrite(&length, sizeof(uint32_t), 1, file_span_super_reads_length);
                                fwrite(&occ_count, sizeof(uint32_t), 1, file_span_super_reads_occ);

                                if (pair_ended){

                                    for (it_id_occ = 0; it_id_occ < occ_count; it_id_occ++){

                                        fseek(file_span_super_reads_pos_mate, pos_length_occ->id_occ[it_id_occ] * sizeof(int64_t), SEEK_SET);
                                        fwrite(&nb_reads, sizeof(int64_t), 1, file_span_super_reads_pos_mate);

                                        nb_reads++;
                                    }

                                    free(pos_length_occ->id_occ);
                                }

                                rev = (uint32_t) pos_length_occ->nb_mismatches;
                                fwrite(&rev, sizeof(uint32_t), 1, file_span_super_reads_nb_mis);

                                qsort(pos_length_occ->list_mismatches, pos_length_occ->nb_mismatches, sizeof(Mismatch), mismatch_cmp);

                                for (it_id_occ = 0, position = 0; it_id_occ < pos_length_occ->nb_mismatches; it_id_occ++){

                                    position = ((uint32_t) pos_length_occ->list_mismatches[it_id_occ].position) - position;
                                    length = (uint32_t) pos_length_occ->list_mismatches[it_id_occ].id_read;

                                    fwrite(&position, sizeof(uint32_t), 1, file_span_super_reads_pos_mis);
                                    fwrite(&length, sizeof(uint32_t), 1, file_span_super_reads_id_mis);
                                    fwrite(&(pos_length_occ->list_mismatches[it_id_occ].mismatch_char), sizeof(char), 1, file_span_super_reads_char_mis);

                                    position = (uint32_t) pos_length_occ->list_mismatches[it_id_occ].position;
                                }

                                position = (uint32_t) pos_length_occ->position;

                                free(pos_length_occ->list_mismatches);
                                free(pos_length_occ);
                            }
                        }
                    }
                }

                pos_start_line = ftell(file_bin);
            }

            fclose(file_bin);
        }
    }

    for (it_buffer = 0; it_buffer < pow(4, size_minimizer); it_buffer++){

        sprintf(buffer_int, "_%" PRIu64, it_buffer);
        strcpy(&bin_name[length_prefix_bin_name], buffer_int);

        if ((file_bin = fopen(bin_name, "r")) != NULL){

            while (getline(&buffer_read, &size_buffer_read, file_bin) != -1){

                if (buffer_read[0] != read_used){

                    first_pos_length_occ = true;
                    position = 0;
                    delimiter_read = strchr(buffer_read, ' ') + 1;

                    delimiter_meta = strchr(delimiter_read, ' ');
                    delimiter_meta[0] = eol;
                    delimiter_meta++;

                    if (delimiter_meta[strlen(delimiter_meta)-1] == nl) delimiter_meta[strlen(delimiter_meta)-1] = eol;

                    fwrite(delimiter_read, sizeof(char), strlen(delimiter_read), file_span_super_reads);
                    fwrite(&nl, sizeof(char), 1, file_span_super_reads);

                    while (delimiter_meta[0] != eol){

                        pos_length_occ = parse_pos_length_occ(&delimiter_meta, pair_ended);

                        position_tmp = (uint32_t) pos_length_occ->position;
                        position = position_tmp - position;

                        length = ((uint32_t) pos_length_occ->length) << 1;
                        rev = (uint32_t) pos_length_occ->rev;
                        occ_count = (uint32_t) pos_length_occ->occ_count;

                        if (first_pos_length_occ){
                            length |= 0x1;
                            first_pos_length_occ = false;
                        }

                        read_info = (read_info << 1) | ((uint8_t) rev);
                        length_read_info++;

                        if (length_read_info == SIZE_BITS_UINT_8T){
                            fwrite(&read_info, sizeof(uint8_t), 1, file_span_super_reads_rev);
                            read_info = 0;
                            length_read_info = 0;
                        }

                        fwrite(&position, sizeof(uint32_t), 1, file_span_super_reads_pos);
                        fwrite(&length, sizeof(uint32_t), 1, file_span_super_reads_length);
                        fwrite(&occ_count, sizeof(uint32_t), 1, file_span_super_reads_occ);

                        if (pair_ended){

                            for (it_id_occ = 0; it_id_occ < occ_count; it_id_occ++){

                                fseek(file_span_super_reads_pos_mate, pos_length_occ->id_occ[it_id_occ] * sizeof(int64_t), SEEK_SET);
                                fwrite(&nb_reads, sizeof(int64_t), 1, file_span_super_reads_pos_mate);

                                nb_reads++;
                            }

                            free(pos_length_occ->id_occ);
                        }

                        rev = (uint32_t) pos_length_occ->nb_mismatches;
                        fwrite(&rev, sizeof(uint32_t), 1, file_span_super_reads_nb_mis);

                        for (it_id_occ = 0, position = 0; it_id_occ < pos_length_occ->nb_mismatches; it_id_occ++){

                            position = ((uint32_t) pos_length_occ->list_mismatches[it_id_occ].position) - position;
                            length = (uint32_t) pos_length_occ->list_mismatches[it_id_occ].id_read;

                            fwrite(&position, sizeof(uint32_t), 1, file_span_super_reads_pos_mis);
                            fwrite(&length, sizeof(uint32_t), 1, file_span_super_reads_id_mis);
                            fwrite(&(pos_length_occ->list_mismatches[it_id_occ].mismatch_char), sizeof(char), 1, file_span_super_reads_char_mis);

                            position = (uint32_t) pos_length_occ->list_mismatches[it_id_occ].position;
                        }

                        position = position_tmp;

                        free(pos_length_occ->list_mismatches);
                        free(pos_length_occ);
                    }
                }
            }

            fclose(file_bin);

            if (remove(bin_name)) printf("Warning: could not delete temporary file.\n");
        }
    }

    if (length_read_info){
        read_info <<= SIZE_BITS_UINT_8T - length_read_info;
        fwrite(&read_info, sizeof(uint8_t), 1, file_span_super_reads_rev);
    }

    fseek(file_span_super_reads, 0 - ((int) sizeof(char)), SEEK_CUR);
    fwrite(&eol, sizeof(char), 1, file_span_super_reads);

    tot_length_seq = ftell(file_span_super_reads);

    fclose(file_span_super_reads_pos);
    fclose(file_span_super_reads_length);
    fclose(file_span_super_reads_rev);
    fclose(file_span_super_reads_occ);
    fclose(file_span_super_reads_pos_mis);
    fclose(file_span_super_reads_char_mis);
    fclose(file_span_super_reads_id_mis);
    fclose(file_span_super_reads_nb_mis);

    fclose(file_span_super_reads);

    if (pair_ended){

        fclose(file_span_super_reads_pos_mate);

        strcpy(bin_name, prefix_bin_name);
        strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_mate");

        file_span_super_reads_mate = fopen(bin_name, "w+b");
        ASSERT_NULL_PTR(file_span_super_reads_mate, "create_spanning_super_reads()\n")

        strcpy(bin_name, prefix_bin_name);
        strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_mate_info");

        file_span_super_reads_mate_info = fopen(bin_name, "w+b");
        ASSERT_NULL_PTR(file_span_super_reads_mate_info, "create_spanning_super_reads()\n")

        strcpy(bin_name, prefix_bin_name);
        strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_tmp");

        file_span_super_reads_pos_mate = fopen(bin_name, "rb");
        ASSERT_NULL_PTR(file_span_super_reads_pos_mate, "create_spanning_super_reads()\n")

        file_span_super_reads_pos_mate_tmp = fopen(bin_name, "rb");
        ASSERT_NULL_PTR(file_span_super_reads_pos_mate_tmp, "create_spanning_super_reads()\n")

        fseek(file_span_super_reads_pos_mate_tmp, nb_reads_per_file * sizeof(int64_t), SEEK_SET);

        ftruncate(fileno(file_span_super_reads_mate_info), CEIL(nb_reads_per_file * 2, SIZE_BITS_UINT_8T));

        for (nb_reads = 0; nb_reads < nb_reads_per_file; nb_reads++){

            fread(&id_occ, sizeof(int64_t), 1, file_span_super_reads_pos_mate);
            fread(&nb_reads_tmp, sizeof(int64_t), 1, file_span_super_reads_pos_mate_tmp);

            fseek(file_span_super_reads_mate_info, (id_occ / SIZE_BITS_UINT_8T) * sizeof(uint8_t), SEEK_SET);
            fread(&mate_info, sizeof(uint8_t), 1, file_span_super_reads_mate_info);

            mate_info |= MASK_POWER_8[id_occ%SIZE_BITS_UINT_8T];

            fseek(file_span_super_reads_mate_info, 0 - sizeof(uint8_t), SEEK_CUR);
            fwrite(&mate_info, sizeof(uint8_t), 1, file_span_super_reads_mate_info);

            nb_reads_tmp -= id_occ;

            fseek(file_span_super_reads_mate, id_occ * sizeof(int64_t), SEEK_SET);
            fwrite(&nb_reads_tmp, sizeof(int64_t), 1, file_span_super_reads_mate);
        }

        fclose(file_span_super_reads_pos_mate);
        fclose(file_span_super_reads_pos_mate_tmp);

        rewind(file_span_super_reads_mate);
        rewind(file_span_super_reads_mate_info);

        id_occ = 0;
        nb_reads = 0;

        while (fread(&mate_info, sizeof(uint8_t), 1, file_span_super_reads_mate_info) == 1){

            for (i = 0; mate_info && (i < SIZE_BITS_UINT_8T); i++, mate_info >>= 1){

                if (mate_info & 0x1){

                    fseek(file_span_super_reads_mate, (nb_reads + i) * sizeof(int64_t), SEEK_SET);
                    fread(&nb_reads_tmp, sizeof(int64_t), 1, file_span_super_reads_mate);

                    fseek(file_span_super_reads_mate, id_occ * sizeof(int64_t), SEEK_SET);
                    fwrite(&nb_reads_tmp, sizeof(int64_t), 1, file_span_super_reads_mate);

                    id_occ++;
                }
            }

            nb_reads += SIZE_BITS_UINT_8T;
        }

        ftruncate(fileno(file_span_super_reads_mate), id_occ * sizeof(int64_t));

        fclose(file_span_super_reads_mate);
        fclose(file_span_super_reads_mate_info);

        if (remove(bin_name)) printf("Warning: could not delete temporary file.\n");
    }

    free(bin_name);
    free(id_minimizer);

    free(array_pos_length_occ);

    free(buffer_read);
    free(buffer_read_tmp);
    free(buffer_read_rev);

    return tot_length_seq;
}

uint64_t create_spanning_super_reads2(char* prefix_bin_name, bool pair_ended, int64_t nb_reads_per_file, int size_minimizer){

    //printf("create_spanning_super_reads()\n");

    FILE* file_bin;
    FILE* file_new_bin;
    FILE* file_span_super_reads;

    FILE* file_span_super_reads_mate;
    FILE* file_span_super_reads_mate_info;
    FILE* file_span_super_reads_pos_mate;
    FILE* file_span_super_reads_pos_mate_tmp;
    FILE* file_span_super_reads_pos;
    FILE* file_span_super_reads_length;
    FILE* file_span_super_reads_occ;
    FILE* file_span_super_reads_rev;
    FILE* file_span_super_reads_pos_mis;
    FILE* file_span_super_reads_char_mis;
    FILE* file_span_super_reads_id_mis;
    FILE* file_span_super_reads_nb_mis;

    List* list_pos_length_occ;

    ListNode* cur = NULL;

    Pos_length_occ* pos_length_occ;

    int match;
    int length_spanning;
    int list_count;

    int it_id_occ;
    int length_read;
    int length_read_tmp;
    int length_after_minimizer;
    int pos_minimizer;
    int pos_minimizer_tmp;
    int diff_pos_minimizers;
    int pos_delimiter_read_tmp;

    int nb_mis_max = 5;

    int length_read_info = 0;
    int minimizer_shift = 0;
    int length_prefix_bin_name = strlen(prefix_bin_name);
    int size_minimizer_bytes = CEIL(size_minimizer * 2, SIZE_BITS_UINT_8T);

    long int pos_start_line;
    long int pos_curr_line;

    long int pos_start_line_tmp;
    long int pos_curr_line_tmp;

    bool continue_loop;
    bool first_pos_length_occ;
    bool is_rev;

    size_t size_buffer_read = SIZE_BUFFER;
    size_t size_buffer_read_tmp = SIZE_BUFFER;
    size_t size_buffer_read_rev = SIZE_BUFFER;

    size_t size_array_pos_length_occ = 1;

    uint32_t i;
    uint32_t position, position_tmp;
    uint32_t length;
    uint32_t occ_count;
    uint32_t rev;

    int64_t nb_reads = 0;
    int64_t nb_reads_tmp = 0;

    int64_t id_occ;

    uint64_t opened_file_minimizer;
    uint64_t minimizer;
    uint64_t it_buffer;
    uint64_t tot_length_seq;

    uint8_t read_info = 0;
    uint8_t mate_info = 0;

    char saved_char;

    char nl = '\n';
    char eol = '\0';
    char read_used = 'X';
    char buffer_int[16] = {0};

    char* delimiter_read;
    char* delimiter_read_tmp;

    char* delimiter_meta;
    char* delimiter_meta_tmp;

    char* read;

    char* curr_minimizer;

    uint8_t* id_minimizer = malloc(size_minimizer_bytes * sizeof(uint8_t));
    ASSERT_NULL_PTR(id_minimizer, "create_spanning_super_reads()\n")

    char* bin_name = malloc((length_prefix_bin_name + size_minimizer + 100) * sizeof(char));
    ASSERT_NULL_PTR(bin_name, "create_spanning_super_reads()\n")

    char* buffer_read = malloc(size_buffer_read * sizeof(char));
    ASSERT_NULL_PTR(buffer_read, "create_spanning_super_reads()\n")

    char* buffer_read_tmp = malloc(size_buffer_read_tmp * sizeof(char));
    ASSERT_NULL_PTR(buffer_read_tmp, "create_spanning_super_reads()\n")

    char* buffer_read_rev = malloc(size_buffer_read_rev * sizeof(char));
    ASSERT_NULL_PTR(buffer_read_rev, "create_spanning_super_reads()\n")

    Pos_length_occ** array_pos_length_occ = malloc(size_array_pos_length_occ * sizeof(Pos_length_occ*));
    ASSERT_NULL_PTR(array_pos_length_occ, "create_spanning_super_reads()\n")

    Mismatch* mis = malloc(nb_mis_max * sizeof(Mismatch));
    ASSERT_NULL_PTR(mis, "binning_reads()\n")

    strcpy(bin_name, prefix_bin_name);

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads");

    file_span_super_reads = fopen(bin_name, "w");
    ASSERT_NULL_PTR(file_span_super_reads, "create_spanning_super_reads()\n")

    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_pos");

    file_span_super_reads_pos = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_pos, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_length");

    file_span_super_reads_length = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_length, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_occ");

    file_span_super_reads_occ = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_occ, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_rev");

    file_span_super_reads_rev = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_rev, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_pos_mis");

    file_span_super_reads_pos_mis = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_pos_mis, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_char_mis");

    file_span_super_reads_char_mis = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_char_mis, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_id_mis");

    file_span_super_reads_id_mis = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_id_mis, "create_spanning_super_reads()\n")

    strcpy(bin_name, prefix_bin_name);
    strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_nb_mis");

    file_span_super_reads_nb_mis = fopen(bin_name, "wb");
    ASSERT_NULL_PTR(file_span_super_reads_nb_mis, "create_spanning_super_reads()\n")

    if (pair_ended){

        strcpy(bin_name, prefix_bin_name);
        strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_tmp");

        file_span_super_reads_pos_mate = fopen(bin_name, "wb");
        ASSERT_NULL_PTR(file_span_super_reads_pos_mate, "create_spanning_super_reads()\n")
    }

    if (size_minimizer % 4) minimizer_shift = SIZE_BITS_UINT_8T - ((size_minimizer % 4) * 2);

    for (it_buffer = 0; it_buffer < pow(4, size_minimizer); it_buffer++){

        sprintf(buffer_int, "_%" PRIu64, it_buffer);
        strcpy(&bin_name[length_prefix_bin_name], buffer_int);

        opened_file_minimizer = it_buffer;

        if ((file_bin = fopen(bin_name, "r+")) != NULL){

            pos_start_line = ftell(file_bin);

            while (getline(&buffer_read, &size_buffer_read, file_bin) != -1){

                if (buffer_read[0] != read_used){

                    length_spanning = 0;

                    pos_minimizer = custom_atoi(buffer_read, &delimiter_read, ' ') + size_minimizer;
                    delimiter_read[0] = eol;
                    delimiter_read++;

                    delimiter_meta = strchr(delimiter_read, ' ');
                    delimiter_meta[0] = eol;
                    delimiter_meta++;

                    if (delimiter_meta[strlen(delimiter_meta)-1] == nl) delimiter_meta[strlen(delimiter_meta)-1] = eol;

                    read = delimiter_read;
                    length_read = strlen(delimiter_read);

                    is_rev = false;

                    curr_minimizer = &read[pos_minimizer];

                    //What if the read AND its rc contains only minimizers with IUPAC characters ?
                    for (i = custom_atoi(buffer_read, &delimiter_read_tmp, eol) + 1; i <= length_read - size_minimizer; i++){

                        saved_char = read[i + size_minimizer];
                        read[i + size_minimizer] = eol;

                        if (is_substring_IUPAC(&read[i]) == 0){

                            if (memcmp(&read[i], curr_minimizer, size_minimizer * sizeof(char)) < 0){

                                if ((strstr(&read[i], "AAA") == NULL) && (strstr(&read[i], "CCC") == NULL)
                                    && (strstr(&read[i], "GGG") == NULL) && (strstr(&read[i], "TTT") == NULL)){

                                    curr_minimizer = &read[i];
                                    pos_minimizer = i;
                                }
                            }
                        }

                        read[i + size_minimizer] = saved_char;
                    }

                    if (size_buffer_read_rev < size_buffer_read){

                        size_buffer_read_rev = size_buffer_read;

                        buffer_read_rev = realloc(buffer_read_rev, size_buffer_read_rev * sizeof(char));
                        ASSERT_NULL_PTR(buffer_read_rev, "create_spanning_super_reads()\n")
                    }

                    //if (i > length_read - size_minimizer) continue;

                    reverse_complement(delimiter_read, buffer_read_rev, length_read);

                    for (i = 0; i < length_read - size_minimizer; i++){

                        saved_char = buffer_read_rev[i + size_minimizer];
                        buffer_read_rev[i + size_minimizer] = eol;

                        if (is_substring_IUPAC(&buffer_read_rev[i]) == 0){

                            if (memcmp(&buffer_read_rev[i], curr_minimizer, size_minimizer * sizeof(char)) < 0){

                                if ((strstr(&buffer_read_rev[i], "AAA") == NULL) && (strstr(&buffer_read_rev[i], "CCC") == NULL)
                                    && (strstr(&buffer_read_rev[i], "GGG") == NULL) && (strstr(&buffer_read_rev[i], "TTT") == NULL)){

                                    curr_minimizer = &buffer_read_rev[i];
                                    pos_minimizer = i;
                                    read = buffer_read_rev;
                                    is_rev = true;
                                }
                            }
                        }

                        buffer_read_rev[i + size_minimizer] = saved_char;
                    }

                    memset(id_minimizer, 0, size_minimizer_bytes * sizeof(uint8_t));
                    parseKmerCount(curr_minimizer, size_minimizer, id_minimizer, 0);

                    for (i = 0, minimizer = 0; i < size_minimizer_bytes; i++)
                        minimizer = (minimizer << SIZE_BITS_UINT_8T) | reverse_word_8(id_minimizer[i]);

                    minimizer >>= minimizer_shift;

                    if (minimizer != opened_file_minimizer){

                        sprintf(buffer_int, "_%" PRIu64, minimizer);
                        strcpy(&bin_name[length_prefix_bin_name], buffer_int);

                        while ((file_new_bin = fopen(bin_name, "r+")) != NULL){

                            continue_loop = false;

                            pos_start_line_tmp = ftell(file_new_bin);

                            while (getline(&buffer_read_tmp, &size_buffer_read_tmp, file_new_bin) != -1){

                                if (buffer_read_tmp[0] != read_used){

                                    pos_minimizer_tmp = custom_atoi(buffer_read_tmp, &delimiter_read_tmp, ' ');
                                    delimiter_read_tmp[0] = eol;
                                    delimiter_read_tmp++;

                                    delimiter_meta_tmp = strchr(delimiter_read_tmp, ' ');
                                    delimiter_meta_tmp[0] = eol;
                                    delimiter_meta_tmp++;

                                    if (delimiter_meta_tmp[strlen(delimiter_meta_tmp)-1] == nl) delimiter_meta_tmp[strlen(delimiter_meta_tmp)-1] = eol;

                                    length_read_tmp = strlen(delimiter_read_tmp);

                                    if (pos_minimizer_tmp <= pos_minimizer){

                                        //Compute matches

                                        match = memcmp_mismatch((uint8_t*) &read[pos_minimizer - pos_minimizer_tmp], (uint8_t*) delimiter_read_tmp,
                                                                (pos_minimizer_tmp + MIN(length_read - pos_minimizer, length_read_tmp - pos_minimizer_tmp)) * sizeof(char),
                                                                nb_mis_max, mis, INT_MAX - 5);

                                        if ((match >= 0) && (match <= nb_mis_max)){

                                            length_spanning++;

                                            if (length_spanning == 1){

                                                if (!is_rev) memmove(buffer_read, read, (length_read + 1) * sizeof(char));
                                                else strcpy(buffer_read, read);

                                                length_after_minimizer = length_read - pos_minimizer - size_minimizer;

                                                pos_curr_line = ftell(file_bin);

                                                fseek(file_bin, pos_start_line, SEEK_SET);
                                                fwrite(&read_used, sizeof(char), 1, file_bin);
                                                fseek(file_bin, pos_curr_line, SEEK_SET);

                                                list_pos_length_occ = List_create();

                                                while (delimiter_meta[0] != eol){

                                                    pos_length_occ = parse_pos_length_occ(&delimiter_meta, pair_ended);
                                                    List_push(list_pos_length_occ, pos_length_occ);
                                                }

                                                for (cur = list_pos_length_occ->first; (cur != NULL) && is_rev; cur = cur->next){

                                                    pos_length_occ = cur->value;
                                                    pos_length_occ->rev = pos_length_occ->rev ? 0 : 1;
                                                    pos_length_occ->position = length_read - pos_length_occ->position - pos_length_occ->length;

                                                    for (it_id_occ = 0; it_id_occ < pos_length_occ->nb_mismatches; it_id_occ++){

                                                        pos_length_occ->list_mismatches[it_id_occ].position = pos_length_occ->length - pos_length_occ->list_mismatches[it_id_occ].position - 1;
                                                        pos_length_occ->list_mismatches[it_id_occ].mismatch_char = reverse_complement_char(pos_length_occ->list_mismatches[it_id_occ].mismatch_char);

                                                        if (pos_length_occ->list_mismatches[it_id_occ].id_read < pos_length_occ->occ_count){
                                                            pos_length_occ->list_mismatches[it_id_occ].id_read = (pos_length_occ->occ_count - 1)
                                                                                                                - pos_length_occ->list_mismatches[it_id_occ].id_read;
                                                        }
                                                    }
                                                }
                                            }

                                            pos_curr_line_tmp = ftell(file_new_bin);

                                            fseek(file_new_bin, pos_start_line_tmp, SEEK_SET);
                                            fwrite(&read_used, sizeof(char), 1, file_new_bin);
                                            fseek(file_new_bin, pos_curr_line_tmp, SEEK_SET);

                                            diff_pos_minimizers = pos_minimizer - pos_minimizer_tmp;

                                            while (delimiter_meta_tmp[0] != eol){

                                                pos_length_occ = parse_pos_length_occ(&delimiter_meta_tmp, pair_ended);

                                                for (it_id_occ = 0; it_id_occ < match; it_id_occ++){

                                                    if ((mis[it_id_occ].position >= pos_length_occ->position) && (mis[it_id_occ].position < pos_length_occ->position + pos_length_occ->length)){

                                                        pos_length_occ->list_mismatches = realloc(pos_length_occ->list_mismatches, (pos_length_occ->nb_mismatches + 1) * sizeof(Mismatch));
                                                        ASSERT_NULL_PTR(pos_length_occ->list_mismatches, "create_spanning_super_reads()\n")

                                                        memcpy(&(pos_length_occ->list_mismatches[pos_length_occ->nb_mismatches]), &mis[it_id_occ], sizeof(Mismatch));

                                                        pos_length_occ->list_mismatches[pos_length_occ->nb_mismatches].position -= diff_pos_minimizers;
                                                        pos_length_occ->list_mismatches[pos_length_occ->nb_mismatches].id_read = pos_length_occ->occ_count + 1;

                                                        pos_length_occ->nb_mismatches += 1;
                                                    }
                                                }

                                                pos_length_occ->position += diff_pos_minimizers;

                                                List_push(list_pos_length_occ, pos_length_occ);
                                            }

                                            pos_delimiter_read_tmp = pos_minimizer_tmp + size_minimizer + length_after_minimizer;

                                            if (pos_delimiter_read_tmp < length_read_tmp){

                                                if (size_buffer_read < length_read + length_read_tmp - pos_delimiter_read_tmp + 1){

                                                    size_buffer_read = length_read + length_read_tmp - pos_delimiter_read_tmp + 1;

                                                    buffer_read = realloc(buffer_read, size_buffer_read * sizeof(char));
                                                    ASSERT_NULL_PTR(buffer_read, "create_spanning_super_reads()\n")
                                                }

                                                strcpy(&buffer_read[length_read], &delimiter_read_tmp[pos_delimiter_read_tmp]);
                                                length_read += length_read_tmp - pos_delimiter_read_tmp;
                                            }

                                            is_rev = false;
                                            read = buffer_read;
                                            pos_minimizer_tmp = pos_minimizer + 1;
                                            curr_minimizer = &buffer_read[pos_minimizer];

                                            //What if the read AND its rc contains only minimizers with IUPAC characters ?
                                            for (i = pos_minimizer_tmp; i <= length_read - size_minimizer; i++){

                                                saved_char = buffer_read[i + size_minimizer];
                                                buffer_read[i + size_minimizer] = eol;

                                                if (is_substring_IUPAC(&buffer_read[i]) == 0){

                                                    if (memcmp(&buffer_read[i], curr_minimizer, size_minimizer * sizeof(char)) < 0){

                                                        if ((strstr(&buffer_read[i], "AAA") == NULL) && (strstr(&buffer_read[i], "CCC") == NULL)
                                                            && (strstr(&buffer_read[i], "GGG") == NULL) && (strstr(&buffer_read[i], "TTT") == NULL)){

                                                            curr_minimizer = &buffer_read[i];
                                                            pos_minimizer = i;
                                                        }
                                                    }
                                                }

                                                buffer_read[i + size_minimizer] = saved_char;
                                            }

                                            if (size_buffer_read_rev < size_buffer_read){

                                                size_buffer_read_rev = size_buffer_read;

                                                buffer_read_rev = realloc(buffer_read_rev, size_buffer_read_rev * sizeof(char));
                                                ASSERT_NULL_PTR(buffer_read_rev, "create_spanning_super_reads()\n")
                                            }

                                            //if (i > length_read - size_minimizer) continue;

                                            reverse_complement(buffer_read, buffer_read_rev, length_read);

                                            for (i = 0; i <= length_read - size_minimizer; i++){

                                                saved_char = buffer_read_rev[i + size_minimizer];
                                                buffer_read_rev[i + size_minimizer] = eol;

                                                if (is_substring_IUPAC(&buffer_read_rev[i]) == 0){

                                                    if (memcmp(&buffer_read_rev[i], curr_minimizer, size_minimizer * sizeof(char)) < 0){

                                                        if ((strstr(&buffer_read_rev[i], "AAA") == NULL) && (strstr(&buffer_read_rev[i], "CCC") == NULL)
                                                            && (strstr(&buffer_read_rev[i], "GGG") == NULL) && (strstr(&buffer_read_rev[i], "TTT") == NULL)){

                                                            curr_minimizer = &buffer_read_rev[i];
                                                            pos_minimizer = i;
                                                            read = buffer_read_rev;
                                                            is_rev = true;
                                                        }
                                                    }
                                                }

                                                buffer_read_rev[i + size_minimizer] = saved_char;
                                            }

                                            if (is_rev){

                                                strcpy(buffer_read, read);

                                                for (cur = list_pos_length_occ->first; cur != NULL; cur = cur->next){

                                                    pos_length_occ = cur->value;
                                                    pos_length_occ->rev = (pos_length_occ->rev ? 0 : 1);
                                                    pos_length_occ->position = length_read - pos_length_occ->position - pos_length_occ->length;

                                                    for (it_id_occ = 0; it_id_occ < pos_length_occ->nb_mismatches; it_id_occ++){

                                                        pos_length_occ->list_mismatches[it_id_occ].position = pos_length_occ->length - pos_length_occ->list_mismatches[it_id_occ].position - 1;
                                                        pos_length_occ->list_mismatches[it_id_occ].mismatch_char = reverse_complement_char(pos_length_occ->list_mismatches[it_id_occ].mismatch_char);

                                                        if (pos_length_occ->list_mismatches[it_id_occ].id_read < pos_length_occ->occ_count)
                                                            pos_length_occ->list_mismatches[it_id_occ].id_read = (pos_length_occ->occ_count - 1) - pos_length_occ->list_mismatches[it_id_occ].id_read;
                                                    }
                                                }
                                            }

                                            length_after_minimizer = length_read - pos_minimizer - size_minimizer;

                                            memset(id_minimizer, 0, size_minimizer_bytes * sizeof(uint8_t));
                                            parseKmerCount(curr_minimizer, size_minimizer, id_minimizer, 0);

                                            for (i = 0, minimizer = 0; i < size_minimizer_bytes; i++)
                                                minimizer = (minimizer << SIZE_BITS_UINT_8T) | reverse_word_8(id_minimizer[i]);

                                            minimizer >>= minimizer_shift;

                                            sprintf(buffer_int, "_%" PRIu64, minimizer);
                                            strcpy(&bin_name[length_prefix_bin_name], buffer_int);

                                            if (minimizer != opened_file_minimizer) continue_loop = true;

                                            break;
                                        }
                                    }
                                }

                                pos_start_line_tmp = ftell(file_new_bin);
                            }

                            fclose(file_new_bin);

                            if (!continue_loop) break;
                        }

                        if (length_spanning){

                            fwrite(read, sizeof(char), length_read, file_span_super_reads);
                            fwrite(&nl, sizeof(char), 1, file_span_super_reads);

                            //Sort read position
                            list_count = list_pos_length_occ->count;

                            if (list_count > size_array_pos_length_occ){

                                size_array_pos_length_occ = list_count;

                                array_pos_length_occ = realloc(array_pos_length_occ, size_array_pos_length_occ * sizeof(Pos_length_occ*));
                                ASSERT_NULL_PTR(array_pos_length_occ, "binning_reads()\n")
                            }

                            for (cur = list_pos_length_occ->first, i = 0; cur != NULL; cur = cur->next, i++)
                                array_pos_length_occ[i] = (Pos_length_occ*) cur->value;

                            List_destroy(list_pos_length_occ);

                            qsort(array_pos_length_occ, list_count, sizeof(Pos_length_occ*), pos_length_occ_cmp);

                            first_pos_length_occ = true;
                            position = 0;

                            for (i = 0; i < list_count; i++){

                                pos_length_occ = array_pos_length_occ[i];

                                position = (uint32_t) pos_length_occ->position - position;
                                occ_count = (uint32_t) pos_length_occ->occ_count;
                                rev = (uint32_t) pos_length_occ->rev;
                                length = ((uint32_t) pos_length_occ->length) << 1;

                                if (first_pos_length_occ){
                                    length |= 0x1;
                                    first_pos_length_occ = false;
                                }

                                read_info = (read_info << 1) | ((uint8_t) rev);
                                length_read_info++;

                                if (length_read_info == SIZE_BITS_UINT_8T){
                                    fwrite(&read_info, sizeof(uint8_t), 1, file_span_super_reads_rev);
                                    read_info = 0;
                                    length_read_info = 0;
                                }

                                fwrite(&position, sizeof(uint32_t), 1, file_span_super_reads_pos);
                                fwrite(&length, sizeof(uint32_t), 1, file_span_super_reads_length);
                                fwrite(&occ_count, sizeof(uint32_t), 1, file_span_super_reads_occ);

                                if (pair_ended){

                                    for (it_id_occ = 0; it_id_occ < occ_count; it_id_occ++){

                                        fseek(file_span_super_reads_pos_mate, pos_length_occ->id_occ[it_id_occ] * sizeof(int64_t), SEEK_SET);
                                        fwrite(&nb_reads, sizeof(int64_t), 1, file_span_super_reads_pos_mate);

                                        nb_reads++;
                                    }

                                    free(pos_length_occ->id_occ);
                                }

                                rev = (uint32_t) pos_length_occ->nb_mismatches;
                                fwrite(&rev, sizeof(uint32_t), 1, file_span_super_reads_nb_mis);

                                qsort(pos_length_occ->list_mismatches, pos_length_occ->nb_mismatches, sizeof(Mismatch), mismatch_cmp);

                                for (it_id_occ = 0, position = 0; it_id_occ < pos_length_occ->nb_mismatches; it_id_occ++){

                                    position = ((uint32_t) pos_length_occ->list_mismatches[it_id_occ].position) - position;
                                    length = (uint32_t) pos_length_occ->list_mismatches[it_id_occ].id_read;

                                    fwrite(&position, sizeof(uint32_t), 1, file_span_super_reads_pos_mis);
                                    fwrite(&length, sizeof(uint32_t), 1, file_span_super_reads_id_mis);
                                    fwrite(&(pos_length_occ->list_mismatches[it_id_occ].mismatch_char), sizeof(char), 1, file_span_super_reads_char_mis);

                                    position = (uint32_t) pos_length_occ->list_mismatches[it_id_occ].position;
                                }

                                position = (uint32_t) pos_length_occ->position;

                                free(pos_length_occ->list_mismatches);
                                free(pos_length_occ);
                            }
                        }
                    }
                }

                pos_start_line = ftell(file_bin);
            }

            fclose(file_bin);
        }
    }

    for (it_buffer = 0; it_buffer < pow(4, size_minimizer); it_buffer++){

        sprintf(buffer_int, "_%" PRIu64, it_buffer);
        strcpy(&bin_name[length_prefix_bin_name], buffer_int);

        if ((file_bin = fopen(bin_name, "r")) != NULL){

            while (getline(&buffer_read, &size_buffer_read, file_bin) != -1){

                if (buffer_read[0] != read_used){

                    first_pos_length_occ = true;
                    position = 0;
                    delimiter_read = strchr(buffer_read, ' ') + 1;

                    delimiter_meta = strchr(delimiter_read, ' ');
                    delimiter_meta[0] = eol;
                    delimiter_meta++;

                    if (delimiter_meta[strlen(delimiter_meta)-1] == nl) delimiter_meta[strlen(delimiter_meta)-1] = eol;

                    fwrite(delimiter_read, sizeof(char), strlen(delimiter_read), file_span_super_reads);
                    fwrite(&nl, sizeof(char), 1, file_span_super_reads);

                    while (delimiter_meta[0] != eol){

                        pos_length_occ = parse_pos_length_occ(&delimiter_meta, pair_ended);

                        position_tmp = (uint32_t) pos_length_occ->position;
                        position = position_tmp - position;

                        length = ((uint32_t) pos_length_occ->length) << 1;
                        rev = (uint32_t) pos_length_occ->rev;
                        occ_count = (uint32_t) pos_length_occ->occ_count;

                        if (first_pos_length_occ){
                            length |= 0x1;
                            first_pos_length_occ = false;
                        }

                        read_info = (read_info << 1) | ((uint8_t) rev);
                        length_read_info++;

                        if (length_read_info == SIZE_BITS_UINT_8T){
                            fwrite(&read_info, sizeof(uint8_t), 1, file_span_super_reads_rev);
                            read_info = 0;
                            length_read_info = 0;
                        }

                        fwrite(&position, sizeof(uint32_t), 1, file_span_super_reads_pos);
                        fwrite(&length, sizeof(uint32_t), 1, file_span_super_reads_length);
                        fwrite(&occ_count, sizeof(uint32_t), 1, file_span_super_reads_occ);

                        if (pair_ended){

                            for (it_id_occ = 0; it_id_occ < occ_count; it_id_occ++){

                                fseek(file_span_super_reads_pos_mate, pos_length_occ->id_occ[it_id_occ] * sizeof(int64_t), SEEK_SET);
                                fwrite(&nb_reads, sizeof(int64_t), 1, file_span_super_reads_pos_mate);

                                nb_reads++;
                            }

                            free(pos_length_occ->id_occ);
                        }

                        rev = (uint32_t) pos_length_occ->nb_mismatches;
                        fwrite(&rev, sizeof(uint32_t), 1, file_span_super_reads_nb_mis);

                        for (it_id_occ = 0, position = 0; it_id_occ < pos_length_occ->nb_mismatches; it_id_occ++){

                            position = ((uint32_t) pos_length_occ->list_mismatches[it_id_occ].position) - position;
                            length = (uint32_t) pos_length_occ->list_mismatches[it_id_occ].id_read;

                            fwrite(&position, sizeof(uint32_t), 1, file_span_super_reads_pos_mis);
                            fwrite(&length, sizeof(uint32_t), 1, file_span_super_reads_id_mis);
                            fwrite(&(pos_length_occ->list_mismatches[it_id_occ].mismatch_char), sizeof(char), 1, file_span_super_reads_char_mis);

                            position = (uint32_t) pos_length_occ->list_mismatches[it_id_occ].position;
                        }

                        position = position_tmp;

                        free(pos_length_occ->list_mismatches);
                        free(pos_length_occ);
                    }
                }
            }

            fclose(file_bin);

            if (remove(bin_name)) printf("Warning: could not delete temporary file.\n");
        }
    }

    if (length_read_info){
        read_info <<= SIZE_BITS_UINT_8T - length_read_info;
        fwrite(&read_info, sizeof(uint8_t), 1, file_span_super_reads_rev);
    }

    fseek(file_span_super_reads, 0 - ((int) sizeof(char)), SEEK_CUR);
    fwrite(&eol, sizeof(char), 1, file_span_super_reads);

    tot_length_seq = ftell(file_span_super_reads);

    fclose(file_span_super_reads_pos);
    fclose(file_span_super_reads_length);
    fclose(file_span_super_reads_rev);
    fclose(file_span_super_reads_occ);
    fclose(file_span_super_reads_pos_mis);
    fclose(file_span_super_reads_char_mis);
    fclose(file_span_super_reads_id_mis);
    fclose(file_span_super_reads_nb_mis);

    fclose(file_span_super_reads);

    if (pair_ended){

        fclose(file_span_super_reads_pos_mate);

        strcpy(bin_name, prefix_bin_name);
        strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_mate");

        file_span_super_reads_mate = fopen(bin_name, "w+b");
        ASSERT_NULL_PTR(file_span_super_reads_mate, "create_spanning_super_reads()\n")

        strcpy(bin_name, prefix_bin_name);
        strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_mate_info");

        file_span_super_reads_mate_info = fopen(bin_name, "w+b");
        ASSERT_NULL_PTR(file_span_super_reads_mate_info, "create_spanning_super_reads()\n")

        strcpy(bin_name, prefix_bin_name);
        strcpy(&bin_name[length_prefix_bin_name], "_span_sup_reads_tmp");

        file_span_super_reads_pos_mate = fopen(bin_name, "rb");
        ASSERT_NULL_PTR(file_span_super_reads_pos_mate, "create_spanning_super_reads()\n")

        file_span_super_reads_pos_mate_tmp = fopen(bin_name, "rb");
        ASSERT_NULL_PTR(file_span_super_reads_pos_mate_tmp, "create_spanning_super_reads()\n")

        fseek(file_span_super_reads_pos_mate_tmp, nb_reads_per_file * sizeof(int64_t), SEEK_SET);

        ftruncate(fileno(file_span_super_reads_mate_info), CEIL(nb_reads_per_file * 2, SIZE_BITS_UINT_8T));

        for (nb_reads = 0; nb_reads < nb_reads_per_file; nb_reads++){

            fread(&id_occ, sizeof(int64_t), 1, file_span_super_reads_pos_mate);
            fread(&nb_reads_tmp, sizeof(int64_t), 1, file_span_super_reads_pos_mate_tmp);

            fseek(file_span_super_reads_mate_info, (id_occ / SIZE_BITS_UINT_8T) * sizeof(uint8_t), SEEK_SET);
            fread(&mate_info, sizeof(uint8_t), 1, file_span_super_reads_mate_info);

            mate_info |= MASK_POWER_8[id_occ%SIZE_BITS_UINT_8T];

            fseek(file_span_super_reads_mate_info, 0 - sizeof(uint8_t), SEEK_CUR);
            fwrite(&mate_info, sizeof(uint8_t), 1, file_span_super_reads_mate_info);

            nb_reads_tmp -= id_occ;

            fseek(file_span_super_reads_mate, id_occ * sizeof(int64_t), SEEK_SET);
            fwrite(&nb_reads_tmp, sizeof(int64_t), 1, file_span_super_reads_mate);
        }

        fclose(file_span_super_reads_pos_mate);
        fclose(file_span_super_reads_pos_mate_tmp);

        rewind(file_span_super_reads_mate);
        rewind(file_span_super_reads_mate_info);

        id_occ = 0;
        nb_reads = 0;

        while (fread(&mate_info, sizeof(uint8_t), 1, file_span_super_reads_mate_info) == 1){

            for (i = 0; mate_info && (i < SIZE_BITS_UINT_8T); i++, mate_info >>= 1){

                if (mate_info & 0x1){

                    fseek(file_span_super_reads_mate, (nb_reads + i) * sizeof(int64_t), SEEK_SET);
                    fread(&nb_reads_tmp, sizeof(int64_t), 1, file_span_super_reads_mate);

                    fseek(file_span_super_reads_mate, id_occ * sizeof(int64_t), SEEK_SET);
                    fwrite(&nb_reads_tmp, sizeof(int64_t), 1, file_span_super_reads_mate);

                    id_occ++;
                }
            }

            nb_reads += SIZE_BITS_UINT_8T;
        }

        ftruncate(fileno(file_span_super_reads_mate), id_occ * sizeof(int64_t));

        fclose(file_span_super_reads_mate);
        fclose(file_span_super_reads_mate_info);

        if (remove(bin_name)) printf("Warning: could not delete temporary file.\n");
    }

    free(bin_name);
    free(id_minimizer);

    free(mis);

    free(array_pos_length_occ);

    free(buffer_read);
    free(buffer_read_tmp);
    free(buffer_read_rev);

    return tot_length_seq;
}

Pos_length_occ* parse_pos_length_occ(char** meta_to_parse, bool pair_ended){

    ASSERT_NULL_PTR(*meta_to_parse, "parse_pos_length_occ()\n")

    Pos_length_occ* pos_length_occ = create_pos_length_occ();

    pos_length_occ->position = custom_atoi(*meta_to_parse, meta_to_parse, ',');
    (*meta_to_parse)++;

    pos_length_occ->length = custom_atoi(*meta_to_parse, meta_to_parse, ',');
    (*meta_to_parse)++;

    pos_length_occ->rev = custom_atoi(*meta_to_parse, meta_to_parse, ',');
    (*meta_to_parse)++;

    pos_length_occ->occ_count = custom_atoi(*meta_to_parse, meta_to_parse, ',');
    (*meta_to_parse)++;

    pos_length_occ->nb_mismatches = custom_atoi(*meta_to_parse, meta_to_parse, ',');
    (*meta_to_parse)++;

    if (pair_ended){

        pos_length_occ->id_occ = malloc(pos_length_occ->occ_count * sizeof(int64_t));
        ASSERT_NULL_PTR(pos_length_occ->id_occ, "parse_pos_length_occ()\n")

        for (int i = 0; i < pos_length_occ->occ_count; i++){
            pos_length_occ->id_occ[i] = (int64_t) strtoll(*meta_to_parse, meta_to_parse, 10);
            (*meta_to_parse)++;
        }
    }

    if (pos_length_occ->nb_mismatches){

        pos_length_occ->list_mismatches = malloc(pos_length_occ->nb_mismatches * sizeof(Mismatch));
        ASSERT_NULL_PTR(pos_length_occ->list_mismatches, "parse_pos_length_occ()\n")

        for (int i = 0; i < pos_length_occ->nb_mismatches; i++){

            pos_length_occ->list_mismatches[i].position = custom_atoi(*meta_to_parse, meta_to_parse, ',');
            (*meta_to_parse)++;

            pos_length_occ->list_mismatches[i].mismatch_char = **meta_to_parse;
            (*meta_to_parse) += 2;

            pos_length_occ->list_mismatches[i].id_read = custom_atoi(*meta_to_parse, meta_to_parse, ',');
            (*meta_to_parse)++;
        }
    }

    return pos_length_occ;
}

int string_cmp(const void *a, const void *b){

    const char **ia = (const char **)a;
    const char **ib = (const char **)b;

    return strcmp(*ia, *ib);
}

int pos_length_occ_cmp(const void *a, const void *b){

    Pos_length_occ** pos_len_occ_a = (Pos_length_occ**) a;
    Pos_length_occ** pos_len_occ_b = (Pos_length_occ**) b;

    if ((*pos_len_occ_a)->position > (*pos_len_occ_b)->position) return 1;
    if ((*pos_len_occ_a)->position < (*pos_len_occ_b)->position)  return -1;

    return 0;
}

int pos_length_occ_cmp_rev(const void *a, const void *b){

    Pos_length_occ** pos_len_occ_a = (Pos_length_occ**) a;
    Pos_length_occ** pos_len_occ_b = (Pos_length_occ**) b;

    if ((*pos_len_occ_a)->position > (*pos_len_occ_b)->position) return -1;
    if ((*pos_len_occ_a)->position < (*pos_len_occ_b)->position)  return 1;

    return 0;
}

int mismatch_cmp(const void *a, const void *b){

    Mismatch* mis_a = (Mismatch*) a;
    Mismatch* mis_b = (Mismatch*) b;

    if (mis_a->position > mis_b->position) return 1;
    if (mis_a->position < mis_b->position)  return -1;

    return 0;
}

int mismatch_cmp_rev(const void *a, const void *b){

    Mismatch* mis_a = (Mismatch*) a;
    Mismatch* mis_b = (Mismatch*) b;

    if (mis_a->position > mis_b->position) return -1;
    if (mis_a->position < mis_b->position)  return 1;

    return 0;
}

int mismatch_pos_cmp(const void *a, const void *b){

    Mismatch* mis_a = (Mismatch*) a;
    Mismatch* mis_b = (Mismatch*) b;

    if (mis_a->id_read > mis_b->id_read) return 1;
    if (mis_a->id_read < mis_b->id_read)  return -1;

    if (mis_a->position > mis_b->position) return 1;
    if (mis_a->position < mis_b->position)  return -1;

    return 0;
}

int custom_atoi(const char *old_str, char** new_str, char end) {
    int x = 0;
    while (*old_str != end) {
        x = (x*10) + (*old_str++ - '0');
    }
    *new_str = old_str;
    return x;
};

void reorder_reads(char* filename_src, char* filename_dest, int size_minimizer, int size_window)
{
    ASSERT_NULL_PTR(filename_src,"reorder_reads()\n")
    ASSERT_NULL_PTR(filename_dest,"reorder_reads()\n")

    bool write_everything = false;

    int nb_bytes_minimizer = CEIL(size_minimizer * 4, SIZE_BITS_UINT_8T);

    uint64_t length_read;

    uint64_t nb_read = 0;
    uint64_t nb_read_processed = 0xffffffffffffffff;

    uint64_t size_window_tmp = size_window;
    uint64_t size_processed_reads = SIZE_BUFFER * SIZE_BITS_UINT_8T;
    uint64_t size_array_minimizers = SIZE_BUFFER;

    double first_reading = 0;

    uint8_t* processed_reads = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(processed_reads, "reorder_reads()\n")

    uint8_t* minimizers = calloc(SIZE_BUFFER, sizeof(char));
    ASSERT_NULL_PTR(minimizers, "reorder_reads()\n")

    char* curr_minimizer;

    char* buffer_read = malloc(SIZE_BUFFER * sizeof(char));
    ASSERT_NULL_PTR(buffer_read, "reorder_reads()\n")

    char* buffer_read_rev = malloc(SIZE_BUFFER * sizeof(char));
    ASSERT_NULL_PTR(buffer_read_rev, "reorder_reads()\n")

    List* list_minimizers = List_create();
    ASSERT_NULL_PTR(list_minimizers, "reorder_reads()\n")

    Pvoid_t PJArray = (PWord_t)NULL;
    PWord_t PValue;
    Word_t Rc_word;
    int Rc_int;

    FILE* file_reads_src;

    FILE* file_reads_dest = fopen(filename_dest, "w");
    ASSERT_NULL_PTR(file_reads_dest, "reorder_reads()\n")

    while (nb_read_processed != nb_read){

        nb_read = 0;
        first_reading++;

        file_reads_src = fopen(filename_src, "r");
        ASSERT_NULL_PTR(file_reads_src, "reorder_reads()\n")

        if (first_reading == 1){

            while (fgets(buffer_read, SIZE_BUFFER, file_reads_src) != NULL){

                if (nb_read+1 > size_processed_reads){
                    processed_reads = realloc(processed_reads, (size_processed_reads/SIZE_BITS_UINT_8T + SIZE_BUFFER) * sizeof(uint8_t));
                    ASSERT_NULL_PTR(processed_reads, "reorder_reads()\n")

                    memset(&processed_reads[size_processed_reads/SIZE_BITS_UINT_8T], 0 , SIZE_BUFFER * sizeof(uint8_t));

                    size_processed_reads += SIZE_BUFFER * SIZE_BITS_UINT_8T;
                }

                if ((nb_read+1) * nb_bytes_minimizer > size_array_minimizers){

                    minimizers = realloc(minimizers, (size_array_minimizers + SIZE_BUFFER) * sizeof(uint8_t));
                    ASSERT_NULL_PTR(minimizers, "reorder_reads()\n")

                    memset(&minimizers[size_array_minimizers], 0, SIZE_BUFFER * sizeof(uint8_t));

                    size_array_minimizers += SIZE_BUFFER;
                }

                length_read = strlen(buffer_read) - 1;

                curr_minimizer = buffer_read;

                reverse_complement(buffer_read, buffer_read_rev, length_read);

                for (uint32_t i = 0; i < length_read - size_minimizer; i++){
                    if (memcmp(&buffer_read[i], curr_minimizer, size_minimizer * sizeof(char)) < 0)
                        curr_minimizer = &buffer_read[i];
                    if (memcmp(&buffer_read_rev[i], curr_minimizer, size_minimizer * sizeof(char)) < 0)
                        curr_minimizer = &buffer_read_rev[i];
                }

                parseKmerCount_IUPAC(curr_minimizer, size_minimizer, &minimizers[nb_read * nb_bytes_minimizer], 0);

                nb_read++;
            }
        }
        else {

            nb_read_processed = 0;

            while (fgets(buffer_read, SIZE_BUFFER, file_reads_src) != NULL){

                if (!(processed_reads[nb_read/SIZE_BITS_UINT_8T] & MASK_POWER_8[nb_read%SIZE_BITS_UINT_8T])){

                    if (!write_everything){

                        length_read = strlen(buffer_read) - 1;
                        curr_minimizer = (char*) &minimizers[nb_read * nb_bytes_minimizer];

                        JHSG(PValue, PJArray, curr_minimizer, nb_bytes_minimizer);

                        if (PValue == NULL){

                            if (list_minimizers->count == size_window_tmp){

                                JHSD(Rc_int, PJArray, list_minimizers->first->value, nb_bytes_minimizer);
                                List_first_become_last(list_minimizers);
                                list_minimizers->last->value = curr_minimizer;
                            }
                            else List_push(list_minimizers, curr_minimizer);

                            JHSI(PValue, PJArray, curr_minimizer, nb_bytes_minimizer);

                            fwrite(buffer_read, sizeof(char), length_read+1, file_reads_dest);

                            nb_read_processed++;
                            processed_reads[nb_read/SIZE_BITS_UINT_8T] |= MASK_POWER_8[nb_read%SIZE_BITS_UINT_8T];
                        }
                    }
                    else{
                        fwrite(buffer_read, sizeof(char), length_read+1, file_reads_dest);
                        nb_read_processed++;
                    }
                }
                else nb_read_processed++;

                nb_read++;
            }

            size_window_tmp = MAX(1, size_window_tmp / acosh(first_reading+1));

            if (size_window_tmp == 1) write_everything = true;

            while (list_minimizers->count > size_window_tmp) JHSD(Rc_int, PJArray, List_pop_first(list_minimizers), nb_bytes_minimizer);
        }

        fclose(file_reads_src);
    }

    fclose(file_reads_dest);

    List_destroy(list_minimizers);

    free(minimizers);
    free(processed_reads);
    free(buffer_read);
    free(buffer_read_rev);

    JHSFA(Rc_word, PJArray);

    return;
}

void min_simReads_max_simkmers2(char* filename_read_1, char* filename_read_2, bool pair_ended,
                               char* output_prefix, int size_kmer, int size_seed, int size_minimizer,
                               bool compress_shift, BFT_Root* graph_no_iupac, BFT_Root* graph_iupac)
{
    ASSERT_NULL_PTR(filename_read_1,"min_simReads_max_simkmers()")
    ASSERT_NULL_PTR(output_prefix,"min_simReads_max_simkmers()")

    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    const char eol = '\n';

    FILE* file_reads;
    FILE* file_kmers;
    FILE* file_runs;
    FILE* file_shift;
    FILE* file_info_read;

    BloomFilter* inserted_kmers;

    memory_Used* mem;

    int pos = 0;
    int shift = 0;
    int nb_read = 0;
    int nb_kmer = 0;
    int nb_kmer_present = 0;
    int nb_kmer_to_insert = 0;
    int best_shift = 0;
    int best_nb_kmer_present = 0;
    int size_buffer_kmer = 0;
    int size_buffer_info_read = 0;
    int prev_size_buffer_shift = 0;
    int size_buffer_shift = 0;

    size_t size_buffer_rev_comp = SIZE_BUFFER;
    size_t size_buffer_read = SIZE_BUFFER;
    size_t size_all_hash_val = SIZE_BUFFER / sizeof(uint64_t);

    int len_output_prefix = strlen(output_prefix);
    int nb_bytes_iupac = CEIL(size_kmer * 4, SIZE_BITS_UINT_8T);

    int remaining = size_kmer - size_seed;

    int i;
    int length_seq;
    int length_nuc_left_read;

    int64_t nb_read_per_file;

    uint64_t hash_v;
    uint64_t hash_v2;
    uint64_t hash_tmp;

    uint64_t nb_kmer_file = 0;
    uint64_t max_nb_kmer_file = 0;
    uint64_t nb_autoloop_kmers = 0;

    uint32_t coef_threshold = 1;

    double fp_rate_bf = 0.01;

    bool is_rev = false;

    char* buffer_rev_comp;
    char* buffer_read;
    char* buffer_kmer;
    char* buffer_shift;
    char* output_prefix_buffer;

    uint8_t* buffer_info_read;
    uint8_t* kmer;
    uint8_t* curr_block;

    uint64_t* hash_val;
    uint64_t* hash_val_rev;
    uint64_t* curr_hash_val;

    //printf("\nRe-ordering reads\n\n");

    nb_read_per_file = binning_reads(filename_read_1, filename_read_2, pair_ended, output_prefix, size_minimizer);
    create_super_reads(output_prefix, pair_ended, size_minimizer);
    max_nb_kmer_file = create_spanning_super_reads(output_prefix, pair_ended, nb_read_per_file, size_minimizer) / remaining;
    //max_nb_kmer_file = create_spanning_super_reads2(output_prefix, pair_ended, nb_read_per_file, size_minimizer) / remaining;

    if ((graph_no_iupac != NULL) && (graph_iupac != NULL)){

        mem = printMemoryUsedFromNode(&(graph_no_iupac->node), (graph_no_iupac->k / NB_CHAR_SUF_PREF) - 1,
                                      graph_no_iupac->k, graph_no_iupac->info_per_lvl);

        max_nb_kmer_file += mem->nb_kmers_in_UCptr;

        free(mem);

        mem = printMemoryUsedFromNode(&(graph_iupac->node), (graph_iupac->k / NB_CHAR_SUF_PREF) - 1,
                                                   graph_iupac->k, graph_iupac->info_per_lvl);

        max_nb_kmer_file += mem->nb_kmers_in_UCptr;

        free(mem);
    }

    inserted_kmers = create_BloomFilter(max_nb_kmer_file, fp_rate_bf);

    if ((graph_no_iupac != NULL) && (graph_iupac != NULL)){
        iterate_over_kmers(graph_iupac, insert_kmer_into_bf_from_graph, true, inserted_kmers);
        iterate_over_kmers(graph_no_iupac, insert_kmer_into_bf_from_graph, false, inserted_kmers);
    }

    output_prefix_buffer = malloc((len_output_prefix + 50) * sizeof(char));
    ASSERT_NULL_PTR(output_prefix_buffer, "min_simReads_max_simkmers() 2\n")

    strcpy(output_prefix_buffer, output_prefix);

    buffer_read = malloc(size_buffer_read * sizeof(char));
    ASSERT_NULL_PTR(buffer_read, "min_simReads_max_simkmers() 7\n")

    buffer_rev_comp = calloc(size_buffer_rev_comp, sizeof(char));
    ASSERT_NULL_PTR(buffer_rev_comp, "min_simReads_max_simkmers() 1\n")

    buffer_kmer = malloc(SIZE_BUFFER * sizeof(char));
    ASSERT_NULL_PTR(buffer_kmer, "min_simReads_max_simkmers() 12\n")

    buffer_shift = malloc(SIZE_BUFFER * sizeof(char));
    ASSERT_NULL_PTR(buffer_shift, "min_simReads_max_simkmers() 13\n")

    buffer_info_read = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(buffer_info_read, "min_simReads_max_simkmers() 16\n")

    kmer = calloc(nb_bytes_iupac, sizeof(uint8_t));
    ASSERT_NULL_PTR(kmer, "min_simReads_max_simkmers() 16\n")

    hash_val = malloc(size_all_hash_val * sizeof(uint64_t));
    ASSERT_NULL_PTR(hash_val,"min_simReads_max_simkmers() 17\n")

    hash_val_rev = malloc(size_all_hash_val * sizeof(uint64_t));
    ASSERT_NULL_PTR(hash_val_rev,"min_simReads_max_simkmers() 18\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads");
    file_reads = fopen(output_prefix_buffer, "r");
    ASSERT_NULL_PTR(file_reads,"min_simReads_max_simkmers() 23\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_kmers_tmp");
    file_kmers = fopen(output_prefix_buffer, "wb");
    ASSERT_NULL_PTR(file_kmers,"min_simReads_max_simkmers() 20\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_runs");
    file_runs = fopen(output_prefix_buffer, "wb");
    ASSERT_NULL_PTR(file_runs,"min_simReads_max_simkmers() 20\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_shifts");
    file_shift = fopen(output_prefix_buffer, "wb");
    ASSERT_NULL_PTR(file_shift,"min_simReads_max_simkmers() 21\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_read_info");
    file_info_read = fopen(output_prefix_buffer, "wb");
    ASSERT_NULL_PTR(file_info_read,"min_simReads_max_simkmers() 22\n")

    if (compress_shift) fwrite((uint8_t []){0x1}, sizeof(uint8_t), 1, file_shift);
    else fwrite((uint8_t []){0x0}, sizeof(uint8_t), 1, file_shift);

    while (getline(&buffer_read, &size_buffer_read, file_reads) != -1){

        length_seq = strlen(buffer_read);
        length_seq -= (buffer_read[length_seq - 1] == '\n');

        is_rev = false;
        best_shift = 0;
        nb_kmer_to_insert = 0;
        best_nb_kmer_present = 0;

        if (length_seq * 2 > size_all_hash_val){

            size_all_hash_val = length_seq * 2;

            hash_val = realloc(hash_val, size_all_hash_val * sizeof(uint64_t));
            ASSERT_NULL_PTR(hash_val,"min_simReads_max_simkmers() 17\n")

            hash_val_rev = realloc(hash_val_rev, size_all_hash_val * sizeof(uint64_t));
            ASSERT_NULL_PTR(hash_val_rev,"min_simReads_max_simkmers() 18\n")
        }

        nb_kmer_to_insert = length_seq / remaining;
        if (nb_kmer_to_insert * remaining + size_seed >= length_seq) nb_kmer_to_insert--;

        length_nuc_left_read = length_seq - size_seed - nb_kmer_to_insert * remaining;

        for (shift = 0; shift <= length_nuc_left_read; shift++){

            nb_kmer_present = 0;

            for (pos = shift, nb_kmer = nb_kmer_to_insert; nb_kmer > 0; pos += remaining, nb_kmer--){

                nb_kmer_present++;

                hash_val[pos * 2] = (hash_v = XXH64(&buffer_read[pos], size_kmer, inserted_kmers->seed_hash1));
                hash_val[pos * 2 + 1] = (hash_v2 = XXH64(&buffer_read[pos], size_kmer, inserted_kmers->seed_hash2));

                curr_block = &(inserted_kmers->bf[(hash_v % inserted_kmers->nb_blocks) * inserted_kmers->size_block_bytes]);
                __builtin_prefetch(curr_block, 0, 3);

                for (i = 1; i <= inserted_kmers->nb_hash; i++){
                    hash_tmp = (hash_v + i * hash_v2) % inserted_kmers->size_block_bits;
                    if (!(curr_block[hash_tmp/SIZE_BITS_UINT_8T] & MASK_POWER_8[hash_tmp%SIZE_BITS_UINT_8T])){
                        nb_kmer_present--;
                        break;
                    }
                }

                if (nb_kmer_present + nb_kmer * coef_threshold <= best_nb_kmer_present) break;
            }

            if (nb_kmer_present > best_nb_kmer_present){
                best_nb_kmer_present = nb_kmer_present;
                best_shift = shift;
                if (best_nb_kmer_present >= nb_kmer_to_insert * coef_threshold) break;
            }
        }

        if (best_nb_kmer_present < nb_kmer_to_insert * coef_threshold){

            if (size_buffer_read > size_buffer_rev_comp){

                size_buffer_rev_comp = size_buffer_read;

                buffer_rev_comp = realloc(buffer_rev_comp, size_buffer_rev_comp * sizeof(char));
                ASSERT_NULL_PTR(buffer_rev_comp, "min_simReads_max_simkmers() 20\n")
            }

            reverse_complement(buffer_read, buffer_rev_comp, length_seq);

            for (shift = 0; shift <= length_nuc_left_read; shift++){

                nb_kmer_present = 0;

                for (pos = shift, nb_kmer = nb_kmer_to_insert; nb_kmer > 0; pos += remaining, nb_kmer--){

                    nb_kmer_present++;

                    hash_val_rev[pos * 2] = (hash_v = XXH64(&buffer_rev_comp[pos], size_kmer, inserted_kmers->seed_hash1));
                    hash_val_rev[pos * 2 + 1] = (hash_v2 = XXH64(&buffer_rev_comp[pos], size_kmer, inserted_kmers->seed_hash2));

                    curr_block = &(inserted_kmers->bf[(hash_v % inserted_kmers->nb_blocks) * inserted_kmers->size_block_bytes]);
                    __builtin_prefetch(curr_block, 0, 3);

                    for (i = 1; i <= inserted_kmers->nb_hash; i++){
                        hash_tmp = (hash_v + i * hash_v2) % inserted_kmers->size_block_bits;
                        if (!(curr_block[hash_tmp/SIZE_BITS_UINT_8T] & MASK_POWER_8[hash_tmp%SIZE_BITS_UINT_8T])){
                            nb_kmer_present--;
                            break;
                        }
                    }

                    if (nb_kmer_present + nb_kmer * coef_threshold <= best_nb_kmer_present) break;
                }

                if (nb_kmer_present > best_nb_kmer_present){
                    best_nb_kmer_present = nb_kmer_present;
                    best_shift = shift;
                    is_rev = true;
                    if (best_nb_kmer_present >= nb_kmer_to_insert * coef_threshold) break;
                }
            }

            if (is_rev){
                curr_hash_val = hash_val_rev;
                memcpy(buffer_read, buffer_rev_comp, length_seq * sizeof(char));
            }
            else curr_hash_val = hash_val;
        }
        else curr_hash_val = hash_val;

        length_nuc_left_read = best_shift + size_seed;

        if (best_shift || ((nb_kmer_to_insert + 1) * remaining + size_seed != length_seq)) length_nuc_left_read += nb_kmer_to_insert * remaining;
        else length_nuc_left_read += (nb_kmer_to_insert + 1) * remaining;

        if (size_buffer_shift + (size_kmer + 1) * 2 > SIZE_BUFFER){
            fwrite(buffer_shift, sizeof(char), size_buffer_shift, file_shift);
            size_buffer_shift = 0;
        }

        if (compress_shift){

            memset(kmer, 0, nb_bytes_iupac * sizeof(uint8_t));

            prev_size_buffer_shift = size_buffer_shift;
            buffer_shift[size_buffer_shift] = best_shift + size_seed;

            if (parseKmerCount_xIUPAC(buffer_read, best_shift + size_seed, kmer, kmer, 0, 0)){
                buffer_shift[size_buffer_shift] |= 0x80;
                memcpy(&buffer_shift[size_buffer_shift+1], kmer, CEIL((best_shift + size_seed) * 4, SIZE_BITS_UINT_8T) * sizeof(char));
                size_buffer_shift += CEIL((best_shift + size_seed) * 4, SIZE_BITS_UINT_8T) + 1;
            }
            else{
                memcpy(&(buffer_shift[size_buffer_shift+1]), kmer, CEIL((best_shift + size_seed) * 2, SIZE_BITS_UINT_8T) * sizeof(char));
                size_buffer_shift += CEIL((best_shift + size_seed) * 2, SIZE_BITS_UINT_8T) + 1;
            }

            memset(kmer, 0, nb_bytes_iupac * sizeof(uint8_t));

            if (parseKmerCount_xIUPAC(&buffer_read[length_nuc_left_read], length_seq - length_nuc_left_read, kmer, kmer, 0, 0)){
                buffer_shift[prev_size_buffer_shift] |= 0x40;
                prev_size_buffer_shift = CEIL((length_seq - length_nuc_left_read) * 4, SIZE_BITS_UINT_8T);
            }
            else prev_size_buffer_shift = CEIL((length_seq - length_nuc_left_read) * 2, SIZE_BITS_UINT_8T);

            memcpy(&buffer_shift[size_buffer_shift], kmer, prev_size_buffer_shift * sizeof(char));
            size_buffer_shift += prev_size_buffer_shift;
        }
        else{
            memcpy(&buffer_shift[size_buffer_shift], buffer_read, (best_shift + size_seed) * sizeof(char));
            buffer_shift[size_buffer_shift + (best_shift + size_seed)] = eol;
            size_buffer_shift += (best_shift + size_seed) + 1;

            memcpy(&buffer_shift[size_buffer_shift], &buffer_read[length_nuc_left_read], (length_seq - length_nuc_left_read) * sizeof(char));
            buffer_shift[size_buffer_shift + (length_seq - length_nuc_left_read)] = eol;
            size_buffer_shift += (length_seq - length_nuc_left_read) + 1;
        }

        if (size_buffer_info_read + 1 > SIZE_BUFFER){
            fwrite(buffer_info_read, sizeof(uint8_t), size_buffer_info_read, file_info_read);
            size_buffer_info_read = 0;
        }

        buffer_info_read[size_buffer_info_read] <<= 1;
        if (is_rev) buffer_info_read[size_buffer_info_read] |= 0x1;

        nb_read++;
        if (nb_read%SIZE_BITS_UINT_8T == 0) size_buffer_info_read++;

        for (pos = best_shift; pos <= length_seq - size_kmer; pos += remaining){

            if (size_buffer_kmer + size_kmer + 1 >= SIZE_BUFFER){

                for (i = 0; i < size_buffer_kmer; i += size_kmer + 1){

                    if (memcmp(&buffer_kmer[i], &buffer_kmer[i + size_kmer - size_seed], size_seed * sizeof(char)) == 0){

                        nb_autoloop_kmers++;

                        fwrite(&buffer_kmer[i], sizeof(char), (size_kmer + 1), file_runs);
                        fwrite(&nb_kmer_file, sizeof(uint64_t), 1, file_runs);
                    }
                    else {

                        memset(kmer, 0, nb_bytes_iupac * sizeof(uint8_t));
                        parseKmerCount_IUPAC(&buffer_kmer[i], size_kmer, kmer, 0);

                        fwrite(kmer, sizeof(uint8_t), nb_bytes_iupac, file_kmers);
                    }

                    nb_kmer_file++;
                }

                size_buffer_kmer = 0;
            }

            memcpy(&buffer_kmer[size_buffer_kmer], &buffer_read[pos], size_kmer * sizeof(char));
            size_buffer_kmer += size_kmer;
            buffer_kmer[size_buffer_kmer] = eol;
            size_buffer_kmer++;

            hash_v = curr_hash_val[pos * 2];
            hash_v2 = curr_hash_val[pos * 2 + 1];

            curr_block = &(inserted_kmers->bf[(hash_v % inserted_kmers->nb_blocks) * inserted_kmers->size_block_bytes]);
            __builtin_prefetch(curr_block, 1, 3);

            for (i = 1; i <= inserted_kmers->nb_hash; i++){
                hash_tmp = (hash_v + i * hash_v2) % inserted_kmers->size_block_bits;
                curr_block[hash_tmp/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash_tmp%SIZE_BITS_UINT_8T];
            }
        }

        if (nb_read%1000000 == 0) printf("# reads = %d\n", nb_read);
    }

    fclose(file_reads);

    if (size_buffer_kmer){

        for (i = 0; i < size_buffer_kmer; i += size_kmer + 1){

            if (memcmp(&(buffer_kmer[i]), &(buffer_kmer[i + size_kmer - size_seed]), size_seed * sizeof(char)) == 0){

                nb_autoloop_kmers++;

                fwrite(&(buffer_kmer[i]), sizeof(char), (size_kmer + 1), file_runs);
                fwrite(&nb_kmer_file, sizeof(uint64_t), 1, file_runs);
            }
            else {

                memset(kmer, 0, nb_bytes_iupac * sizeof(uint8_t));
                parseKmerCount_IUPAC(&buffer_kmer[i], size_kmer, kmer, 0);

                fwrite(kmer, sizeof(uint8_t), nb_bytes_iupac, file_kmers);
            }

            nb_kmer_file++;
        }
    }

    fwrite(buffer_shift, sizeof(uint8_t), size_buffer_shift, file_shift);

    nb_read %= SIZE_BUFFER * SIZE_BITS_UINT_8T;

    if (nb_read % SIZE_BITS_UINT_8T){

        buffer_info_read[size_buffer_info_read] <<= SIZE_BITS_UINT_8T - (nb_read % SIZE_BITS_UINT_8T);
        fwrite(buffer_info_read, sizeof(char), size_buffer_info_read + 1, file_info_read);
    }
    else fwrite(buffer_info_read, sizeof(char), size_buffer_info_read, file_info_read);

    printf("Number of autoloops k-mers: %" PRIu64 "\n", nb_autoloop_kmers);

    fclose(file_info_read);
    fclose(file_kmers);
    fclose(file_runs);
    fclose(file_shift);

    free_BloomFilter(inserted_kmers);

    free(buffer_read);
    free(buffer_kmer);
    free(buffer_shift);
    free(buffer_info_read);
    free(buffer_rev_comp);
    free(kmer);

    free(hash_val);
    free(hash_val_rev);

    free(output_prefix_buffer);

    gettimeofday(&tval_after, NULL);
    time_spent(&tval_before, &tval_after, &tval_result);
    printf("\nElapsed time: %ld.%06ld s\n\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

    return;
}

uint32_t compressKmers_from_KmerFiles_bis2(char* output_prefix, int size_seed, int size_kmer, uint32_t prev_nb_parts){

    ASSERT_NULL_PTR(output_prefix,"compressKmers_from_KmerFiles_bis2()")

    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    Pvoid_t prefix_judy = (PWord_t)NULL;
    Pvoid_t suffix_judy = (PWord_t)NULL;

    Pvoid_t pref_pos = (PWord_t)NULL;
    Pvoid_t suf_pos = (PWord_t)NULL;

    PWord_t PValue;
    Word_t Rc_word;
    int Rc_int;

    FILE* file_kmers;
    FILE* file_parts;
    FILE* file_parts_count;

    List* list_part = List_create();

    ListNode* curr_node = NULL;

    int i = 0, j = 0, l = 0, z = 0;
    int shift = 0;

    int len_output_prefix = strlen(output_prefix);

    int nb_bytes_seed = CEIL(size_seed * 2, SIZE_BITS_UINT_8T);
    int nb_bytes_seed_iupac = CEIL(size_seed * 4, SIZE_BITS_UINT_8T);

    int nb_bytes_kmer_iupac = CEIL(size_kmer * 4, SIZE_BITS_UINT_8T);

    int nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer_iupac;

    int pos_start_comp_suffix = (size_kmer - size_seed) / 2;
    int size_after_pos_start_comp_suffix = size_kmer - pos_start_comp_suffix * 2;

    int nb_combi = pow(4,size_seed);
    int size_annot_seeds = CEIL(nb_combi,SIZE_BITS_UINT_8T);

    uint32_t prev_suffix;
    uint32_t suffix;
    uint32_t prefix;

    uint64_t min;

    uint64_t pref_div;
    uint64_t suf_div;

    uint64_t pref_iupac_tmp = 0;
    uint64_t suf_iupac_tmp = 0;
    uint64_t nb_recycled_part = 0;

    uint64_t kmers_read = 0;
    uint32_t current_partition = prev_nb_parts;
    uint32_t recycled_partition = 0;

    uint8_t pref_mod;
    uint8_t suf_mod;

    bool prev_suf_iupac;
    bool suf_iupac;
    bool pref_iupac;

    bool same_part_pref = true;
    bool same_part = true;

    bool first_kmers_part = true;
    bool last_kmers_part = true;

    uint64_t* pref_part = calloc(1, sizeof(uint64_t));
    uint64_t* suf_part = calloc(1, sizeof(uint64_t));

    uint64_t* union_pref_part = NULL;

    uint64_t** pref_suf;

    uint8_t* kmers_iupac = malloc(nb_kmer_in_buf * nb_bytes_kmer_iupac * sizeof(uint8_t));
    ASSERT_NULL_PTR(kmers_iupac, "compressKmers_from_KmerFiles()")

    uint8_t* prefixes_iupac = calloc(nb_bytes_seed_iupac, sizeof(uint8_t));
    ASSERT_NULL_PTR(prefixes_iupac, "compressKmers_from_KmerFiles()")

    uint8_t* suffixes_iupac = calloc(nb_bytes_seed_iupac, sizeof(uint8_t));
    ASSERT_NULL_PTR(suffixes_iupac, "compressKmers_from_KmerFiles()")

    uint8_t* prev_suffix_iupac = calloc(nb_bytes_seed_iupac, sizeof(uint8_t));
    ASSERT_NULL_PTR(prev_suffix_iupac, "compressKmers_from_KmerFiles()")

    uint8_t* prefixes = calloc(nb_bytes_seed, sizeof(uint8_t));
    ASSERT_NULL_PTR(prefixes, "compressKmers_from_KmerFiles()")

    uint8_t* suffixes = calloc(nb_bytes_seed, sizeof(uint8_t));
    ASSERT_NULL_PTR(suffixes, "compressKmers_from_KmerFiles()")

    uint8_t* prefixes_in_partition = calloc(size_annot_seeds, sizeof(uint8_t));
    ASSERT_NULL_PTR(prefixes_in_partition, "compressKmers_from_KmerFiles()")

    uint8_t* suffixes_in_partition = calloc(size_annot_seeds, sizeof(uint8_t));
    ASSERT_NULL_PTR(suffixes_in_partition, "compressKmers_from_KmerFiles()")

    char* line = calloc(size_kmer + 1, sizeof(char));
    ASSERT_NULL_PTR(line,"compressKmers_from_KmerFiles()")

    char* seed_ASCII = calloc(size_seed + 1, sizeof(char));
    ASSERT_NULL_PTR(seed_ASCII, "compressKmers_from_KmerFiles()")

    char* output_prefix_buffer = malloc((len_output_prefix + 50) * sizeof(char));
    ASSERT_NULL_PTR(output_prefix_buffer, "min_simReads_max_simkmers()")

    strcpy(output_prefix_buffer, output_prefix);

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_kmers_tmp");
    file_kmers = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_kmers,"compressKmers_from_KmerFiles()")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_partitions");
    file_parts = fopen(output_prefix_buffer, "wb");
    ASSERT_NULL_PTR(file_parts,"compressKmers_from_KmerFiles()")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_parts_count_tmp");
    file_parts_count = fopen(output_prefix_buffer, "wb");
    ASSERT_NULL_PTR(file_parts_count,"compressKmers_from_KmerFiles()")

    if (size_seed % 4) shift = SIZE_BITS_UINT_8T - ((size_seed % 4) * 2);

    while ((j = fread(kmers_iupac, nb_bytes_kmer_iupac * sizeof(uint8_t), nb_kmer_in_buf, file_kmers))){

        for (i = 0; i < j; i++, kmers_read++){

            memset(prefixes, 0, nb_bytes_seed * sizeof(uint8_t));
            memset(suffixes, 0, nb_bytes_seed * sizeof(uint8_t));

            memset(prefixes_iupac, 0, nb_bytes_seed_iupac * sizeof(uint8_t));
            memset(suffixes_iupac, 0, nb_bytes_seed_iupac * sizeof(uint8_t));

            // prefix
            kmer_iupac_comp_to_ascii(&kmers_iupac[i * nb_bytes_kmer_iupac], size_seed, seed_ASCII);

            if (is_substring_IUPAC(seed_ASCII)) pref_iupac = true;
            else pref_iupac = false;

            parseKmerCount(seed_ASCII, size_seed, prefixes, 0);
            parseKmerCount_IUPAC(seed_ASCII, size_seed, prefixes_iupac, 0);

            // suffix
            kmer_iupac_comp_to_ascii(&kmers_iupac[i * nb_bytes_kmer_iupac + pos_start_comp_suffix],
                                     size_after_pos_start_comp_suffix, line);

            strncpy(seed_ASCII, &line[size_after_pos_start_comp_suffix - size_seed], size_seed);

            if (is_substring_IUPAC(seed_ASCII)) suf_iupac = true;
            else suf_iupac = false;

            parseKmerCount(seed_ASCII, size_seed, suffixes, 0);
            parseKmerCount_IUPAC(seed_ASCII, size_seed, suffixes_iupac, 0);

            same_part_pref = true;
            same_part = true;

            if (!pref_iupac){

                for (l = 0, prefix = 0; l < nb_bytes_seed; l++)
                    prefix = (prefix << SIZE_BITS_UINT_8T) | reverse_word_8(prefixes[l]);

                prefix >>= shift;

                pref_div = prefix/SIZE_BITS_UINT_8T;
                pref_mod = prefix%SIZE_BITS_UINT_8T;

                if (prefixes_in_partition[pref_div] & MASK_POWER_8[pref_mod]){
                    same_part = false;
                    same_part_pref = false;
                }
                else if (suffixes_in_partition[pref_div] & MASK_POWER_8[pref_mod]) same_part = false;
            }
            else {

                JHSG(PValue, prefix_judy, prefixes_iupac, nb_bytes_seed_iupac);

                if (PValue != NULL){
                    same_part = false;
                    same_part_pref = false;
                }
                else{

                    JHSG(PValue, suffix_judy, prefixes_iupac, nb_bytes_seed_iupac);
                    if (PValue != NULL) same_part = false;
                }
            }

            if (!suf_iupac){

                for (l = 0, suffix = 0; l < nb_bytes_seed; l++)
                    suffix = (suffix << SIZE_BITS_UINT_8T) | reverse_word_8(suffixes[l]);

                suffix >>= shift;

                suf_div = suffix/SIZE_BITS_UINT_8T;
                suf_mod = suffix%SIZE_BITS_UINT_8T;

                if (same_part){

                    if (prefixes_in_partition[suf_div] & MASK_POWER_8[suf_mod]) same_part = false;
                    else if (suffixes_in_partition[suf_div] & MASK_POWER_8[suf_mod]) same_part = false;
                    else if (!first_kmers_part && (memcmp(suffixes_iupac, prev_suffix_iupac, nb_bytes_seed_iupac) == 0)) same_part = false;
                }
            }
            else if (same_part){

                JHSG(PValue, prefix_judy, suffixes_iupac, nb_bytes_seed_iupac);

                if (PValue != NULL) same_part = false;
                else{

                    JHSG(PValue, suffix_judy, suffixes_iupac, nb_bytes_seed_iupac);

                    if (PValue != NULL) same_part = false;
                    else if (!first_kmers_part && (memcmp(suffixes_iupac, prev_suffix_iupac, nb_bytes_seed_iupac) == 0)) same_part = false;
                }
            }

            LAST_PART:

            if (!same_part){

                Rc_word = 0;
                J1F(Rc_int, pref_pos, Rc_word);

                while (Rc_int){
                    prefixes_in_partition[Rc_word] = 0;
                    J1N(Rc_int, pref_pos, Rc_word);
                }

                Rc_word = 0;
                J1F(Rc_int, suf_pos, Rc_word);

                while (Rc_int){
                    suffixes_in_partition[Rc_word] = 0;
                    J1N(Rc_int, suf_pos, Rc_word);
                }

                J1FA(Rc_word, pref_pos);
                J1FA(Rc_word, suf_pos);

                JHSFA(Rc_word, prefix_judy)
                JHSFA(Rc_word, suffix_judy)

                qsort(&suf_part[1], suf_part[0], sizeof(uint64_t), comp_uint64);

                if (same_part_pref){

                    for (z = 0, pref_iupac_tmp = 0; z < nb_bytes_seed_iupac; z++)
                        pref_iupac_tmp = (pref_iupac_tmp << 8) | prefixes_iupac[z];

                    pref_part[0]++;

                    pref_part = realloc(pref_part, (pref_part[0]+1) * sizeof(uint64_t));
                    ASSERT_NULL_PTR(pref_part, "compressKmers_from_KmerFiles_bis2()\n")

                    pref_part[pref_part[0]] = pref_iupac_tmp;
                }

                qsort(&pref_part[1], pref_part[0], sizeof(uint64_t), comp_uint64);

                if (list_part->count > 3){

                    min = 0xfffffffffffffff;

                    for (z = 0, curr_node = list_part->first; z < list_part->count - 1; z++){

                        pref_suf = curr_node->value;

                        if (is_intersecting_uint64_SIMD(pref_part, pref_suf[1], 0) <= 0){
                            if (is_intersecting_uint64_SIMD(pref_part, pref_suf[0], 0) <= 0){
                                if (is_intersecting_uint64_SIMD(suf_part, pref_suf[0], 0) <= 0){
                                    if (is_intersecting_uint64_SIMD(suf_part, pref_suf[1], 0) <= 0){
                                        min = 0;
                                        break;
                                    }
                                }
                            }
                        }

                        curr_node = curr_node->next;
                    }

                    if (min <= 0){

                        recycled_partition = current_partition - list_part->count + z;

                        fwrite(&recycled_partition, sizeof(uint32_t), 1, file_parts);
                        fwrite(suf_part, sizeof(uint32_t), 1, file_parts_count);

                        pref_suf[0] = realloc(pref_suf[0], (pref_part[0] + pref_suf[0][0] + 1) * sizeof(uint64_t));
                        ASSERT_NULL_PTR(pref_suf[0], "compressKmers_from_KmerFiles_bis2()\n")

                        memcpy(&(pref_suf[0][pref_suf[0][0]+1]), &pref_part[1], pref_part[0] * sizeof(uint64_t));

                        pref_suf[0][0] += pref_part[0];

                        free(pref_part);

                        pref_suf[1] = realloc(pref_suf[1], (suf_part[0] + pref_suf[1][0] + 1) * sizeof(uint64_t));
                        ASSERT_NULL_PTR(pref_suf[1], "compressKmers_from_KmerFiles_bis2()\n")

                        memcpy(&(pref_suf[1][pref_suf[1][0]+1]), &suf_part[1], suf_part[0] * sizeof(uint64_t));

                        pref_suf[1][0] += suf_part[0];

                        free(suf_part);

                        qsort(&(pref_suf[0][1]), pref_suf[0][0], sizeof(uint64_t), comp_uint64);
                        qsort(&(pref_suf[1][1]), pref_suf[1][0], sizeof(uint64_t), comp_uint64);

                        nb_recycled_part++;
                    }
                    else{
                        pref_suf = malloc(2 * sizeof(uint64_t*));
                        ASSERT_NULL_PTR(pref_suf, "compressKmers_from_KmerFiles_bis2()\n")

                        pref_suf[0] = pref_part;
                        pref_suf[1] = suf_part;

                        List_push(list_part, pref_suf);

                        if (list_part->count > MAX_SIZE_LIST_PART){
                            pref_suf = List_pop_first(list_part);
                            free(pref_suf[0]);
                            free(pref_suf[1]);
                            free(pref_suf);
                        }

                        fwrite(&current_partition, sizeof(uint32_t), 1, file_parts);
                        fwrite(suf_part, sizeof(uint32_t), 1, file_parts_count);

                        current_partition++;
                    }
                }
                else{
                    pref_suf = malloc(2 * sizeof(uint64_t*));
                    ASSERT_NULL_PTR(pref_suf, "compressKmers_from_KmerFiles_bis2()\n")

                    pref_suf[0] = pref_part;
                    pref_suf[1] = suf_part;

                    List_push(list_part, pref_suf);

                    if (list_part->count > MAX_SIZE_LIST_PART){
                        pref_suf = List_pop_first(list_part);
                        free(pref_suf[0]);
                        free(pref_suf[1]);
                        free(pref_suf);
                    }

                    fwrite(&current_partition, sizeof(uint32_t), 1, file_parts);
                    fwrite(suf_part, sizeof(uint32_t), 1, file_parts_count);

                    current_partition++;
                }

                first_kmers_part = true;
                pref_part = calloc(1, sizeof(uint64_t));
                suf_part = calloc(1, sizeof(uint64_t));
            }

            if (pref_iupac) JHSI(PValue, prefix_judy, prefixes_iupac, nb_bytes_seed_iupac)
            else {
                prefixes_in_partition[pref_div] |= MASK_POWER_8[pref_mod];
                J1S(Rc_int, pref_pos, pref_div);
            }

            if (!first_kmers_part){

                if (prev_suf_iupac) JHSI(PValue, suffix_judy, prev_suffix_iupac, nb_bytes_seed_iupac)
                else {
                    suffixes_in_partition[prev_suffix/SIZE_BITS_UINT_8T] |= MASK_POWER_8[prev_suffix%SIZE_BITS_UINT_8T];
                    J1S(Rc_int, suf_pos, prev_suffix/SIZE_BITS_UINT_8T);
                }
            }

            memcpy(prev_suffix_iupac, suffixes_iupac, nb_bytes_seed_iupac * sizeof(uint8_t));
            prev_suffix = suffix;
            prev_suf_iupac = suf_iupac;

            pref_part[0]++;

            pref_part = realloc(pref_part, (pref_part[0]+1) * sizeof(uint64_t));
            ASSERT_NULL_PTR(pref_part, "compressKmers_from_KmerFiles_bis2()\n")

            for (z = 0, pref_iupac_tmp = 0; z < nb_bytes_seed_iupac; z++)
                pref_iupac_tmp = (pref_iupac_tmp << 8) | prefixes_iupac[z];

            pref_part[pref_part[0]] = pref_iupac_tmp;

            suf_part[0]++;

            suf_part = realloc(suf_part, (suf_part[0]+1) * sizeof(uint64_t));
            ASSERT_NULL_PTR(suf_part, "compressKmers_from_KmerFiles_bis2()\n")

            for (z = 0, suf_iupac_tmp = 0; z < nb_bytes_seed_iupac; z++)
                suf_iupac_tmp = (suf_iupac_tmp << 8) | suffixes_iupac[z];

            suf_part[suf_part[0]] = suf_iupac_tmp;

            first_kmers_part = false;

            if (kmers_read % PRINT_EVERY_X_KMERS == 0){
                printf("%" PRIu64 " kmers read, %" PRIu64 " partition recycled, current_partition = %" PRIu32 "\n",
                       kmers_read, nb_recycled_part, current_partition);
            }
        }
    }

    if (last_kmers_part){
        same_part = false;
        last_kmers_part = false;
        goto LAST_PART;
    }

    fclose(file_kmers);
    fclose(file_parts);
    fclose(file_parts_count);

    J1FA(Rc_word, pref_pos);
    J1FA(Rc_word, suf_pos);

    JHSFA(Rc_word, prefix_judy)
    JHSFA(Rc_word, suffix_judy)

    while (list_part->count){
        pref_suf = List_pop(list_part);
        free(pref_suf[0]);
        free(pref_suf[1]);
        free(pref_suf);
    }

    List_destroy(list_part);

    free(prefixes_in_partition);
    free(suffixes_in_partition);

    free(seed_ASCII);

    free(kmers_iupac);

    free(suffixes_iupac);
    free(prefixes_iupac);
    free(prev_suffix_iupac);

    free(prefixes);
    free(suffixes);

    free(pref_part);
    free(suf_part);

    free(output_prefix_buffer);

    free(line);

    gettimeofday(&tval_after, NULL);
    time_spent(&tval_before, &tval_after, &tval_result);
    printf("\nElapsed time: %ld.%06ld s\n\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

    return current_partition;
}

uint32_t compressKmers_from_KmerFiles_bis3(char* output_prefix, int size_seed, int size_kmer, uint32_t prev_nb_parts){

    ASSERT_NULL_PTR(output_prefix,"compressKmers_from_KmerFiles_bis2()")

    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    Pvoid_t partition_iupac = (PWord_t)NULL;
    Pvoid_t partition_iupac_tmp;
    Pvoid_t pref_suf_pos = (PWord_t)NULL;

    PWord_t PValue;
    PWord_t PValue2;
    Word_t Rc_word;
    int Rc_int;

    FILE* file_kmers;
    FILE* file_parts;
    FILE* file_parts_count;

    List* list_part = List_create();

    ListNode* curr_node = NULL;

    int i = 0, j = 0, z = 0, it_list = 0;
    int shift_pref_suf = 0;

    int len_output_prefix = strlen(output_prefix);

    int nb_bytes_seed = CEIL(size_seed * 2, SIZE_BITS_UINT_8T);

    int nb_bytes_kmer_iupac = CEIL(size_kmer * 4, SIZE_BITS_UINT_8T);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer_iupac;

    int pos_start_comp_suffix = (size_kmer - size_seed) / 2;
    int size_after_pos_start_comp_suffix = size_kmer - pos_start_comp_suffix * 2;

    int nb_combi = pow(4,size_seed);
    int size_annot_seeds = CEIL(nb_combi, SIZE_BITS_UINT_8T);

    uint64_t pref_div;
    uint64_t suf_div;

    uint64_t nb_recycled_part = 0;

    uint64_t kmers_read = 0;

    uint32_t current_partition = prev_nb_parts;
    uint32_t recycled_partition = 0;
    uint32_t nb_kmer_part = 0;

    uint32_t prev_suffix;
    uint32_t suffix;
    uint32_t prefix;

    uint8_t pref_mod;
    uint8_t suf_mod;

    bool overlap;
    bool prev_suf_iupac;
    bool suf_iupac;
    bool pref_iupac;

    bool same_part = true;

    bool first_kmers_part = true;
    bool last_kmers_part = true;

    void** partition_all = NULL;

    uint8_t* partition_tmp;

    uint8_t* kmers_iupac = malloc(nb_kmer_in_buf * nb_bytes_kmer_iupac * sizeof(uint8_t));
    ASSERT_NULL_PTR(kmers_iupac, "compressKmers_from_KmerFiles()")

    uint8_t* prefixes = calloc(nb_bytes_seed, sizeof(uint8_t));
    ASSERT_NULL_PTR(prefixes, "compressKmers_from_KmerFiles()")

    uint8_t* suffixes = calloc(nb_bytes_seed, sizeof(uint8_t));
    ASSERT_NULL_PTR(suffixes, "compressKmers_from_KmerFiles()")

    uint8_t* partition = calloc(size_annot_seeds, sizeof(uint8_t));
    ASSERT_NULL_PTR(partition, "compressKmers_from_KmerFiles()")

    char* prefixes_iupac = calloc(size_seed + 1, sizeof(char));
    ASSERT_NULL_PTR(prefixes_iupac, "compressKmers_from_KmerFiles()")

    char* suffixes_iupac = calloc(size_seed + 1, sizeof(char));
    ASSERT_NULL_PTR(suffixes_iupac, "compressKmers_from_KmerFiles()")

    char* prev_suffix_iupac = calloc(size_seed + 1, sizeof(char));
    ASSERT_NULL_PTR(prev_suffix_iupac, "compressKmers_from_KmerFiles()")

    char* it_iupac = calloc(size_seed + 1, sizeof(char));
    ASSERT_NULL_PTR(it_iupac, "compressKmers_from_KmerFiles()")

    char* line = calloc(size_kmer + 1, sizeof(char));
    ASSERT_NULL_PTR(line,"compressKmers_from_KmerFiles()")

    char* output_prefix_buffer = malloc((len_output_prefix + 50) * sizeof(char));
    ASSERT_NULL_PTR(output_prefix_buffer, "min_simReads_max_simkmers()")

    strcpy(output_prefix_buffer, output_prefix);

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_kmers_tmp");
    file_kmers = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_kmers,"compressKmers_from_KmerFiles()")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_partitions");
    file_parts = fopen(output_prefix_buffer, "wb");
    ASSERT_NULL_PTR(file_parts,"compressKmers_from_KmerFiles()")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_parts_count_tmp");
    file_parts_count = fopen(output_prefix_buffer, "wb");
    ASSERT_NULL_PTR(file_parts_count,"compressKmers_from_KmerFiles()")

    if (size_seed % 4) shift_pref_suf = SIZE_BITS_UINT_8T - ((size_seed % 4) * 2);

    while ((j = fread(kmers_iupac, nb_bytes_kmer_iupac * sizeof(uint8_t), nb_kmer_in_buf, file_kmers))){

        for (i = 0; i < j; i++, kmers_read++){

            memset(prefixes, 0, nb_bytes_seed * sizeof(uint8_t));
            memset(suffixes, 0, nb_bytes_seed * sizeof(uint8_t));

            // prefix
            kmer_iupac_comp_to_ascii(&kmers_iupac[i * nb_bytes_kmer_iupac], size_seed, prefixes_iupac);

            if (is_substring_IUPAC(prefixes_iupac)) pref_iupac = true;
            else pref_iupac = false;

            // suffix
            kmer_iupac_comp_to_ascii(&kmers_iupac[i * nb_bytes_kmer_iupac + pos_start_comp_suffix],
                                     size_after_pos_start_comp_suffix, line);

            strcpy(suffixes_iupac, &line[size_after_pos_start_comp_suffix - size_seed]);

            if (is_substring_IUPAC(suffixes_iupac)) suf_iupac = true;
            else suf_iupac = false;

            same_part = true;

            if (!pref_iupac){

                parseKmerCount(prefixes_iupac, size_seed, prefixes, 0);

                for (z = 0, prefix = 0; z < nb_bytes_seed; z++)
                    prefix = (prefix << SIZE_BITS_UINT_8T) | reverse_word_8(prefixes[z]);

                prefix >>= shift_pref_suf;

                pref_div = prefix/SIZE_BITS_UINT_8T;
                pref_mod = prefix%SIZE_BITS_UINT_8T;

                if (partition[pref_div] & MASK_POWER_8[pref_mod]) same_part = false;
            }
            else {

                JSLG(PValue, partition_iupac, (uint8_t*) prefixes_iupac);

                if (PValue != NULL) same_part = false;
            }

            if (!suf_iupac){

                parseKmerCount(suffixes_iupac, size_seed, suffixes, 0);

                for (z = 0, suffix = 0; z < nb_bytes_seed; z++)
                    suffix = (suffix << SIZE_BITS_UINT_8T) | reverse_word_8(suffixes[z]);

                suffix >>= shift_pref_suf;

                suf_div = suffix/SIZE_BITS_UINT_8T;
                suf_mod = suffix%SIZE_BITS_UINT_8T;

                if (same_part){

                    if (partition[suf_div] & MASK_POWER_8[suf_mod]) same_part = false;
                    else if (!first_kmers_part && (strcmp(suffixes_iupac, prev_suffix_iupac) == 0)) same_part = false;
                }
            }
            else if (same_part){

                JSLG(PValue, partition_iupac, (uint8_t*) suffixes_iupac);

                if (PValue != NULL) same_part = false;
                else if (!first_kmers_part && (strcmp(suffixes_iupac, prev_suffix_iupac) == 0)) same_part = false;
            }

            LAST_PART:

            if (!same_part){

                if (prev_suf_iupac) JSLI(PValue, partition_iupac, (uint8_t*) prev_suffix_iupac)
                else {
                    partition[prev_suffix/SIZE_BITS_UINT_8T] |= MASK_POWER_8[prev_suffix%SIZE_BITS_UINT_8T];
                    J1S(Rc_int, pref_suf_pos, prev_suffix/SIZE_BITS_UINT_8T);
                }

                if (pref_iupac) JSLI(PValue, partition_iupac, (uint8_t*) prefixes_iupac)
                else {
                    partition[pref_div] |= MASK_POWER_8[pref_mod];
                    J1S(Rc_int, pref_suf_pos, pref_div);
                }

                if (list_part->count > 3){

                    recycled_partition = current_partition - list_part->count;

                    for (curr_node = list_part->first, it_list = 0; it_list < list_part->count - 1; curr_node = curr_node->next, it_list++){

                        overlap = false;

                        partition_tmp = (uint8_t*) ((void**) curr_node->value)[0];

                        Rc_word = 0;
                        J1F(Rc_int, pref_suf_pos, Rc_word);

                        while (Rc_int){

                            if (partition_tmp[Rc_word] & partition[Rc_word]){
                                overlap = true;
                                break;
                            }

                            J1N(Rc_int, pref_suf_pos, Rc_word);
                        }

                        if (!overlap){

                            partition_iupac_tmp = (Pvoid_t) ((void**) curr_node->value)[1];

                            memset(it_iupac, 0, size_seed * sizeof(uint8_t));
                            it_iupac[size_seed] = '\0';

                            JSLF(PValue, partition_iupac, (uint8_t*) it_iupac);

                            while (PValue != NULL){

                                JSLG(PValue2, partition_iupac_tmp, (uint8_t*) it_iupac);

                                if (PValue2 != NULL){
                                    overlap = true;
                                    break;
                                }

                                JSLN(PValue, partition_iupac, (uint8_t*) it_iupac);
                            }

                            if (!overlap) break;
                        }
                    }

                    if (!overlap){

                        recycled_partition = current_partition - list_part->count + it_list;

                        fwrite(&recycled_partition, sizeof(uint32_t), 1, file_parts);
                        fwrite(&nb_kmer_part, sizeof(uint32_t), 1, file_parts_count);

                        nb_recycled_part++;

                        Rc_word = 0;
                        J1F(Rc_int, pref_suf_pos, Rc_word);

                        while (Rc_int){

                            partition_tmp[Rc_word] |= partition[Rc_word];
                            J1N(Rc_int, pref_suf_pos, Rc_word);
                        }

                        memset(it_iupac, 0, size_seed * sizeof(uint8_t));
                        it_iupac[size_seed] = '\0';

                        JSLF(PValue, partition_iupac, (uint8_t*) it_iupac);

                        while (PValue != NULL){

                            JSLI(PValue2, partition_iupac_tmp, (uint8_t*) it_iupac);
                            JSLN(PValue, partition_iupac, (uint8_t*) it_iupac);
                        }

                        if (pref_iupac) JSLI(PValue2, partition_iupac_tmp, (uint8_t*) prefixes_iupac)
                        else partition_tmp[pref_div] |= MASK_POWER_8[pref_mod];

                        ((void**) curr_node->value)[1] = (void*) partition_iupac_tmp;

                        free(partition);

                        JSLFA(Rc_word, partition_iupac);
                    }
                    else{

                        partition_all = malloc(2 * sizeof(void*));
                        ASSERT_NULL_PTR(partition_all, "compressKmers_from_KmerFiles()")

                        partition_all[0] = (void*) partition;
                        partition_all[1] = (void*) partition_iupac;

                        List_push(list_part, partition_all);

                        if (list_part->count > MAX_SIZE_LIST_PART){

                            partition_all = List_pop_first(list_part);

                            free(partition_all[0]);

                            partition_iupac_tmp = (Pvoid_t) partition_all[1];
                            JSLFA(Rc_word, partition_iupac_tmp);

                            free(partition_all);
                        }

                        fwrite(&current_partition, sizeof(uint32_t), 1, file_parts);
                        fwrite(&nb_kmer_part, sizeof(uint32_t), 1, file_parts_count);

                        current_partition++;
                    }
                }
                else{

                    partition_all = malloc(2 * sizeof(void*));
                    ASSERT_NULL_PTR(partition_all, "compressKmers_from_KmerFiles()")

                    partition_all[0] = (void*) partition;
                    partition_all[1] = (void*) partition_iupac;

                    List_push(list_part, partition_all);

                    if (list_part->count > MAX_SIZE_LIST_PART){

                        partition_all = List_pop_first(list_part);

                        free(partition_all[0]);

                        partition_iupac_tmp = (Pvoid_t) partition_all[1];
                        JSLFA(Rc_word, partition_iupac_tmp);

                        free(partition_all);
                    }

                    fwrite(&current_partition, sizeof(uint32_t), 1, file_parts);
                    fwrite(&nb_kmer_part, sizeof(uint32_t), 1, file_parts_count);

                    current_partition++;
                }

                J1FA(Rc_word, pref_suf_pos);

                partition = calloc(size_annot_seeds, sizeof(uint8_t));
                ASSERT_NULL_PTR(partition, "compressKmers_from_KmerFiles()")

                partition_iupac = NULL;

                first_kmers_part = true;

                nb_kmer_part = 0;
            }

            if (pref_iupac) JSLI(PValue, partition_iupac, (uint8_t*) prefixes_iupac)
            else {
                partition[pref_div] |= MASK_POWER_8[pref_mod];
                J1S(Rc_int, pref_suf_pos, pref_div);
            }

            if (!first_kmers_part){

                if (prev_suf_iupac) JSLI(PValue, partition_iupac, (uint8_t*) prev_suffix_iupac)
                else {
                    partition[prev_suffix/SIZE_BITS_UINT_8T] |= MASK_POWER_8[prev_suffix%SIZE_BITS_UINT_8T];
                    J1S(Rc_int, pref_suf_pos, prev_suffix/SIZE_BITS_UINT_8T);
                }
            }

            strcpy(prev_suffix_iupac, suffixes_iupac);
            prev_suffix = suffix;
            prev_suf_iupac = suf_iupac;

            first_kmers_part = false;

            nb_kmer_part++;

            if (kmers_read % PRINT_EVERY_X_KMERS == 0){
                printf("%" PRIu64 " kmers read, %" PRIu64 " partition recycled, current_partition = %" PRIu32 "\n",
                       kmers_read, nb_recycled_part, current_partition);
            }
        }
    }

    if (last_kmers_part){
        same_part = false;
        last_kmers_part = false;
        goto LAST_PART;
    }

    fclose(file_kmers);
    fclose(file_parts);
    fclose(file_parts_count);

    J1FA(Rc_word, pref_suf_pos);

    JSLFA(Rc_word, partition_iupac);

    while (list_part->count){

        partition_all = List_pop(list_part);

        free(partition_all[0]);

        partition_iupac_tmp = (Pvoid_t) partition_all[1];
        JSLFA(Rc_word, partition_iupac_tmp);

        free(partition_all);
    }

    List_destroy(list_part);

    free(partition);

    free(kmers_iupac);

    free(suffixes_iupac);
    free(prefixes_iupac);
    free(prev_suffix_iupac);
    free(it_iupac);

    free(prefixes);
    free(suffixes);

    free(output_prefix_buffer);

    free(line);

    gettimeofday(&tval_after, NULL);
    time_spent(&tval_before, &tval_after, &tval_result);
    printf("\nElapsed time: %ld.%06ld s\n\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

    return current_partition;
}

void insert_kmers_list(List* list_kmers, BFT_annotation* bft_annot, BFT_Root* root_iupac, int lvl_root_iupac,
                       BFT_Root* root_no_iupac, int lvl_root_no_iupac, uint32_t part, uint32_t size_part_bytes){

    resultPresence* res;
    UC* uc;

    uint8_t* kmer;

    while ((kmer = List_pop_first(list_kmers)) != NULL){

        if (kmer[0]){

            res = isKmerPresent(&(root_iupac->node), root_iupac, lvl_root_iupac, &(kmer[1]), root_iupac->k);

            if (res->posFilter2 != 0) uc = (UC*)res->container;
            else uc = &(((UC*)((CC*)res->container)->children)[res->bucket]);

            get_annot(uc, &bft_annot->annot, &bft_annot->annot_ext, &bft_annot->annot_cplx, &bft_annot->size_annot,
                      &bft_annot->size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);

            insert_partition(res, bft_annot, root_iupac, uc, part, size_part_bytes);
        }
        else{

            res = isKmerPresent(&(root_no_iupac->node), root_no_iupac, lvl_root_no_iupac, &(kmer[1]), root_no_iupac->k);

            if (res->posFilter2 != 0) uc = (UC*)res->container;
            else uc = &(((UC*)((CC*)res->container)->children)[res->bucket]);

            get_annot(uc, &bft_annot->annot, &bft_annot->annot_ext, &bft_annot->annot_cplx, &bft_annot->size_annot,
                      &bft_annot->size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);

            insert_partition(res, bft_annot, root_no_iupac, uc, part, size_part_bytes);
        }

        free(res);
        free(kmer);
    }

    return;
}

void insert_partition(resultPresence* res, BFT_annotation* bft_annot, BFT* root, UC* uc, uint32_t part, uint32_t size_part_bytes){

    uint32_t* list_id_genomes = get_list_id_genomes(bft_annot, root);

    uint32_t pow2_part = 0;
    uint32_t k = 1;

    uint32_t current_part;

    bool part_to_insert = true;

    if (part < list_id_genomes[list_id_genomes[0]]){

        if (bft_annot->annot[0] & 0x2){

            uint32_t i = 0, pos = 0;

            for (; k <= list_id_genomes[0]; k++)
                if (part < list_id_genomes[k]) break;

            while ((i < bft_annot->size_annot) && (bft_annot->annot[i] & 0x2)){

                if (pos == k-1) break;

                i++;
                pos++;
                while ((i < bft_annot->size_annot) && (bft_annot->annot[i] & 0x1)) i++;
            }

            memset(&bft_annot->annot[i], 0, (bft_annot->size_annot - i) * sizeof(uint8_t));
            bft_annot->size_annot = i;
        }
        else{
            memset(bft_annot->annot, 0, bft_annot->size_annot * sizeof(uint8_t));
            bft_annot->size_annot = 0;
        }

        memset(bft_annot->annot_cplx, 0, bft_annot->size_annot_cplx * sizeof(uint8_t));
        bft_annot->size_annot_cplx = 0;

        if (bft_annot->annot_ext != NULL)
            delete_extend_annots(uc, res->posFilter2, res->posFilter3, res->pos_sub_bucket, res->pos_sub_bucket, 0, 0, 1);

        for (; k <= list_id_genomes[0]; k++){

            get_annot(uc, &bft_annot->annot, &bft_annot->annot_ext, &bft_annot->annot_cplx, &bft_annot->size_annot,
                      &bft_annot->size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);

            if (part_to_insert && (part < list_id_genomes[k])){
                current_part = part;
                part_to_insert = false;
                k--;
            }
            else current_part = list_id_genomes[k];

            if (current_part >= pow2_part){
                pow2_part = round_up_next_highest_power2(current_part);
                size_part_bytes = get_nb_bytes_power2_annot_bis(current_part, pow2_part);
            }

            compute_best_mode(root->ann_inf, root->comp_set_colors, bft_annot->annot, bft_annot->size_annot, bft_annot->annot_ext,
                              1, current_part, size_part_bytes);

            if (root->ann_inf->last_added != current_part){

                if (root->ann_inf->min_size > uc->size_annot){
                    if ((bft_annot->annot_ext == NULL) || (root->ann_inf->min_size > uc->size_annot+1)){
                        bft_annot->annot_ext = realloc_annotation(uc, res->posFilter2, res->posFilter3, root->ann_inf->min_size,
                                                                  0, res->pos_sub_bucket);
                    }
                }

                modify_mode_annotation(root->ann_inf, &(uc->suffixes[res->pos_sub_bucket * (res->posFilter2 + uc->size_annot) + res->posFilter2]),
                                       uc->size_annot, bft_annot->annot_ext, 1, current_part, size_part_bytes);

                if ((bft_annot->annot_ext != NULL) && (bft_annot->annot_ext[0] == 0))
                    delete_extend_annots(uc, res->posFilter2, res->posFilter3, res->pos_sub_bucket, res->pos_sub_bucket, 0, 0, 1);
            }

            reinit_annotation_inform(root->ann_inf);
        }
    }
    else{

        compute_best_mode(root->ann_inf, root->comp_set_colors, bft_annot->annot, bft_annot->size_annot, bft_annot->annot_ext, 1,
                          part, size_part_bytes);

        if (root->ann_inf->last_added != part){

            if (root->ann_inf->min_size > uc->size_annot){
                if ((bft_annot->annot_ext == NULL) || (root->ann_inf->min_size > uc->size_annot+1)){
                    bft_annot->annot_ext = realloc_annotation(uc, res->posFilter2, res->posFilter3, root->ann_inf->min_size,
                                                              0, res->pos_sub_bucket);
                }
            }

            modify_mode_annotation(root->ann_inf, &(uc->suffixes[res->pos_sub_bucket * (res->posFilter2 + uc->size_annot) + res->posFilter2]),
                                   uc->size_annot, bft_annot->annot_ext, 1, part, size_part_bytes);

            if ((bft_annot->annot_ext != NULL) && (bft_annot->annot_ext[0] == 0))
                delete_extend_annots(uc, res->posFilter2, res->posFilter3, res->pos_sub_bucket, res->pos_sub_bucket, 0, 0, 1);
        }

        reinit_annotation_inform(root->ann_inf);
    }

    free(list_id_genomes);

    return;
}

uint32_t* recycle_sub_paths(char* seed, uint32_t recycled_part, Pvoid_t* recycl_subpaths){

    ASSERT_NULL_PTR(seed, "recycle_sub_paths()\n")
    ASSERT_NULL_PTR(recycl_subpaths, "recycle_sub_paths()\n")

    uint32_t* array_pref_part_recycl;

    PWord_t PValue;

    JSLG(PValue, *recycl_subpaths, (uint8_t*) seed);

    if (PValue == PJERR) ERROR("recycle_sub_paths(): Something went wrong with the JudyArray... \n")

    if (PValue == NULL){

        array_pref_part_recycl = malloc(2 * sizeof(uint32_t));
        ASSERT_NULL_PTR(array_pref_part_recycl, "load_part_recycl()\n")

        array_pref_part_recycl[0] = 1;
        array_pref_part_recycl[1] = recycled_part;

        JSLI(PValue, *recycl_subpaths, (uint8_t*) seed);
    }
    else{

        array_pref_part_recycl = (uint32_t*) *PValue;

        array_pref_part_recycl[0]++;

        array_pref_part_recycl = realloc(array_pref_part_recycl, (array_pref_part_recycl[0] + 1) * sizeof(uint32_t));
        ASSERT_NULL_PTR(array_pref_part_recycl, "load_part_recycl()\n")

        array_pref_part_recycl[array_pref_part_recycl[0]] = recycled_part;
    }

    *PValue = (Word_t) array_pref_part_recycl;

    return &array_pref_part_recycl[array_pref_part_recycl[0]];
}

void serialize_subpaths_recycling(char* filename_subpaths_recycling, Pvoid_t* recycl_subpaths, int size_seed){

    ASSERT_NULL_PTR(filename_subpaths_recycling, "serialize_subpaths_recycling()\n")
    ASSERT_NULL_PTR(recycl_subpaths, "serialize_subpaths_recycling()\n")

    PWord_t PValue;

    int i;

    int seed_shift = 0;
    int seed_shift_iupac = 0;

    int nb_bytes_seed = CEIL(size_seed * 2, SIZE_BITS_UINT_8T);
    int nb_bytes_seed_iupac = CEIL(size_seed * 4, SIZE_BITS_UINT_8T);

    long int pos_after_no_iupac_seeds;

    uint32_t seed;

    uint32_t size_seed_present = CEIL((uint32_t) pow(4, size_seed), SIZE_BITS_UINT_8T);

    uint32_t* recycled_part;

    uint64_t seed_iupac;
    uint64_t seed_iupac_tmp;
    uint64_t prev_seed_iupac = 0;

    uint8_t* seed_comp = malloc(nb_bytes_seed_iupac * sizeof(uint8_t));
    ASSERT_NULL_PTR(seed_comp, "serialize_subpaths_recycling()\n")

    uint8_t* seed_ascii = malloc((size_seed + 1) * sizeof(uint8_t));
    ASSERT_NULL_PTR(seed_ascii, "serialize_subpaths_recycling()\n")

    uint8_t* seed_present = calloc(size_seed_present, sizeof(uint8_t));
    ASSERT_NULL_PTR(seed_present, "serialize_subpaths_recycling()\n")

    FILE* file_subpaths_recycling = fopen(filename_subpaths_recycling, "wb");
    ASSERT_NULL_PTR(file_subpaths_recycling, "serialize_subpaths_recycling()\n")

    if (size_seed % 4) seed_shift = SIZE_BITS_UINT_8T - ((size_seed % 4) * 2);
    if (size_seed % 2) seed_shift_iupac = 4;

    fseek(file_subpaths_recycling, size_seed_present, SEEK_SET);

    seed_ascii[0] = '\0';

    JSLF(PValue, *recycl_subpaths, seed_ascii);

    while (PValue != NULL){

        if (is_substring_IUPAC((char*) seed_ascii) == 0){

            memset(seed_comp, 0, nb_bytes_seed * sizeof(uint8_t));
            parseKmerCount((char*)seed_ascii, size_seed, seed_comp, 0);

            for (i = 0, seed = 0; i < nb_bytes_seed; i++)
                seed = (seed << SIZE_BITS_UINT_8T) | reverse_word_8(seed_comp[i]);

            seed >>= seed_shift;

            seed_present[seed/SIZE_BITS_UINT_8T] |= MASK_POWER_8[seed%SIZE_BITS_UINT_8T];

            recycled_part = (uint32_t*) *PValue;

            fwrite(recycled_part, sizeof(uint32_t), recycled_part[0]+1, file_subpaths_recycling);
        }

        JSLN(PValue, *recycl_subpaths, seed_ascii);
    }

    pos_after_no_iupac_seeds = ftell(file_subpaths_recycling);

    rewind(file_subpaths_recycling);

    fwrite(seed_present, sizeof(uint8_t), size_seed_present, file_subpaths_recycling);

    fseek(file_subpaths_recycling, pos_after_no_iupac_seeds, SEEK_SET);

    free(seed_present);

    seed_ascii[0] = '\0';

    JSLF(PValue, *recycl_subpaths, seed_ascii);

    while (PValue != NULL){

        if (is_substring_IUPAC((char*) seed_ascii)){

            memset(seed_comp, 0, nb_bytes_seed_iupac * sizeof(uint8_t));
            parseKmerCount_IUPAC((char*)seed_ascii, size_seed, seed_comp, 0);

            for (i = 0, seed_iupac = 0; i < nb_bytes_seed_iupac; i++)
                seed_iupac = (seed_iupac << SIZE_BITS_UINT_8T) | reverse_word_8(seed_comp[i]);

            seed_iupac >>= seed_shift_iupac;
            seed_iupac_tmp = seed_iupac - prev_seed_iupac;
            prev_seed_iupac = seed_iupac;

            fwrite(&seed_iupac_tmp, sizeof(uint64_t), 1, file_subpaths_recycling);

            recycled_part = (uint32_t*) *PValue;

            fwrite(recycled_part, sizeof(uint32_t), recycled_part[0]+1, file_subpaths_recycling);
        }

        JSLN(PValue, *recycl_subpaths, seed_ascii);
    }

    free(seed_comp);
    free(seed_ascii);

    fclose(file_subpaths_recycling);

    return;
}

void insert_kmers_partitions2(char* output_prefix, int size_seed, int size_kmer, bool recycle_paths,
                              BFT_Root* root_no_IUPAC, BFT_Root* root_IUPAC){

    ASSERT_NULL_PTR(output_prefix,"compressKmers_from_KmerFiles()")

    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    Pvoid_t recycled_parts = (Pvoid_t) NULL;
    Pvoid_t recycl_subpaths = (Pvoid_t) NULL;
    PWord_t PValue;
    Word_t Rc_word;
    int Rc_int;

    FILE* file_kmers;
    FILE* file_parts;
    FILE* file_parts_count;
    FILE* file_recycl_parts;

    resultPresence* res;

    UC* uc;

    BFT_Root* root;

    int i = 0, j = 0, k = 0, m = 0;

    int min_count_l_recycl_kmers = 1;

    int len_output_prefix = strlen(output_prefix);

    int nb_bytes_kmer = CEIL(size_kmer * 2, SIZE_BITS_UINT_8T);
    int nb_bytes_kmer_iupac = CEIL(size_kmer * 4, SIZE_BITS_UINT_8T);

    int nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer_iupac;

    int lvl_root;
    int lvl_root_no_iupac = (size_kmer / NB_CHAR_SUF_PREF) - 1;
    int lvl_root_iupac = (size_kmer * 2 / NB_CHAR_SUF_PREF) - 1;

    int nb_cell_to_delete;
    int size_suffix;
    int current_level;
    int nb_cell;

    uint64_t kmers_read = 0;

    bool kmer_is_iupac;

    char* seed = malloc((size_seed + 1) * sizeof(char));
    ASSERT_NULL_PTR(seed, "compressKmers_from_KmerFiles()")

    char* seed_tmp = malloc((size_seed + 1) * sizeof(char));
    ASSERT_NULL_PTR(seed_tmp, "compressKmers_from_KmerFiles()")

    char* ascii_kmer = malloc((size_kmer + 1) * sizeof(char));
    ASSERT_NULL_PTR(ascii_kmer, "compressKmers_from_KmerFiles()")

    char* line = calloc(size_kmer + 1, sizeof(char));
    ASSERT_NULL_PTR(line, "compressKmers_from_KmerFiles()")

    uint8_t* kmers_iupac = malloc(nb_kmer_in_buf * nb_bytes_kmer_iupac * sizeof(uint8_t));
    ASSERT_NULL_PTR(kmers_iupac, "compressKmers_from_KmerFiles()")

    uint8_t* kmer_no_iupac = calloc(nb_bytes_kmer, sizeof(uint8_t));
    ASSERT_NULL_PTR(kmer_no_iupac, "compressKmers_from_KmerFiles()")

    uint8_t* substring = calloc(nb_bytes_kmer_iupac, sizeof(uint8_t));
    ASSERT_NULL_PTR(substring, "compressKmers_from_KmerFiles()")

    char* output_prefix_buffer = malloc((len_output_prefix + 50) * sizeof(char));
    ASSERT_NULL_PTR(output_prefix_buffer, "min_simReads_max_simkmers()")

    //uint8_t mask_seed_iupac = SIZE_BITS_UINT_8T - ((size_seed % 2) * 4);
    //uint8_t mask_seed_no_iupac = SIZE_BITS_UINT_8T - ((size_seed % 4) * 2);

    uint8_t* kmer;
    uint8_t* kmer_tmp;

    uint32_t* inter_part = NULL;
    uint32_t* inter_part_tmp = NULL;
    uint32_t* inter_part_new_node = NULL;

    uint32_t tmp;

    uint32_t nb_kmer_curr_part = 0;
    uint32_t count_nb_kmer_curr_part = 0xfffffffe;

    uint32_t parts_to_insert;

    uint32_t delta_part = 0xffffffff;
    uint32_t delta_part_tmp = 0;
    uint32_t prev_delta_part = 0;

    uint32_t prev_part_to_insert = 0xffffffff;
    uint32_t part_to_insert = 0xffffffff;
    uint32_t prev_size_part_to_insert = 0xffffffff;
    uint32_t size_part_to_insert = 0xffffffff;

    uint32_t nb_kmers_recycled;

    List* l_recycl_kmers = List_create();

    if ((root_no_IUPAC == NULL) && (root_IUPAC == NULL)){
        root_no_IUPAC = createBFT_Root(size_kmer, 1, 1);
        root_IUPAC = createBFT_Root(size_kmer * 2, 1, 1);

        root_no_IUPAC->ann_inf = create_annotation_inform(-1);
        root_IUPAC->ann_inf = create_annotation_inform(-1);

        for (i = 0; i < size_kmer / NB_CHAR_SUF_PREF; i++){
            root_no_IUPAC->info_per_lvl[i].nb_ucs_skp = 8;
            root_no_IUPAC->info_per_lvl[i].nb_bits_per_cell_skip_filter2 = root_no_IUPAC->info_per_lvl[i].nb_ucs_skp;
            root_no_IUPAC->info_per_lvl[i].nb_bits_per_cell_skip_filter3 = root_no_IUPAC->info_per_lvl[i].nb_ucs_skp;
            root_no_IUPAC->info_per_lvl[i].nb_bytes_per_cell_skip_filter2 = CEIL(root_no_IUPAC->info_per_lvl[i].nb_ucs_skp, SIZE_BITS_UINT_8T);
            root_no_IUPAC->info_per_lvl[i].nb_bytes_per_cell_skip_filter3 = CEIL(root_no_IUPAC->info_per_lvl[i].nb_ucs_skp, SIZE_BITS_UINT_8T);
        }

        for (i = 0; i < (size_kmer * 2) / NB_CHAR_SUF_PREF; i++){
            root_IUPAC->info_per_lvl[i].nb_ucs_skp = 8;
            root_IUPAC->info_per_lvl[i].nb_bits_per_cell_skip_filter2 = root_IUPAC->info_per_lvl[i].nb_ucs_skp;
            root_IUPAC->info_per_lvl[i].nb_bits_per_cell_skip_filter3 = root_IUPAC->info_per_lvl[i].nb_ucs_skp;
            root_IUPAC->info_per_lvl[i].nb_bytes_per_cell_skip_filter2 = CEIL(root_IUPAC->info_per_lvl[i].nb_ucs_skp, SIZE_BITS_UINT_8T);
            root_IUPAC->info_per_lvl[i].nb_bytes_per_cell_skip_filter3 = CEIL(root_IUPAC->info_per_lvl[i].nb_ucs_skp, SIZE_BITS_UINT_8T);
        }
    }

    BFT_annotation* bft_annot = create_BFT_annotation();
    bft_annot->from_BFT = 1;

    strcpy(output_prefix_buffer, output_prefix);

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_kmers_tmp");
    file_kmers = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_kmers,"compressKmers_from_KmerFiles()")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_partitions");
    file_parts = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_parts,"compressKmers_from_KmerFiles()")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_parts_count_tmp");
    file_parts_count = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_parts_count,"compressKmers_from_KmerFiles()")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_recycl_parts1");
    file_recycl_parts = fopen(output_prefix_buffer, "wb");
    ASSERT_NULL_PTR(file_recycl_parts,"compressKmers_from_KmerFiles()")

    //mask_seed_iupac = 0xff >> (mask_seed_iupac == SIZE_BITS_UINT_8T ? 0 : mask_seed_iupac);
    //mask_seed_no_iupac = 0xff >> (mask_seed_no_iupac == SIZE_BITS_UINT_8T ? 0 : mask_seed_no_iupac);

    while ((j = fread(kmers_iupac, nb_bytes_kmer_iupac * sizeof(uint8_t), nb_kmer_in_buf, file_kmers))){

        for (i = 0; i < j; i++){

            count_nb_kmer_curr_part++;

            kmer_iupac_comp_to_ascii(&kmers_iupac[i * nb_bytes_kmer_iupac], size_kmer, line);

            if (is_substring_IUPAC(line)){

                kmer = &kmers_iupac[i * nb_bytes_kmer_iupac];

                kmer_is_iupac = true;

                root = root_IUPAC;
                lvl_root = lvl_root_iupac;
            }
            else {

                memset(kmer_no_iupac, 0, nb_bytes_kmer * sizeof(uint8_t));
                parseKmerCount(line, size_kmer, kmer_no_iupac, 0);

                kmer = kmer_no_iupac;

                kmer_is_iupac = false;

                root = root_no_IUPAC;
                lvl_root = lvl_root_no_iupac;
            }

            res = isKmerPresent(&(root->node), root, lvl_root, kmer, root->k);

            if (res->link_child != NULL){

                if (res->posFilter2 != 0) uc = (UC*)res->container;
                else uc = &(((UC*)((CC*)res->container)->children)[res->bucket]);

                get_annot(uc, &bft_annot->annot, &bft_annot->annot_ext, &bft_annot->annot_cplx, &bft_annot->size_annot,
                          &bft_annot->size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);

                if (count_nb_kmer_curr_part > nb_kmer_curr_part){

                    count_nb_kmer_curr_part = 1;

                    if (fread(&parts_to_insert, sizeof(uint32_t), 1, file_parts) != 1)
                        ERROR("k-mer file and partition file don't have the same nb. of lines.\n");

                    if (fread(&nb_kmer_curr_part, sizeof(uint32_t), 1, file_parts_count) != 1)
                        ERROR("k-mer file and partition file don't have the same nb. of lines.\n");

                    prev_size_part_to_insert = size_part_to_insert;
                    size_part_to_insert = get_nb_bytes_power2_annot(parts_to_insert);

                    insert_partition(res, bft_annot, root, uc, parts_to_insert, size_part_to_insert);

                    if (recycle_paths){

                        if (l_recycl_kmers->count < min_count_l_recycl_kmers){

                            insert_kmers_list(l_recycl_kmers, bft_annot, root_IUPAC, lvl_root_iupac, root_no_IUPAC,
                                              lvl_root_no_iupac, part_to_insert, prev_size_part_to_insert);
                        }
                        else if (inter_part != NULL){

                            nb_kmers_recycled = (uint32_t) l_recycl_kmers->count;

                            for (k = 1; k <= inter_part[0]; k++){
                                J1T(Rc_int, recycled_parts, inter_part[k]);
                                if (Rc_int == 1) break;
                            }

                            if (k > inter_part[0]){
                                k = 1;
                                J1S(Rc_int, recycled_parts, inter_part[1]);
                            }

                            seed_tmp[0] = '\0';
                            seed[size_seed] = '\0';

                            while ((kmer_tmp = List_pop_first(l_recycl_kmers)) != NULL){

                                if (kmer_tmp[0]) kmer_iupac_comp_to_ascii(&kmer_tmp[1], size_kmer, ascii_kmer);
                                else kmer_comp_to_ascii(&kmer_tmp[1], size_kmer, ascii_kmer);

                                memcpy(seed, ascii_kmer, size_seed * sizeof(char));

                                if (strcmp(seed, seed_tmp)) recycle_sub_paths(seed, inter_part[k], &recycl_subpaths);

                                strcpy(seed_tmp, &ascii_kmer[size_kmer - size_seed]);

                                free(kmer_tmp);
                            }

                            if ((delta_part - prev_delta_part) || (ftell(file_recycl_parts) == 0)){
                                delta_part_tmp = ((delta_part - prev_delta_part) << 1) | 0x1;
                                fwrite(&delta_part_tmp, sizeof(uint32_t), 1, file_recycl_parts);
                            }

                            prev_delta_part = delta_part;

                            nb_kmers_recycled <<= 1;
                            fwrite(&nb_kmers_recycled, sizeof(uint32_t), 1, file_recycl_parts);
                        }

                        if (inter_part != NULL){
                            free(inter_part);
                            inter_part = NULL;
                        }

                        delta_part++;
                    }

                    part_to_insert = parts_to_insert;
                }
                else if (recycle_paths){

                    if (inter_part == NULL) inter_part = get_list_id_genomes(bft_annot, root);
                    else{

                        inter_part_new_node = get_list_id_genomes(bft_annot, root);
                        inter_part_tmp = intersection_uint32_SIMD(inter_part, inter_part_new_node);

                        if (inter_part_tmp[0] == 0){

                            free(inter_part_tmp);

                            if (l_recycl_kmers->count < min_count_l_recycl_kmers){

                                insert_kmers_list(l_recycl_kmers, bft_annot, root_IUPAC, lvl_root_iupac, root_no_IUPAC,
                                                  lvl_root_no_iupac, part_to_insert, size_part_to_insert);
                            }
                            else{

                                nb_kmers_recycled = (uint32_t) l_recycl_kmers->count;

                                for (k = 1; k <= inter_part[0]; k++){
                                    J1T(Rc_int, recycled_parts, inter_part[k]);
                                    if (Rc_int == 1) break;
                                }

                                if (k > inter_part[0]){
                                    k = 1;
                                    J1S(Rc_int, recycled_parts, inter_part[1]);
                                }

                                seed_tmp[0] = '\0';
                                seed[size_seed] = '\0';

                                while ((kmer_tmp = List_pop_first(l_recycl_kmers)) != NULL){

                                    if (kmer_tmp[0]) kmer_iupac_comp_to_ascii(&kmer_tmp[1], size_kmer, ascii_kmer);
                                    else kmer_comp_to_ascii(&kmer_tmp[1], size_kmer, ascii_kmer);

                                    memcpy(seed, ascii_kmer, size_seed * sizeof(char));

                                    if (strcmp(seed, seed_tmp)) recycle_sub_paths(seed, inter_part[k], &recycl_subpaths);

                                    strcpy(seed_tmp, &ascii_kmer[size_kmer - size_seed]);

                                    free(kmer_tmp);
                                }

                                if ((delta_part - prev_delta_part) || (ftell(file_recycl_parts) == 0)){
                                    delta_part_tmp = ((delta_part - prev_delta_part) << 1) | 0x1;
                                    fwrite(&delta_part_tmp, sizeof(uint32_t), 1, file_recycl_parts);
                                }

                                prev_delta_part = delta_part;

                                nb_kmers_recycled <<= 1;
                                fwrite(&nb_kmers_recycled, sizeof(uint32_t), 1, file_recycl_parts);
                            }

                            free(inter_part);

                            inter_part = inter_part_new_node;
                        }
                        else {

                            free(inter_part);
                            free(inter_part_new_node);

                            inter_part = inter_part_tmp;
                        }
                    }

                    if (inter_part != NULL){

                        if (kmer_is_iupac){
                            kmer_tmp = malloc((nb_bytes_kmer_iupac + 1) * sizeof(uint8_t));
                            memcpy(&kmer_tmp[1], kmer, nb_bytes_kmer_iupac * sizeof(uint8_t));
                            kmer_tmp[0] = 0xff;
                        }
                        else{
                            kmer_tmp = malloc((nb_bytes_kmer + 1) * sizeof(uint8_t));
                            memcpy(&kmer_tmp[1], kmer, nb_bytes_kmer * sizeof(uint8_t));
                            kmer_tmp[0] = 0x0;
                        }

                        List_push(l_recycl_kmers, kmer_tmp);
                    }
                }
                else insert_partition(res, bft_annot, root, uc, part_to_insert, size_part_to_insert);
            }
            else{

                memset(substring, 0, nb_bytes_kmer_iupac * sizeof(uint8_t));

                if (kmer_is_iupac) memcpy(substring, kmer, nb_bytes_kmer_iupac * sizeof(uint8_t));
                else memcpy(substring, kmer, nb_bytes_kmer * sizeof(uint8_t));

                size_suffix = root->k;
                current_level = root->k/NB_CHAR_SUF_PREF - 1;
                nb_cell = root->info_per_lvl[current_level].size_kmer_in_bytes;

                while (current_level > res->level_node){

                    nb_cell_to_delete = 2 + ((size_suffix == 45) || (size_suffix == 81) || (size_suffix == 117));

                    for (m=0; m < nb_cell - nb_cell_to_delete; m++){
                        substring[m] = substring[m+2] >> 2;
                        if (m+3 < nb_cell) substring[m] |= substring[m+3] << 6;
                    }

                    substring[m-1] &= root->info_per_lvl[current_level].mask_shift_kmer;

                    current_level--;
                    size_suffix -= NB_CHAR_SUF_PREF;
                    nb_cell = root->info_per_lvl[current_level].size_kmer_in_bytes;
                }

                delta_part_tmp = delta_part;
                prev_part_to_insert = part_to_insert;
                prev_size_part_to_insert = size_part_to_insert;

                if (count_nb_kmer_curr_part > nb_kmer_curr_part){

                    count_nb_kmer_curr_part = 1;
                    delta_part++;

                    if (fread(&part_to_insert, sizeof(uint32_t), 1, file_parts) != 1)
                        ERROR("k-mer file and partition file don't have the same nb. of lines.\n");

                    if (fread(&nb_kmer_curr_part, sizeof(uint32_t), 1, file_parts_count) != 1)
                        ERROR("k-mer file and partition file don't have the same nb. of lines.\n");

                    size_part_to_insert = get_nb_bytes_power2_annot(part_to_insert);
                }

                insertKmer_Node(res->node, root, res->level_node, substring, size_suffix, kmer,
                                part_to_insert, size_part_to_insert, res->pos_container);

                if (recycle_paths){

                    if (l_recycl_kmers->count < min_count_l_recycl_kmers){

                        insert_kmers_list(l_recycl_kmers, bft_annot, root_IUPAC, lvl_root_iupac, root_no_IUPAC,
                                          lvl_root_no_iupac, prev_part_to_insert, prev_size_part_to_insert);
                    }
                    else if (inter_part != NULL){

                        nb_kmers_recycled = (uint32_t) l_recycl_kmers->count;

                        for (k = 1; k <= inter_part[0]; k++){
                            J1T(Rc_int, recycled_parts, inter_part[k]);
                            if (Rc_int == 1) break;
                        }

                        if (k > inter_part[0]){
                            k = 1;
                            J1S(Rc_int, recycled_parts, inter_part[1]);
                        }

                        seed_tmp[0] = '\0';
                        seed[size_seed] = '\0';

                        while ((kmer_tmp = List_pop_first(l_recycl_kmers)) != NULL){

                            if (kmer_tmp[0]) kmer_iupac_comp_to_ascii(&kmer_tmp[1], size_kmer, ascii_kmer);
                            else kmer_comp_to_ascii(&kmer_tmp[1], size_kmer, ascii_kmer);

                            memcpy(seed, ascii_kmer, size_seed * sizeof(char));

                            if (strcmp(seed, seed_tmp)) recycle_sub_paths(seed, inter_part[k], &recycl_subpaths);

                            strcpy(seed_tmp, &ascii_kmer[size_kmer - size_seed]);

                            free(kmer_tmp);
                        }

                        if ((delta_part_tmp - prev_delta_part) || (ftell(file_recycl_parts) == 0)){
                            tmp = ((delta_part_tmp - prev_delta_part) << 1) | 0x1;
                            fwrite(&tmp, sizeof(uint32_t), 1, file_recycl_parts);
                        }

                        prev_delta_part = delta_part_tmp;

                        nb_kmers_recycled <<= 1;
                        fwrite(&nb_kmers_recycled, sizeof(uint32_t), 1, file_recycl_parts);
                    }

                    if (inter_part != NULL){
                        free(inter_part);
                        inter_part = NULL;
                    }
                }
            }

            free(res);
        }

        if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+j)%PRINT_EVERY_X_KMERS)){
            printf("%" PRIu64 " kmers read\n", kmers_read+j);
        }

        kmers_read += j;
    }

    if (recycle_paths){

        if (l_recycl_kmers->count < min_count_l_recycl_kmers){

            insert_kmers_list(l_recycl_kmers, bft_annot, root_IUPAC, lvl_root_iupac, root_no_IUPAC,
                              lvl_root_no_iupac, part_to_insert, prev_size_part_to_insert);
        }
        else if (inter_part != NULL){

            nb_kmers_recycled = (uint32_t) l_recycl_kmers->count;

            for (k = 1; k <= inter_part[0]; k++){
                J1T(Rc_int, recycled_parts, inter_part[k]);
                if (Rc_int == 1) break;
            }

            if (k > inter_part[0]){
                k = 1;
                J1S(Rc_int, recycled_parts, inter_part[1]);
            }

            seed_tmp[0] = '\0';
            seed[size_seed] = '\0';

            while ((kmer_tmp = List_pop_first(l_recycl_kmers)) != NULL){

                if (kmer_tmp[0]) kmer_iupac_comp_to_ascii(&kmer_tmp[1], size_kmer, ascii_kmer);
                else kmer_comp_to_ascii(&kmer_tmp[1], size_kmer, ascii_kmer);

                memcpy(seed, ascii_kmer, size_seed * sizeof(char));

                if (strcmp(seed, seed_tmp)) recycle_sub_paths(seed, inter_part[k], &recycl_subpaths);

                strcpy(seed_tmp, &ascii_kmer[size_kmer - size_seed]);

                free(kmer_tmp);
            }

            if ((delta_part - prev_delta_part) || (ftell(file_recycl_parts) == 0)){
                delta_part_tmp = ((delta_part - prev_delta_part) << 1) | 0x1;
                fwrite(&delta_part_tmp, sizeof(uint32_t), 1, file_recycl_parts);
            }

            nb_kmers_recycled <<= 1;
            fwrite(&nb_kmers_recycled, sizeof(uint32_t), 1, file_recycl_parts);
        }
    }

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_recycl_parts2");
    serialize_subpaths_recycling(output_prefix_buffer, &recycl_subpaths, size_seed);

    fclose(file_kmers);
    fclose(file_parts);
    fclose(file_parts_count);
    fclose(file_recycl_parts);

    J1FA(Rc_word, recycled_parts);

    ascii_kmer[0] = '\0';

    JSLF(PValue, recycl_subpaths, (uint8_t*) ascii_kmer);

    while (PValue != NULL){

        free((uint32_t*) *PValue);
        JSLN(PValue, recycl_subpaths, (uint8_t*) ascii_kmer);
    }

    JSLFA(Rc_word, recycl_subpaths);

    List_clear_destroy(l_recycl_kmers);

    if (inter_part != NULL) free(inter_part);

    free(kmer_no_iupac);
    free(kmers_iupac);

    free(line);
    free(ascii_kmer);

    free(seed);
    free(seed_tmp);

    free(substring);

    free(bft_annot);

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_graph_no_iupac");
    compress_annotations_BFT_disk(root_no_IUPAC, output_prefix_buffer);

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_graph_iupac");
    compress_annotations_BFT_disk(root_IUPAC, output_prefix_buffer);

    free(output_prefix_buffer);

    gettimeofday(&tval_after, NULL);
    time_spent(&tval_before, &tval_after, &tval_result);
    printf("\nElapsed time: %ld.%06ld s\n\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

    return;
}

void compress_annotations_BFT_disk(BFT_Root* bft, char* filename_bft){

    ASSERT_NULL_PTR(bft, "compress_annotations_BFT_disk()\n")
    ASSERT_NULL_PTR(filename_bft, "compress_annotations_BFT_disk()\n")

    FILE* file;

    Pvoid_t comp_annots = (PWord_t)NULL;
    Word_t Rc_word;

    memory_Used* mem;

    int len_longest_annot;
    int lvl_bft = (bft->k / NB_CHAR_SUF_PREF) - 1;

    char* filename_bft_tmp = malloc((strlen(filename_bft) + 50) * sizeof(char));
    ASSERT_NULL_PTR(filename_bft_tmp, "compress_annotations_BFT_disk()\n")

    strcpy(filename_bft_tmp, filename_bft);
    strcpy(&filename_bft_tmp[strlen(filename_bft)], "_tmp");

    mem = printMemoryUsedFromNode(&bft->node, lvl_bft, bft->k, bft->info_per_lvl);
    len_longest_annot = (int) MAX(mem->size_biggest_annot+1, getMaxSize_annotation_array_elem(bft->comp_set_colors));
    free(mem);

    bft->compressed = 0;
    write_BFT_Root_sparse(bft, filename_bft, false);
    bft->compressed = 1;

    load_annotation_from_Node2(&bft->node, lvl_bft, bft->k, len_longest_annot,
                              bft->info_per_lvl, &comp_annots, bft->comp_set_colors, bft->ann_inf);

    freeNode(&bft->node, lvl_bft, bft->info_per_lvl);
    free_annotation_array_elem(&bft->comp_set_colors, &bft->length_comp_set_colors);

    sort_annotations3(&comp_annots, len_longest_annot);
    write_partial_comp_set_colors(filename_bft_tmp, &comp_annots, len_longest_annot);

    write_BFT_Root_sparse2(bft, filename_bft_tmp, true);
    freeBFT_Root(bft);

    bft = read_BFT_Root_sparse2(filename_bft, filename_bft_tmp, &comp_annots, len_longest_annot);
    bft->compressed = 1;

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

    rename(filename_bft_tmp, filename_bft);

    file = fopen(filename_bft, "ab");
    ASSERT_NULL_PTR(file, "compress_annotations_BFT_disk()\n")

    write_Node_sparse(&bft->node, bft, lvl_bft, file, bft->k);

    fclose(file);

    freeBFT_Root(bft);

    free(filename_bft_tmp);

    return;
}

void create_subgraph_decomp(char* output_prefix, FILE* file_partitions, Pvoid_t* prefix_part_recycl,
                            BFT_Root** graph_no_iupac, BFT_Root** graph_iupac, int size_seed){

    ASSERT_NULL_PTR(output_prefix, "create_subgraph_decomp()\n")
    ASSERT_NULL_PTR(file_partitions, "create_subgraph_decomp()\n")
    ASSERT_NULL_PTR(prefix_part_recycl, "create_subgraph_decomp()\n")
    ASSERT_NULL_PTR(graph_no_iupac, "create_subgraph_decomp()\n")
    ASSERT_NULL_PTR(graph_iupac, "create_subgraph_decomp()\n")
    ASSERT_NULL_PTR(*graph_no_iupac, "create_subgraph_decomp()\n")
    ASSERT_NULL_PTR(*graph_iupac, "create_subgraph_decomp()\n")

    Pvoid_t parts = (Pvoid_t) NULL;
    Word_t Rc_word;
    int Rc_int;

    memory_Used* mem;

    BFT_Root* new_graph_no_iupac;
    BFT_Root* new_graph_iupac;

    int len_output_prefix = strlen(output_prefix);

    int lvl_root_no_iupac = (*graph_no_iupac)->k / NB_CHAR_SUF_PREF - 1;
    int lvl_root_iupac = (*graph_iupac)->k / NB_CHAR_SUF_PREF - 1;

    uint32_t part_tmp;

    uint32_t part_max = 0;
    uint32_t part_min = 0xffffffff;

    uint32_t nb_kmers = 0;

    char* seed = calloc(size_seed * 2 + 1, sizeof(char));
    ASSERT_NULL_PTR(seed, "create_subgraph_decomp()\n")

    char* output_prefix_buffer = malloc((len_output_prefix + 50) * sizeof(char));
    ASSERT_NULL_PTR(output_prefix_buffer, "create_subgraph_decomp()\n")

    strcpy(output_prefix_buffer, output_prefix);

    new_graph_no_iupac = createBFT_Root((*graph_no_iupac)->k, 1, 1);
    new_graph_no_iupac->ann_inf = create_annotation_inform(-1);

    for (int i = 0; i <= lvl_root_no_iupac; i++){
        new_graph_no_iupac->info_per_lvl[i].nb_ucs_skp = 8;
        new_graph_no_iupac->info_per_lvl[i].nb_bits_per_cell_skip_filter2 = new_graph_no_iupac->info_per_lvl[i].nb_ucs_skp;
        new_graph_no_iupac->info_per_lvl[i].nb_bits_per_cell_skip_filter3 = new_graph_no_iupac->info_per_lvl[i].nb_ucs_skp;
        new_graph_no_iupac->info_per_lvl[i].nb_bytes_per_cell_skip_filter2 = CEIL(new_graph_no_iupac->info_per_lvl[i].nb_ucs_skp, SIZE_BITS_UINT_8T);
        new_graph_no_iupac->info_per_lvl[i].nb_bytes_per_cell_skip_filter3 = CEIL(new_graph_no_iupac->info_per_lvl[i].nb_ucs_skp, SIZE_BITS_UINT_8T);
    }

    new_graph_iupac = createBFT_Root((*graph_iupac)->k, 1, 1);
    new_graph_iupac->ann_inf = create_annotation_inform(-1);

    for (int i = 0; i <= lvl_root_iupac; i++){
        new_graph_iupac->info_per_lvl[i].nb_ucs_skp = 8;
        new_graph_iupac->info_per_lvl[i].nb_bits_per_cell_skip_filter2 = new_graph_iupac->info_per_lvl[i].nb_ucs_skp;
        new_graph_iupac->info_per_lvl[i].nb_bits_per_cell_skip_filter3 = new_graph_iupac->info_per_lvl[i].nb_ucs_skp;
        new_graph_iupac->info_per_lvl[i].nb_bytes_per_cell_skip_filter2 = CEIL(new_graph_iupac->info_per_lvl[i].nb_ucs_skp, SIZE_BITS_UINT_8T);
        new_graph_iupac->info_per_lvl[i].nb_bytes_per_cell_skip_filter3 = CEIL(new_graph_iupac->info_per_lvl[i].nb_ucs_skp, SIZE_BITS_UINT_8T);
    }

    while (fread(&part_tmp, sizeof(uint32_t), 1, file_partitions) == 1){
        part_min = MIN(part_min, part_tmp);
        part_max = MAX(part_max, part_tmp);
    }

    for (uint32_t i = part_min; i <= part_max; i++) J1S(Rc_int, parts, i);

    iterate_over_kmers(*graph_iupac, load_subgraph, *graph_no_iupac, new_graph_no_iupac, new_graph_iupac, prefix_part_recycl, &parts,
                       true, size_seed, seed, &nb_kmers);

    iterate_over_kmers(*graph_no_iupac, load_subgraph, *graph_iupac, new_graph_no_iupac, new_graph_iupac, prefix_part_recycl, &parts,
                       false, size_seed, seed, &nb_kmers);

    freeBFT_Root(*graph_iupac);
    freeBFT_Root(*graph_no_iupac);

    J1FA(Rc_word, parts);

    free(seed);

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_graph_no_iupac");
    compress_annotations_BFT_disk(new_graph_no_iupac, output_prefix_buffer);

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_graph_iupac");
    compress_annotations_BFT_disk(new_graph_iupac, output_prefix_buffer);

    *graph_iupac = read_BFT_Root_sparse(output_prefix_buffer);
    remove(output_prefix_buffer);

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_graph_no_iupac");
    *graph_no_iupac = read_BFT_Root_sparse(output_prefix_buffer);
    remove(output_prefix_buffer);

    free(output_prefix_buffer);

    mem = printMemoryUsedFromNode(&((*graph_no_iupac)->node), ((*graph_no_iupac)->k / NB_CHAR_SUF_PREF) - 1,
                                  (*graph_no_iupac)->k, (*graph_no_iupac)->info_per_lvl);
    printMemory(mem);
    free(mem);

    mem = printMemoryUsedFromNode(&((*graph_iupac)->node), ((*graph_iupac)->k / NB_CHAR_SUF_PREF) - 1,
                                  (*graph_iupac)->k, (*graph_iupac)->info_per_lvl);
    printMemory(mem);
    free(mem);

    return;
}

size_t load_subgraph(BFT_kmer* kmer, BFT* graph, va_list args){

    UC* uc;
    resultPresence* res;
    BFT_annotation* bft_annot;
    BFT_Root* new_graph;

    Word_t* PValue;

    int Rc_int;
    int nb_bytes_kmer;
    int lvl_new_graph;

    int m;
    int nb_cell_to_delete;

    int lvl_suffix;
    int size_suffix;
    int nb_cell;

    uint32_t i, j;
    uint32_t current_part;

    uint32_t pow2_part = 0;
    uint32_t size_part_bytes = 1;

    uint32_t* list_part_recycl;
    uint32_t* list_pref_part_recycl;
    uint32_t* list_kmer_partitions;

    uint8_t* substring;

    BFT_Root* accessory_graph = va_arg(args, BFT_Root*);
    BFT_Root* new_graph_no_iupac = va_arg(args, BFT_Root*);
    BFT_Root* new_graph_iupac = va_arg(args, BFT_Root*);

    Pvoid_t* prefix_part_recycl = va_arg(args, Pvoid_t*);
    Pvoid_t* part_recycl = va_arg(args, Pvoid_t*);

    bool graph_is_iupac = va_arg(args, int);

    int size_seed = va_arg(args, int);

    char* seed = va_arg(args, char*);

    uint32_t* nb_kmers = va_arg(args, int*);

    bft_annot = get_annotation(kmer);

    list_kmer_partitions = get_list_id_genomes(bft_annot, graph);

    if (graph_is_iupac){
        kmer_iupac_comp_to_ascii(kmer->kmer_comp, graph->k/2, kmer->kmer);
        new_graph = new_graph_iupac;
    }
    else new_graph = new_graph_no_iupac;

    memcpy(seed, kmer->kmer, size_seed * sizeof(char));
    seed[size_seed] = '\0';

    JSLG(PValue, *prefix_part_recycl, (uint8_t*)seed);

    list_part_recycl = malloc((list_kmer_partitions[0] + 1) * sizeof(uint32_t));
    ASSERT_NULL_PTR(list_part_recycl, "load_subgraph()\n")

    list_part_recycl[0] = 0;

    if (PValue == NULL){

        for (i = 1; i <= list_kmer_partitions[0]; i++){

            J1T(Rc_int, *part_recycl, list_kmer_partitions[i]);

            if (Rc_int){
                list_part_recycl[0]++;
                list_part_recycl[list_part_recycl[0]] = list_kmer_partitions[i] << 1;
            }
        }
    }
    else {

        list_pref_part_recycl = (uint32_t*) *PValue;

        for (i = 1, j = 1; i <= list_kmer_partitions[0]; i++){

            J1T(Rc_int, *part_recycl, list_kmer_partitions[i]);

            if (Rc_int){
                list_part_recycl[0]++;
                list_part_recycl[list_part_recycl[0]] = list_kmer_partitions[i] << 1;
            }

            while((j <= list_pref_part_recycl[0]) && (list_pref_part_recycl[j] < list_kmer_partitions[i])) j++;

            if ((j <= list_pref_part_recycl[0]) && (list_pref_part_recycl[j] == list_kmer_partitions[i])){

                if (Rc_int) list_part_recycl[list_part_recycl[0]] |= 0x1;
                else{
                    list_part_recycl[0]++;
                    list_part_recycl[list_part_recycl[0]] = (list_kmer_partitions[i] << 1) | 0x1;
                }

                j++;
            }
        }
    }

    free(list_kmer_partitions);

    if (list_part_recycl[0]){

        lvl_new_graph = (new_graph->k / NB_CHAR_SUF_PREF) - 1;

        res = isKmerPresent(&(new_graph->node), new_graph, lvl_new_graph, kmer->kmer_comp, new_graph->k);

        if (res->link_child == NULL){

            current_part = list_part_recycl[1] >> 1;

            if (current_part >= pow2_part){
                pow2_part = round_up_next_highest_power2(current_part);
                size_part_bytes = get_nb_bytes_power2_annot_bis(current_part, pow2_part);
            }

            lvl_suffix = lvl_new_graph;
            size_suffix = new_graph->k;
            nb_cell = new_graph->info_per_lvl[lvl_suffix].size_kmer_in_bytes;
            nb_bytes_kmer = CEIL(new_graph->k * 2, SIZE_BITS_UINT_8T);

            substring = malloc(nb_bytes_kmer * sizeof(uint8_t));
            ASSERT_NULL_PTR(substring, "kmer_has_partition()\n")

            memcpy(substring, kmer->kmer_comp, nb_bytes_kmer * sizeof(uint8_t));

            while (lvl_suffix > res->level_node){

                nb_cell_to_delete = 2 + ((size_suffix == 45) || (size_suffix == 81) || (size_suffix == 117));

                for (m=0; m < nb_cell - nb_cell_to_delete; m++){
                    substring[m] = substring[m+2] >> 2;
                    if (m+3 < nb_cell) substring[m] |= substring[m+3] << 6;
                }

                substring[m-1] &= new_graph->info_per_lvl[lvl_suffix].mask_shift_kmer;

                lvl_suffix--;
                size_suffix -= NB_CHAR_SUF_PREF;
                nb_cell = new_graph->info_per_lvl[lvl_suffix].size_kmer_in_bytes;
            }

            insertKmer_Node(res->node, new_graph, res->level_node, substring, size_suffix, kmer->kmer_comp,
                            current_part, size_part_bytes, res->pos_container);

            free(substring);


            if (list_part_recycl[0] > 1){

                res = isKmerPresent(&(new_graph->node), new_graph, lvl_new_graph, kmer->kmer_comp, new_graph->k);

                if (res->posFilter2 != 0) uc = (UC*)res->container;
                else uc = &(((UC*)((CC*)res->container)->children)[res->bucket]);

                for (uint32_t k = 2; k <= list_part_recycl[0]; k++){

                    current_part = list_part_recycl[k] >> 1;

                    if (current_part >= pow2_part){
                        pow2_part = round_up_next_highest_power2(current_part);
                        size_part_bytes = get_nb_bytes_power2_annot_bis(current_part, pow2_part);
                    }

                    get_annot(uc, &bft_annot->annot, &bft_annot->annot_ext, &bft_annot->annot_cplx, &bft_annot->size_annot,
                              &bft_annot->size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);

                    compute_best_mode(new_graph->ann_inf, new_graph->comp_set_colors, bft_annot->annot, bft_annot->size_annot, bft_annot->annot_ext,
                                      1, current_part, size_part_bytes);

                    if (new_graph->ann_inf->last_added != current_part){

                        if (new_graph->ann_inf->min_size > uc->size_annot){
                            if ((bft_annot->annot_ext == NULL) || (new_graph->ann_inf->min_size > uc->size_annot+1)){
                                bft_annot->annot_ext = realloc_annotation(uc, res->posFilter2, res->posFilter3,
                                                                          new_graph->ann_inf->min_size, 0, res->pos_sub_bucket);
                            }
                        }

                        modify_mode_annotation(new_graph->ann_inf, &(uc->suffixes[res->pos_sub_bucket * (res->posFilter2 + uc->size_annot) + res->posFilter2]),
                                               uc->size_annot, bft_annot->annot_ext, 1, current_part, size_part_bytes);

                        if ((bft_annot->annot_ext != NULL) && (bft_annot->annot_ext[0] == 0))
                            delete_extend_annots(uc, res->posFilter2, res->posFilter3, res->pos_sub_bucket, res->pos_sub_bucket, 0, 0, 1);
                    }

                    reinit_annotation_inform(new_graph->ann_inf);
                }
            }
        }
        else {

            if (res->posFilter2 != 0) uc = (UC*)res->container;
            else uc = &(((UC*)((CC*)res->container)->children)[res->bucket]);

            for (uint32_t k = 1; k <= list_part_recycl[0]; k++){

                current_part = list_part_recycl[k] >> 1;

                if (current_part >= pow2_part){
                    pow2_part = round_up_next_highest_power2(current_part);
                    size_part_bytes = get_nb_bytes_power2_annot_bis(current_part, pow2_part);
                }

                get_annot(uc, &bft_annot->annot, &bft_annot->annot_ext, &bft_annot->annot_cplx, &bft_annot->size_annot,
                          &bft_annot->size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);

                insert_partition(res, bft_annot, new_graph, uc, current_part, size_part_bytes);
            }
        }

        free(res);

        pow2_part = 0;
        size_part_bytes = 1;

        if (graph_is_iupac){

            for (uint32_t k = 1; k <= list_part_recycl[0]; k++){

                if (list_part_recycl[k] & 0x1){

                    list_part_recycl[k] >>= 1;

                    if (list_part_recycl[k] >= pow2_part){
                        pow2_part = round_up_next_highest_power2(list_part_recycl[k]);
                        size_part_bytes = get_nb_bytes_power2_annot_bis(list_part_recycl[k], pow2_part);
                    }

                    strcpy(seed, &(kmer->kmer[graph->k / 2 - size_seed]));

                    cdbg_traverse_overlap_insert(accessory_graph, graph, new_graph_no_iupac, new_graph_iupac,
                                                 seed, list_part_recycl[k], size_part_bytes);
                }
            }
        }
        else {

            for (uint32_t k = 1; k <= list_part_recycl[0]; k++){

                if (list_part_recycl[k] & 0x1){

                    list_part_recycl[k] >>= 1;

                    if (list_part_recycl[k] >= pow2_part){
                        pow2_part = round_up_next_highest_power2(list_part_recycl[k]);
                        size_part_bytes = get_nb_bytes_power2_annot_bis(list_part_recycl[k], pow2_part);
                    }

                    strcpy(seed, &(kmer->kmer[graph->k - size_seed]));

                    cdbg_traverse_overlap_insert(graph, accessory_graph, new_graph_no_iupac, new_graph_iupac,
                                                 seed, list_part_recycl[k], size_part_bytes);
                }
            }
        }
    }

    free_BFT_annotation(bft_annot);
    free(list_part_recycl);

    *nb_kmers += 1;

    if (*nb_kmers % 1000000 == 0) printf("%" PRIu32 "k-mers analyzed\n", *nb_kmers);

    return 1;
}

void cdbg_traverse_overlap_insert(BFT* bft_no_iupac, BFT* bft_iupac, BFT* insert_bft_no_iupac, BFT* insert_bft_iupac,
                                  char* overlap, uint32_t partition, uint32_t size_partition){

    ASSERT_NULL_PTR(bft_no_iupac, "cdbg_traverse_overlap_insert()\n")
    ASSERT_NULL_PTR(bft_iupac, "cdbg_traverse_overlap_insert()\n")
    ASSERT_NULL_PTR(insert_bft_no_iupac, "cdbg_traverse_overlap_insert()\n")
    ASSERT_NULL_PTR(insert_bft_iupac, "cdbg_traverse_overlap_insert()\n")
    ASSERT_NULL_PTR(overlap, "cdbg_traverse_overlap_insert()\n")

    bool kmer_has_partition = true;

    int length_overlap = strlen(overlap);
    int size_overlap_comp = CEIL(length_overlap * 4, SIZE_BITS_UINT_8T);

    uint8_t* overlap_comp;

    while (kmer_has_partition){

        kmer_has_partition = false;

        if (is_substring_IUPAC(overlap)){

            overlap_comp = calloc(size_overlap_comp, sizeof(uint8_t));
            ASSERT_NULL_PTR(overlap_comp, "cdbg_traverse_overlap_insert()\n")

            parseKmerCount_IUPAC(overlap, length_overlap, overlap_comp, 0);
            kmer_comp_to_ascii(overlap_comp, length_overlap * 2, overlap);

            prefix_matching(bft_iupac, overlap, insert_kmer_with_partition, overlap,
                            insert_bft_iupac, partition, size_partition, &kmer_has_partition);

            if (kmer_has_partition){

                memset(overlap_comp, 0, size_overlap_comp * sizeof(uint8_t));
                parseKmerCount(overlap, length_overlap * 2, overlap_comp, 0);
                kmer_iupac_comp_to_ascii(overlap_comp, length_overlap, overlap);
            }

            free(overlap_comp);
        }
        else {

            prefix_matching(bft_no_iupac, overlap, insert_kmer_with_partition, overlap,
                            insert_bft_no_iupac, partition, size_partition, &kmer_has_partition);

            if (!kmer_has_partition){

                overlap_comp = calloc(size_overlap_comp, sizeof(uint8_t));
                ASSERT_NULL_PTR(overlap_comp, "cdbg_traverse_overlap_insert()\n")

                parseKmerCount_IUPAC(overlap, length_overlap, overlap_comp, 0);
                kmer_comp_to_ascii(overlap_comp, length_overlap * 2, overlap);

                prefix_matching(bft_iupac, overlap, insert_kmer_with_partition, overlap,
                                insert_bft_iupac, partition, size_partition, &kmer_has_partition);

                if (kmer_has_partition){

                    memset(overlap_comp, 0, size_overlap_comp * sizeof(uint8_t));
                    parseKmerCount(overlap, length_overlap * 2, overlap_comp, 0);
                    kmer_iupac_comp_to_ascii(overlap_comp, length_overlap, overlap);
                }

                free(overlap_comp);
            }
        }
    }

    return;
}

bool cdbg_get_successor_overlap(BFT* bft, bool graph_is_iupac, char** overlap, uint32_t partition){

    ASSERT_NULL_PTR(bft, "cdbg_get_successor_overlap()\n")
    ASSERT_NULL_PTR(overlap, "cdbg_get_successor_overlap()\n")
    ASSERT_NULL_PTR(*overlap, "cdbg_get_successor_overlap()\n")

    bool kmer_has_partition_bool = false;

    if (graph_is_iupac){

        int length_overlap = strlen(*overlap);

        char* overlap_no_iupac = malloc((bft->k + 1) * sizeof(char));
        ASSERT_NULL_PTR(overlap_no_iupac, "cdbg_get_successor_overlap()\n")

        uint8_t* overlap_comp = calloc(CEIL(bft->k * 2, SIZE_BITS_UINT_8T), sizeof(uint8_t));
        ASSERT_NULL_PTR(overlap_comp, "cdbg_get_successor_overlap()\n")

        parseKmerCount_IUPAC(*overlap, length_overlap, overlap_comp, 0);

        kmer_comp_to_ascii(overlap_comp, length_overlap * 2, overlap_no_iupac);

        prefix_matching(bft, overlap_no_iupac, kmer_has_partition, overlap_no_iupac, partition, &kmer_has_partition_bool);

        if (kmer_has_partition_bool){

            *overlap = realloc(*overlap, (bft->k / 2 + 1) * sizeof(char));
            ASSERT_NULL_PTR(*overlap, "cdbg_get_successor_overlap()\n")

            parseKmerCount(overlap_no_iupac, bft->k, overlap_comp, 0);

            kmer_iupac_comp_to_ascii(overlap_comp, bft->k / 2, *overlap);
        }

        free(overlap_no_iupac);
        free(overlap_comp);
    }
    else{

        *overlap = realloc(*overlap, (bft->k + 1) * sizeof(char));
        ASSERT_NULL_PTR(*overlap, "cdbg_get_successor_overlap()\n")

        prefix_matching(bft, *overlap, kmer_has_partition, *overlap, partition, &kmer_has_partition_bool);
    }

    return kmer_has_partition_bool;
}

bool cdbg_get_successor_overlap_custom(BFT* bft, bool graph_is_iupac, char** overlap, resultPresence** junction_overlap, uint32_t partition){

    ASSERT_NULL_PTR(bft, "cdbg_get_successor_overlap()\n")
    ASSERT_NULL_PTR(overlap, "cdbg_get_successor_overlap()\n")
    ASSERT_NULL_PTR(*overlap, "cdbg_get_successor_overlap()\n")

    bool kmer_has_partition_bool = false;

    if (graph_is_iupac){

        int length_overlap = strlen(*overlap);

        char* overlap_no_iupac = malloc((bft->k + 1) * sizeof(char));
        ASSERT_NULL_PTR(overlap_no_iupac, "cdbg_get_successor_overlap()\n")

        uint8_t* overlap_comp = calloc(CEIL(bft->k * 2, SIZE_BITS_UINT_8T), sizeof(uint8_t));
        ASSERT_NULL_PTR(overlap_comp, "cdbg_get_successor_overlap()\n")

        parseKmerCount_IUPAC(*overlap, length_overlap, overlap_comp, 0);

        kmer_comp_to_ascii(overlap_comp, length_overlap * 2, overlap_no_iupac);

        prefix_matching(bft, overlap_no_iupac, kmer_has_partition, overlap_no_iupac, partition, &kmer_has_partition_bool);

        if (kmer_has_partition_bool){

            *overlap = realloc(*overlap, (bft->k / 2 + 1) * sizeof(char));
            ASSERT_NULL_PTR(*overlap, "cdbg_get_successor_overlap()\n")

            parseKmerCount(overlap_no_iupac, bft->k, overlap_comp, 0);

            kmer_iupac_comp_to_ascii(overlap_comp, bft->k / 2, *overlap);
        }

        free(overlap_no_iupac);
        free(overlap_comp);
    }
    else{

        *overlap = realloc(*overlap, (bft->k + 1) * sizeof(char));
        ASSERT_NULL_PTR(*overlap, "cdbg_get_successor_overlap()\n")

        prefix_matching_custom(bft, *overlap, junction_overlap, kmer_has_partition, *overlap, partition, &kmer_has_partition_bool);
    }

    return kmer_has_partition_bool;
}

size_t insert_kmer_with_partition(BFT_kmer* kmer, BFT* graph, va_list args){

    char* overlap = va_arg(args, char*);

    BFT* bft_to_insert = va_arg(args, BFT*);

    uint32_t partition = va_arg(args, uint32_t);
    uint32_t size_partition = va_arg(args, uint32_t);

    bool* kmer_has_partition_bool = va_arg(args, int*);

    BFT_annotation* bft_annot = get_annotation(kmer);

    *kmer_has_partition_bool = is_genome_present(graph->ann_inf, graph->comp_set_colors, bft_annot->annot, bft_annot->size_annot,
                                                 bft_annot->annot_ext, bft_annot->annot_ext != NULL, partition);

    if (*kmer_has_partition_bool){

        int lvl_bft_to_insert = bft_to_insert->k / NB_CHAR_SUF_PREF - 1;

        strcpy(overlap, &(kmer->kmer[graph->k - strlen(overlap)]));

        resultPresence* res = isKmerPresent(&(bft_to_insert->node), bft_to_insert, lvl_bft_to_insert,
                                            kmer->kmer_comp, bft_to_insert->k);

        if (res->link_child != NULL){

            UC* uc;

            if (res->posFilter2 != 0) uc = (UC*)res->container;
            else uc = &(((UC*)((CC*)res->container)->children)[res->bucket]);

            get_annot(uc, &bft_annot->annot, &bft_annot->annot_ext, &bft_annot->annot_cplx, &bft_annot->size_annot,
                      &bft_annot->size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);

            insert_partition(res, bft_annot, bft_to_insert, uc, partition, size_partition);
        }
        else{

            int m;
            int nb_cell_to_delete;

            int size_suffix = bft_to_insert->k;
            int nb_cell = bft_to_insert->info_per_lvl[lvl_bft_to_insert].size_kmer_in_bytes;
            int nb_bytes_kmer = CEIL(bft_to_insert->k * 2, SIZE_BITS_UINT_8T);

            uint8_t* substring = malloc(nb_bytes_kmer * sizeof(uint8_t));
            ASSERT_NULL_PTR(substring, "insert_kmer_with_partition()\n")

            memcpy(substring, kmer->kmer_comp, nb_bytes_kmer * sizeof(uint8_t));

            while (lvl_bft_to_insert > res->level_node){

                nb_cell_to_delete = 2 + ((size_suffix == 45) || (size_suffix == 81) || (size_suffix == 117));

                for (m=0; m < nb_cell - nb_cell_to_delete; m++){
                    substring[m] = substring[m+2] >> 2;
                    if (m+3 < nb_cell) substring[m] |= substring[m+3] << 6;
                }

                substring[m-1] &= bft_to_insert->info_per_lvl[lvl_bft_to_insert].mask_shift_kmer;

                lvl_bft_to_insert--;
                size_suffix -= NB_CHAR_SUF_PREF;
                nb_cell = bft_to_insert->info_per_lvl[lvl_bft_to_insert].size_kmer_in_bytes;
            }

            insertKmer_Node(res->node, bft_to_insert, res->level_node, substring, size_suffix, kmer->kmer_comp,
                            partition, size_partition, res->pos_container);

            free(substring);
        }

        free(res);
        free_BFT_annotation(bft_annot);

        return 0;
    }

    free_BFT_annotation(bft_annot);

    return 1;
}

size_t kmer_has_partition(BFT_kmer* kmer, BFT* graph, va_list args){

    char* overlap = va_arg(args, char*);

    uint32_t partition = va_arg(args, uint32_t);

    bool* kmer_has_partition_bool = va_arg(args, int*);

    BFT_annotation* bft_annot = get_annotation(kmer);

    *kmer_has_partition_bool = is_genome_present(graph->ann_inf, graph->comp_set_colors, bft_annot->annot, bft_annot->size_annot,
                                                 bft_annot->annot_ext, bft_annot->annot_ext != NULL, partition);

    if (*kmer_has_partition_bool) strcpy(overlap, kmer->kmer);

    free_BFT_annotation(bft_annot);

    if (*kmer_has_partition_bool) return 0;
    return 1;
}

void erase_partitions(uint8_t* partition_no_iupac, Pvoid_t* partition_iupac, Pvoid_t* pref_suf_pos){

    Word_t Rc_word;
    int Rc_int;

    Rc_word = 0;
    J1F(Rc_int, *pref_suf_pos, Rc_word);

    while (Rc_int){
        partition_no_iupac[Rc_word] = 0;
        J1N(Rc_int, *pref_suf_pos, Rc_word);
    }

    J1FA(Rc_word, *pref_suf_pos);
    JSLFA(Rc_word, *partition_iupac);

    return;
}

void write_partition(char* filename, uint8_t* partition_no_iupac, Pvoid_t* partition_iupac, int size_seed){

    ASSERT_NULL_PTR(filename, "write_partition()\n")
    ASSERT_NULL_PTR(partition_no_iupac, "write_partition()\n")
    ASSERT_NULL_PTR(partition_iupac, "write_partition()\n")

    PWord_t PValue;

    int nb_bytes_seed_iupac = CEIL(size_seed * 4, SIZE_BITS_UINT_8T);

    uint32_t size_seed_present = CEIL((uint32_t) pow(4, size_seed), SIZE_BITS_UINT_8T);

    uint8_t* seed_comp = malloc(nb_bytes_seed_iupac * sizeof(uint8_t));
    ASSERT_NULL_PTR(seed_comp, "write_partition()\n")

    uint8_t* seed_ascii = malloc((size_seed + 1) * sizeof(uint8_t));
    ASSERT_NULL_PTR(seed_ascii, "write_partition()\n")

    FILE* file_part = fopen(filename, "wb");
    ASSERT_NULL_PTR(file_part, "write_partition()\n")

    fwrite(partition_no_iupac, sizeof(uint8_t), size_seed_present, file_part);

    seed_ascii[0] = '\0';

    JSLF(PValue, *partition_iupac, seed_ascii);

    while (PValue != NULL){

        memset(seed_comp, 0, nb_bytes_seed_iupac * sizeof(uint8_t));
        parseKmerCount_IUPAC((char*)seed_ascii, size_seed, seed_comp, 0);

        fwrite(seed_comp, sizeof(uint8_t), nb_bytes_seed_iupac, file_part);

        JSLN(PValue, *partition_iupac, seed_ascii);
    }

    free(seed_comp);
    free(seed_ascii);

    fclose(file_part);

    return;
}

void read_partition(char* filename, uint8_t* partition_no_iupac, Pvoid_t* partition_iupac, Pvoid_t* pref_suf_pos, int size_seed){

    ASSERT_NULL_PTR(filename, "write_partition()\n")
    ASSERT_NULL_PTR(partition_no_iupac, "write_partition()\n")
    ASSERT_NULL_PTR(partition_iupac, "write_partition()\n")
    ASSERT_NULL_PTR(pref_suf_pos, "write_partition()\n")

    PWord_t PValue;

    int Rc_int;
    int nb_bytes_seed_iupac = CEIL(size_seed * 4, SIZE_BITS_UINT_8T);

    uint32_t size_seed_present = CEIL((uint32_t) pow(4, size_seed), SIZE_BITS_UINT_8T);

    uint8_t* seed_comp = malloc(nb_bytes_seed_iupac * sizeof(uint8_t));
    ASSERT_NULL_PTR(seed_comp, "write_partition()\n")

    uint8_t* seed_ascii = malloc((size_seed + 1) * sizeof(uint8_t));
    ASSERT_NULL_PTR(seed_ascii, "write_partition()\n")

    FILE* file_part = fopen(filename, "rb");

    if (file_part != NULL){

        fread(partition_no_iupac, sizeof(uint8_t), size_seed_present, file_part);

        for (uint32_t i = 0; i < size_seed_present; i++){
            if (partition_no_iupac[i]) J1S(Rc_int, *pref_suf_pos, i);
        }

        while (fread(seed_comp, sizeof(uint8_t), nb_bytes_seed_iupac, file_part) == nb_bytes_seed_iupac){

            kmer_iupac_comp_to_ascii(seed_comp, size_seed, (char*) seed_ascii);
            JSLI(PValue, *partition_iupac, (uint8_t*) seed_ascii);
        }

        fclose(file_part);
    }

    free(seed_comp);
    free(seed_ascii);

    return;
}

uint8_t find_kmer(BFT* graph_iupac, BFT* graph_no_iupac,
                  uint8_t* partition_no_iupac, Pvoid_t* partition_iupac, Pvoid_t* pref_suf_pos,
                  char** seed, int length_seed, uint32_t traversed_part, bool verify_prefix){

    ASSERT_NULL_PTR(graph_iupac, "find_kmer()\n")
    ASSERT_NULL_PTR(graph_no_iupac, "find_kmer()\n")
    ASSERT_NULL_PTR(seed, "find_kmer()\n")
    ASSERT_NULL_PTR(*seed, "find_kmer()\n")
    ASSERT_NULL_PTR(partition_no_iupac, "find_kmer()\n")
    ASSERT_NULL_PTR(partition_iupac, "find_kmer()\n")
    ASSERT_NULL_PTR(pref_suf_pos, "find_kmer()\n")

    PWord_t PValue_partition_iupac;

    bool prefix_is_IUPAC;
    bool kmer_found = true;

    int z;
    int nb_bytes_seed = CEIL(length_seed * 2, SIZE_BITS_UINT_8T);
    int shift_pref_suf = 0;

    uint64_t prefix;
    uint64_t suffix;

    uint8_t ret = 3;

    uint8_t* seed_comp;

    if (length_seed % 4) shift_pref_suf = SIZE_BITS_UINT_8T - ((length_seed % 4) * 2);

    // Check if the prefix is present in the partition

    prefix_is_IUPAC = is_substring_IUPAC(*seed);

    if (verify_prefix){

        if (prefix_is_IUPAC){

            JSLG(PValue_partition_iupac, *partition_iupac, (uint8_t*) *seed);
            if (PValue_partition_iupac != NULL) ret = 1;
        }
        else{

            seed_comp = calloc(nb_bytes_seed,sizeof(char));
            ASSERT_NULL_PTR(seed_comp, "decompress() 1")

            parseKmerCount(*seed, length_seed, seed_comp, 0);

            for (z = 0, prefix = 0; z < nb_bytes_seed; z++)
                prefix = (prefix << SIZE_BITS_UINT_8T) | reverse_word_8(seed_comp[z]);

            prefix >>= shift_pref_suf;

            if (partition_no_iupac[prefix/SIZE_BITS_UINT_8T] & MASK_POWER_8[prefix%SIZE_BITS_UINT_8T]) ret = 1;

            free(seed_comp);
        }

        //if (ret == 1) printf("Prefix found in partition\n");
    }

    // If the prefix is not present in the partition, look for a k-mer with the corresponding
    // prefix and partition

    if (ret == 3){

        kmer_found = cdbg_get_successor_overlap(graph_iupac, true, seed, traversed_part);

        if (!kmer_found && !prefix_is_IUPAC) kmer_found = cdbg_get_successor_overlap(graph_no_iupac, false, seed, traversed_part);

        if (kmer_found){ // Such a k-mer is found

            // Check if the suffix is the partition

            if (is_substring_IUPAC(&((*seed)[graph_no_iupac->k - length_seed]))){

                JSLG(PValue_partition_iupac, *partition_iupac, (uint8_t*) &((*seed)[graph_no_iupac->k - length_seed]));
                if (PValue_partition_iupac != NULL) ret = 2;
            }
            else{

                seed_comp = calloc(nb_bytes_seed,sizeof(char));
                ASSERT_NULL_PTR(seed_comp, "decompress() 1")

                parseKmerCount(&((*seed)[graph_no_iupac->k - length_seed]), length_seed, seed_comp, 0);

                for (z = 0, suffix = 0; z < nb_bytes_seed; z++)
                    suffix = (suffix << SIZE_BITS_UINT_8T) | reverse_word_8(seed_comp[z]);

                suffix >>= shift_pref_suf;

                if (partition_no_iupac[suffix/SIZE_BITS_UINT_8T] & MASK_POWER_8[suffix%SIZE_BITS_UINT_8T]) ret = 2;

                free(seed_comp);
            }

            //if (ret == 2) printf("Suffix found in partition\n");
        }
        else ret = 0;
    }

    return ret;
}

uint8_t find_kmer_custom(BFT* graph_iupac, BFT* graph_no_iupac,
                        uint8_t* partition_no_iupac, Pvoid_t* partition_iupac, Pvoid_t* pref_suf_pos,
                        char** seed, int length_seed, resultPresence** junction_seed, uint32_t traversed_part, bool verify_prefix){

    ASSERT_NULL_PTR(graph_iupac, "find_kmer()\n")
    ASSERT_NULL_PTR(graph_no_iupac, "find_kmer()\n")
    ASSERT_NULL_PTR(seed, "find_kmer()\n")
    ASSERT_NULL_PTR(*seed, "find_kmer()\n")
    ASSERT_NULL_PTR(partition_no_iupac, "find_kmer()\n")
    ASSERT_NULL_PTR(partition_iupac, "find_kmer()\n")
    ASSERT_NULL_PTR(pref_suf_pos, "find_kmer()\n")

    PWord_t PValue_partition_iupac;

    bool prefix_is_IUPAC;
    bool kmer_found;

    int z;
    int nb_bytes_seed = CEIL(length_seed * 2, SIZE_BITS_UINT_8T);
    int shift_pref_suf = 0;

    uint64_t prefix;
    uint64_t suffix;

    uint8_t ret = 3;

    uint8_t* seed_comp;

    if (length_seed % 4) shift_pref_suf = SIZE_BITS_UINT_8T - ((length_seed % 4) * 2);

    // Check if the prefix is present in the partition

    prefix_is_IUPAC = is_substring_IUPAC(*seed);

    if (verify_prefix){

        if (prefix_is_IUPAC){

            JSLG(PValue_partition_iupac, *partition_iupac, (uint8_t*) *seed);
            if (PValue_partition_iupac != NULL) ret = 1;
        }
        else{

            seed_comp = calloc(nb_bytes_seed,sizeof(char));
            ASSERT_NULL_PTR(seed_comp, "decompress() 1")

            parseKmerCount(*seed, length_seed, seed_comp, 0);

            for (z = 0, prefix = 0; z < nb_bytes_seed; z++)
                prefix = (prefix << SIZE_BITS_UINT_8T) | reverse_word_8(seed_comp[z]);

            prefix >>= shift_pref_suf;

            if (partition_no_iupac[prefix/SIZE_BITS_UINT_8T] & MASK_POWER_8[prefix%SIZE_BITS_UINT_8T]) ret = 1;

            free(seed_comp);
        }
    }

    // If the prefix is not present in the partition, look for a k-mer with the corresponding
    // prefix and partition

    if (ret == 3){

        if (*junction_seed == NULL) kmer_found = cdbg_get_successor_overlap(graph_iupac, true, seed, traversed_part);

        if (!kmer_found && !prefix_is_IUPAC)
            kmer_found = cdbg_get_successor_overlap_custom(graph_no_iupac, false, seed, junction_seed, traversed_part);

        if (kmer_found){ // Such a k-mer is found

            // Check if the suffix is the partition

            if (is_substring_IUPAC(&((*seed)[graph_no_iupac->k - length_seed]))){

                JSLG(PValue_partition_iupac, *partition_iupac, (uint8_t*) &((*seed)[graph_no_iupac->k - length_seed]));
                if (PValue_partition_iupac != NULL) ret = 2;
            }
            else{

                seed_comp = calloc(nb_bytes_seed,sizeof(char));
                ASSERT_NULL_PTR(seed_comp, "decompress() 1")

                parseKmerCount(&((*seed)[graph_no_iupac->k - length_seed]), length_seed, seed_comp, 0);

                for (z = 0, suffix = 0; z < nb_bytes_seed; z++)
                    suffix = (suffix << SIZE_BITS_UINT_8T) | reverse_word_8(seed_comp[z]);

                suffix >>= shift_pref_suf;

                if (partition_no_iupac[suffix/SIZE_BITS_UINT_8T] & MASK_POWER_8[suffix%SIZE_BITS_UINT_8T]) ret = 2;

                free(seed_comp);
            }
        }
        else ret = 0;
    }

    return ret;
}

void insert_overlap(char* overlap, int length_overlap, uint8_t* partition_no_iupac, Pvoid_t* partition_iupac, Pvoid_t* pref_suf_pos){

    ASSERT_NULL_PTR(overlap, "insert_overlap()\n")
    ASSERT_NULL_PTR(partition_no_iupac, "insert_overlap()\n")
    ASSERT_NULL_PTR(partition_iupac, "insert_overlap()\n")
    ASSERT_NULL_PTR(pref_suf_pos, "insert_overlap()\n")

    PWord_t PValue_partition_iupac;

    int z;
    int Rc_int;
    int shift_pref_suf = 0;
    int nb_bytes_overlap = CEIL(length_overlap * 2, SIZE_BITS_UINT_8T);

    char saved_char;

    uint8_t* overlap_comp;

    uint64_t overlap_int;

    if (length_overlap % 4) shift_pref_suf = SIZE_BITS_UINT_8T - ((length_overlap % 4) * 2);

    saved_char = overlap[length_overlap];
    overlap[length_overlap] = '\0';

    if (is_substring_IUPAC(overlap)) JSLI(PValue_partition_iupac, *partition_iupac, (uint8_t*) overlap)
    else{

        overlap_comp = calloc(nb_bytes_overlap, sizeof(uint8_t));
        ASSERT_NULL_PTR(overlap_comp, "insert_overlap()\n")

        parseKmerCount(overlap, length_overlap, overlap_comp, 0);

        for (z = 0, overlap_int = 0; z < nb_bytes_overlap; z++)
            overlap_int = (overlap_int << SIZE_BITS_UINT_8T) | reverse_word_8(overlap_comp[z]);

        overlap_int >>= shift_pref_suf;

        partition_no_iupac[overlap_int/SIZE_BITS_UINT_8T] |= MASK_POWER_8[overlap_int%SIZE_BITS_UINT_8T];
        J1S(Rc_int, *pref_suf_pos, overlap_int/SIZE_BITS_UINT_8T);

        free(overlap_comp);
    }

    overlap[length_overlap] = saved_char;

    return;
}

void decompress2(char* output_prefix, bool pair_ended, int size_kmer, int size_seed){

    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    FILE* file_runs;
    FILE* file_parts;
    FILE* file_recycl1;
    FILE* file_recycl2;
    FILE* file_shifts;
    FILE* file_reads;
    FILE* file_read_info;
    FILE* file_super_read_pos;
    FILE* file_super_read_length;

    BFT_Root* graph_no_iupac;
    BFT_Root* graph_iupac;

    Pvoid_t prefix_part_recycl = (Pvoid_t) NULL;
    Pvoid_t partition_iupac = (Pvoid_t) NULL;
    Pvoid_t pref_suf_pos = (Pvoid_t) NULL;
    PWord_t PValue_prefix_part_recycl;
    Word_t Rc_word;

    bool comp_shifts;
    bool recycled_part;
    bool is_kmer_run = false;
    bool verify_prefix = true;

    int nb_kmers_in_read;
    int curr_read_length;

    int size_annot_seeds = CEIL((int) pow(4,size_seed), SIZE_BITS_UINT_8T);

    int len_output_prefix = strlen(output_prefix);
    int remaining = size_kmer - size_seed;

    size_t size_read_buffer;

    uint8_t header_shift;
    uint8_t kmer_found;

    uint64_t pos_kmer_run;

    uint64_t nb_kmers = 0;
    uint64_t nb_kmers_in_part = 0;
    uint64_t nb_kmers_in_subpart = 0;

    uint32_t nb_kmers_recycled = 0;
    uint32_t nb_reads = 0;

    uint32_t pos_read;
    uint32_t pos_read_tmp;
    uint32_t read_length;
    uint32_t read_length_tmp;

    uint32_t traversed_part;
    uint32_t smallest_part = 0xffffffff;
    uint32_t largest_part = 0;

    uint32_t it_buf_read_info = 0;
    uint32_t it_part_loaded = 1;
    uint32_t it_part_recycled = 1;

    uint32_t* loaded_recycl_part;
    uint32_t* loaded_part;

    uint8_t* buf_read_info;
    uint8_t* partition_no_iupac;

    char eol = '\0';
    char nl = '\n';

    char* extracted_subseq;
    char* read;
    char* read_rc;

    char* suffix;
    char* kmer_run;
    char* traversed_kmer;
    char* traversed_kmer_tmp;

    char buff_int[64];

    char* output_prefix_buffer = malloc((len_output_prefix + 50) * sizeof(char));
    ASSERT_NULL_PTR(output_prefix_buffer, "decompress() 3")

    strcpy(output_prefix_buffer, output_prefix);

    //Open files necessary to extract the subgraph
    strcpy(&(output_prefix_buffer[len_output_prefix]), "_graph_no_iupac");
    graph_no_iupac = read_BFT_Root_sparse(output_prefix_buffer);
    if (remove(output_prefix_buffer)) printf("Warning: Could not remove temporary file.\n");

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_graph_iupac");
    graph_iupac = read_BFT_Root_sparse(output_prefix_buffer);
    if (remove(output_prefix_buffer)) printf("Warning: Could not remove temporary file.\n");

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_partitions");
    file_parts = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_parts,"decompress() 4")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_recycl_parts2");
    file_recycl2 = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_recycl2,"decompress() 5")

    traversed_kmer_tmp = calloc(size_kmer + 1, sizeof(char));
    ASSERT_NULL_PTR(traversed_kmer_tmp, "decompress() 1")

    //Extract the subgraph
    /*unserialize_subpaths_recycling(file_recycl2, &prefix_part_recycl, size_seed,
                                   1, 0xffffffff, false);

    sort_subpaths_recycling(&prefix_part_recycl, size_seed, true);

    create_subgraph_decomp(output_prefix, file_parts, &prefix_part_recycl,
                           &graph_no_iupac, &graph_iupac, size_seed);

    rewind(file_parts);

    traversed_kmer_tmp[0] = '\0';
    JSLF(PValue_prefix_part_recycl, prefix_part_recycl, (uint8_t*)traversed_kmer_tmp);

    while (PValue_prefix_part_recycl != NULL){
        free((uint32_t*)*PValue_prefix_part_recycl);
        JSLN(PValue_prefix_part_recycl, prefix_part_recycl, (uint8_t*)traversed_kmer_tmp);
    }

    JSLFA(Rc_word, prefix_part_recycl);*/

    unserialize_subpaths_recycling(file_recycl2, &prefix_part_recycl, size_seed,
                                   1, 0xffffffff, false);

    partition_no_iupac = calloc(size_annot_seeds, sizeof(uint8_t));
    ASSERT_NULL_PTR(partition_no_iupac, "compressKmers_from_KmerFiles()")

    buf_read_info = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(buf_read_info, "decompress() 2.2")

    kmer_run = malloc((size_kmer + 1) * sizeof(char));
    ASSERT_NULL_PTR(kmer_run, "decompress() 1")

    traversed_kmer = calloc(size_kmer + 1, sizeof(char));
    ASSERT_NULL_PTR(traversed_kmer, "decompress() 1")

    suffix = malloc((size_seed + 1) * sizeof(char));
    ASSERT_NULL_PTR(suffix, "decompress() 3")

    loaded_part = malloc((DECOMP_NB_MAX_LOAD_PART + 1) * sizeof(uint32_t));
    ASSERT_NULL_PTR(loaded_part, "decompress() 2")

    loaded_part[0] = DECOMP_NB_MAX_LOAD_PART;

    //Open remaining files
    strcpy(&(output_prefix_buffer[len_output_prefix]), ".super_reads");
    file_reads = fopen(output_prefix_buffer, "w");
    ASSERT_NULL_PTR(file_reads,"decompress() 6")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_runs");
    file_runs = fopen(output_prefix_buffer, "r");
    ASSERT_NULL_PTR(file_runs,"decompress() 4.5")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_recycl_parts1");
    file_recycl1 = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_recycl1,"decompress() 5")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_read_info");
    file_read_info = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_read_info,"decompress() 6")

    fread(buf_read_info, sizeof(uint8_t), SIZE_BUFFER, file_read_info);

    if (fread(kmer_run, sizeof(char), size_kmer+1, file_runs) == size_kmer+1){
        if (fread(&pos_kmer_run, sizeof(uint64_t), 1, file_runs) != 1)
            ERROR("decompress(): truncated autoloop k-mer file.\n")
    }
    else pos_kmer_run = 0xffffffffffffffff;

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_shifts");
    comp_shifts = init_extract_subseq_shifted(output_prefix_buffer, &file_shifts);

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_pos");
    file_super_read_pos = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_super_read_pos,"decompress() 7")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_length");
    file_super_read_length = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_super_read_length,"decompress() 8")

    if (fread(&read_length_tmp, sizeof(uint32_t), 1, file_super_read_length) != 1)
        ERROR("decompress(): truncated spanning super read length file.\n")

    if (fread(&pos_read_tmp, sizeof(uint32_t), 1, file_super_read_pos) != 1)
        ERROR("decompress(): truncated spanning super read length file.\n")

    read_length_tmp &= 0xfffffffe;

    read_length = read_length_tmp;
    pos_read = pos_read_tmp;

    while ((read_length_tmp & 0x1) == 0){

        if (fread(&read_length_tmp, sizeof(uint32_t), 1, file_super_read_length) == 1){

            if (fread(&pos_read_tmp, sizeof(uint32_t), 1, file_super_read_pos) != 1)
                ERROR("decompress(): truncated spanning super read position file.\n")

            if ((read_length_tmp & 0x1) == 0){
                read_length = read_length_tmp;
                pos_read += pos_read_tmp;
            }
        }
        else break;
    }

    fseek(file_super_read_length, 0 - ((int) sizeof(uint32_t)), SEEK_CUR);
    fseek(file_super_read_pos, 0 - ((int) sizeof(uint32_t)), SEEK_CUR);

    read_length = pos_read + (read_length >> 1);
    size_read_buffer = read_length;

    //printf("read_length = %" PRIu32 "\n", read_length);

    read = malloc((size_read_buffer + 2) * sizeof(char));
    ASSERT_NULL_PTR(read, "decompress() 3")

    read_rc = malloc((size_read_buffer + 2) * sizeof(char));
    ASSERT_NULL_PTR(read_rc, "decompress() 3")

    header_shift = get_header_next_subseq_shifted(file_shifts);

    extracted_subseq = extract_next_subseq_shifted(file_shifts, comp_shifts,
                                                   (bool) ((header_shift & 0x80) >> 7),
                                                   header_shift & 0x3f);

    curr_read_length = strlen(extracted_subseq);
    //printf("extracted_subseq: %s\n", extracted_subseq);

    memcpy(traversed_kmer, &extracted_subseq[curr_read_length - size_seed], size_seed * sizeof(uint8_t));
    traversed_kmer[size_kmer] = eol;

    memcpy(read, extracted_subseq, curr_read_length * sizeof(uint8_t));

    nb_kmers_in_read = read_length / remaining;

    if (nb_kmers_in_read * remaining + size_seed >= read_length){
        if ((curr_read_length != size_seed) || (nb_kmers_in_read * remaining + size_seed != read_length)) nb_kmers_in_read--;
    }

    //printf("nb_kmers_in_read = %d\n", nb_kmers_in_read);

    free(extracted_subseq);

    suffix[0] = '\0';

    resultPresence* junction_prefix = NULL;

    while (!feof(file_parts)){

        //printf("New loading\n");

        loaded_part[0] = fread(&loaded_part[1], sizeof(uint32_t), loaded_part[0], file_parts);

        it_part_loaded = 1;
        traversed_part = loaded_part[it_part_loaded];

        sprintf(buff_int, "_part%" PRIu32, traversed_part);
        strcpy(&(output_prefix_buffer[len_output_prefix]), buff_int);

        read_partition(output_prefix_buffer, partition_no_iupac, &partition_iupac, &pref_suf_pos, size_seed);

        smallest_part = MIN(smallest_part, traversed_part);
        largest_part = MAX(largest_part, traversed_part);

        recycled_part = false;

        while ((it_part_loaded <= loaded_part[0]) || recycled_part){

            nb_kmers_in_subpart = 0;

            strcpy(traversed_kmer_tmp, traversed_kmer);
            memset(&traversed_kmer_tmp[size_seed], 0, remaining * sizeof(char));

            if (recycled_part){

                JSLG(PValue_prefix_part_recycl, prefix_part_recycl, (uint8_t*) traversed_kmer_tmp);

                if (PValue_prefix_part_recycl == PJERR) ERROR("decompress(): Something went wrong with the Judy Array...\n")
                if (PValue_prefix_part_recycl == NULL) ERROR("decompress(): Could not find recycled partition.\n")

                loaded_recycl_part = (uint32_t*) *PValue_prefix_part_recycl;
                traversed_part = loaded_recycl_part[1];

                loaded_recycl_part[0]--;

                if (loaded_recycl_part[0]) memmove(&loaded_recycl_part[1], &loaded_recycl_part[2], loaded_recycl_part[0] * sizeof(uint32_t));
                else {
                    free(loaded_recycl_part);
                    JSLD(PValue_prefix_part_recycl, prefix_part_recycl, (uint8_t*) traversed_kmer_tmp);
                }
            }

            /*while ((kmer_found = find_kmer_custom(graph_iupac, graph_no_iupac, partition_no_iupac, &partition_iupac,
                                           &pref_suf_pos, &traversed_kmer_tmp, size_seed, &junction_prefix, traversed_part, verify_prefix)) == 3){

                free(junction_prefix);
                junction_prefix = NULL;*/

            while ((kmer_found = find_kmer(graph_iupac, graph_no_iupac, partition_no_iupac, &partition_iupac,
                                           &pref_suf_pos, &traversed_kmer_tmp, size_seed, traversed_part, verify_prefix)) == 3){

                if (strlen(suffix)){
                    insert_overlap(suffix, size_seed, partition_no_iupac, &partition_iupac, &pref_suf_pos);
                    suffix[0] = '\0';
                }

                insert_overlap(traversed_kmer_tmp, size_seed, partition_no_iupac, &partition_iupac, &pref_suf_pos);

                nb_kmers_in_subpart++;

                NO_AUTOLOOP_KMER: memcpy(&read[curr_read_length], &traversed_kmer_tmp[size_seed], remaining * sizeof(char));
                curr_read_length += remaining;

                //printf("recycled_part: %s, ", recycled_part ? "true" : "false");
                //printf("Good: %s, partition %" PRIu32 ", nb_kmers_in_subpart: %" PRIu64 "\n", traversed_kmer_tmp, traversed_part, nb_kmers_in_subpart);

                nb_kmers_in_read--;
                nb_kmers++;

                if (nb_kmers_in_read){

                    if ((nb_kmers == pos_kmer_run) && (!is_kmer_run)) strcpy(suffix, &traversed_kmer_tmp[remaining]);

                    memmove(traversed_kmer_tmp, &traversed_kmer_tmp[remaining], size_seed * sizeof(char));
                    memset(&traversed_kmer_tmp[size_seed], 0, remaining * sizeof(char));
                    memcpy(traversed_kmer, traversed_kmer_tmp, size_seed * sizeof(char));
                }
                else{

                    if (!is_kmer_run) strcpy(suffix, &traversed_kmer_tmp[remaining]);

                    extracted_subseq = extract_next_subseq_shifted(file_shifts, comp_shifts,
                                                                   (bool) ((header_shift & 0x40) >> 6),
                                                                   read_length - curr_read_length);

                    //printf("extracted_subseq: %s, strlen(extracted_subseq): %d\n", extracted_subseq, strlen(extracted_subseq));
                    //printf("length to extract = %d\n", read_length - curr_read_length);

                    if (extracted_subseq != NULL){

                        strcpy(&read[curr_read_length], extracted_subseq);
                        free(extracted_subseq);

                        read[read_length] = nl;
                        read[read_length+1] = eol;
                    }
                    else break;

                    if ((buf_read_info[it_buf_read_info/SIZE_BITS_UINT_8T] >> (SIZE_BITS_UINT_8T - 1 - (it_buf_read_info % SIZE_BITS_UINT_8T))) & 0x1){

                        reverse_complement(read, read_rc, read_length);
                        read_rc[read_length] = nl;
                        read_rc[read_length+1] = eol;

                        //printf("read rc: %s", read_rc);

                        if (fwrite(read_rc, sizeof(char), read_length+1, file_reads) != read_length+1)
                            ERROR("decompress(): Cannot write reads in output file.\n")
                    }
                    else {

                        //printf("read: %s", read);

                        if (fwrite(read, sizeof(char), read_length+1, file_reads) != read_length+1)
                            ERROR("decompress(): Cannot write reads in output file.\n")
                    }

                    header_shift = get_header_next_subseq_shifted(file_shifts);

                    extracted_subseq = extract_next_subseq_shifted(file_shifts, comp_shifts, (bool) ((header_shift & 0x80) >> 7),
                                                                   header_shift & 0x3f);

                    if ((extracted_subseq != NULL) && (header_shift & 0x3f)){

                        //printf("extracted_subseq: %s, strlen(extracted_subseq): %d\n", extracted_subseq, strlen(extracted_subseq));
                        //printf("recycled_part: %s\n", recycled_part ? "true" : "false");
                        //if (recycled_part) printf("nb_kmers_in_subpart: %d, nb_kmers_recycled: %d\n", nb_kmers_in_subpart, nb_kmers_recycled);

                        nb_reads++;
                        it_buf_read_info++;

                        //printf("nb_reads: %" PRIu32 "\n", nb_reads);

                        if (it_buf_read_info >= SIZE_BUFFER * SIZE_BITS_UINT_8T){
                            fread(buf_read_info, sizeof(uint8_t), SIZE_BUFFER, file_read_info);
                            it_buf_read_info = 0;
                        }

                        if (fread(&read_length_tmp, sizeof(uint32_t), 1, file_super_read_length) != 1)
                            ERROR("decompress(): truncated spanning super read length file. 1\n")

                        if (fread(&pos_read_tmp, sizeof(uint32_t), 1, file_super_read_pos) != 1)
                            ERROR("decompress(): truncated spanning super read position file.\n")

                        read_length_tmp &= 0xfffffffe;

                        read_length = read_length_tmp;
                        pos_read = pos_read_tmp;

                        while ((read_length_tmp & 0x1) == 0){

                            if (fread(&read_length_tmp, sizeof(uint32_t), 1, file_super_read_length) == 1){

                                if (fread(&pos_read_tmp, sizeof(uint32_t), 1, file_super_read_pos) != 1)
                                    ERROR("decompress(): truncated spanning super read position file.\n")

                                if ((read_length_tmp & 0x1) == 0){
                                    read_length = read_length_tmp;
                                    pos_read += pos_read_tmp;
                                }
                            }
                            else break;
                        }

                        fseek(file_super_read_length, 0 - ((int) sizeof(uint32_t)), SEEK_CUR);
                        fseek(file_super_read_pos, 0 - ((int) sizeof(uint32_t)), SEEK_CUR);

                        read_length = pos_read + (read_length >> 1);

                        if (read_length > size_read_buffer){

                            read = realloc(read, (read_length + 2) * sizeof(char));
                            ASSERT_NULL_PTR(read, "decompress() 3")

                            read_rc = realloc(read_rc, (read_length + 2) * sizeof(char));
                            ASSERT_NULL_PTR(read_rc, "decompress() 3")

                            size_read_buffer = read_length;
                        }

                        curr_read_length = strlen(extracted_subseq);
                        memcpy(read, extracted_subseq, curr_read_length * sizeof(uint8_t));

                        nb_kmers_in_read = read_length / remaining;

                        if (nb_kmers_in_read * remaining + size_seed >= read_length){
                            if ((curr_read_length != size_seed) || (nb_kmers_in_read * remaining + size_seed != read_length)) nb_kmers_in_read--;
                        }

                        //printf("read_length = %" PRIu32 "\n", read_length);
                        //printf("nb_kmers_in_read = %d\n", nb_kmers_in_read);

                        memcpy(traversed_kmer_tmp, &extracted_subseq[curr_read_length - size_seed], size_seed * sizeof(char));
                        memset(&traversed_kmer_tmp[size_seed], 0, remaining * sizeof(char));
                        memcpy(traversed_kmer, traversed_kmer_tmp, size_seed * sizeof(char));

                        free(extracted_subseq);
                    }
                    else {

                        free(extracted_subseq);
                        break;
                    }
                }

                if (nb_kmers == pos_kmer_run){

                    //printf("Insert this kmer, partition = %" PRIu32 "\n", traversed_part);

                    memcpy(traversed_kmer_tmp, kmer_run, size_kmer * sizeof(char));
                    traversed_kmer_tmp[size_kmer] = eol;

                    if (fread(kmer_run, sizeof(char), size_kmer+1, file_runs) == size_kmer+1){
                        if (fread(&pos_kmer_run, sizeof(uint64_t), 1, file_runs) != 1)
                            ERROR("decompress(): truncated autoloop k-mer file.\n")
                    }
                    else pos_kmer_run = 0xffffffffffffffff;

                    is_kmer_run = true;

                    goto NO_AUTOLOOP_KMER;
                }

                is_kmer_run = false;
                verify_prefix = true;

                if (recycled_part){

                    if ((nb_kmers_in_subpart == nb_kmers_recycled)) break;
                    else if (strlen(suffix) && memcmp(traversed_kmer_tmp, suffix, size_seed * sizeof(char))){

                        JSLG(PValue_prefix_part_recycl, prefix_part_recycl, (uint8_t*) traversed_kmer_tmp);

                        if (PValue_prefix_part_recycl == PJERR) ERROR("decompress(): Something went wrong with the Judy Array...\n")
                        if (PValue_prefix_part_recycl == NULL) ERROR("decompress(): Could not find recycled partition.\n")

                        loaded_recycl_part = (uint32_t*) *PValue_prefix_part_recycl;
                        traversed_part = loaded_recycl_part[1];

                        loaded_recycl_part[0]--;

                        if (loaded_recycl_part[0]) memmove(&loaded_recycl_part[1], &loaded_recycl_part[2], loaded_recycl_part[0] * sizeof(uint32_t));
                        else {
                            free(loaded_recycl_part);
                            JSLD(PValue_prefix_part_recycl, prefix_part_recycl, (uint8_t*) traversed_kmer_tmp);
                        }

                        //printf("Trigger recycling 1\n");
                    }
                }
            }

            nb_kmers_in_part += nb_kmers_in_subpart;

            if (recycled_part){

                if (nb_kmers_in_subpart == 0) ERROR("No kmers in this subpartition\n")
                if (nb_kmers_in_subpart != nb_kmers_recycled) ERROR("Could not find enough recycled k-mers\n")

                recycled_part = false;
                traversed_part = loaded_part[it_part_loaded];

                //printf("Trigger non-recycling\n");
            }
            else {

                if (fread(&nb_kmers_recycled, sizeof(uint32_t), 1, file_recycl1) == 1){

                    if (nb_kmers_recycled & 0x1) it_part_recycled += nb_kmers_recycled >> 1;
                    else fseek(file_recycl1, 0 - ((int) sizeof(uint32_t)), SEEK_CUR);
                }
                else it_part_recycled = 0xffffffff;

                if (it_part_loaded == it_part_recycled){

                    if (fread(&nb_kmers_recycled, sizeof(uint32_t), 1, file_recycl1) != 1) ERROR("Could not read length of recycled subpath.\n");

                    nb_kmers_recycled >>= 1;
                    recycled_part = true;
                    verify_prefix = false;

                    //printf("Trigger recycling 2, nb_kmers_recycled: %d\n", nb_kmers_recycled);
                }
                else {

                    //printf("New part\n");

                    if (nb_kmers_in_part == 0) ERROR("No kmers in this subpartition\n")

                    if (strlen(suffix)){
                        insert_overlap(suffix, size_seed, partition_no_iupac, &partition_iupac, &pref_suf_pos);
                        suffix[0] = '\0';
                    }

                    insert_overlap(traversed_kmer_tmp, size_seed, partition_no_iupac, &partition_iupac, &pref_suf_pos);

                    sprintf(buff_int, "_part%" PRIu32, traversed_part);
                    strcpy(&(output_prefix_buffer[len_output_prefix]), buff_int);

                    write_partition(output_prefix_buffer, partition_no_iupac, &partition_iupac, size_seed);
                    erase_partitions(partition_no_iupac, &partition_iupac, &pref_suf_pos);

                    it_part_loaded++;

                    if (it_part_loaded <= loaded_part[0]){

                        traversed_part = loaded_part[it_part_loaded];

                        sprintf(buff_int, "_part%" PRIu32, traversed_part);
                        strcpy(&(output_prefix_buffer[len_output_prefix]), buff_int);

                        read_partition(output_prefix_buffer, partition_no_iupac, &partition_iupac, &pref_suf_pos, size_seed);

                        smallest_part = MIN(smallest_part, traversed_part);
                        largest_part = MAX(largest_part, traversed_part);
                    }
                    else it_part_recycled -= loaded_part[0];
                }
            }
        }

        //printf("Number of k-mers extracted = %" PRIu64 "\n", nb_kmers);

        gettimeofday(&tval_after, NULL);
        time_spent(&tval_before, &tval_after, &tval_result);
        //printf("\nElapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
    }

    fseek(file_reads, 0 - ((int) sizeof(char)), SEEK_CUR);
    if (fwrite(&eol, sizeof(char), 1, file_reads) != 1) ERROR("decompress(): Cannot write reads in output file.\n")

    fclose(file_reads);

    fclose(file_read_info);
    fclose(file_runs);
    fclose(file_parts);
    fclose(file_recycl1);
    fclose(file_recycl2);
    fclose(file_super_read_pos);
    fclose(file_super_read_length);

    deinit_extract_subseq_shifted(&file_shifts);

    freeBFT_Root(graph_iupac);
    freeBFT_Root(graph_no_iupac);

    traversed_kmer_tmp[0] = '\0';
    JSLF(PValue_prefix_part_recycl, prefix_part_recycl, (uint8_t*)traversed_kmer_tmp);

    while (PValue_prefix_part_recycl != NULL){
        free((uint32_t*)*PValue_prefix_part_recycl);
        JSLN(PValue_prefix_part_recycl, prefix_part_recycl, (uint8_t*)traversed_kmer_tmp);
    }

    JSLFA(Rc_word, prefix_part_recycl);

    J1FA(Rc_word, pref_suf_pos);
    JSLFA(Rc_word, partition_iupac);

    free(read);
    free(suffix);
    free(read_rc);
    free(kmer_run);
    free(buf_read_info);
    free(traversed_kmer);
    free(traversed_kmer_tmp);
    free(partition_no_iupac);

    free(loaded_part);

    free(junction_prefix);

    for (; smallest_part <= largest_part; smallest_part++){

        sprintf(buff_int, "_part%" PRIu32, smallest_part);
        strcpy(&(output_prefix_buffer[len_output_prefix]), buff_int);

        remove(output_prefix_buffer);
    }

    free(output_prefix_buffer);

    gettimeofday(&tval_after, NULL);
    time_spent(&tval_before, &tval_after, &tval_result);
    //printf("\nElapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
}

uint32_t unserialize_subpaths_recycling(FILE* file_part_recycl2, Pvoid_t* recycled_parts, int size_seed,
                                    int64_t load_from_pos_start, int64_t nb_parts_to_load, bool seeds_already_loaded){

    ASSERT_NULL_PTR(file_part_recycl2, "serialize_subpaths_recycling()\n")
    ASSERT_NULL_PTR(recycled_parts, "serialize_subpaths_recycling()\n")

    PWord_t PValue;

    int l;
    int seed_shift = 0;
    int seed_shift_iupac = 0;
    int nb_bytes_seed = CEIL(size_seed * 2, SIZE_BITS_UINT_8T);
    int nb_bytes_seed_iupac = CEIL(size_seed * 4, SIZE_BITS_UINT_8T);

    uint32_t j;
    uint32_t seed;
    uint32_t recycled_part_size;
    uint32_t recycled_part_size_to_load;

    uint32_t max_size_array = 0;

    uint32_t size_seed_present = CEIL((uint32_t) pow(4, size_seed), SIZE_BITS_UINT_8T);

    uint32_t* recycled_part;

    uint64_t seed_iupac;
    uint64_t prev_seed_iupac = 0;

    uint8_t tmp;

    uint8_t* seed_ascii = malloc((size_seed + 1) * sizeof(uint8_t));
    ASSERT_NULL_PTR(seed_ascii, "serialize_subpaths_recycling()\n")

    uint8_t* seed_present;
    uint8_t* seed_comp;

    rewind(file_part_recycl2);

    if (seeds_already_loaded){

        fseek(file_part_recycl2, size_seed_present, SEEK_SET);

        seed_ascii[0] = '\0';

        JSLF(PValue, *recycled_parts, seed_ascii);

        while (PValue != NULL){

            if (is_substring_IUPAC((char*) seed_ascii) == 0){

                recycled_part = (uint32_t*) *PValue;

                fread(&recycled_part_size, sizeof(uint32_t), 1, file_part_recycl2);

                max_size_array = MAX(max_size_array, recycled_part_size);

                recycled_part_size_to_load = MIN(nb_parts_to_load, MAX(0, ((int64_t) recycled_part_size) - load_from_pos_start + 1));

                recycled_part = realloc(recycled_part, (recycled_part[0] + recycled_part_size_to_load + 1) * sizeof(uint32_t));
                ASSERT_NULL_PTR(recycled_part, "serialize_subpaths_recycling()\n")

                fseek(file_part_recycl2, MIN(load_from_pos_start-1, recycled_part_size) * sizeof(uint32_t), SEEK_CUR);

                fread(&recycled_part[recycled_part[0]+1], sizeof(uint32_t), recycled_part_size_to_load, file_part_recycl2);

                fseek(file_part_recycl2,
                      MAX(0, ((int64_t) recycled_part_size) - load_from_pos_start - nb_parts_to_load + 1) * sizeof(uint32_t),
                      SEEK_CUR);

                recycled_part[0] += recycled_part_size_to_load;

                *PValue = (Word_t) recycled_part;
            }

            JSLN(PValue, *recycled_parts, seed_ascii);
        }

        seed_ascii[0] = '\0';

        JSLF(PValue, *recycled_parts, seed_ascii);

        while (PValue != NULL){

            if (is_substring_IUPAC((char*) seed_ascii)){

                recycled_part = (uint32_t*) *PValue;

                fread(&seed_iupac, sizeof(uint64_t), 1, file_part_recycl2);
                fread(&recycled_part_size, sizeof(uint32_t), 1, file_part_recycl2);

                max_size_array = MAX(max_size_array, recycled_part_size);

                recycled_part_size_to_load = MIN(nb_parts_to_load, MAX(0, ((int64_t) recycled_part_size) - load_from_pos_start + 1));

                recycled_part = realloc(recycled_part, (recycled_part[0] + recycled_part_size_to_load + 1) * sizeof(uint32_t));
                ASSERT_NULL_PTR(recycled_part, "serialize_subpaths_recycling()\n")

                fseek(file_part_recycl2, MIN(load_from_pos_start-1, recycled_part_size) * sizeof(uint32_t), SEEK_CUR);

                fread(&recycled_part[recycled_part[0]+1], sizeof(uint32_t), recycled_part_size_to_load, file_part_recycl2);

                fseek(file_part_recycl2,
                      MAX(0, ((int64_t) recycled_part_size) - load_from_pos_start - nb_parts_to_load + 1) * sizeof(uint32_t),
                      SEEK_CUR);

                recycled_part[0] += recycled_part_size_to_load;

                *PValue = (Word_t) recycled_part;
            }

            JSLN(PValue, *recycled_parts, seed_ascii);
        }
    }
    else{

        seed_present = calloc(size_seed_present, sizeof(uint8_t));
        ASSERT_NULL_PTR(seed_present, "serialize_subpaths_recycling()\n")

        seed_comp = malloc(nb_bytes_seed_iupac * sizeof(uint8_t));
        ASSERT_NULL_PTR(seed_comp, "serialize_subpaths_recycling()\n")

        if (size_seed % 4) seed_shift = SIZE_BITS_UINT_8T - ((size_seed % 4) * 2);
        if (size_seed % 2) seed_shift_iupac = 4;

        fread(seed_present, sizeof(uint8_t), size_seed_present, file_part_recycl2);

        for (uint32_t i = 0; i < size_seed_present; i++){

            if (seed_present[i]){

                for (j = 0, tmp = seed_present[i]; j < SIZE_BITS_UINT_8T; j++, tmp >>= 1){

                    if (tmp & 0x1){

                        seed = (i * SIZE_BITS_UINT_8T + j) << seed_shift;

                        for (l = nb_bytes_seed-1; l >= 0; l--, seed >>= SIZE_BITS_UINT_8T)
                            seed_comp[l] = reverse_word_8(seed & 0xff);

                        kmer_comp_to_ascii(seed_comp, size_seed, (char*) seed_ascii);

                        JSLI(PValue, *recycled_parts, seed_ascii);

                        if (PValue == PJERR)
                            ERROR("unserialize_subpaths_recycling(): Something went wrong with the JudyArray...\n")

                        if (PValue == NULL)
                            ERROR("unserialize_subpaths_recycling(): Seed already seen, therefore there is a bug in the code. Do not panic.\n")

                        fread(&recycled_part_size, sizeof(uint32_t), 1, file_part_recycl2);

                        max_size_array = MAX(max_size_array, recycled_part_size);

                        recycled_part_size_to_load = MIN(nb_parts_to_load, MAX(0, ((int64_t) recycled_part_size) - load_from_pos_start + 1));

                        recycled_part = malloc((recycled_part_size_to_load + 1) * sizeof(uint32_t));
                        ASSERT_NULL_PTR(recycled_part, "unserialize_subpaths_recycling()")

                        fseek(file_part_recycl2, MIN(load_from_pos_start-1, recycled_part_size) * sizeof(uint32_t), SEEK_CUR);

                        fread(&recycled_part[1], sizeof(uint32_t), recycled_part_size_to_load, file_part_recycl2);

                        fseek(file_part_recycl2,
                              MAX(0, ((int64_t) recycled_part_size) - load_from_pos_start - nb_parts_to_load + 1) * sizeof(uint32_t),
                              SEEK_CUR);

                        recycled_part[0] = recycled_part_size_to_load;

                        *PValue = (Word_t) recycled_part;
                    }
                }
            }
        }

        free(seed_present);

        while (fread(&seed_iupac, sizeof(uint64_t), 1, file_part_recycl2) == 1){

            seed_iupac += prev_seed_iupac;
            prev_seed_iupac = seed_iupac;

            seed_iupac <<= seed_shift_iupac;

            for (l = nb_bytes_seed_iupac-1; l >= 0; l--, seed_iupac >>= SIZE_BITS_UINT_8T)
                seed_comp[l] = reverse_word_8(seed_iupac & 0xff);

            kmer_iupac_comp_to_ascii(seed_comp, size_seed, (char*) seed_ascii);

            JSLI(PValue, *recycled_parts, seed_ascii);

            if (PValue == PJERR) ERROR("unserialize_subpaths_recycling(): Something went wrong with the JudyArray...\n")
            if (PValue == NULL) ERROR("unserialize_subpaths_recycling(): Seed already seen, therefore there is a bug in the code. Do not panic.\n")

            fread(&recycled_part_size, sizeof(uint32_t), 1, file_part_recycl2);

            max_size_array = MAX(max_size_array, recycled_part_size);

            recycled_part_size_to_load = MIN(nb_parts_to_load, MAX(0, ((int64_t) recycled_part_size) - load_from_pos_start + 1));

            recycled_part = malloc((recycled_part_size_to_load + 1) * sizeof(uint32_t));
            ASSERT_NULL_PTR(recycled_part, "unserialize_subpaths_recycling()")

            fseek(file_part_recycl2, MIN(load_from_pos_start-1, recycled_part_size) * sizeof(uint32_t), SEEK_CUR);

            fread(&recycled_part[1], sizeof(uint32_t), recycled_part_size_to_load, file_part_recycl2);

            fseek(file_part_recycl2,
                  MAX(0, ((int64_t) recycled_part_size) - load_from_pos_start - nb_parts_to_load + 1) * sizeof(uint32_t),
                  SEEK_CUR);

            recycled_part[0] = recycled_part_size_to_load;

            *PValue = (Word_t) recycled_part;
        }

        free(seed_comp);
    }

    free(seed_ascii);

    return max_size_array;
}

void sort_subpaths_recycling(Pvoid_t* recycled_parts, int size_seed, bool replace_by_sorted_parts){

    ASSERT_NULL_PTR(recycled_parts, "serialize_subpaths_recycling()\n")

    PWord_t PValue_prefix_part_recycl;

    uint32_t size_recycled_part_sorted = 1;

    uint32_t it_recycled_part;
    uint32_t it_recycled_part_uniq;

    uint32_t* recycled_part_sorted = calloc(1, sizeof(uint32_t));
    ASSERT_NULL_PTR(recycled_part_sorted, "decompress()")

    uint32_t* recycled_part_not_sorted;

    uint8_t* seed = calloc(size_seed + 1, sizeof(uint8_t));
    ASSERT_NULL_PTR(seed, "decompress()")

    seed[0] = '\0';

    JSLF(PValue_prefix_part_recycl, *recycled_parts, (uint8_t*) seed);

    while (PValue_prefix_part_recycl != NULL){

        recycled_part_not_sorted = (uint32_t*) *PValue_prefix_part_recycl;

        if (recycled_part_not_sorted[0] > 1){

            if (recycled_part_not_sorted[0]+1 > size_recycled_part_sorted){

                size_recycled_part_sorted = recycled_part_not_sorted[0]+1;

                recycled_part_sorted = realloc(recycled_part_sorted, size_recycled_part_sorted * sizeof(uint32_t));
                ASSERT_NULL_PTR(recycled_part_sorted, "decompress()");
            }

            recycled_part_sorted[0] = recycled_part_not_sorted[0];

            memcpy(&recycled_part_sorted[1], &recycled_part_not_sorted[1], recycled_part_sorted[0] * sizeof(uint32_t));

            qsort(&recycled_part_sorted[1], recycled_part_sorted[0], sizeof(uint32_t), comp_uint32);

            it_recycled_part = 2;
            it_recycled_part_uniq = 2;

            while (it_recycled_part <= recycled_part_sorted[0]){

                if (recycled_part_sorted[it_recycled_part] != recycled_part_sorted[it_recycled_part-1]){
                    recycled_part_sorted[it_recycled_part_uniq] = recycled_part_sorted[it_recycled_part];
                    it_recycled_part_uniq++;
                }

                it_recycled_part++;
            }

            recycled_part_sorted[0] = it_recycled_part_uniq-1;

            if (replace_by_sorted_parts){

                recycled_part_not_sorted = realloc(recycled_part_not_sorted, (recycled_part_sorted[0] + 1) * sizeof(uint32_t));
                ASSERT_NULL_PTR(recycled_part_not_sorted, "unserialize_subpaths_recycling()")

                memcpy(recycled_part_not_sorted, recycled_part_sorted, (recycled_part_sorted[0] + 1) * sizeof(uint32_t));
            }
            else{

                recycled_part_not_sorted = realloc(recycled_part_not_sorted,
                                                   (recycled_part_not_sorted[0] + recycled_part_sorted[0] + 2) * sizeof(uint32_t));

                ASSERT_NULL_PTR(recycled_part_not_sorted, "unserialize_subpaths_recycling()")

                memcpy(&recycled_part_not_sorted[recycled_part_not_sorted[0]+1],
                       recycled_part_sorted,
                       (recycled_part_sorted[0] + 1) * sizeof(uint32_t));
            }

            *PValue_prefix_part_recycl = (Word_t) recycled_part_not_sorted;
        }

        JSLN(PValue_prefix_part_recycl, *recycled_parts, (uint8_t*) seed);
    }

    free(recycled_part_sorted);
    free(seed);

    return;
}

bool init_extract_subseq_shifted(char* filename_shifts, FILE** file_shift){

    ASSERT_NULL_PTR(filename_shifts, "init_extract_subseq_shifted()\n")
    ASSERT_NULL_PTR(file_shift, "init_extract_subseq_shifted()\n")

    uint8_t comp;

    *file_shift = fopen(filename_shifts, "r");
    ASSERT_NULL_PTR(file_shift, "init_extract_subseq_shifted()\n")

    if (fread(&comp, sizeof(uint8_t), 1, *file_shift) != 1)
        ERROR("init_extract_subseq_shifted(): Cannot read skipped subsequences file.\n")

    if (comp == 0) return false;
    return true;
}

uint8_t get_header_next_subseq_shifted(FILE* file_shift){

    ASSERT_NULL_PTR(file_shift, "get_header_next_subseq_shifted()\n")

    uint8_t header = 0;

    fread(&header, sizeof(uint8_t), 1, file_shift);

    return header;
}

char* extract_next_subseq_shifted(FILE* file_shift, bool compress_shift, bool is_iupac, int length){

    ASSERT_NULL_PTR(file_shift, "extract_next_subseq_shifted()\n")

    char* kmer = NULL;

    kmer = calloc(length + 1, sizeof(char));
    ASSERT_NULL_PTR(kmer, "extract_next_subseq_shifted()\n")

    int info;

    if (compress_shift){

        info = CEIL(length * 4, SIZE_BITS_UINT_8T);

        uint8_t* kmer_comp = calloc(info, sizeof(uint8_t));
        ASSERT_NULL_PTR(kmer_comp, "extract_next_subseq_shifted()\n")

        if (is_iupac){ //IUPAC

            if (fread(kmer_comp, sizeof(uint8_t), info, file_shift) == info) kmer_iupac_comp_to_ascii(kmer_comp, length, kmer);
            else {
                free(kmer);
                kmer = NULL;
            }
        }
        else{

            info = CEIL(length * 2, SIZE_BITS_UINT_8T);

            if (fread(kmer_comp, sizeof(uint8_t), info, file_shift) == info) kmer_comp_to_ascii(kmer_comp, length, kmer);
            else {
                free(kmer);
                kmer = NULL;
            }
        }

        free(kmer_comp);
    }
    else{

        int pos_kmer = 0;
        int c = fgetc(file_shift);

        while ((c != EOF) && (c != '\0') && (c != '\n')){
            kmer[pos_kmer] = (char) c;
            pos_kmer++;
            c = fgetc(file_shift);
        }

        if (c == EOF){
            free(kmer);
            return NULL;
        }

        kmer[pos_kmer] = '\0';
    }

    return kmer;
}

void deinit_extract_subseq_shifted(FILE** file_shift){

    ASSERT_NULL_PTR(file_shift, "deinit_extract_subseq_shifted()\n")

    fclose(*file_shift);

    return;
}

void extract_reads(char* output_prefix, bool paired_end){

    FILE* file_super_read_pos_mis;
    FILE* file_super_read_id_mis;
    FILE* file_super_read_char_mis;
    FILE* file_super_read_nb_mis;
    FILE* file_super_read_length;
    FILE* file_super_read_pos;
    FILE* file_super_read_occ;
    FILE* file_super_read_rev;
    FILE* file_super_reads;
    FILE* file_reads;
    FILE* file_reads_pos;

    bool first_read;

    int len_output_prefix = strlen(output_prefix);

    long int read_pos = 0;

    uint32_t i, it_mis;
    uint32_t position_mis, position_mis_tmp, id_mis;
    uint32_t read_length;
    uint32_t pos_read, pos_read_tmp;
    uint32_t occ_count;
    uint32_t nb_mismatches;

    int64_t nb_reads = 0;

    size_t size_buffer_read = SIZE_BUFFER;
    size_t size_read = 1;
    size_t size_mis = 1;

    char nl = '\n';
    char eol = '\0';

    Pos_length_occ* pos_length_occ;

    List* l = List_create();

    Mismatch* mis = malloc(size_mis * sizeof(Mismatch));
    ASSERT_NULL_PTR(mis, "binning_reads()\n")

    char* buffer_read = malloc(size_buffer_read * sizeof(char));
    ASSERT_NULL_PTR(buffer_read, "binning_reads()\n")

    char* buffer_rev = malloc(SIZE_BUFFER * sizeof(char));
    ASSERT_NULL_PTR(buffer_rev, "binning_reads()\n")

    char* read = malloc(size_read * sizeof(char));
    ASSERT_NULL_PTR(read, "binning_reads()\n")

    char* read_tmp = malloc(size_read * sizeof(char));
    ASSERT_NULL_PTR(read_tmp, "binning_reads()\n")

    char* read_rev = malloc(size_read * sizeof(char));
    ASSERT_NULL_PTR(read_rev, "binning_reads()\n")

    char* output_prefix_buffer = malloc((len_output_prefix + 50) * sizeof(char));
    ASSERT_NULL_PTR(output_prefix_buffer, "extract_reads()\n")

    strcpy(output_prefix_buffer, output_prefix);

    strcpy(&(output_prefix_buffer[len_output_prefix]), ".reads");
    file_reads = fopen(output_prefix_buffer, "w");
    ASSERT_NULL_PTR(file_reads,"extract_reads()\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), ".super_reads");
    file_super_reads = fopen(output_prefix_buffer, "r");
    ASSERT_NULL_PTR(file_super_reads,"extract_reads()\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_pos");
    file_super_read_pos = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_super_read_pos,"extract_reads()\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_length");
    file_super_read_length = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_super_read_length,"extract_reads()\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_occ");
    file_super_read_occ = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_super_read_occ,"extract_reads()\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_rev");
    file_super_read_rev = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_super_read_rev,"extract_reads()\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_nb_mis");
    file_super_read_nb_mis = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_super_read_nb_mis,"extract_reads()\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_pos_mis");
    file_super_read_pos_mis = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_super_read_pos_mis,"extract_reads()\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_id_mis");
    file_super_read_id_mis = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_super_read_id_mis,"extract_reads()\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_char_mis");
    file_super_read_char_mis = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_super_read_char_mis,"extract_reads()\n")

    if (paired_end){

        strcpy(&(output_prefix_buffer[len_output_prefix]), "_reads_pos");
        file_reads_pos = fopen(output_prefix_buffer, "w");
        ASSERT_NULL_PTR(file_reads_pos,"extract_reads()\n")
    }

    while (getline(&buffer_read, &size_buffer_read, file_super_reads) != -1){

        first_read = true;
        pos_read = 0;

        while (true){

            if (fread(&read_length, sizeof(uint32_t), 1, file_super_read_length) != 1) break;

            if (!first_read && (read_length & 0x1)) break;
            else {

                if (fread(&pos_read_tmp, sizeof(uint32_t), 1, file_super_read_pos) != 1)
                    ERROR("extract_reads(): truncated spanning super read position file.\n")

                if (fread(&occ_count, sizeof(uint32_t), 1, file_super_read_occ) != 1)
                    ERROR("extract_reads(): truncated spanning super read occurence file.\n")

                if (fread(&nb_mismatches, sizeof(uint32_t), 1, file_super_read_nb_mis) != 1)
                    ERROR("extract_reads(): truncated spanning super read nb mismatches file.\n")

                if ((nb_reads % (SIZE_BUFFER * SIZE_BITS_UINT_8T)) == 0){
                    nb_reads = 0;
                    fread(buffer_rev, sizeof(uint8_t), SIZE_BUFFER, file_super_read_rev);
                }

                first_read = false;
                position_mis = 0;
                read_length >>= 1;
                pos_read += pos_read_tmp;

                if ((read_length + 1) > size_read){

                    size_read = read_length + 1;

                    read = realloc(read, size_read * sizeof(char));
                    ASSERT_NULL_PTR(read, "extract_reads()\n")

                    read_tmp = realloc(read_tmp, size_read * sizeof(char));
                    ASSERT_NULL_PTR(read_tmp, "extract_reads()\n")

                    read_rev = realloc(read_rev, size_read * sizeof(char));
                    ASSERT_NULL_PTR(read_rev, "extract_reads()\n")
                }

                if (nb_mismatches){

                    mis = malloc(nb_mismatches * sizeof(Mismatch));
                    ASSERT_NULL_PTR(mis, "extract_reads()\n")

                    for (i = 0; i < nb_mismatches; i++){

                        fread(&position_mis_tmp, sizeof(uint32_t), 1, file_super_read_pos_mis);
                        fread(&id_mis, sizeof(uint32_t), 1, file_super_read_id_mis);
                        fread(&(mis[i].mismatch_char), sizeof(char), 1, file_super_read_char_mis);

                        position_mis += position_mis_tmp;

                        mis[i].position = (int) position_mis;
                        mis[i].id_read = (int) id_mis;
                    }

                    qsort(mis, nb_mismatches, sizeof(Mismatch), mismatch_pos_cmp);
                }
                else mis = NULL;

                pos_length_occ = create_set_Pos_length_occ(pos_read, read_length, occ_count,
                                                           (buffer_rev[nb_reads/SIZE_BITS_UINT_8T] >> (SIZE_BITS_UINT_8T - 1 - (nb_reads%SIZE_BITS_UINT_8T))) & 0x1,
                                                           nb_mismatches, NULL, mis);

                List_push(l, pos_length_occ);

                nb_reads++;
            }
        }

        fseek(file_super_read_length, 0 - ((int) sizeof(uint32_t)), SEEK_CUR);

        while ((pos_length_occ = List_pop_first(l)) != NULL){

            memcpy(read, &buffer_read[pos_length_occ->position], pos_length_occ->length * sizeof(char));
            read[pos_length_occ->length] = eol;

            it_mis = 0;

            while ((it_mis < pos_length_occ->nb_mismatches) && (pos_length_occ->list_mismatches[it_mis].id_read < pos_length_occ->occ_count))
                it_mis++;

            while ((it_mis < pos_length_occ->nb_mismatches) && (pos_length_occ->list_mismatches[it_mis].id_read >= pos_length_occ->occ_count)){
                read[pos_length_occ->list_mismatches[it_mis].position] = pos_length_occ->list_mismatches[it_mis].mismatch_char;
                it_mis++;
            }

            for (i = 0, it_mis = 0; i < pos_length_occ->occ_count; i++){

                strcpy(read_tmp, read);

                while ((it_mis < pos_length_occ->nb_mismatches) && (pos_length_occ->list_mismatches[it_mis].id_read < i)) it_mis++;

                while ((it_mis < pos_length_occ->nb_mismatches) && (pos_length_occ->list_mismatches[it_mis].id_read == i)){
                    read_tmp[pos_length_occ->list_mismatches[it_mis].position] = pos_length_occ->list_mismatches[it_mis].mismatch_char;
                    it_mis++;
                }

                if (paired_end){
                    read_pos = ftell(file_reads);
                    fwrite(&read_pos, sizeof(long int), 1, file_reads_pos);
                }

                if (pos_length_occ->rev){
                    reverse_complement(read_tmp, read_rev, pos_length_occ->length);
                    fwrite(read_rev, sizeof(char), pos_length_occ->length, file_reads);
                }
                else fwrite(read_tmp, sizeof(char), pos_length_occ->length, file_reads);

                fwrite(&nl, sizeof(char), 1, file_reads);
            }

            free(pos_length_occ->list_mismatches);
            free(pos_length_occ);
        }
    }

    fseek(file_reads, 0 - ((int) sizeof(char)), SEEK_CUR);
    fwrite(&eol, sizeof(char), 1, file_reads);

    fclose(file_super_read_pos_mis);
    fclose(file_super_read_id_mis);
    fclose(file_super_read_char_mis);
    fclose(file_super_read_nb_mis);
    fclose(file_super_read_length);
    fclose(file_super_read_occ);
    fclose(file_super_read_pos);
    fclose(file_super_read_rev);
    fclose(file_super_reads);
    fclose(file_reads);

    if (paired_end) fclose(file_reads_pos);

    List_destroy(l);

    free(read);
    free(read_tmp);
    free(read_rev);
    free(buffer_rev);
    free(buffer_read);
    free(output_prefix_buffer);
}

void reorder_paired_reads(char* output_prefix){

    FILE* file_reads;
    FILE* file_reads_mate1;
    FILE* file_reads_mate2;
    FILE* file_reads_pos;
    FILE* file_mate_pos;
    FILE* file_mate_info;

    int len_output_prefix = strlen(output_prefix);

    long int mate2_pos = 0;
    long int mate1_pos = 0;

    int64_t mate_pos;

    int64_t nb_reads = 0;
    int64_t nb_reads_buff = 0;

    size_t size_buffer_read = SIZE_BUFFER;

    char* buffer_read = malloc(size_buffer_read * sizeof(char));
    ASSERT_NULL_PTR(buffer_read, "reorder_paired_reads() 1\n")

    char* output_prefix_buffer = malloc((len_output_prefix + 50) * sizeof(char));
    ASSERT_NULL_PTR(output_prefix_buffer, "reorder_paired_reads() 2\n")

    uint8_t* buffer_mate_info = malloc(SIZE_BUFFER * sizeof(uint8_t));
    ASSERT_NULL_PTR(buffer_mate_info, "reorder_paired_reads() 3\n")

    strcpy(output_prefix_buffer, output_prefix);

    strcpy(&(output_prefix_buffer[len_output_prefix]), ".reads");
    file_reads = fopen(output_prefix_buffer, "r");
    ASSERT_NULL_PTR(file_reads,"reorder_paired_reads() 4\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), ".reads_1");
    file_reads_mate1 = fopen(output_prefix_buffer, "w");
    ASSERT_NULL_PTR(file_reads_mate1,"reorder_paired_reads() 5\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), ".reads_2");
    file_reads_mate2 = fopen(output_prefix_buffer, "w");
    ASSERT_NULL_PTR(file_reads_mate2,"reorder_paired_reads() 5\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_reads_pos");
    file_reads_pos = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_reads_pos,"reorder_paired_reads() 6\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_mate");
    file_mate_pos = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_mate_pos,"reorder_paired_reads() 7\n")

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_span_sup_reads_mate_info");
    file_mate_info = fopen(output_prefix_buffer, "rb");
    ASSERT_NULL_PTR(file_mate_info,"reorder_paired_reads() 8\n")

    while (getline(&buffer_read, &size_buffer_read, file_reads) != -1){

        if ((nb_reads % (SIZE_BUFFER * SIZE_BITS_UINT_8T)) == 0){

            fread(buffer_mate_info, sizeof(uint8_t), SIZE_BUFFER, file_mate_info);
            nb_reads_buff = 0;
        }

        if (buffer_mate_info[nb_reads_buff/SIZE_BITS_UINT_8T] & MASK_POWER_8[nb_reads_buff%SIZE_BITS_UINT_8T]){ //If current read is left mate

            fwrite(buffer_read, sizeof(char), strlen(buffer_read), file_reads_mate1); //Write left mate
            mate1_pos = ftell(file_reads); //Save current position in read file

            fread(&mate_pos, sizeof(int64_t), 1, file_mate_pos); //Read location right mate in read file

            mate_pos += nb_reads;

            fseek(file_reads_pos, mate_pos * sizeof(long int), SEEK_SET); //Read position right mate in read file
            fread(&mate2_pos, sizeof(long int), 1, file_reads_pos);

            fseek(file_reads, mate2_pos * sizeof(char), SEEK_SET); //Move to position right mate in read file

            getline(&buffer_read, &size_buffer_read, file_reads); //Read right mate

            fwrite(buffer_read, sizeof(char), strlen(buffer_read), file_reads_mate2); //Write right mate

            fseek(file_reads, mate1_pos * sizeof(char), SEEK_SET); //Move back to the original position in read file
        }

        nb_reads++;
        nb_reads_buff++;
    }

    fclose(file_reads);
    fclose(file_reads_mate1);
    fclose(file_reads_mate2);
    fclose(file_reads_pos);
    fclose(file_mate_pos);
    fclose(file_mate_info);

    free(buffer_read);
    free(buffer_mate_info);
    free(output_prefix_buffer);

    strcpy(&(output_prefix_buffer[len_output_prefix]), ".reads");
    if (remove(output_prefix_buffer)) printf("Warning: Could not remove temporary file.\n");

    strcpy(&(output_prefix_buffer[len_output_prefix]), "_reads_pos");
    if (remove(output_prefix_buffer)) printf("Warning: Could not remove temporary file.\n");

    return;
}
