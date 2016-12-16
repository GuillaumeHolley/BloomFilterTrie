#include "bft.h"

extern inline uint8_t intersection_annots(const uint8_t a, const uint8_t b);
extern inline uint8_t union_annots(const uint8_t a, const uint8_t b);
extern inline uint8_t sym_difference_annots(const uint8_t a, const uint8_t b);

/** Function creating a colored de Bruijn graph stored in a BFT.
* @param k is the length of k-mers.
* @param treshold_compression indicates when the color compression should be triggered (every treshold_compression genome inserted).
* @return a BFT pointer.
*/
BFT* create_cdbg(int k, int treshold_compression){
    return createBFT_Root(k, treshold_compression, 0);
}

/** Free an allocated colored de Bruijn graph stored in a BFT.
* @param bft is an allocated BFT.
*/
void free_cdbg(BFT* bft){
    freeBFT_Root(bft);
    return;
}

/** Function inserting genomes (k-mer file) in a BFT.
* @param nb_files is the number of files to insert.
* @param paths is an nb_files size array of strings (char*). Each string is the name of a file (+ eventually its path) to insert.
* @param bft is a BFT where the genomes are inserted.
* @param prefix_bft_filename is a prefix filename (including path) where temporary data can be written to. The prefix must be
            unique in its directory.
*/
void insert_genomes_from_files(int nb_files, char** paths, BFT* bft, char* prefix_bft_filename){

    if (bft->marked == 0) insert_Genomes_from_KmerFiles(bft, nb_files, paths, 0, prefix_bft_filename);
    else ERROR("insert_genomes_from_files(): graph is locked for marking or/and neighbors traversal. Unlock it before insertion.\n");

    return;
}

/** Function inserting k-mers of a new genome in a BFT.
* @param nb_kmers is the number of k-mers to insert.
* @param kmers is a pointer to an array of strings (char *) that are the k-mers to insert. The array is of length nb_kmers.
* @param genome_name is the name of the new genome to which the inserted k-mers come from.
* @param bft is a colored de Bruijn graph stored as a BFT.
*/
void insert_kmers_new_genome(int nb_kmers, char** kmers, char* genome_name, BFT* bft){

    int size_id_genome;
    int nb_bytes_kmer_comp;

    uint8_t* kmer_comp;

    if (bft->marked == 0){

        add_genomes_BFT_Root(1, &genome_name, bft);

        if (nb_kmers > 0){

            nb_bytes_kmer_comp = CEIL(bft->k * 2, SIZE_BITS_UINT_8T);

            kmer_comp = malloc(nb_bytes_kmer_comp * sizeof(uint8_t));
            ASSERT_NULL_PTR(kmer_comp,"insert_kmers_new_genome()\n")

            size_id_genome = get_nb_bytes_power2_annot(bft->nb_genomes-1);

            for (int i = 0; i < nb_kmers; i++){

                memset(kmer_comp, 0, nb_bytes_kmer_comp * sizeof(uint8_t));

                if (parseKmerCount(kmers[i], bft->k, kmer_comp, 0)) insertKmers(bft, kmer_comp, 1, bft->nb_genomes-1, size_id_genome);
                else ERROR("insert_kmers_new_genome(): could not insert k-mer in graph, it probably contains unvalid characters.\n")
            }

            free(kmer_comp);
        }
    }
    else ERROR("insert_kmers_new_genome(): graph is locked for marking or/and neighbors traversal. Unlock it before insertion.\n");

    return;
}

/** Function inserting k-mers of the last inserted genome in a BFT.
* @param nb_kmers is the number of k-mers to insert.
* @param kmers is a pointer to an array of strings (char *) that are the k-mers to insert. The arrayis of length nb_kmers.
* @param bft is a colored de Bruijn graph stored as a BFT.
*/
void insert_kmers_last_genome(int nb_kmers, char** kmers, BFT* bft){

    int size_id_genome;
    int nb_bytes_kmer_comp;

    uint8_t* kmer_comp;

    if (bft->marked == 0){

        if (bft->nb_genomes > 0){

            if (nb_kmers > 0){

                nb_bytes_kmer_comp = CEIL(bft->k * 2, SIZE_BITS_UINT_8T);

                kmer_comp = malloc(nb_bytes_kmer_comp * sizeof(uint8_t));
                ASSERT_NULL_PTR(kmer_comp,"insert_kmers_last_genome()\n")

                size_id_genome = get_nb_bytes_power2_annot(bft->nb_genomes-1);

                for (int i = 0; i < nb_kmers; i++){

                    memset(kmer_comp, 0, nb_bytes_kmer_comp * sizeof(uint8_t));

                    if (parseKmerCount(kmers[i], bft->k, kmer_comp, 0)) insertKmers(bft, kmer_comp, 1, bft->nb_genomes-1, size_id_genome);
                    else ERROR("insert_kmers_last_genome(): could not insert k-mer in graph, it probably contains unvalid characters.\n")
                }

                free(kmer_comp);
            }
        }
        else ERROR("insert_kmers_last_genome(): the graph is empty, there is no last genome.\n");
    }
    else ERROR("insert_kmers_last_genome(): graph is locked for marking or/and neighbors traversal. Unlock it before insertion.\n");

    return;
}

/** Function creating a BFT_kmer object from a k-mer encoded as an ASCII string (char*).
* @param kmer is an an ASCII encoded k-mer string (char*).
* @param k is the k-mer length.
* @return a BFT_kmer pointer.
*/
BFT_kmer* create_kmer(const char* kmer, int k){

    ASSERT_NULL_PTR(kmer, "create_kmer()\n")

    if (strlen(kmer) != k) ERROR("create_kmer(): k-mer length is not the one used in the graph.\n")

    BFT_kmer* bft_kmer = malloc(sizeof(BFT_kmer));
    ASSERT_NULL_PTR(bft_kmer, "create_kmer()\n")

    bft_kmer->kmer = malloc((k + 1) * sizeof(char));
    ASSERT_NULL_PTR(bft_kmer->kmer, "create_kmer()\n")

    bft_kmer->kmer_comp = calloc(CEIL(k * 2, SIZE_BITS_UINT_8T), sizeof(uint8_t));
    ASSERT_NULL_PTR(bft_kmer->kmer_comp, "create_kmer()\n")

    memcpy(bft_kmer->kmer, kmer, k * sizeof(char));
    bft_kmer->kmer[k] = '\0';

    if (parseKmerCount(bft_kmer->kmer, k, bft_kmer->kmer_comp, 0)){

        bft_kmer->res = malloc(sizeof(resultPresence));
        ASSERT_NULL_PTR(bft_kmer->res, "create_kmer()\n");

        initialize_resultPresence(bft_kmer->res);

        return bft_kmer;
    }

    ERROR("create_kmer(): Unexpected character encountered in k-mer.\n")
}

/** Function creating an empty BFT_kmer object (all its components are NULL).
* @return a BFT_kmer pointer.
*/
BFT_kmer* create_empty_kmer(){

    BFT_kmer* bft_kmer = malloc(sizeof(BFT_kmer));
    ASSERT_NULL_PTR(bft_kmer, "create_empty_kmer()\n")

    bft_kmer->kmer = NULL;
    bft_kmer->kmer_comp = NULL;
    bft_kmer->res = NULL;

    return bft_kmer;
}

/** Function freeing allocated BFT_kmers.
* @param bft_kmer is a pointer to an array of at least one BFT_kmer.
* @param nb_bft_kmer is the number of BFT_kmer in bft_kmer.
*/
void free_BFT_kmer(BFT_kmer* bft_kmer, int nb_bft_kmer){

    ASSERT_NULL_PTR(bft_kmer, "free_BFT_kmer()\n")

    for (int i = 0; i < nb_bft_kmer; i++){
        if (bft_kmer[i].kmer != NULL) free(bft_kmer[i].kmer);
        if (bft_kmer[i].kmer_comp != NULL) free(bft_kmer[i].kmer_comp);
        if (bft_kmer[i].res != NULL) free(bft_kmer[i].res);
    }

    free(bft_kmer);

    return;
}

/** Function freeing the content of allocated BFT_kmers.
* @param bft_kmer is a pointer to an array of at least one BFT_kmer.
* @param nb_bft_kmer is the number of BFT_kmer in bft_kmer.
*/
void free_BFT_kmer_content(BFT_kmer* bft_kmer, int nb_bft_kmer){

    ASSERT_NULL_PTR(bft_kmer, "free_BFT_kmer_content()\n")

    for (int i = 0; i < nb_bft_kmer; i++){
        if (bft_kmer[i].kmer != NULL) free(bft_kmer[i].kmer);
        if (bft_kmer[i].kmer_comp != NULL) free(bft_kmer[i].kmer_comp);
        if (bft_kmer[i].res != NULL) free(bft_kmer[i].res);
    }

    return;
}

/** Function searching for a k-mer in a BFT.
* @param kmer is an an ASCII encoded k-mer string (char*) to search for in the BFT.
* @param bft is a BFT in which k-mer is searched
* @return a BFT_kmer pointer.
*/
BFT_kmer* get_kmer(const char* kmer, BFT* bft){

    ASSERT_NULL_PTR(kmer, "get_kmer()\n")
    ASSERT_NULL_PTR(bft, "get_kmer()\n")

    BFT_kmer* bft_kmer = malloc(sizeof(BFT_kmer));
    ASSERT_NULL_PTR(bft_kmer, "isKmer_in_cdbg()\n")

    bft_kmer->kmer = malloc((bft->k + 1) * sizeof(char));
    ASSERT_NULL_PTR(bft_kmer->kmer, "get_kmer()\n")

    bft_kmer->kmer_comp = calloc(CEIL(bft->k * 2, SIZE_BITS_UINT_8T), sizeof(uint8_t));
    ASSERT_NULL_PTR(bft_kmer->kmer_comp, "get_kmer()\n")

    memcpy(bft_kmer->kmer, kmer, bft->k * sizeof(char));
    bft_kmer->kmer[bft->k] = '\0';

    if (parseKmerCount(bft_kmer->kmer, bft->k, bft_kmer->kmer_comp, 0)){

        bft_kmer->res = isKmerPresent(&(bft->node), bft, (bft->k / NB_CHAR_SUF_PREF) - 1, bft_kmer->kmer_comp, bft->k);
        return bft_kmer;
    }

    ERROR("get_kmer(): Unexpected character encountered in k-mer.\n")
}

/** Function testing if a k-mer is in a BFT.
* @param bft_kmer is a k-mer obtained via search or iteration over a BFT (via get_kmer() for example).
* @return a boolean indicating the presence (true) or absence (false) of the k-mer in a BFT.
*/
bool is_kmer_in_cdbg(BFT_kmer* bft_kmer){
    return (bft_kmer->res->link_child == NULL ? false : true);
}

/** Function extracting the k-mers of a BFT in a file.
* @param bft is a BFT containing the k-mers to iterate over.
* @param filename_output is the name of a file to which the k-mers are written. File is overwritten if it already exists.
* @param compressed_output is a boolean indicating if the k-mers should be written in their 2 bits form (true) or ASCII form (false).
*/
void extract_kmers_to_disk(BFT* bft, char* filename_output, bool compressed_output){

    ASSERT_NULL_PTR(bft,"extract_kmers_to_disk()\n")
    ASSERT_NULL_PTR(filename_output,"extract_kmers_to_disk()\n")

    char int_to_string[20];

    int length_string;
    int lvl_root = (bft->k / NB_CHAR_SUF_PREF) - 1;
    int nb_bytes_kmer_comp = CEIL(bft->k * 2, SIZE_BITS_UINT_8T);

    FILE* file_extract_kmers;

    if ((file_extract_kmers = fopen(filename_output, "w")) == NULL)
        ERROR("extract_kmers_to_disk(): failed to create/open output file.\n")

    if (compressed_output){

        length_string = sprintf(int_to_string, "%d", bft->k);
        int_to_string[length_string] = '\n';
        fwrite(int_to_string, sizeof(char), length_string+1, file_extract_kmers);

        memory_Used* mem = printMemoryUsedFromNode(&(bft->node), lvl_root, bft->k, bft->info_per_lvl);

        length_string = sprintf(int_to_string, "%d", (int)mem->nb_kmers_in_UCptr);
        int_to_string[length_string] = '\n';
        fwrite(int_to_string, sizeof(char), length_string+1, file_extract_kmers);

        free(mem);

        iterate_over_kmers(bft, write_kmer_comp_to_disk, nb_bytes_kmer_comp, file_extract_kmers);
    }
    else iterate_over_kmers(bft, write_kmer_ascii_to_disk, file_extract_kmers);

    fclose(file_extract_kmers);
}

/** Function writing an ASCII k-mer in a file.
*   This function is of type BFT_func_ptr and is intended to be a parameter of iterate_over_kmers() or v_iterate_over_kmers().
* @param bft_kmer is a k-mer to write to disk.
* @param bft is a BFT from which bft_kmer was extracted.
* @param args is a variable list of arguments. It contains a pointer to a file where to write bft_kmer.
*/
size_t write_kmer_ascii_to_disk(BFT_kmer* bft_kmer, BFT* bft, va_list args){

    FILE* file = va_arg(args, FILE*);

    bft_kmer->kmer[bft->k] = '\n';

    fwrite(bft_kmer->kmer, sizeof(char), bft->k+1, file);

    return 1;
}

/** Function writing an 2 bits encoded k-mer in a file.
*   This function is of type BFT_func_ptr and is intended to be a parameter of iterate_over_kmers() or v_iterate_over_kmers().
* @param bft_kmer is a k-mer to write to disk.
* @param bft is a BFT from which bft_kmer was extracted.
* @param args is a variable list of arguments. It contains a pointer to a file where to write bft_kmer.
*/
size_t write_kmer_comp_to_disk(BFT_kmer* bft_kmer, BFT* bft, va_list args){

    int nb_bytes_kmer_comp = va_arg(args, int);
    FILE* file = va_arg(args, FILE*);

    fwrite(bft_kmer->kmer_comp, sizeof(uint8_t), nb_bytes_kmer_comp, file);

    return 1;
}

 /** Function creating an empty BFT_annotation. */
inline BFT_annotation* create_BFT_annotation(){

    BFT_annotation* bft_annot = malloc(sizeof(BFT_annotation));
    ASSERT_NULL_PTR(bft_annot, "create_BFT_annotation()\n")

    bft_annot->annot = NULL;
    bft_annot->annot_ext = NULL;
    bft_annot->annot_cplx = NULL;

    bft_annot->size_annot = -1;
    bft_annot->size_annot_cplx = -1;

    bft_annot->from_BFT = 0;

    return bft_annot;
}

/** Function freeing a BFT_annotation.
* @param bft_annot is a pointer to the BFT_annotation to free.
*/
void free_BFT_annotation(BFT_annotation* bft_annot){

    ASSERT_NULL_PTR(bft_annot, "free_BFT_annotation()\n")

    if (!bft_annot->from_BFT){
        if (bft_annot->annot != NULL) free(bft_annot->annot);
        if (bft_annot->annot_ext != NULL) free(bft_annot->annot_ext);
        if (bft_annot->annot_cplx != NULL) free(bft_annot->annot_cplx);
    }

    free(bft_annot);
}

/** Function extracting the annotation (set of colors) associated with a k-mer of a BFT.
* @param bft_kmer is a k-mer obtained via search or iteration over a BFT (via get_kmer() for example).
* @return a BFT_annotation pointer.
*/
BFT_annotation* get_annotation(BFT_kmer* bft_kmer){

    ASSERT_NULL_PTR(bft_kmer, "get_annotation()\n")

    if (is_kmer_in_cdbg(bft_kmer)){

        BFT_annotation* bft_annot = create_BFT_annotation();

        if (bft_kmer->res->posFilter2 != 0){
            get_annot((UC*)bft_kmer->res->container, &(bft_annot->annot), &(bft_annot->annot_ext), &(bft_annot->annot_cplx),
                      &(bft_annot->size_annot), &(bft_annot->size_annot_cplx), bft_kmer->res->posFilter2, bft_kmer->res->posFilter3,
                      bft_kmer->res->pos_sub_bucket);
        }
        else{
            get_annot(&(((UC*)((CC*)bft_kmer->res->container)->children)[bft_kmer->res->bucket]), &(bft_annot->annot),
                      &(bft_annot->annot_ext), &(bft_annot->annot_cplx), &(bft_annot->size_annot), &(bft_annot->size_annot_cplx),
                      bft_kmer->res->posFilter2, bft_kmer->res->posFilter3, bft_kmer->res->pos_sub_bucket);
        }

        bft_annot->from_BFT = 1;

        return bft_annot;
    }
    else ERROR("get_annotation(): k-mer is not present in the graph.\n")
}

/** Function testing if a k-mer occured in a genome.
* @param id_genome is the genome identifier.
* @param bft_annot is the annotation of the k-mer to test the presence in genome.
* @param bft is a BFT in which the k-mer is is stored.
* @return a boolean indicating the presence (true) or absence (false) of the k-mer in a the genome.
*/
bool presence_genome(uint32_t id_genome, BFT_annotation* bft_annot, BFT* bft){

    ASSERT_NULL_PTR(bft_annot, "is_genome_present()\n")
    ASSERT_NULL_PTR(bft, "is_genome_present()\n")

    if (id_genome < bft->nb_genomes){

        if (id_genome < bft->nb_genomes/2){

            return is_genome_present(bft->ann_inf, bft->comp_set_colors, bft_annot->annot, bft_annot->size_annot,
                                     bft_annot->annot_ext, bft_annot->annot_ext != NULL, id_genome);
        }

        return is_genome_present_from_end_annot(bft->ann_inf, bft->comp_set_colors, bft_annot->annot, bft_annot->size_annot,
                                                bft_annot->annot_ext, bft_annot->annot_ext != NULL, id_genome);
    }

    return false;
}

/** Function computing the intersection of a set of annotations.
* @param bft is a BFT from which the input annotations are originated.
* @param nb_annotations indicates how many annotations must be included in the intersection.
* @param ... is a list of nb_annotations BFT_annotation pointers of which the intersection is computed.
* @return a BFT_annotation pointer to an annotation which is the intersection of the input annotations.
*/
BFT_annotation* intersection_annotations(BFT* bft, uint32_t nb_annotations, ... ){

    va_list args;
    BFT_annotation* bft_annot;
    annotation_array_elem* ann_arr_elem1;
    annotation_array_elem* ann_arr_elem2;

    va_start(args, nb_annotations);

    if (nb_annotations == 0) ERROR("intersection_annotations(): no annotations given as parameters.\n")
    else if (nb_annotations == 1){

        ann_arr_elem1 = cmp_annots(NULL, 0, NULL, 0, NULL, 0, NULL, 0, bft->nb_genomes-1,
                                  intersection_annots, bft->comp_set_colors);
    }
    else{

        BFT_annotation* bft_annot_tmp;

        bft_annot = va_arg(args, BFT_annotation*);
        ASSERT_NULL_PTR(bft_annot, "intersection_annotations()\n");

        bft_annot_tmp = va_arg(args, BFT_annotation*);
        ASSERT_NULL_PTR(bft_annot_tmp, "intersection_annotations()\n");

        int size_annot_ext = bft_annot->annot_ext != NULL;
        int size_annot_ext_tmp = bft_annot_tmp->annot_ext != NULL;

        ann_arr_elem1 = cmp_annots(bft_annot->annot, bft_annot->size_annot, bft_annot->annot_ext, size_annot_ext,
                                   bft_annot_tmp->annot, bft_annot_tmp->size_annot, bft_annot_tmp->annot_ext, size_annot_ext_tmp,
                                   bft->nb_genomes-1, intersection_annots, bft->comp_set_colors);

        for (uint32_t i = 2; i < nb_annotations; i++){

            bft_annot_tmp = va_arg(args, BFT_annotation*);
            ASSERT_NULL_PTR(bft_annot_tmp, "intersection_annotations()\n");

            size_annot_ext_tmp = bft_annot_tmp->annot_ext != NULL;

            ann_arr_elem2 = cmp_annots(ann_arr_elem1->annot_array, ann_arr_elem1->size_annot, NULL, 0,
                                       bft_annot_tmp->annot, bft_annot_tmp->size_annot, bft_annot_tmp->annot_ext, size_annot_ext_tmp,
                                       bft->nb_genomes-1, intersection_annots, bft->comp_set_colors);

            free(ann_arr_elem1->annot_array);
            free(ann_arr_elem1);

            ann_arr_elem1 = ann_arr_elem2;
        }
    }

    bft_annot = create_BFT_annotation();
    bft_annot->annot = ann_arr_elem1->annot_array;
    bft_annot->size_annot = ann_arr_elem1->size_annot;

    free(ann_arr_elem1);

    va_end(args);

    return bft_annot;
}

/** Function computing the union of a set of annotations.
* @param bft is a BFT from which the input annotations are originated.
* @param nb_annotations indicates how many annotations must be included in the union.
* @param ... is a list of nb_annotations BFT_annotation pointers of which the union is computed.
* @return a BFT_annotation pointer to an annotation which is the union of the input annotations.
*/
BFT_annotation* union_annotations(BFT* bft, uint32_t nb_annotations, ... ){

    va_list args;

    BFT_annotation* bft_annot;

    annotation_array_elem* ann_arr_elem1;
    annotation_array_elem* ann_arr_elem2;

    int size_annot_ext;

    va_start(args, nb_annotations);

    if (nb_annotations == 0) ERROR("union_annotations(): no annotations given as parameters.\n")
    else if (nb_annotations == 1){

        bft_annot = va_arg(args, BFT_annotation*);
        ASSERT_NULL_PTR(bft_annot, "union_annotations()\n");

        size_annot_ext = bft_annot->annot_ext != NULL;

        ann_arr_elem1 = cmp_annots(bft_annot->annot, bft_annot->size_annot, bft_annot->annot_ext, size_annot_ext,
                                   NULL, 0, NULL, 0, bft->nb_genomes-1, union_annots, bft->comp_set_colors);
    }
    else{

        BFT_annotation* bft_annot_tmp;

        bft_annot = va_arg(args, BFT_annotation*);
        ASSERT_NULL_PTR(bft_annot, "union_annotations()\n");

        bft_annot_tmp = va_arg(args, BFT_annotation*);
        ASSERT_NULL_PTR(bft_annot_tmp, "union_annotations()\n");

        size_annot_ext = bft_annot->annot_ext != NULL;

        int size_annot_ext_tmp = bft_annot_tmp->annot_ext != NULL;

        ann_arr_elem1 = cmp_annots(bft_annot->annot, bft_annot->size_annot, bft_annot->annot_ext, size_annot_ext,
                                   bft_annot_tmp->annot, bft_annot_tmp->size_annot, bft_annot_tmp->annot_ext, size_annot_ext_tmp,
                                   bft->nb_genomes-1, union_annots, bft->comp_set_colors);

        for (uint32_t i = 2; i < nb_annotations; i++){

            bft_annot_tmp = va_arg(args, BFT_annotation*);
            ASSERT_NULL_PTR(bft_annot_tmp, "union_annotations()\n");

            size_annot_ext_tmp = bft_annot_tmp->annot_ext != NULL;

            ann_arr_elem2 = cmp_annots(ann_arr_elem1->annot_array, ann_arr_elem1->size_annot, NULL, 0,
                                       bft_annot_tmp->annot, bft_annot_tmp->size_annot, bft_annot_tmp->annot_ext, size_annot_ext_tmp,
                                       bft->nb_genomes-1, union_annots, bft->comp_set_colors);

            free(ann_arr_elem1->annot_array);
            free(ann_arr_elem1);

            ann_arr_elem1 = ann_arr_elem2;
        }
    }

    bft_annot = create_BFT_annotation();
    bft_annot->annot = ann_arr_elem1->annot_array;
    bft_annot->size_annot = ann_arr_elem1->size_annot;

    free(ann_arr_elem1);

    va_end(args);

    return bft_annot;
}

/** Function computing the symmetric difference of a set of annotations.
* @param bft is a BFT from which the input annotations are originated.
* @param nb_annotations indicates how many annotations must be included in the symmetric difference.
* @param ... is a list of nb_annotations BFT_annotation pointers of which the symmetric difference is computed.
* @return a BFT_annotation pointer to an annotation which is the symmetric difference of the input annotations.
*/
BFT_annotation* sym_difference_annotations(BFT* bft, uint32_t nb_annotations, ... ){

    va_list args;

    BFT_annotation* bft_annot;

    annotation_array_elem* ann_arr_elem1;

    int size_annot_ext;

    va_start(args, nb_annotations);

    if (nb_annotations == 0) ERROR("sym_difference_annotations(): no annotations given as parameters.\n")
    else if (nb_annotations == 1){

        bft_annot = va_arg(args, BFT_annotation*);
        ASSERT_NULL_PTR(bft_annot, "sym_difference_annotations()\n");

        size_annot_ext = bft_annot->annot_ext != NULL;

        ann_arr_elem1 = cmp_annots(bft_annot->annot, bft_annot->size_annot, bft_annot->annot_ext, size_annot_ext,
                                   NULL, 0, NULL, 0, bft->nb_genomes-1, sym_difference_annots, bft->comp_set_colors);
    }
    else{

        BFT_annotation* bft_annot_tmp;

        bft_annot = intersection_annotations(bft, nb_annotations, args);
        bft_annot_tmp = union_annotations(bft, nb_annotations, args);

        size_annot_ext = bft_annot->annot_ext != NULL;

        int size_annot_ext_tmp = bft_annot_tmp->annot_ext != NULL;

        ann_arr_elem1 = cmp_annots(bft_annot->annot, bft_annot->size_annot, bft_annot->annot_ext, size_annot_ext,
                                   bft_annot_tmp->annot, bft_annot_tmp->size_annot, bft_annot_tmp->annot_ext, size_annot_ext_tmp,
                                   bft->nb_genomes-1, sym_difference_annots, bft->comp_set_colors);
    }

    bft_annot = create_BFT_annotation();
    bft_annot->annot = ann_arr_elem1->annot_array;
    bft_annot->size_annot = ann_arr_elem1->size_annot;

    free(ann_arr_elem1);

    va_end(args);

    return bft_annot;
}

/** Function extracting a list of genome identifiers from an annotation.
* @param bft_annot is an annotation from which the ids must be extracted.
* @param bft is a BFT from which the annotation was extracted.
* @return a pointer to an array of genome identifiers (uint32_t).
    The first element of this array (position 0) indicates how many ids there are in this array.
    Therefore, the length of the array is array[0] + 1.
*/
uint32_t* get_list_id_genomes(BFT_annotation* bft_annot, BFT* bft){

    ASSERT_NULL_PTR(bft_annot, "get_list_id_genomes()\n")
    ASSERT_NULL_PTR(bft, "get_list_id_genomes()\n")

    int size_annot_ext = bft_annot->annot_ext != NULL;

    get_id_genomes_from_annot(bft->ann_inf, bft->comp_set_colors, bft_annot->annot, bft_annot->size_annot,
                              bft_annot->annot_ext, size_annot_ext);

    uint32_t* ids = malloc((bft->ann_inf->nb_id_stored + 1) * sizeof(uint32_t));
    ASSERT_NULL_PTR(ids, "get_list_id_genomes()\n")

    ids[0] = bft->ann_inf->nb_id_stored;
    memcpy(&(ids[1]), bft->ann_inf->id_stored, ids[0] * sizeof(uint32_t));

    reinit_annotation_inform(bft->ann_inf);

    return ids;
}

/** Function counting the number of genome identifiers in an annotation.
* @param bft_annot is an annotation.
* @param bft is a BFT from which the annotation was extracted.
* @return a count of genome identifiers.
*/
uint32_t get_count_id_genomes(BFT_annotation* bft_annot, BFT* bft){

    ASSERT_NULL_PTR(bft_annot, "get_count_id_genomes()\n")
    ASSERT_NULL_PTR(bft, "get_count_id_genomes()\n")

    int size_annot_ext = bft_annot->annot_ext != NULL;

    return get_count_id_genomes_from_annot(bft->ann_inf, bft->comp_set_colors, bft_annot->annot,
                                           bft_annot->size_annot, bft_annot->annot_ext, size_annot_ext);
}

uint32_t* intersection_list_id_genomes(uint32_t* list_a, uint32_t* list_b){

    ASSERT_NULL_PTR(list_a, "intersection_list_id_genomes()\n")
    ASSERT_NULL_PTR(list_b, "intersection_list_id_genomes()\n")

    int i = 1, j = 1;

    uint32_t size_a = list_a[0]+1, size_b = list_b[0]+1;

    uint32_t* list_c = malloc(MIN(size_a, size_b) * sizeof(uint32_t));
    ASSERT_NULL_PTR(list_c, "intersection_list_id_genomes()\n")

    list_c[0] = 0;

    while (i < size_a && j < size_b){
        if (list_a[i] > list_b[j]) j++;
        else if (list_b[j] > list_a[i]) i++;
        else {
            list_c[0]++;
            list_c[list_c[0]] = list_a[i];
            i++; j++;
        }
    }

    return list_c;
}

/** Function locking and preparing the graph for vertices marking (no insertion can happen before unlocking).
* By default, all k-mers of the graph are initialized with a 0 flag value.
* @param bft is a BFT to lock and prepare for vertices marking.
*/
void set_marking(BFT* bft){

    ASSERT_NULL_PTR(bft, "set_marking()\n")

    create_marking_Node_4states(&(bft->node), (bft->k / NB_CHAR_SUF_PREF) - 1, bft->k, bft->info_per_lvl);
    bft->marked |= 0x1;

    return;
}

/** Function unlocking and the graph locked for vertices marking.
* @param bft is a BFT locked for vertices marking.
*/
void unset_marking(BFT* bft){

    ASSERT_NULL_PTR(bft, "unset_marking()\n")

    if (bft->marked & 0x1){
        delete_marking_Node_4states(&(bft->node), (bft->k / NB_CHAR_SUF_PREF) - 1, bft->k, bft->info_per_lvl);
        bft->marked &= 0xfe;
    }

    return;
}

/** Function marking a k-mer of a BFT with a flag.
* @param flag is the mark to add to a k-mer. It can have value 0, 1, 2 or 3.
* @param bft_kmer is a k-mer obtained via search/iteration over a BFT that must be marked.
* @param bft is a BFT locked for vertices marking.
*/
void set_flag_kmer(uint8_t flag, BFT_kmer* bft_kmer, BFT* bft){

    ASSERT_NULL_PTR(bft_kmer, "set_flag_kmer()\n")
    ASSERT_NULL_PTR(bft, "set_flag_kmer()\n")

    if (flag > 3) ERROR("set_flag_kmer(): a flag can only have as value 0, 1, 2 or 3.\n")
    if ((bft->marked & 0x1) == 0) ERROR("set_flag_kmer(): the graph is not initialized for marking.\n")

    if (is_kmer_in_cdbg(bft_kmer)){

        if (bft_kmer->res->posFilter2 != 0){
            mark_UC_4states((UC*)bft_kmer->res->container, bft_kmer->res->posFilter2,
                            bft_kmer->res->posFilter3, bft_kmer->res->pos_sub_bucket, flag);
        }
        else{
            mark_UC_4states(&(((UC*)((CC*)bft_kmer->res->container)->children)[bft_kmer->res->bucket]),
                            bft_kmer->res->posFilter2, bft_kmer->res->posFilter3, bft_kmer->res->pos_sub_bucket, flag);
        }
    }
    else ERROR("set_flag_kmer(): k-mer is not present in the graph.\n")
}

/** Function getting a k-mer of a BFT with a flag.
* @param bft_kmer is a k-mer obtained via search/iteration over a BFT for which the function returns the flag.
* @param bft is a BFT locked for vertices marking.
*/
uint8_t get_flag_kmer(BFT_kmer* bft_kmer, BFT* bft){

    ASSERT_NULL_PTR(bft_kmer, "get_flag_kmer()\n")
    ASSERT_NULL_PTR(bft, "get_flag_kmer()\n")

    if ((bft->marked & 0x1) == 0) ERROR("get_flag_kmer(): the graph is not initialized for marking.\n")

    if (is_kmer_in_cdbg(bft_kmer)){

        if (bft_kmer->res->posFilter2 != 0){
            return get_mark_UC_4states((UC*)bft_kmer->res->container, bft_kmer->res->posFilter2,
                                       bft_kmer->res->posFilter3, bft_kmer->res->pos_sub_bucket);
        }

        return get_mark_UC_4states(&(((UC*)((CC*)bft_kmer->res->container)->children)[bft_kmer->res->bucket]),
                                   bft_kmer->res->posFilter2, bft_kmer->res->posFilter3, bft_kmer->res->pos_sub_bucket);
    }
    else ERROR("get_flag_kmer(): k-mer is not present in the graph.\n")
}

/** Function locking the graph for traversal.
* It is not necessary to lock the graph for traversal (no insertion can happen during the locking)
* but traversing a locked graph is faster than traversing an unlocked graph.
* @param bft is a BFT to lock for traversal.
*/
void set_neighbors_traversal(BFT* bft){

    ASSERT_NULL_PTR(bft, "set_neighbors_traversal()\n")

    bft->skip_sp = build_skip_nodes(&(bft->node));
    bft->marked |= 0x2;

    return;
}

/** Function unlocking a locked graph for traversal.
* @param bft is a locked BFT for traversal that must be unlocked.
*/
void unset_neighbors_traversal(BFT* bft){

    ASSERT_NULL_PTR(bft, "unset_neighbors_traversal()\n")

    if (bft->marked & 0x2){
        free_skip_nodes(&(bft->node), bft->skip_sp);
        bft->marked &= 0xfd;
    }

    return;
}

/** Function extracting the neighbors of a k-mer.
* @param bft_kmer is a k-mer obtained via search/iteration over a BFT.
* @param bft is a BFT from which was extracted bft_kmer
* @return a pointer to an array of 8 BFT_kmer: positions 0 to 3 are
* the possible predecessors and 4 to 7 the possible successors.
*/
BFT_kmer* get_neighbors(BFT_kmer* bft_kmer, BFT* bft){

    ASSERT_NULL_PTR(bft_kmer, "get_neighbors()\n")
    ASSERT_NULL_PTR(bft, "get_neighbors()\n")

    if (is_kmer_in_cdbg(bft_kmer)){

        int lvl_root = (bft->k / NB_CHAR_SUF_PREF) - 1;
        int nb_bytes_kmer_comp = CEIL(bft->k * 2, SIZE_BITS_UINT_8T);

        BFT_kmer* neighbors = malloc(8 * sizeof(BFT_kmer));

        uint8_t* kmer_comp_tmp = malloc(nb_bytes_kmer_comp * sizeof(uint8_t));
        ASSERT_NULL_PTR(kmer_comp_tmp, "get_neighbors()\n")

        memcpy(kmer_comp_tmp, bft_kmer->kmer_comp, nb_bytes_kmer_comp * sizeof(uint8_t));

        resultPresence* neigh_res = getLeftNeighbors(&(bft->node), bft, lvl_root, kmer_comp_tmp, bft->k);

        for (int i = 0; i < 4; i++){

            neighbors[i].res = malloc(sizeof(resultPresence));
            ASSERT_NULL_PTR(neighbors[i].res, "get_neighbors()\n")

            memcpy(neighbors[i].res, &(neigh_res[i]), sizeof(resultPresence));

            neighbors[i].kmer = malloc((bft->k + 1) * sizeof(char));
            ASSERT_NULL_PTR(neighbors[i].kmer, "get_neighbors()\n")

            neighbors[i].kmer_comp = calloc(nb_bytes_kmer_comp, sizeof(uint8_t));
            ASSERT_NULL_PTR(neighbors[i].kmer_comp, "get_neighbors()\n")

            switch(i){
                case 0: neighbors[i].kmer[0] = 'A'; break;
                case 1: neighbors[i].kmer[0] = 'C'; break;
                case 2: neighbors[i].kmer[0] = 'G'; break;
                case 3: neighbors[i].kmer[0] = 'T'; break;
            }

            memcpy(&(neighbors[i].kmer[1]), bft_kmer->kmer, (bft->k - 1) * sizeof(char));
            neighbors[i].kmer[bft->k] = '\0';

            parseKmerCount(neighbors[i].kmer, bft->k, neighbors[i].kmer_comp, 0);
        }

        memcpy(kmer_comp_tmp, bft_kmer->kmer_comp, nb_bytes_kmer_comp * sizeof(uint8_t));

        free(neigh_res);
        neigh_res = getRightNeighbors(&(bft->node), bft, lvl_root, kmer_comp_tmp, bft->k);

        for (int i = 4; i < 8; i++){

            neighbors[i].res = malloc(sizeof(resultPresence));
            ASSERT_NULL_PTR(neighbors[i].res, "get_neighbors()\n")

            memcpy(neighbors[i].res, &(neigh_res[i-4]), sizeof(resultPresence));

            neighbors[i].kmer = malloc((bft->k + 1) * sizeof(char));
            ASSERT_NULL_PTR(neighbors[i].kmer, "get_neighbors()\n")

            neighbors[i].kmer_comp = calloc(nb_bytes_kmer_comp, sizeof(uint8_t));
            ASSERT_NULL_PTR(neighbors[i].kmer_comp, "get_neighbors()\n")

            memcpy(neighbors[i].kmer, &(bft_kmer->kmer[1]), (bft->k - 1) * sizeof(char));

            switch(i){
                case 4: neighbors[i].kmer[bft->k - 1] = 'A'; break;
                case 5: neighbors[i].kmer[bft->k - 1] = 'C'; break;
                case 6: neighbors[i].kmer[bft->k - 1] = 'G'; break;
                case 7: neighbors[i].kmer[bft->k - 1] = 'T'; break;
            }

            neighbors[i].kmer[bft->k] = '\0';

            parseKmerCount(neighbors[i].kmer, bft->k, neighbors[i].kmer_comp, 0);
        }

        free(neigh_res);
        free(kmer_comp_tmp);
        return neighbors;
    }
    else ERROR("get_neighbors(): k-mer is not present in the graph.\n")
}

/** Function extracting the predecessors of a k-mer.
* @param bft_kmer is a k-mer obtained via search/iteration over a BFT.
* @param bft is a BFT from which was extracted bft_kmer
* @return a pointer to an array of 4 BFT_kmer that are the possible predecessors.
*/
BFT_kmer* get_predecessors(BFT_kmer* bft_kmer, BFT* bft){

    ASSERT_NULL_PTR(bft_kmer, "get_predecessors()\n")
    ASSERT_NULL_PTR(bft, "get_predecessors()\n")

    if (is_kmer_in_cdbg(bft_kmer)){

        int lvl_root = (bft->k / NB_CHAR_SUF_PREF) - 1;
        int nb_bytes_kmer_comp = CEIL(bft->k * 2, SIZE_BITS_UINT_8T);

        BFT_kmer* predecessors = malloc(4 * sizeof(BFT_kmer));

        uint8_t* kmer_comp_tmp = malloc(nb_bytes_kmer_comp * sizeof(uint8_t));
        ASSERT_NULL_PTR(kmer_comp_tmp, "get_predecessors()\n")

        memcpy(kmer_comp_tmp, bft_kmer->kmer_comp, nb_bytes_kmer_comp * sizeof(uint8_t));

        resultPresence* pred = getLeftNeighbors(&(bft->node), bft, lvl_root, kmer_comp_tmp, bft->k);

        for (int i = 0; i < 4; i++){

            predecessors[i].res = malloc(sizeof(resultPresence));
            ASSERT_NULL_PTR(predecessors[i].res, "get_predecessors()\n")

            memcpy(predecessors[i].res, &(pred[i]), sizeof(resultPresence));

            predecessors[i].kmer = malloc((bft->k + 1) * sizeof(char));
            ASSERT_NULL_PTR(predecessors[i].kmer, "get_predecessors()\n")

            predecessors[i].kmer_comp = calloc(nb_bytes_kmer_comp, sizeof(uint8_t));
            ASSERT_NULL_PTR(predecessors[i].kmer_comp, "get_predecessors()\n")

            switch(i){
                case 0: predecessors[i].kmer[0] = 'A'; break;
                case 1: predecessors[i].kmer[0] = 'C'; break;
                case 2: predecessors[i].kmer[0] = 'G'; break;
                case 3: predecessors[i].kmer[0] = 'T'; break;
            }

            memcpy(&(predecessors[i].kmer[1]), bft_kmer->kmer, (bft->k - 1) * sizeof(char));
            predecessors[i].kmer[bft->k] = '\0';

            parseKmerCount(predecessors[i].kmer, bft->k, predecessors[i].kmer_comp, 0);
        }

        free(pred);
        free(kmer_comp_tmp);

        return predecessors;
    }
    else ERROR("get_predecessors(): k-mer is not present in the graph.\n")
}

/** Function extracting the successors of a k-mer.
* @param bft_kmer is a k-mer obtained via search/iteration over a BFT.
* @param bft is a BFT from which was extracted bft_kmer
* @return a pointer to an array of 4 BFT_kmer that are the possible successors.
*/
BFT_kmer* get_successors(BFT_kmer* bft_kmer, BFT* bft){

    ASSERT_NULL_PTR(bft_kmer, "get_successors()\n")
    ASSERT_NULL_PTR(bft, "get_successors()\n")

    if (is_kmer_in_cdbg(bft_kmer)){

        int lvl_root = (bft->k / NB_CHAR_SUF_PREF) - 1;
        int nb_bytes_kmer_comp = CEIL(bft->k * 2, SIZE_BITS_UINT_8T);

        BFT_kmer* successors = malloc(4 * sizeof(BFT_kmer));

        uint8_t* kmer_comp_tmp = malloc(nb_bytes_kmer_comp * sizeof(uint8_t));
        ASSERT_NULL_PTR(kmer_comp_tmp, "get_successors()\n")

        memcpy(kmer_comp_tmp, bft_kmer->kmer_comp, nb_bytes_kmer_comp * sizeof(uint8_t));

        resultPresence* succ = getRightNeighbors(&(bft->node), bft, lvl_root, kmer_comp_tmp, bft->k);

        for (int i = 0; i < 4; i++){

            successors[i].res = malloc(sizeof(resultPresence));
            ASSERT_NULL_PTR(successors[i].res, "get_successors()\n")

            memcpy(successors[i].res, &(succ[i]), sizeof(resultPresence));

            successors[i].kmer = malloc((bft->k + 1) * sizeof(char));
            ASSERT_NULL_PTR(successors[i].kmer, "get_successors()\n")

            successors[i].kmer_comp = calloc(nb_bytes_kmer_comp, sizeof(uint8_t));
            ASSERT_NULL_PTR(successors[i].kmer_comp, "get_successors()\n")

            memcpy(successors[i].kmer, &(bft_kmer->kmer[1]), (bft->k - 1) * sizeof(char));

            switch(i){
                case 0: successors[i].kmer[bft->k - 1] = 'A'; break;
                case 1: successors[i].kmer[bft->k - 1] = 'C'; break;
                case 2: successors[i].kmer[bft->k - 1] = 'G'; break;
                case 3: successors[i].kmer[bft->k - 1] = 'T'; break;
            }

            successors[i].kmer[bft->k] = '\0';

            parseKmerCount(successors[i].kmer, bft->k, successors[i].kmer_comp, 0);
        }

        free(succ);
        free(kmer_comp_tmp);

        return successors;
    }
    else ERROR("get_successors(): k-mer is not present in the graph.\n")
}

/** Function iterating over the k-mers of a BFT.
* This function should be used only when called from a function with a variable number of arguments.
* If not, you must use iterate_over_kmers().
* @param bft is a BFT containing the k-mers to iterate over.
* @param f is a pointer on function that will be called on each k-mer.
*        If f returns 0, the calling function returns.
* @param args should contain all additional arguments to pass to f.
*        They can be extracted in f via its parameter of type va_list.
*/
void v_iterate_over_kmers(BFT* bft, BFT_func_ptr f, va_list args){

    ASSERT_NULL_PTR(bft,"v_iterate_over_kmers()\n")

    int nb_bytes_kmer_comp = CEIL(bft->k * 2, SIZE_BITS_UINT_8T);
    int lvl_root = (bft->k / NB_CHAR_SUF_PREF) - 1;

    uint8_t* kmer = calloc(nb_bytes_kmer_comp, sizeof(uint8_t));
    ASSERT_NULL_PTR(kmer,"v_iterate_over_kmers()\n")

    BFT_kmer* bft_kmer = create_empty_kmer();

    bft_kmer->kmer = malloc((bft->k + 1) * sizeof(char));
    ASSERT_NULL_PTR(bft_kmer->kmer, "v_iterate_over_kmers()\n")

    bft_kmer->kmer_comp = calloc(nb_bytes_kmer_comp, sizeof(uint8_t));
    ASSERT_NULL_PTR(bft_kmer->kmer_comp, "v_iterate_over_kmers()\n")

    bft_kmer->res = create_resultPresence();

    bft->marked |= 0x2;

    iterate_over_kmers_from_node(&(bft->node), bft, lvl_root, kmer, bft_kmer, bft->k, 0, 0, f, args);
    free_BFT_kmer(bft_kmer, 1);

    bft->marked &= 0xfd;

    free(kmer);
}

/** Function iterating over the k-mers of a BFT.
* @param bft is a BFT containing the k-mers to iterate over.
* @param f is a pointer on function that will be called on each k-mer.
*        If f returns 0, the calling function returns.
* @param ... are the additional arguments that must be transmitted to f.
*        They can be extracted in f via its parameter of type va_list.
*/
void iterate_over_kmers(BFT* bft, BFT_func_ptr f, ... ){

    ASSERT_NULL_PTR(bft,"iterate_over_kmers()\n")

    va_list args;

    va_start(args, f);

    int nb_bytes_kmer_comp = CEIL(bft->k * 2, SIZE_BITS_UINT_8T);
    int lvl_root = (bft->k / NB_CHAR_SUF_PREF) - 1;

    uint8_t* kmer = calloc(nb_bytes_kmer_comp, sizeof(uint8_t));
    ASSERT_NULL_PTR(kmer,"iterate_over_kmers()\n")

    BFT_kmer* bft_kmer = create_empty_kmer();

    bft_kmer->kmer = malloc((bft->k + 1) * sizeof(char));
    ASSERT_NULL_PTR(bft_kmer->kmer, "iterate_over_kmers()\n")

    bft_kmer->kmer_comp = calloc(nb_bytes_kmer_comp, sizeof(uint8_t));
    ASSERT_NULL_PTR(bft_kmer->kmer_comp, "iterate_over_kmers()\n")

    bft_kmer->res = create_resultPresence();

    bft->marked |= 0x2;

    iterate_over_kmers_from_node(&(bft->node), bft, lvl_root, kmer, bft_kmer, bft->k, 0, 0, f, args);
    free_BFT_kmer(bft_kmer, 1);

    va_end(args);

    bft->marked &= 0xfd;

    free(kmer);
}

/** Function for prefix matching over the k-mers of a BFT.
* @param bft is a BFT containing the k-mers to match.
* @param prefix is string containing a prefix the k-mers must match.
* @param f is a pointer on function that will be called on each k-mer matching the prefix.
*        If f returns 0, the calling function returns.
* @param ... are the additional arguments that must be transmitted to f.
*        They can be extracted in f via its parameter of type va_list.
* @return A boolean indicating whether at least one k-mer matched the prefix (true) or not (false).
*/
bool prefix_matching(BFT* bft, char* prefix, BFT_func_ptr f, ...){

    ASSERT_NULL_PTR(bft, "prefix_matching()\n")
    ASSERT_NULL_PTR(prefix, "prefix_matching()\n")

    va_list args;

    va_start(args, f);

    int length_prefix = strlen(prefix);

    if (length_prefix > bft->k) ERROR("prefix_matching(): Prefix length is larger than k-mer length.\n")
    if (length_prefix == 0) ERROR("prefix_matching(): Prefix length is 0.\n")

    int size_kmer_comp = CEIL(bft->k * 2, SIZE_BITS_UINT_8T);

    size_t return_value;

    uint8_t* prefix_comp = calloc(size_kmer_comp, sizeof(uint8_t));
    ASSERT_NULL_PTR(prefix_comp, "prefix_matching()\n")

    uint8_t* shifted_prefix_comp = malloc(size_kmer_comp * sizeof(uint8_t));
    ASSERT_NULL_PTR(shifted_prefix_comp, "prefix_matching()\n")

    BFT_kmer* bft_kmer = create_empty_kmer();

    bft_kmer->kmer = malloc((bft->k + 1) * sizeof(char));
    ASSERT_NULL_PTR(bft_kmer->kmer, "prefix_matching()\n")

    bft_kmer->kmer_comp = calloc(CEIL(bft->k * 2, SIZE_BITS_UINT_8T), sizeof(uint8_t));
    ASSERT_NULL_PTR(bft_kmer->kmer_comp, "prefix_matching()\n")

    bft_kmer->res = create_resultPresence();

    if (parseKmerCount(prefix, length_prefix, prefix_comp, 0) == 0)
        ERROR("prefix_matching(): Non-ACGT char. encountered in prefix.\n")

    memcpy(shifted_prefix_comp, prefix_comp, size_kmer_comp * sizeof(uint8_t));

    return_value = v_prefix_matching(&(bft->node), bft, bft->k / NB_CHAR_SUF_PREF - 1,
                                     shifted_prefix_comp, prefix_comp, length_prefix, length_prefix,
                                     bft_kmer, f, args);

    free_BFT_kmer(bft_kmer, 1);
    free(prefix_comp);
    free(shifted_prefix_comp);

    va_end(args);

    if (return_value > 1) return false;
    return true;
}

bool prefix_matching_custom(BFT* bft, char* prefix, resultPresence** junction_prefix, BFT_func_ptr f, ...){

    ASSERT_NULL_PTR(bft, "prefix_matching()\n")
    ASSERT_NULL_PTR(prefix, "prefix_matching()\n")

    va_list args;

    va_start(args, f);

    int length_prefix = strlen(prefix);

    if (length_prefix > bft->k) ERROR("prefix_matching(): Prefix length is larger than k-mer length.\n")
    if (length_prefix == 0) ERROR("prefix_matching(): Prefix length is 0.\n")

    int size_kmer_comp = CEIL(bft->k * 2, SIZE_BITS_UINT_8T);

    size_t return_value;

    uint8_t* prefix_comp = calloc(size_kmer_comp, sizeof(uint8_t));
    ASSERT_NULL_PTR(prefix_comp, "prefix_matching()\n")

    uint8_t* shifted_prefix_comp = malloc(size_kmer_comp * sizeof(uint8_t));
    ASSERT_NULL_PTR(shifted_prefix_comp, "prefix_matching()\n")

    BFT_kmer* bft_kmer = create_empty_kmer();

    bft_kmer->kmer = malloc((bft->k + 1) * sizeof(char));
    ASSERT_NULL_PTR(bft_kmer->kmer, "prefix_matching()\n")

    bft_kmer->kmer_comp = calloc(CEIL(bft->k * 2, SIZE_BITS_UINT_8T), sizeof(uint8_t));
    ASSERT_NULL_PTR(bft_kmer->kmer_comp, "prefix_matching()\n")

    bft_kmer->res = create_resultPresence();

    if (parseKmerCount(prefix, length_prefix, prefix_comp, 0) == 0)
        ERROR("prefix_matching(): Non-ACGT char. encountered in prefix.\n")

    memcpy(shifted_prefix_comp, prefix_comp, size_kmer_comp * sizeof(uint8_t));

    return_value = v_prefix_matching_custom(&(bft->node), bft, bft->k / NB_CHAR_SUF_PREF - 1,
                                            shifted_prefix_comp, prefix_comp, length_prefix, length_prefix,
                                            bft_kmer, junction_prefix, f, args);

    free_BFT_kmer(bft_kmer, 1);
    free(prefix_comp);
    free(shifted_prefix_comp);

    va_end(args);

    if (return_value > 1) return false;
    return true;
}

/** Function writing a BFT to disk.
* @param bft is the BFT to write on disk.
* @param filename is the name of the file in which bft will be written.
* @param compress_annotations is a boolean indicating if the annotations of the BFT
            must be compressed before writing to disk.
*/
void write_BFT(BFT* bft, char* filename, bool compress_annotations){

    ASSERT_NULL_PTR(bft, "write_BFT()\n")
    ASSERT_NULL_PTR(filename, "write_BFT()\n")

    if (compress_annotations) compress_annotations_disk(bft, filename);
    else write_BFT_Root(bft, filename, false);

    return;
}

/** Function loading a BFT from disk.
* @param filename_and_path is the path and name of the file in which the BFT to load is be written.
*/
BFT* load_BFT(char* filename_and_path){

    ASSERT_NULL_PTR(filename_and_path, "load_BFT()\n")

    return read_BFT_Root(filename_and_path);
}

/** Function querying a BFT for a sequence.
* @param bft is a BFT to be queried.
* @param sequence is a string to query.
* @param threshold is a float (0 < threshold <= 1) indicating the minimum percentage of k-mers from the queried
*           sequence that must be present in a genome to have the queried sequence reported present in this genome.
* @param canonical_search is a boolean indicating if the searched k-mers of the queried sequence must be canonical
            (lexicographically smaller one between a k-mer and its reverse-complement) or not.
* @return a pointer to a sorted array of genome identifiers in which the queried sequence occurs (according to
*           parameter threshold) or NULL if the queried sequence is not present in at least one genome
*           (according to parameter threshold). The first element of the array (position 0) indicates how many
*           ids are in this array.
*/
uint32_t* query_sequence(BFT* bft, char* sequence, double threshold, bool canonical_search){

    ASSERT_NULL_PTR(bft, "query_sequence()\n")
    ASSERT_NULL_PTR(sequence, "query_sequence()\n")

    if (threshold <= 0) ERROR("query_sequence(): the threshold must be superior to 0.\n");
    if (threshold > 1) ERROR("query_sequence(): the threshold must be inferior or equal to 1.\n");

    BFT_kmer* bft_kmer;
    BFT_annotation* bft_annot;

    int64_t it_kmers_query;
    int64_t nb_kmers_found;
    int64_t nb_kmers_colors_set;
    int64_t nb_kmers_query_min;

    uint32_t it_annot, it_colors;
    uint32_t nb_colors = 0;

    uint32_t* colors_set;

    uint32_t* count_colors = calloc(bft->nb_genomes, sizeof(uint32_t));
    ASSERT_NULL_PTR(count_colors, "query_sequence()\n");

    char* kmer;

    char* kmer_not_rev_comp = malloc((bft->k + 1) * sizeof(char));
    ASSERT_NULL_PTR(kmer_not_rev_comp, "query_sequence()\n");

    char* kmer_rev_comp = malloc((bft->k + 1) * sizeof(char));
    ASSERT_NULL_PTR(kmer_rev_comp, "query_sequence()\n");

    kmer_not_rev_comp[bft->k] = '\0';
    kmer_rev_comp[bft->k] = '\0';

    nb_kmers_colors_set = strlen(sequence) - bft->k + 1;
    if (nb_kmers_colors_set < 0) printf("query_sequence(): query %s is too small and must be at least of length k.\n", sequence);

    nb_kmers_query_min = (int64_t) ceil(nb_kmers_colors_set * threshold);

    kmer = kmer_not_rev_comp;

    for (it_kmers_query = 0, nb_kmers_found = 0; it_kmers_query < nb_kmers_colors_set; it_kmers_query++){

        memcpy(kmer_not_rev_comp, &sequence[it_kmers_query], bft->k * sizeof(char));

        if (canonical_search){

            reverse_complement(kmer_not_rev_comp, kmer_rev_comp, bft->k);

            if (strcmp(kmer_not_rev_comp, kmer_rev_comp) >= 0) kmer = kmer_rev_comp;
            else kmer = kmer_not_rev_comp;
        }

        if (!is_substring_IUPAC(kmer)){

            bft_kmer = get_kmer(kmer, bft);

            if (is_kmer_in_cdbg(bft_kmer)){

                nb_kmers_found++;

                bft_annot = get_annotation(bft_kmer);
                colors_set = get_list_id_genomes(bft_annot, bft);
                free_BFT_annotation(bft_annot);

                for (it_annot = 1; it_annot <= colors_set[0]; it_annot++){

                    nb_colors += (count_colors[colors_set[it_annot]] == 0);
                    count_colors[colors_set[it_annot]]++;
                }

                free(colors_set);
            }

            free_BFT_kmer(bft_kmer, 1);
        }

        if (nb_kmers_found + nb_kmers_colors_set - it_kmers_query < nb_kmers_query_min) break;
    }

    colors_set = malloc((nb_colors + 1) * sizeof(uint32_t));
    ASSERT_NULL_PTR(colors_set, "query_sequence()\n");

    it_colors = 0;

    for (it_annot = 0; (it_annot < bft->nb_genomes) && nb_colors; it_annot++){

        if (count_colors[it_annot]){

            nb_colors--;

            if (count_colors[it_annot] >= nb_kmers_query_min){

                it_colors++;
                colors_set[it_colors] = it_annot;
            }
        }
    }

    colors_set[0] = it_colors;

    colors_set = realloc(colors_set, (it_colors + 1) * sizeof(uint32_t));
    ASSERT_NULL_PTR(colors_set, "query_sequence()\n");

    free(kmer_not_rev_comp);
    free(kmer_rev_comp);
    free(count_colors);

    return colors_set;
}
