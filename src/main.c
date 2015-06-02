#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <libgen.h>

#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>

#include <jemalloc/jemalloc.h>
const char* malloc_conf = "narenas:1,tcache:false,lg_dirty_mult:8,lg_chunk:22";

#include <Judy.h>

#include "./../lib/UC_annotation.h"
#include "./../lib/insertNode.h"
#include "./../lib/branchingNode.h"
#include "./../lib/deleteColorsNode.h"
#include "./../lib/fasta.h"
#include "./../lib/printMemory.h"
#include "./../lib/replaceAnnotation.h"
#include "./../lib/write_to_disk.h"

#define PRINT_EVERY_X_KMERS 1000000

#define SIZE_BUFFER 4096 //Size of the buffer (in bytes) used to store kmers read from input files
#define TRESH_DEL_ANNOT 43

#include "./../lib/kseq.h"
//Initialize parser for FASTx files
KSEQ_INIT(int, read, SIZE_BUFFER)

void insertKmers(Root* restrict root,
                 uint8_t* restrict array_kmers,
                 int size_kmers,
                 int nb_kmers,
                 int id_genome,
                 ptrs_on_func* restrict func_on_types,
                 annotation_inform* ann_inf,
                 resultPresence* res);

void insert_Genomes_from_KmerFiles(Root* root, char** filenames, int binary_files,
                                   int size_kmer, int id_start_genome, ptrs_on_func* func_on_types);
void insert_Genomes_from_FASTxFiles(Root* root, char** filenames, int size_kmer,
                                    int id_start_genome, ptrs_on_func* func_on_types);

int queryBFT_kmerPresences_from_KmerFiles(Root* root, char* query_filename,
                                          int binary_file, char* output_filename, int size_kmer);
int queryBFT_kmerBranching_from_KmerFiles(Root* root, char* query_filename,
                                          int binary_file, int size_kmer);

int get_nb_cplx_nodes_from_KmerCounting(Root* tree, char* name_file, int size_kmer, uint16_t** skip_node_root,
                                        ptrs_on_func* func_on_types, annotation_inform* ann_inf);
Root* get_nb_cplx_nodes_from_FASTx(Root* tree, int size_kmer, uint16_t** skip_node_root, ptrs_on_func* func_on_types,
                                   annotation_inform* ann_inf, resultPresence* res);

extern void freeRoot(Root* root);
extern void time_spent(struct timeval *start_time, struct timeval *end_time, struct timeval *resulting_time);

int main(int argc, char *argv[])
{

    Root* root = NULL;

    int i = 0;
    int j = 0;
    int cpt = 0;
    int size_kmer = 27;
    int binary_files = 0;
    int nb_files_2_read = 0;

    const char csv_ext[5] = ".csv\0";

    char buffer[2048];

    char* dot;
    char* str_tmp;
    char* filename_output;

    char** filenames = NULL;
    char** paths_and_names = NULL;

    FILE* file_input = NULL;
    FILE* file_tmp = NULL;

    ptrs_on_func* func_on_types = NULL;

    if (argc == 1){
        ERROR("\nUsage:\n"
              "./bft build k {fastx|kmers|kmers_comp} list_genome_files output_file [Options]\n"
              "./bft load file_bft [-add_genomes {fastx|kmers|kmers_comp} list_genome_files output_file] [Options]\n"
              "\nOptions:\n"
              "[-query_kmers {kmers|kmers_comp} list_kmer_files]\n"
              "[-query_branching {kmers|kmers_comp} list_kmer_files]\n\n")
    }
    else {

        if (strcmp("build", argv[1]) == 0){

            size_kmer = atoi(argv[2]); //Length k argument reading

            //Test if k is valid
            if (size_kmer < 0) ERROR("Provided length k (for k-mers) is < 0\n")
            else if (size_kmer == 0) ERROR("Provided length k (for k-mers) either 0 or not a number\n")
            else if (size_kmer > 63) ERROR("Length k (for k-mers) cannot be superior to 63\n")
            else if (size_kmer%SIZE_SEED != 0) ERROR("Length k (for k-mers) must be a multiple of 9\n")

            //Open the file containing the input files
            if ((file_input = fopen(argv[4], "r")) == NULL) ERROR("Invalid list_genome_files.\n")

            while (fgets(buffer, 2048, file_input)){ //Test if the input files can be opened and read

                buffer[strcspn(buffer, "\r\n")] = 0;

                if ((file_tmp = fopen(buffer, "r")) == NULL){
                    fprintf(stderr, "Invalid input file at line %d of list_genome_files.\n", nb_files_2_read);
                    exit(EXIT_FAILURE);
                }

                nb_files_2_read++;

                fclose(file_tmp);
            }

            fclose(file_input);

            for (i=6; i<argc; i+=3){ //Test if we can open the files for querying the k-mers/branching vertices

                if ((strcmp("-query_kmers", argv[i]) == 0) || (strcmp("-query_branching", argv[i]) == 0)){ //User wants to query the BFT for k-mers

                    if ((strcmp("kmers_comp", argv[i+1]) != 0) && (strcmp("kmers", argv[i+1]) != 0)){
                        if (strcmp("-query_kmers", argv[i]) == 0){
                            ERROR("Unrecognized type of input files for -query_kmers.\n"
                                  "Choice must be 'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")
                        }
                        else{
                            ERROR("Unrecognized type of input files for -query_branching.\n"
                                  "Choice must be 'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")
                        }
                    }

                    if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid k-mer queries file.\n")

                    while (fgets(buffer, 2048, file_input)){ //Test if the input files can be opened and read

                        buffer[strcspn(buffer, "\r\n")] = 0;

                        if ((file_tmp = fopen(buffer, "r")) == NULL){
                            fprintf(stderr, "Invalid input file at line %d of the list of k-mer queries files.\n", cpt);
                            exit(EXIT_FAILURE);
                        }

                        cpt++;

                        fclose(file_tmp);
                    }

                    fclose(file_input);
                }
                else{
                    fprintf(stderr, "Unrecognized command %s.\n", argv[i]);
                    exit(EXIT_FAILURE);
                }

                cpt = 0;
            }

            if ((file_input = fopen(argv[4], "r")) == NULL) ERROR("Invalid list_genome_files.\n")

            filenames = malloc(nb_files_2_read*sizeof(char*)); //Allocate the array of filenames
            ASSERT_NULL_PTR(filenames,"main()")
            paths_and_names = malloc(nb_files_2_read*sizeof(char*)); //Allocate the array of paths + filenames
            ASSERT_NULL_PTR(paths_and_names,"main()")

            i = 0;

            while (fgets(buffer, 2048, file_input)){ //Copy the filenames in an array, the paths and filenames in another array

                str_tmp = basename(buffer);

                filenames[i] = malloc((strlen(str_tmp)+1)*sizeof(char));
                ASSERT_NULL_PTR(filenames[i],"main()")

                strcpy(filenames[i], str_tmp);

                paths_and_names[i] = malloc((strlen(buffer)+1)*sizeof(char));
                ASSERT_NULL_PTR(paths_and_names[i],"main()")

                strcpy(paths_and_names[i], buffer);

                paths_and_names[i][strcspn(paths_and_names[i], "\r\n")] = 0;

                i++;
            }

            fclose(file_input);

            root = createRoot(filenames, nb_files_2_read, size_kmer); //Create a BFT
            func_on_types = create_ptrs_on_func(SIZE_SEED, root->k);

            //Read and test if the type of input files is valid
            //Insert k-mers of the input files in the BFT
            if (strcmp("kmers_comp", argv[3]) == 0){
                insert_Genomes_from_KmerFiles(root, paths_and_names, 1, root->k, 0, func_on_types);
            }
            else if (strcmp("kmers", argv[3]) == 0){
                insert_Genomes_from_KmerFiles(root, paths_and_names, 0, root->k, 0, func_on_types);
            }
            else if (strcmp("fastx", argv[3]) == 0){
                insert_Genomes_from_FASTxFiles(root, paths_and_names, root->k, 0, func_on_types);
            }
            else
                ERROR("Unrecognized type of input files.\nChoice must be 'fastx' for FASTA/FASTQ files, "
                      "'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")

            write_Root(root, argv[5], func_on_types);

            for (i=6; i<argc; i+=3){

                binary_files = 0;

                if (strcmp("-query_kmers", argv[i]) == 0){

                    if (strcmp("kmers_comp", argv[i+1]) == 0) binary_files = 1;
                    if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid k-mer queries files list.\n")

                    while (fgets(buffer, 2048, file_input)){

                        buffer[strcspn(buffer, "\r\n")] = 0;

                        filename_output = malloc((strlen(basename(buffer))+4) * sizeof(char));
                        ASSERT_NULL_PTR(filename_output, "main()")

                        strcpy(filename_output, basename(buffer));

                        if ((dot = strrchr(filename_output, '.')) != NULL) strcpy(dot, csv_ext);
                        else strcpy(&(filename_output[strlen(filename_output)]), csv_ext);

                        //Query the BFT for presence of k-mer queries
                        printf("\nNb k-mers present = %d\n", queryBFT_kmerPresences_from_KmerFiles(root, buffer, binary_files, filename_output, root->k));

                        free(filename_output);
                    }

                    fclose(file_input);
                }
                else if (strcmp("-query_branching", argv[i]) == 0){

                    if (strcmp("kmers_comp", argv[i+1]) == 0) binary_files = 1;
                    if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid branching k-mer queries files list.\n")

                    while (fgets(buffer, 2048, file_input)){

                        buffer[strcspn(buffer, "\r\n")] = 0;
                        //Query the BFT for the number of k-mer queries that are branching vertices
                        printf("\nNb branching k-mers = %d\n", queryBFT_kmerBranching_from_KmerFiles(root, buffer, binary_files, root->k));
                    }

                    fclose(file_input);
                }
                else{
                    fprintf(stderr, "Unrecognized command %s.\n", argv[i]);
                    exit(EXIT_FAILURE);
                }
            }

            if (paths_and_names != NULL){
                int i = 0;
                for (i=0; i<root->nb_genomes; i++) free(paths_and_names[i]);
                free(paths_and_names);
            }

            freeRoot(root);
            free(func_on_types);
        }
        else if (strcmp("load", argv[1]) == 0){

            root = read_Root(argv[2]);

            for (i=3; i<argc; i+=3){ //Test if we can open the files for querying the k-mers/branching vertices

                if ((strcmp("-query_kmers", argv[i]) == 0) || (strcmp("-query_branching", argv[i]) == 0)){ //User wants to query the BFT for k-mers

                    if ((strcmp("kmers_comp", argv[i+1]) != 0) && (strcmp("kmers", argv[i+1]) != 0)){
                        if (strcmp("-query_kmers", argv[i]) == 0){
                            ERROR("Unrecognized type of input files for -query_kmers.\n"
                                  "Choice must be 'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")
                        }
                        else{
                            ERROR("Unrecognized type of input files for -query_branching.\n"
                                  "Choice must be 'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")
                        }
                    }

                    if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid k-mer queries file.\n")

                    while (fgets(buffer, 2048, file_input)){ //Test if the input files can be opened and read

                        buffer[strcspn(buffer, "\r\n")] = 0;

                        if ((file_tmp = fopen(buffer, "r")) == NULL){
                            fprintf(stderr, "Invalid input file at line %d of the list of k-mer queries files.\n", cpt);
                            exit(EXIT_FAILURE);
                        }

                        cpt++;

                        fclose(file_tmp);
                    }

                    fclose(file_input);
                }
                else if (strcmp("-add_genomes", argv[i]) == 0){

                    if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid list_genome_files.\n")

                    while (fgets(buffer, 2048, file_input)){ //Test if the input files can be opened and read

                        buffer[strcspn(buffer, "\r\n")] = 0;

                        if ((file_tmp = fopen(buffer, "r")) == NULL){
                            fprintf(stderr, "Invalid input file at line %d of list_genome_files.\n", nb_files_2_read);
                            exit(EXIT_FAILURE);
                        }

                        nb_files_2_read++;

                        fclose(file_tmp);
                    }

                    fclose(file_input);

                    i++;
                }
                else{
                    fprintf(stderr, "Unrecognized command %s.\n", argv[i]);
                    exit(EXIT_FAILURE);
                }

                cpt = 0;
            }

            for (i=3; i<argc; i+=3){

                binary_files = 0;

                if (strcmp("-query_kmers", argv[i]) == 0){

                    if (strcmp("kmers_comp", argv[i+1]) == 0) binary_files = 1;
                    if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid k-mer query files list.\n")

                    while (fgets(buffer, 2048, file_input)){

                        buffer[strcspn(buffer, "\r\n")] = 0;

                        filename_output = malloc((strlen(basename(buffer))+5) * sizeof(char));
                        ASSERT_NULL_PTR(filename_output, "main()")

                        strcpy(filename_output, basename(buffer));

                        if ((dot = strrchr(filename_output, '.')) != NULL) strcpy(dot, csv_ext);
                        else strcpy(&(filename_output[strlen(filename_output)]), csv_ext);

                        //Query the BFT for presence of k-mer queries
                        printf("\nNb k-mers present = %d\n", queryBFT_kmerPresences_from_KmerFiles(root, buffer, binary_files, filename_output, root->k));

                        free(filename_output);
                    }

                    fclose(file_input);
                }
                else if (strcmp("-query_branching", argv[i]) == 0){

                    if (strcmp("kmers_comp", argv[i+1]) == 0) binary_files = 1;
                    if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid branching k-mer query file list.\n")

                    while (fgets(buffer, 2048, file_input)){

                        buffer[strcspn(buffer, "\r\n")] = 0;
                        //Query the BFT for the number of k-mer queries that are branching vertices
                        printf("\nNb branching k-mers = %d\n", queryBFT_kmerBranching_from_KmerFiles(root, buffer, binary_files, root->k));
                    }

                    fclose(file_input);
                }
                else if (strcmp("-add_genomes", argv[i]) == 0){

                    if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid list_genome_files.\n")

                    paths_and_names = malloc(nb_files_2_read*sizeof(char*)); //Allocate the array of paths + filenames
                    ASSERT_NULL_PTR(paths_and_names,"main()")
                    root->filenames = realloc(root->filenames, (root->nb_genomes + nb_files_2_read) * sizeof(char*)); //Allocate the array of paths + filenames
                    ASSERT_NULL_PTR(root->filenames,"main()")

                    j = 0;

                    while (fgets(buffer, 2048, file_input)){ //Copy the filenames in an array, the paths and filenames in another array

                        str_tmp = basename(buffer);

                        root->filenames[j + root->nb_genomes] = malloc((strlen(str_tmp)+1)*sizeof(char));
                        ASSERT_NULL_PTR(root->filenames[j + root->nb_genomes],"main()")

                        strcpy(root->filenames[j + root->nb_genomes], str_tmp);

                        paths_and_names[j] = malloc((strlen(buffer)+1)*sizeof(char));
                        ASSERT_NULL_PTR(paths_and_names[j],"main()")

                        strcpy(paths_and_names[j], buffer);

                        paths_and_names[j][strcspn(paths_and_names[j], "\r\n")] = 0;

                        j++;
                    }

                    fclose(file_input);

                    func_on_types = create_ptrs_on_func(SIZE_SEED, root->k);

                    root->nb_genomes += nb_files_2_read;

                    //Read and test if the type of input files is valid
                    //Insert k-mers of the input files in the BFT
                    if (strcmp("kmers_comp", argv[i+1]) == 0){
                        insert_Genomes_from_KmerFiles(root, paths_and_names, 1, root->k, root->nb_genomes - nb_files_2_read, func_on_types);
                    }
                    else if (strcmp("kmers", argv[i+1]) == 0){
                        insert_Genomes_from_KmerFiles(root, paths_and_names, 0, root->k, root->nb_genomes - nb_files_2_read, func_on_types);
                    }
                    else if (strcmp("fastx", argv[i+1]) == 0){
                        insert_Genomes_from_FASTxFiles(root, paths_and_names, root->k, root->nb_genomes - nb_files_2_read, func_on_types);
                    }
                    else
                        ERROR("Unrecognized type of input files.\nChoice must be 'fastx' for FASTA/FASTQ files, "
                              "'kmers' for k-mers files or 'kmers_comp' for compressed k-mer files.\n")

                    write_Root(root, argv[i+3], func_on_types);

                    free(func_on_types);

                    i++;
                }
                else{
                    fprintf(stderr, "Unrecognized command %s.\n", argv[i]);
                    exit(EXIT_FAILURE);
                }
            }

            if (paths_and_names != NULL){
                int i = 0;
                for (i=0; i<nb_files_2_read; i++) free(paths_and_names[i]);
                free(paths_and_names);
            }

            freeRoot(root);
        }
        else{
            fprintf(stderr, "Unrecognized command %s.\n", argv[1]);
            exit(EXIT_FAILURE);
        }
    }

    return EXIT_SUCCESS;
}

/* ---------------------------------------------------------------------------------------------------------------
*  insertKmers(root, array_kmers, size_kmers, nb_kmers, id_genome, func_on_types, ann_inf, res, annot_sorted)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert kmers into the BFT
*  ---------------------------------------------------------------------------------------------------------------
*  root: ptr to the root of a BFT
*  array_kmers: array of kmers
*  size_kmers: length k of kmers in array_kmers
*  nb_kmers: number of kmers in array_kmers
*  id_genome: genome id to which belongs the kmers in array_kmers
*  func_on_types: ptr to ptrs_on_func structure, used to manipulate compressed container field children_type
*  ann_inf: ptr to annotation_inform structure, used to make the transition between reading and modifying an annot
*  res: ptr to resultPresence structure, used to store information about the presence of a suffix prefix in a CC
*  ---------------------------------------------------------------------------------------------------------------
*/
void insertKmers(Root* restrict root,
                 uint8_t* restrict array_kmers,
                 int size_kmers,
                 int nb_kmers,
                 int id_genome,
                 ptrs_on_func* restrict func_on_types,
                 annotation_inform* ann_inf,
                 resultPresence* res){

    ASSERT_NULL_PTR(root,"insertKmers()")
    ASSERT_NULL_PTR(func_on_types,"insertKmers()")
    ASSERT_NULL_PTR(ann_inf,"insertKmers()")
    ASSERT_NULL_PTR(res,"insertKmers()")

    int i = 0;
    int nb_bytes = CEIL(size_kmers*2, SIZE_CELL); //Nb bytes used to store a k-mer

    uint8_t kmer[nb_bytes];

    for (i=0; i<nb_kmers; i++){

        memcpy(kmer, &(array_kmers[i*nb_bytes]), nb_bytes * sizeof(uint8_t));

        insertKmer_Node(&(root->node), &(root->node), &(array_kmers[i*nb_bytes]), size_kmers,
                        kmer, size_kmers, id_genome, func_on_types, ann_inf, res, root->comp_set_colors);
    }
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
void insert_Genomes_from_KmerFiles(Root* root, char** filenames, int binary_files, int size_kmer, int id_start_genome, ptrs_on_func* func_on_types){

    ASSERT_NULL_PTR(root,"insert_Genomes_from_KmerFiles()")
    ASSERT_NULL_PTR(filenames,"insert_Genomes_from_KmerFiles()")

    //Structure for time measurement
    struct timeval tval_before, tval_after, tval_last, tval_result;
    gettimeofday(&tval_before, NULL);
    tval_last = tval_before;

    FILE* file;

    //JudyArray structure, needed for compressing the colors sets of the BFT
    Pvoid_t PJArray = (PWord_t)NULL;
    Word_t Rc_word;

    annotation_array_elem* comp_set_colors_tmp = NULL;

    annotation_inform* ann_inf = calloc(1,sizeof(annotation_inform));
    ASSERT_NULL_PTR(ann_inf,"insert_Genomes_from_FASTx()")

    resultPresence* res = create_resultPresence();

    //Buffer used to store k-mers read in the input files
    uint8_t* array_kmers = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(array_kmers,"insert_Genomes_from_KmerFiles()")

    uint8_t* kmer;

    char* line = calloc(100, sizeof(char));
    ASSERT_NULL_PTR(line,"insert_Genomes_from_KmerFiles()")

    uint16_t** skip_node_root = NULL;

    int i = 0;
    int j = 0;
    int k = 0;
    int nb_bytes_kmer = CEIL(size_kmer*2, SIZE_CELL);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer;

    size_t return_fread;

    int length_comp_set_colors_tmp = 0;

    double count = 0;

    uint64_t kmers_read;

    for (i=id_start_genome; i<root->nb_genomes; i++){ //For each file in input

        kmers_read = 0;
        k = 0;
        j = 0;

        file = fopen(filenames[i-id_start_genome], "r");
        ASSERT_NULL_PTR(file,"insert_Genomes_from_KmerFiles()")

        printf("\nFile %d: %s\n\n", i, filenames[i-id_start_genome]);

        if (binary_files){

            if (fgets(line, 100, file) != NULL) k = atoi(line);
            else ERROR("Cannot read header of the file")

            if (fgets(line, 100, file) != NULL) printf("%d %d-mers in the file\n\n", atoi(line), k);
            else ERROR("Cannot read header of the file")

            while ((!ferror(file)) && (!feof(file))){

                return_fread = fread(array_kmers, (size_t)nb_bytes_kmer, (size_t)nb_kmer_in_buf, file);

                insertKmers(root, array_kmers, size_kmer, return_fread, i, func_on_types, ann_inf, res);

                memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));

                if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+return_fread)%PRINT_EVERY_X_KMERS)){
                    printf("%" PRIu64 " kmers read\n", kmers_read+return_fread);
                    //break;
                }

                kmers_read += return_fread;
            }
        }
        else {
            while (fgets(line, 100, file) != NULL){

                if (parseKmerCount(line, size_kmer, array_kmers, k) == 1){
                    k += nb_bytes_kmer;
                    j++;

                    if (j == nb_kmer_in_buf){
                        insertKmers(root, array_kmers, size_kmer, nb_kmer_in_buf, i, func_on_types, ann_inf, res);

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

            insertKmers(root, array_kmers, size_kmer, j, i, func_on_types, ann_inf, res);
            kmers_read += j;

            memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));
        }

        fclose(file);

        if ((i > 5) && (i%TRESH_DEL_ANNOT == 0)){

            load_annotation_from_Node(&(root->node), size_kmer, func_on_types, &PJArray, root->comp_set_colors);

            comp_set_colors_tmp = root->comp_set_colors;
            length_comp_set_colors_tmp = root->length_comp_set_colors;

            root->comp_set_colors = sort_annotations(&PJArray, &(root->length_comp_set_colors));

            compress_annotation_from_Node(&(root->node), size_kmer, func_on_types, &PJArray, comp_set_colors_tmp);

            free_annotation_array_elem(comp_set_colors_tmp, length_comp_set_colors_tmp);

            JSLFA(Rc_word, PJArray);
        }

        //Delete unecessary annotations
        /*skip_node_root = build_skip_nodes(&(root->node), size_kmer, func_on_types);

        printf("\nCreating marking structure\n");

        create_marking_Node_4states(&(root->node), size_kmer, func_on_types);

        printf("\nGetting complex nodes\n");

        count = get_nb_cplx_nodes_from_KmerCounting(root, filenames[i-id_start_genome], size_kmer, skip_node_root, &annot_sorted, func_on_types, ann_inf);

        if (count/((double)kmers_read) <= TRESH_DEL_ANNOT){

            kmer = calloc(CEIL(size_kmer*2, SIZE_CELL), sizeof(uint8_t));
            ASSERT_NULL_PTR(kmer,"insert_Genomes_from_KmerFiles()")

            count = deleteColors_from_branchingNodes(&(root->node), &(root->node), kmer, size_kmer, 0, 0, size_kmer, i, func_on_types, skip_node_root, ann_inf, annot_sorted);
            printf("\nNumber of annotations that could be erased = %f\n", count);

            if (((double)kmers_read-count)/((double)kmers_read) <= TRESH_DEL_ANNOT){

                count = resize_annotation_Node(&(root->node), size_kmer, func_on_types);
                printf("\n%f unnecessary annoations has been deleted, BFT has been reallocated\n\n", count);
            }
            else{
                printf("\nFor this genome, memory cannot be optimized by deleting unnecessary annotations in graph's paths\n\n");
                delete_marking_Node_4states(&(root->node), size_kmer, func_on_types);
            }

            free(kmer);
        }
        else {
            printf("\nFor this genome, memory cannot be optimized by deleting unnecessary annotations in graph's paths\n\n");
            delete_marking_Node_4states(&(root->node), size_kmer, func_on_types);
        }

        free_skip_nodes(&(root->node), skip_node_root);

        if (annot_sorted != NULL) free(annot_sorted);*/

        gettimeofday(&tval_after, NULL);

        time_spent(&tval_last, &tval_after, &tval_result);
        printf("\nElapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

        time_spent(&tval_before, &tval_after, &tval_result);
        printf("Total elapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

        printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);
        printf("Current memory: %llu mb\n", ((unsigned long long int)getCurrentRSS())/1024);

        tval_last = tval_after;
    }

    memory_Used* mem = printMemoryUsedFromNode(&(root->node), size_kmer, func_on_types);

    if (skip_node_root != NULL) free_skip_nodes(&(root->node), skip_node_root);

    printMemory(mem);

    free(mem);

    free(line);
    free(array_kmers);

    free(ann_inf);
    free(res);

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
void insert_Genomes_from_FASTxFiles(Root* root, char** filenames, int size_kmer, int id_start_genome, ptrs_on_func* func_on_types){

    ASSERT_NULL_PTR(root,"insert_Genomes_from_FASTxFiles()")
    ASSERT_NULL_PTR(filenames,"insert_Genomes_from_FASTxFiles()")

    struct timeval tval_before, tval_after, tval_last, tval_result;
    gettimeofday(&tval_before, NULL);
    tval_last = tval_before;

    int i = 0;
    int size_buf_tmp = 0; //How many characters are stored in buf_tmp
    int nb_kmers_buf = 0;
    int length_comp_set_colors_tmp = 0;
    int nb_cell_kmer = CEIL(size_kmer*2, SIZE_CELL); //Size of kmers in bytes

    annotation_array_elem* comp_set_colors_tmp = NULL;

    annotation_inform* ann_inf = calloc(1,sizeof(annotation_inform)); //Initialize structure to pass information between reading and modifying an annotation
    ASSERT_NULL_PTR(ann_inf,"insert_Genomes_from_FASTxFiles()")

    resultPresence* res = create_resultPresence();

    Pvoid_t PJArray = (PWord_t)NULL;
    Word_t Rc_word;

    char* buf_tmp = calloc((size_kmer-1)*2, sizeof(char)); //Allocate temporary buffer
    ASSERT_NULL_PTR(buf_tmp,"insert_Genomes_from_FASTxFiles()")

    uint8_t* tab_kmers = calloc(SIZE_BUFFER*nb_cell_kmer, sizeof(uint8_t)); //Allocate buffer for kmers
    ASSERT_NULL_PTR(tab_kmers,"insert_Genomes_from_FASTxFiles()")

    uint64_t kmers_read = 0;
    uint64_t tmp_kmers_read = 0;

    for (i=id_start_genome; i<root->nb_genomes; i++){ //For each file in input
        size_buf_tmp = 0;
        kmers_read = 0;
        tmp_kmers_read = 0;
        nb_kmers_buf = 0;

        int fp = open(filenames[i-id_start_genome], O_RDONLY); //Open it
        kseq_t *seq = kseq_init(fp); //initialize the parser for this file
        int size_seq = kseq_read(seq, -1); //Start reading file, seq contains a buffer with a part of a sequence from the file

        printf("\nFile : %s\n\n", filenames[i-id_start_genome]);

        while (size_seq > -1) { //While the end of the file is not reached

            if (size_seq > 0) size_buf_tmp = 0; //New sequence

            int current_buf_length = seq->seq.l - seq->seq.z; //Number of characters put into the seq buffer

            if (current_buf_length > 0){ //If the seq buffer is not empty

                nb_kmers_buf = MAX(current_buf_length-size_kmer+1, 0); //Number of kmers it is possible to read in seq buffer

                if (size_buf_tmp == 0){ //If the number of characters in the temporary buffer is 0
                    if (nb_kmers_buf != 0){ //If there is at least one kmer in the seq buffer
                        memcpy(buf_tmp, &(seq->seq.s[nb_kmers_buf]), size_kmer-1); //Copy the last size_kmer-1 characters of seq-buffer into buf_tmp
                        size_buf_tmp = (size_kmer-1)/* *2 */;
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
                        insertKmers(root, tab_kmers, size_kmer, nb_kmers, i, func_on_types, ann_inf, res); //Insert the kmers into the tree
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
                    insertKmers(root, tab_kmers, size_kmer, nb_kmers_buf, i, func_on_types, ann_inf, res);
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

        if ((i > 5) && (i%TRESH_DEL_ANNOT == 0)){

            load_annotation_from_Node(&(root->node), size_kmer, func_on_types, &PJArray, root->comp_set_colors);

            comp_set_colors_tmp = root->comp_set_colors;
            length_comp_set_colors_tmp = root->length_comp_set_colors;

            root->comp_set_colors = sort_annotations(&PJArray, &(root->length_comp_set_colors));

            compress_annotation_from_Node(&(root->node), size_kmer, func_on_types, &PJArray, comp_set_colors_tmp);

            free_annotation_array_elem(comp_set_colors_tmp, length_comp_set_colors_tmp);

            JSLFA(Rc_word, PJArray);
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

    memory_Used* mem = printMemoryUsedFromNode(&(root->node), size_kmer, func_on_types);

    printMemory(mem);

    free(mem);

    free(ann_inf);
    free(res);

    free(buf_tmp);
    free(tab_kmers);

    return;
}

int queryBFT_kmerPresences_from_KmerFiles(Root* root, char* query_filename, int binary_file, char* output_filename, int size_kmer){

    ASSERT_NULL_PTR(root,"queryBFT_kmerPresences_from_KmerFiles()")
    ASSERT_NULL_PTR(query_filename,"queryBFT_kmerPresences_from_KmerFiles()")

    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    const char comma = ',';

    int annot_present;
    int size_annot;
    int size_annot_cplx;
    int size_annot_res;

    int i = 0;
    int j = 0;
    int k = 0;
    int nb_kmers_present = 0;
    int nb_bytes_kmer = CEIL(size_kmer*2, SIZE_CELL);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer;

    uint64_t kmers_read = 0;

    FILE* file_query;
    FILE* file_output;

    ptrs_on_func* func_on_types = create_ptrs_on_func(SIZE_SEED, size_kmer);

    resultPresence* res;

    size_t return_fread;

    uint8_t* annot;
    uint8_t* annot_ext;
    uint8_t* annot_cplx;

    uint8_t* annot_res = calloc(CEIL(root->nb_genomes+2, SIZE_CELL), sizeof(uint8_t));
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

                res = isKmerPresent(&(root->node), &(array_kmers[k*nb_bytes_kmer]), size_kmer, func_on_types);

                if (res->link_child != NULL){

                    annot_present = get_annotation((UC*)res->container, &annot, &annot_ext, &annot_cplx, &size_annot,
                                                   &size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);

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

            if (parseKmerCount(line, size_kmer, array_kmers, k) == 1){
                k += nb_bytes_kmer;
                j++;

                if (j == nb_kmer_in_buf){

                    for (i=0; i<nb_kmer_in_buf; i++){

                        res = isKmerPresent(&(root->node), &(array_kmers[i*nb_bytes_kmer]), size_kmer, func_on_types);

                        if (res->link_child != NULL){

                            annot_present = get_annotation((UC*)res->container, &annot, &annot_ext, &annot_cplx, &size_annot,
                                                           &size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);

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

            res = isKmerPresent(&(root->node), &(array_kmers[i*nb_bytes_kmer]), size_kmer, func_on_types);

            if (res->link_child != NULL){

                annot_present = get_annotation((UC*)res->container, &annot, &annot_ext, &annot_cplx, &size_annot,
                                               &size_annot_cplx, res->posFilter2, res->posFilter3, res->pos_sub_bucket);

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

    return nb_kmers_present;
}

int queryBFT_kmerBranching_from_KmerFiles(Root* root, char* query_filename, int binary_file, int size_kmer){

    ASSERT_NULL_PTR(root,"queryBFT_kmerBranching_from_KmerFiles()")
    ASSERT_NULL_PTR(query_filename,"queryBFT_kmerBranching_from_KmerFiles()")

    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    int i = 0;
    int j = 0;
    int k = 0;
    int count_branching_node = 0;
    int nb_bytes_kmer = CEIL(size_kmer*2, SIZE_CELL);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_bytes_kmer;

    uint64_t kmers_read = 0;

    FILE* file;

    ptrs_on_func* func_on_types = create_ptrs_on_func(SIZE_SEED, size_kmer);

    size_t return_fread;

    uint8_t* array_kmers = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(array_kmers,"queryBFT_kmerBranching_from_KmerFiles()")

    char* line = calloc(100, sizeof(char));
    ASSERT_NULL_PTR(line,"queryBFT_kmerBranching_from_KmerFiles()")

    uint16_t** skip_node_root = build_skip_nodes(&(root->node), size_kmer, func_on_types);

    file = fopen(query_filename, "r");
    ASSERT_NULL_PTR(file,"queryBFT_kmerBranching_from_KmerFiles()")

    printf("\nQuerying BFT for branching k-mers in %s\n\n", query_filename);

    if (binary_file){

        if (fgets(line, 100, file) == NULL) ERROR("Cannot read header of the file")
        if (fgets(line, 100, file) == NULL) ERROR("Cannot read header of the file")

        while ((!ferror(file)) && (!feof(file))){

            return_fread = fread(array_kmers, (size_t)nb_bytes_kmer, (size_t)nb_kmer_in_buf, file);

            for (k=0; k<(int)return_fread; k++){

                if (isBranchingRight(&(root->node), &(array_kmers[k*nb_bytes_kmer]), size_kmer, func_on_types, skip_node_root) > 1){
                    count_branching_node++;
                }
                else if (isBranchingLeft(&(root->node), &(array_kmers[k*nb_bytes_kmer]), size_kmer, func_on_types, skip_node_root) > 1){
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

            if (parseKmerCount(line, size_kmer, array_kmers, k) == 1){
                k += nb_bytes_kmer;
                j++;

                if (j == nb_kmer_in_buf){

                    for (i=0; i<nb_kmer_in_buf; i++){

                        if (isBranchingRight(&(root->node), &(array_kmers[i*nb_bytes_kmer]), size_kmer, func_on_types, skip_node_root) > 1){
                            count_branching_node++;
                        }
                        else if (isBranchingLeft(&(root->node), &(array_kmers[i*nb_bytes_kmer]), size_kmer, func_on_types, skip_node_root) > 1){
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

            if (isBranchingRight(&(root->node), &(array_kmers[i*nb_bytes_kmer]), size_kmer, func_on_types, skip_node_root) > 1){
                count_branching_node++;
            }
            else if (isBranchingLeft(&(root->node), &(array_kmers[i*nb_bytes_kmer]), size_kmer, func_on_types, skip_node_root) > 1){
                count_branching_node++;
            }
        }

        memset(array_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));
    }

    fclose(file);

    free(array_kmers);
    free(line);

    if (skip_node_root != NULL) free_skip_nodes(&(root->node), skip_node_root);

    gettimeofday(&tval_after, NULL);
    time_spent(&tval_before, &tval_after, &tval_result);
    printf("\nElapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

    return count_branching_node;
}

/* ---------------------------------------------------------------------------------------------------------------
*  Functions get_nb_cplx_nodes_from_* are currently in beta test and should not be used.
*  ---------------------------------------------------------------------------------------------------------------
*/

int get_nb_cplx_nodes_from_KmerCounting(Root* tree, char* name_file, int size_kmer, uint16_t** skip_node_root, ptrs_on_func* func_on_types, annotation_inform* ann_inf){

    ASSERT_NULL_PTR(tree,"get_nb_cplx_nodes_from_KmerCounting()")
    ASSERT_NULL_PTR(name_file,"get_nb_cplx_nodes_from_KmerCounting()")
    ASSERT_NULL_PTR(skip_node_root,"get_nb_cplx_nodes_from_KmerCounting()")
    ASSERT_NULL_PTR(func_on_types,"get_nb_cplx_nodes_from_KmerCounting()")
    ASSERT_NULL_PTR(ann_inf,"get_nb_cplx_nodes_from_KmerCounting()")

    resultPresence* res;

    UC* uc;

    uint8_t* tab_kmers = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(tab_kmers,"get_nb_cplx_nodes_from_KmerCounting()")

    char* line = calloc(100, sizeof(char));
    ASSERT_NULL_PTR(line,"get_nb_cplx_nodes_from_KmerCounting()")

    int j = 0;
    int k = 0;
    int kmers_read = 0;
    int count_branching = 0;
    int count_without_left_n = 0;
    int degree_left = 0;

    int nb_cell = CEIL(size_kmer*2, SIZE_CELL);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_cell;

    FILE* file = fopen(name_file, "r");
    ASSERT_NULL_PTR(file,"get_nb_cplx_nodes_from_KmerCounting()")

    while (fgets(line, 100, file) != NULL){

        if (parseKmerCount(line, size_kmer, tab_kmers, k) == 1){
            k += nb_cell;
            j++;

            if (j == nb_kmer_in_buf){

                for (k = 0; k < nb_kmer_in_buf * nb_cell; k += nb_cell){

                    if (isBranchingRight(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types, skip_node_root) > 1){ //Flag 0

                        count_branching++;
                        res = isKmerPresent(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types);

                        uc = (UC*)res->container;

                        if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 2);
                        else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 2);

                        free(res);
                    }
                    else if ((degree_left = isBranchingLeft(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types, skip_node_root)) != 1){

                        res = isKmerPresent(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types);

                        uc = (UC*)res->container;

                        if (degree_left == 0){

                            count_without_left_n++;

                            if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 3);
                            else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 3);
                        }
                        else{

                            count_branching++;

                            if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 2);
                            else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 2);
                        }

                        free(res);
                    }
                }

                j = 0;
                k = 0;
                memset(tab_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));

                if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+nb_kmer_in_buf)%PRINT_EVERY_X_KMERS)){
                    //break;
                }

                kmers_read += nb_kmer_in_buf;
            }
        }
    }

    for (k = 0; k < j * nb_cell; k += nb_cell){

        if (isBranchingRight(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types, skip_node_root) > 1){ //Flag 0

            count_branching++;
            res = isKmerPresent(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types);

            uc = (UC*)res->container;

            if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 2);
            else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 2);

            free(res);
        }
        else if ((degree_left = isBranchingLeft(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types, skip_node_root)) != 1){

            res = isKmerPresent(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types);
            uc = (UC*)res->container;

            if (degree_left == 0){

                count_without_left_n++;

                if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 3);
                else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 3);
            }
            else{

                count_branching++;

                if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 2);
                else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 2);
            }

            free(res);
        }
    }

    fclose(file);

    printf("\nNumber of branching nodes = %d\n", count_branching);
    printf("Number of nodes with 0 left neighbor and 1 right neighbor = %d\n", count_without_left_n);
    printf("Total = %d\n", count_branching + count_without_left_n);

    free(line);
    free(tab_kmers);

    return count_branching /*+ count_without_left_n*/;
}

Root* get_nb_cplx_nodes_from_FASTx(Root* tree, int size_kmer, uint16_t** skip_node_root, ptrs_on_func* func_on_types,
                                   annotation_inform* ann_inf, resultPresence* res){

    ASSERT_NULL_PTR(tree,"get_nb_cplx_nodes_from_FASTx()")
    ASSERT_NULL_PTR(tree->filenames,"get_nb_cplx_nodes_from_FASTx()")
    ASSERT_NULL_PTR(skip_node_root,"get_nb_cplx_nodes_from_FASTx()")
    ASSERT_NULL_PTR(func_on_types,"get_nb_cplx_nodes_from_FASTx()")
    ASSERT_NULL_PTR(ann_inf,"get_nb_cplx_nodes_from_FASTx()")
    ASSERT_NULL_PTR(res,"get_nb_cplx_nodes_from_FASTx()")

    Root* bft_cplx_nodes = createRoot(NULL, 0, size_kmer);

    int i = 0, j = 0;
    int size_buf_tmp = 0; //How many characters are stored in buf_tmp
    int nb_kmers_buf = 0;
    int nb_cell_kmer = CEIL(size_kmer*2, SIZE_CELL); //Size of kmers in bytes

    char* buf_tmp = calloc((size_kmer-1)*2, sizeof(char)); //Allocate temporary buffer
    ASSERT_NULL_PTR(buf_tmp,"get_nb_cplx_nodes_from_FASTx()")

    uint8_t* tab_kmers = calloc(SIZE_BUFFER*nb_cell_kmer, sizeof(uint8_t)); //Allocate buffer for kmers
    ASSERT_NULL_PTR(tab_kmers,"get_nb_cplx_nodes_from_FASTx()")

    uint64_t kmers_read = 0;
    uint64_t tmp_kmers_read = 0;

    for (i=0; i<tree->nb_genomes; i++){ //For each file in input
        size_buf_tmp = 0;
        kmers_read = 0;
        tmp_kmers_read = 0;
        nb_kmers_buf = 0;

        int fp = open(tree->filenames[i], O_RDONLY); //Open it
        kseq_t *seq = kseq_init(fp); //initialize the parser for this file
        int size_seq = kseq_read(seq, -1); //Start reading file, seq contains a buffer with a part of a sequence from the file

        while (size_seq > -1) { //While the end of the file is not reached

            if (size_seq > 0) size_buf_tmp = 0; //New sequence

            int current_buf_length = seq->seq.l - seq->seq.z; //Number of characters put into the seq buffer

            if (current_buf_length > 0){ //If the seq buffer is not empty

                nb_kmers_buf = MAX(current_buf_length-size_kmer+1, 0); //Number of kmers it is possible to read in seq buffer

                if (size_buf_tmp == 0){ //If the number of characters in the temporary buffer is 0
                    if (nb_kmers_buf != 0){ //If there is at least one kmer in the seq buffer
                        memcpy(buf_tmp, &(seq->seq.s[nb_kmers_buf]), size_kmer-1); //Copy the last size_kmer-1 characters of seq-buffer into buf_tmp
                        size_buf_tmp = size_kmer-1;
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

                        for (j=0; j < nb_kmers * nb_cell_kmer; j += nb_cell_kmer){
                            if ((isBranchingRight(&(tree->node), &(tab_kmers[j]), size_kmer, func_on_types, skip_node_root) > 1) ||
                                (isBranchingLeft(&(tree->node), &(tab_kmers[j]), size_kmer, func_on_types, skip_node_root) != 1))
                                   insertKmers(bft_cplx_nodes, &(tab_kmers[j]), size_kmer, 1, 0, func_on_types, ann_inf, res);
                        }

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

                    for (j=0; j< nb_kmers_buf * nb_cell_kmer; j += nb_cell_kmer){
                        if ((isBranchingRight(&(tree->node), &(tab_kmers[j]), size_kmer, func_on_types, skip_node_root) > 1) ||
                            (isBranchingLeft(&(tree->node), &(tab_kmers[j]), size_kmer, func_on_types, skip_node_root) != 1))
                                insertKmers(bft_cplx_nodes, &(tab_kmers[j]), size_kmer, 1, 0, func_on_types, ann_inf, res);
                    }

                    memset(tab_kmers, 0, nb_kmers_buf*nb_cell_kmer*sizeof(uint8_t));
                    tmp_kmers_read += nb_kmers_buf;
                }

                //Display how many kmers were read
                if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+tmp_kmers_read)%PRINT_EVERY_X_KMERS)){
                    printf("%" PRIu64 " kmers read\n", kmers_read+tmp_kmers_read);
                    //break;
                }

                kmers_read += tmp_kmers_read;
            }

            size_seq = kseq_read(seq, size_seq);
        }

        kseq_destroy(seq);
        close(fp);
    }

    free(func_on_types);
    free(buf_tmp);
    free(tab_kmers);

    return bft_cplx_nodes;
}
