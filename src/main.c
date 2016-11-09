#include <libgen.h>

#include <jemalloc/jemalloc.h>
const char* malloc_conf = "narenas:1,tcache:false,lg_dirty_mult:8,lg_chunk:22";

#include <Judy.h>

#include "insertNode.h"
#include "printMemory.h"
#include "write_to_disk.h"
#include "snippets.h"

int main(int argc, char *argv[])
{
    BFT_Root* root = NULL;

    int i = 0, j = 0, cpt = 0;
    int size_kmer = 27;
    int binary_files = 0;
    int nb_files_2_read = 0;

    const char csv_ext[5] = ".csv\0";

    char buffer[2048];

    char* dot;
    char* filename_output;

    char** paths_and_names = NULL;

    FILE* file_input = NULL;
    FILE* file_tmp = NULL;

    if (argc == 1){
        ERROR("\nUsage:\n"
              "./bft build k treshold_compression {kmers|kmers_comp} list_genome_files output_file [Options]\n"
              "./bft load file_bft [-add_genomes {kmers|kmers_comp} list_genome_files output_file] [Options]\n"
              "\nOptions:\n"
              "[-query_kmers {kmers|kmers_comp} list_kmer_files]\n"
              "[-query_branching {kmers|kmers_comp} list_kmer_files]\n"
              "[-extract_kmers {kmers|kmers_comp} compressed_kmers_file]\n\n")
    }
    else {

        if ((strcmp("--version", argv[1]) == 0) || (strcmp("-v", argv[1]) == 0)){
            printf("0.8\n");
            exit(0);
        }
        else if (strcmp("build", argv[1]) == 0){

            size_kmer = atoi(argv[2]); //Length k argument reading

            //Test if k is valid
            if (size_kmer <= 0) ERROR("Provided length k (for k-mers) is either <= 0 or not a number\n")
            else if (size_kmer > KMER_LENGTH_MAX) ERROR("Length k (for k-mers) cannot be superior to 63\n")
            else if (size_kmer%NB_CHAR_SUF_PREF != 0) ERROR("Length k (for k-mers) must be a multiple of 9\n")

            //Open the file containing the input files
            if ((file_input = fopen(argv[5], "r")) == NULL) ERROR("Invalid list_genome_files.\n")

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

            for (i=7; i<argc; i+=3){ //Test if we can open the files for querying the k-mers/branching vertices

                //User wants to query the BFT for k-mers
                if ((strcmp("-query_kmers", argv[i]) == 0) || (strcmp("-query_branching", argv[i]) == 0) || (strcmp("-extract_kmers", argv[i]) == 0)){

                    if ((strcmp("kmers_comp", argv[i+1]) != 0) && (strcmp("kmers", argv[i+1]) != 0)){
                        if (strcmp("-query_kmers", argv[i]) == 0){
                            ERROR("Unrecognized type of input files for -query_kmers.\n"
                                  "Choice must be 'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")
                        }
                        else if (strcmp("-query_branching", argv[i]) == 0){
                            ERROR("Unrecognized type of input files for -query_branching.\n"
                                  "Choice must be 'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")
                        }
                        else if (strcmp("-extract_kmers", argv[i]) == 0){
                            ERROR("Unrecognized type of output files for -extract_kmers.\n"
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

            if ((file_input = fopen(argv[5], "r")) == NULL) ERROR("Invalid list_genome_files.\n")

            paths_and_names = malloc(nb_files_2_read * sizeof(char*)); //Allocate the array of paths + filenames
            ASSERT_NULL_PTR(paths_and_names,"main()")

            i = 0;

            while (fgets(buffer, 2048, file_input)){ //Copy the filenames in an array, the paths and filenames in another array

                paths_and_names[i] = malloc((strlen(buffer)+1)*sizeof(char));
                ASSERT_NULL_PTR(paths_and_names[i],"main()")

                strcpy(paths_and_names[i], buffer);

                paths_and_names[i][strcspn(paths_and_names[i], "\r\n")] = 0;

                i++;
            }

            fclose(file_input);

            //Read and test if the type of input files is valid
            //Insert k-mers of the input files in the BFT
            if (strcmp("kmers_comp", argv[4]) == 0){

                root = createBFT_Root(size_kmer, atoi(argv[3]), 0);
                insert_Genomes_from_KmerFiles(root, nb_files_2_read, paths_and_names, 1, argv[6]);
            }
            else if (strcmp("kmers", argv[4]) == 0){

                root = createBFT_Root(size_kmer, atoi(argv[3]), 0);
                insert_Genomes_from_KmerFiles(root, nb_files_2_read, paths_and_names, 0, argv[6]);
            }
            else
                ERROR("Unrecognized type of input files.\nChoice must be 'fastx' for FASTA/FASTQ files, "
                      "'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")

            write_BFT_Root(root, argv[6], false);

            memory_Used* mem = printMemoryUsedFromNode(&(root->node), (root->k / NB_CHAR_SUF_PREF) - 1, root->k, root->info_per_lvl);
            printMemory(mem);
            free(mem);

            //printf("Memory used by external colors = %f\n", getTotalSize_annotation_array_elem(root->comp_set_colors, root->length_comp_set_colors));

            for (i=7; i<argc; i+=3){

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
                        printf("\nNb k-mers present = %d\n", queryBFT_kmerPresences_from_KmerFiles(root, buffer, binary_files, filename_output));

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
                        printf("\nNb branching k-mers = %d\n", queryBFT_kmerBranching_from_KmerFiles(root, buffer, binary_files));
                    }

                    fclose(file_input);
                }
                else if (strcmp("-extract_kmers", argv[i]) == 0){

                    write_kmers_2disk(root, argv[i+2], strcmp("kmers_comp", argv[i+1]) == 0);
                }
                else{
                    fprintf(stderr, "Unrecognized command %s.\n", argv[i]);
                    exit(EXIT_FAILURE);
                }
            }

            if (paths_and_names != NULL){
                for (int i=0; i<root->nb_genomes; i++) free(paths_and_names[i]);
                free(paths_and_names);
            }

            freeBFT_Root(root);
        }
        else if (strcmp("load", argv[1]) == 0){

            root = read_BFT_Root(argv[2]);

            memory_Used* mem = printMemoryUsedFromNode(&(root->node), (root->k / NB_CHAR_SUF_PREF) - 1, root->k, root->info_per_lvl);
            printMemory(mem);
            free(mem);

            printf("Memory used by external colors = %f\n", getTotalSize_annotation_array_elem(root->comp_set_colors, root->length_comp_set_colors));

            for (i=3; i<argc; i+=3){ //Test if we can open the files for querying the k-mers/branching vertices

                if ((strcmp("-query_kmers", argv[i]) == 0) || (strcmp("-query_branching", argv[i]) == 0) || (strcmp("-extract_kmers", argv[i]) == 0)){ //User wants to query the BFT for k-mers

                    if ((strcmp("kmers_comp", argv[i+1]) != 0) && (strcmp("kmers", argv[i+1]) != 0)){
                        if (strcmp("-query_kmers", argv[i]) == 0){
                            ERROR("Unrecognized type of input files for -query_kmers.\n"
                                  "Choice must be 'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")

                            if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid k-mer queries file.\n")
                        }
                        else if (strcmp("-query_branching", argv[i]) == 0){
                            ERROR("Unrecognized type of input files for -query_branching.\n"
                                  "Choice must be 'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")

                            if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid k-mer queries file.\n")
                        }
                        else if (strcmp("-extract_kmers", argv[i]) == 0){
                            ERROR("Unrecognized type of output files for -extract_kmers.\n"
                                  "Choice must be 'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")
                        }
                    }


                    if (file_input != NULL){
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
                        printf("\nNb k-mers present = %d\n", queryBFT_kmerPresences_from_KmerFiles(root, buffer, binary_files, filename_output));

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
                        printf("\nNb branching k-mers = %d\n", queryBFT_kmerBranching_from_KmerFiles(root, buffer, binary_files));
                    }

                    fclose(file_input);
                }
                else if (strcmp("-add_genomes", argv[i]) == 0){

                    if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid list_genome_files.\n")

                    paths_and_names = malloc(nb_files_2_read*sizeof(char*)); //Allocate the array of paths + filenames
                    ASSERT_NULL_PTR(paths_and_names,"main()")

                    j = 0;

                    while (fgets(buffer, 2048, file_input)){ //Copy the filenames in an array, the paths and filenames in another array

                        paths_and_names[j] = malloc((strlen(buffer)+1)*sizeof(char));
                        ASSERT_NULL_PTR(paths_and_names[j],"main()")

                        strcpy(paths_and_names[j], buffer);

                        paths_and_names[j][strcspn(paths_and_names[j], "\r\n")] = 0;

                        j++;
                    }

                    fclose(file_input);

                    //Read and test if the type of input files is valid
                    //Insert k-mers of the input files in the BFT
                    if (strcmp("kmers_comp", argv[i+1]) == 0){
                        insert_Genomes_from_KmerFiles(root, nb_files_2_read, paths_and_names, 1, argv[i+3]);
                    }
                    else if (strcmp("kmers", argv[i+1]) == 0){
                        insert_Genomes_from_KmerFiles(root, nb_files_2_read, paths_and_names, 0, argv[i+3]);
                    }
                    /*else if (strcmp("fastx", argv[i+1]) == 0){
                        insert_Genomes_from_FASTxFiles(root, nb_files_2_read, paths_and_names);
                    }*/
                    else
                        ERROR("Unrecognized type of input files.\nChoice must be 'fastx' for FASTA/FASTQ files, "
                              "'kmers' for k-mers files or 'kmers_comp' for compressed k-mer files.\n")

                    write_BFT_Root(root, argv[i+3], false);

                    i++;
                }
                else if (strcmp("-extract_kmers", argv[i]) == 0){
                    write_kmers_2disk(root, argv[i+2], strcmp("kmers_comp", argv[i+1]) == 0);
                }
                else{
                    fprintf(stderr, "Unrecognized command %s.\n", argv[i]);
                    exit(EXIT_FAILURE);
                }
            }

            if (paths_and_names != NULL){
                for (int i=0; i<nb_files_2_read; i++) free(paths_and_names[i]);
                free(paths_and_names);
            }

            freeBFT_Root(root);
        }
        else{
            fprintf(stderr, "Unrecognized command %s.\n", argv[1]);
            exit(EXIT_FAILURE);
        }
    }

    return EXIT_SUCCESS;
}
