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

    double threshold;

    bool canonical;

    const char* csv_ext = ".csv";

    char buffer[2048];

    char* dot;
    char* filename_output;

    char** paths_and_names = NULL;

    FILE* file_input = NULL;
    FILE* file_tmp = NULL;

    if (argc == 1){

        ERROR("\nUsage:\n"
              "bft build k {kmers|kmers_comp} list_genome_files output_file [Options]\n"
              "bft load file_bft [-add_genomes {kmers|kmers_comp} list_genome_files output_file] [Options]\n"
              "\nOptions:\n"
              "[-query_sequences threshold {canonical|non_canonical} list_sequence_files]\n"
              "[-query_kmers {kmers|kmers_comp} list_kmer_files]\n"
              "[-query_branching {kmers|kmers_comp} list_kmer_files]\n"
              "[-extract_kmers {kmers|kmers_comp} compressed_kmers_file]\n\n")
    }
    else {

        if ((strcmp("--version", argv[1]) == 0) || (strcmp("-v", argv[1]) == 0)){

            fprintf(stderr, "0.8\n");
            exit(0);
        }
        else if (strcmp("build", argv[1]) == 0){

            size_kmer = atoi(argv[2]); //Length k argument reading

            //Test if k is valid
            if (size_kmer <= 0) ERROR("Provided length k (for k-mers) is either <= 0 or not a number.\n")
            else if (size_kmer > KMER_LENGTH_MAX) ERROR("Length k (for k-mers) cannot be superior to 63.\n")
            else if (size_kmer%NB_CHAR_SUF_PREF != 0) ERROR("Length k (for k-mers) must be a multiple of 9.\n")

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

            i = 6;
        }
        else if (strcmp("load", argv[1]) == 0){

            if ((file_input = fopen(argv[2], "rb")) == NULL) ERROR("Invalid file_bft.\n")
            fclose(file_input);

            i = 3;
        }

        for (; i < argc; i += 3){ //Test if we can open the files for querying the k-mers/branching vertices

            if ((strcmp("-query_kmers", argv[i]) == 0) || (strcmp("-query_branching", argv[i]) == 0) || (strcmp("-extract_kmers", argv[i]) == 0) ||
                ((strcmp("load", argv[1]) == 0) && (strcmp("-add_genomes", argv[i]) == 0))){

                if ((strcmp("kmers_comp", argv[i+1]) != 0) && (strcmp("kmers", argv[i+1]) != 0)){

                    fprintf(stderr, "Unrecognized type of input files for %s.\n"
                            "Choice must be 'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n", argv[i]);
                    exit(EXIT_FAILURE);
                }
            }
            else if (strcmp("-query_sequences", argv[i]) == 0){

                threshold = atof(argv[i+1]);
                if (threshold == 0) ERROR("Could not parse threshold for command -query_sequences.\n");

                if ((strcmp("canonical", argv[i+2]) != 0) && (strcmp("non_canonical", argv[i+2]) != 0)){

                    fprintf(stderr, "Unrecognized type of k-mers to search for %s.\n"
                            "Choice must be 'canonical' or 'non_canonical' to consider canonical or non-canonical k-mers of the queries.\n", argv[i]);
                    exit(EXIT_FAILURE);
                }
            }
            else {

                fprintf(stderr, "Unrecognized command %s.\n", argv[i]);
                exit(EXIT_FAILURE);
            }

            if (strcmp("-extract_kmers", argv[i])){

                if ((file_input = fopen(argv[i + (strcmp("-query_sequences", argv[i]) == 0 ? 3 : 2)], "r")) == NULL){

                    fprintf(stderr, "Invalid list of input files for command %s.\n", argv[i]);
                    exit(EXIT_FAILURE);
                }

                cpt = 0;

                while (fgets(buffer, 2048, file_input)){ //Test if the input files can be opened and read

                    buffer[strcspn(buffer, "\r\n")] = '\0';

                    if ((file_tmp = fopen(buffer, "r")) == NULL){

                        fprintf(stderr, "Invalid file at line %d of the list of input files for command %s.\n", cpt, argv[i]);
                        exit(EXIT_FAILURE);
                    }

                    cpt++;

                    fclose(file_tmp);
                }

                fclose(file_input);
            }

            if (((strcmp("load", argv[1]) == 0) && (strcmp("-add_genomes", argv[i]) == 0)) || (strcmp("-query_sequences", argv[i]) == 0)) i++;
        }

        if (strcmp("build", argv[1]) == 0){

            if ((file_input = fopen(argv[4], "r")) == NULL) ERROR("Invalid list_genome_files.\n")

            paths_and_names = malloc(nb_files_2_read * sizeof(char*)); //Allocate the array of paths + filenames
            ASSERT_NULL_PTR(paths_and_names,"main()")

            i = 0;

            while (fgets(buffer, 2048, file_input)){ //Copy the filenames in an array, the paths and filenames in another array

                paths_and_names[i] = malloc((strlen(buffer) + 1) * sizeof(char));
                ASSERT_NULL_PTR(paths_and_names[i], "main()")

                strcpy(paths_and_names[i], buffer);

                paths_and_names[i][strcspn(paths_and_names[i], "\r\n")] = '\0';

                i++;
            }

            fclose(file_input);

            if (strcmp("kmers_comp", argv[3]) == 0){

                root = createBFT_Root(size_kmer, 1, 0);
                insert_Genomes_from_KmerFiles(root, nb_files_2_read, paths_and_names, 1, argv[5]);
            }
            else if (strcmp("kmers", argv[3]) == 0){

                root = createBFT_Root(size_kmer, 1, 0);
                insert_Genomes_from_KmerFiles(root, nb_files_2_read, paths_and_names, 0, argv[5]);
            }
            else {

                ERROR("Unrecognized type of input files.\nChoice must be "
                      "'kmers' for k-mers files or 'kmers_comp' for compressed k-mers files.\n")
            }

            memory_Used* mem = printMemoryUsedFromNode(&(root->node), (root->k / NB_CHAR_SUF_PREF) - 1, root->k, root->info_per_lvl);
            printMemory(mem);
            free(mem);

            write_BFT_Root(root, argv[5], false);

            i = 6;
        }
        else{

            root = read_BFT_Root(argv[2]);

            memory_Used* mem = printMemoryUsedFromNode(&(root->node), (root->k / NB_CHAR_SUF_PREF) - 1, root->k, root->info_per_lvl);
            printMemory(mem);
            free(mem);

            printf("Memory used by external colors = %f\n", getTotalSize_annotation_array_elem(root->comp_set_colors, root->length_comp_set_colors));

            i = 3;
        }

        for (; i < argc; i += 3){

            if ((strcmp("load", argv[1]) == 0) && (strcmp("-add_genomes", argv[i]) == 0)){

                if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid list_genome_files.\n")

                paths_and_names = malloc(nb_files_2_read * sizeof(char*)); //Allocate the array of paths + filenames
                ASSERT_NULL_PTR(paths_and_names,"main()")

                j = 0;

                while (fgets(buffer, 2048, file_input)){ //Copy the filenames in an array, the paths and filenames in another array

                    paths_and_names[j] = malloc((strlen(buffer)+1)*sizeof(char));
                    ASSERT_NULL_PTR(paths_and_names[j],"main()")

                    strcpy(paths_and_names[j], buffer);

                    paths_and_names[j][strcspn(paths_and_names[j], "\r\n")] = '\0';

                    j++;
                }

                fclose(file_input);

                if (strcmp("kmers_comp", argv[i+1]) == 0) insert_Genomes_from_KmerFiles(root, nb_files_2_read, paths_and_names, 1, argv[i+3]);
                else insert_Genomes_from_KmerFiles(root, nb_files_2_read, paths_and_names, 0, argv[i+3]);

                write_BFT_Root(root, argv[i+3], false);

                i++;
            }
            else if (strcmp("-query_kmers", argv[i]) == 0){

                if (strcmp("kmers_comp", argv[i+1]) == 0) binary_files = 1;
                else binary_files = 0;

                if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid k-mer queries files list.\n")

                while (fgets(buffer, 2048, file_input)){

                    buffer[strcspn(buffer, "\r\n")] = 0;

                    filename_output = malloc((strlen(basename(buffer))+4) * sizeof(char));
                    ASSERT_NULL_PTR(filename_output, "main()")

                    strcpy(filename_output, basename(buffer));

                    if ((dot = strrchr(filename_output, '.')) != NULL) strcpy(dot, csv_ext);
                    else strcpy(&(filename_output[strlen(filename_output)]), csv_ext);

                    printf("\nNb k-mers present = %d\n", queryBFT_kmerPresences_from_KmerFiles(root, buffer, binary_files, filename_output));

                    free(filename_output);
                }

                fclose(file_input);
            }
            else if (strcmp("-query_sequences", argv[i]) == 0){

                threshold = atof(argv[i+1]);

                canonical = (strcmp("canonical", argv[i+2]) == 0 ? true : false);

                if ((file_input = fopen(argv[i+3], "r")) == NULL) ERROR("Invalid sequence query file list.\n")

                while (fgets(buffer, 2048, file_input)){

                    buffer[strcspn(buffer, "\r\n")] = '\0';

                    filename_output = malloc((strlen(basename(buffer))+4) * sizeof(char));
                    ASSERT_NULL_PTR(filename_output, "main()")

                    strcpy(filename_output, basename(buffer));

                    if ((dot = strrchr(filename_output, '.')) != NULL) strcpy(dot, csv_ext);
                    else strcpy(&(filename_output[strlen(filename_output)]), csv_ext);

                    query_sequences_outputCSV(root, buffer, filename_output, threshold, canonical);

                    free(filename_output);
                }

                fclose(file_input);

                i++;
            }
            else if (strcmp("-query_branching", argv[i]) == 0){

                if (strcmp("kmers_comp", argv[i+1]) == 0) binary_files = 1;
                else binary_files = 0;

                if ((file_input = fopen(argv[i+2], "r")) == NULL) ERROR("Invalid branching k-mer queries files list.\n")

                while (fgets(buffer, 2048, file_input)){

                    buffer[strcspn(buffer, "\r\n")] = 0;
                    printf("\nNb branching k-mers = %d\n", queryBFT_kmerBranching_from_KmerFiles(root, buffer, binary_files));
                }

                fclose(file_input);
            }
            else write_kmers_2disk(root, argv[i+2], strcmp("kmers_comp", argv[i+1]) == 0);
        }

        if (paths_and_names != NULL){

            for (i = 0; i < root->nb_genomes; i++) free(paths_and_names[i]);
            free(paths_and_names);
        }

        freeBFT_Root(root);
    }

    return EXIT_SUCCESS;
}
