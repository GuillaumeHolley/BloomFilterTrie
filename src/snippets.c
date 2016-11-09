#include "snippets.h"

/** Function of type BFT_func_ptr extracting a core k-mer to disk.
*   A core k-mer contains in its annotation all genome ids inserted in the graph.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param args contains all additional parameters given to extract_pangenome_kmers_to_disk():
*   A pointer to a file where to write the k-mer and a pointer to the current number of k-mers written.
*/
size_t extract_core_kmers(BFT_kmer* kmer, BFT* graph, va_list args){

    //Extract from the variable arguments list args a pointer to a file where to extract core k-mers
    FILE* file = va_arg(args, FILE*);
    //Extract from the variable arguments list args a pointer to the current number of core k-mers extracted
    int* nb_core_kmers = va_arg(args, int*);

    BFT_annotation* kmer_annot = get_annotation(kmer); //Get the k-mer annotation

    //If the k-mer annotation has all the genome ids inserted, k-mer is cpre
    if (get_count_id_genomes(kmer_annot, graph) == graph->nb_genomes){
        fwrite(kmer->kmer, sizeof(char), strlen(kmer->kmer)+1, file); //Write the k-mer in the file
        *nb_core_kmers += 1; //Increase the current number of core k-mers extracted
    }

    return true;
}

/** Function of type BFT_func_ptr extracting a dispensable k-mer to disk.
*   A dispensable k-mer contains in its annotation less than all genome ids inserted in the graph.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param args contains all additional parameters given to extract_pangenome_kmers_to_disk():
*   A pointer to a file where to write the k-mer and a pointer to the current number of k-mers written.
*/
size_t extract_dispensable_kmers(BFT_kmer* kmer, BFT* graph, va_list args){

    //Extract from the variable arguments list args a pointer to a file where to dispensable k-mers
    FILE* file = va_arg(args, FILE*);
    //Extract from the variable arguments list args a pointer to the current number of dispensable k-mers extracted
    int* nb_disp_kmers = va_arg(args, int*);

    BFT_annotation* kmer_annot = get_annotation(kmer); //Get the k-mer annotation

    //If the k-mer annotation does not contain the total number of genome ids inserted, k-mer is dispensable
    if (get_count_id_genomes(kmer_annot, graph) < graph->nb_genomes){
        fwrite(kmer->kmer, sizeof(char), strlen(kmer->kmer)+1, file); //Write the k-mer in the file
        *nb_disp_kmers += 1; //Increase the current number of dispensable k-mers extracted
    }

    return true;
}

/** Function of type BFT_func_ptr extracting a singleton k-mer to disk.
*   Singleton k-mer contains in its annotation one genome id inserted in the graph.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param args contains all additional parameters given to extract_pangenome_kmers_to_disk():
*   A pointer to a file where to write the k-mer and a pointer to the current number of k-mers written.
*/
size_t extract_singleton_kmers(BFT_kmer* kmer, BFT* graph, va_list args){

    //Extract from the variable arguments list args a pointer to a file where to singleton k-mers
    FILE* file = va_arg(args, FILE*);
    //Extract from the variable arguments list args a pointer to the current number of singleton k-mers extracted
    int* nb_sing_kmers = va_arg(args, int*);

    BFT_annotation* kmer_annot = get_annotation(kmer); //Get the k-mer annotation

    if (get_count_id_genomes(kmer_annot, graph) == 1){ //If the k-mer annotation contains one genome id
        fwrite(kmer->kmer, sizeof(char), strlen(kmer->kmer)+1, file); //Write the k-mer in the file
        *nb_sing_kmers += 1; //Increase the current number of singleton k-mers extracted
    }

    return true;
}

/** Function extracting pan-genome (core/dispensable/singleton) k-mers from a BFT to disk.
*   A core k-mer contains in its annotation all genome ids inserted in the graph.
*   A dispensable k-mer contains in its annotation less than all genome ids inserted in the graph.
*   A singleton k-mer contains in its annotation one genome id inserted in the graph.
*   @param graph is a BFT from which k-mers must be extracted.
*   @param filename_output is the name of the file where to write the k-mers.
*   @param f is a function of type BFT_func_ptr (like extract_core_kmers(), extract_dispensable_kmers(),
*   and extract_singleton_kmers()) which write a pan-genome k-mer to disk.
*/
void extract_pangenome_kmers_to_disk(BFT* graph, char* filename_output, BFT_func_ptr f){

    ASSERT_NULL_PTR(graph, "extract_pangenome_kmers_to_disk()\n")
    ASSERT_NULL_PTR(filename_output, "extract_pangenome_kmers_to_disk()\n")

    FILE* file;
    int nb_kmers = 0;

    if ((file = fopen(filename_output, "w")) == NULL) //Create a file where to write k-mers
        ERROR("extract_pangenome_kmers_to_disk(): failed to create/open output file.\n")

    //Iterate over the k-mers of the graph. f() will be called on each k-mer
    //with file and &nb_kmers as additional arguments
    iterate_over_kmers(graph, f, file, &nb_kmers);

    fclose(file); //Close the file

    printf("Number of extracted k-mers is %d.\n", nb_kmers);

    return;
}

/** Function of type BFT_func_ptr extracting from a non-branching k-mer the simple (non branching)
*   path in which it is.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param args contains all additional parameters given to extract_simple_paths_to_disk():
*   A pointer to a file where to write paths and a pointer to the current max. size of a path written in the file.
*/
size_t extract_simple_paths(BFT_kmer* kmer, BFT* graph, va_list args){

    int l_simple_path = 0;
    int nb_pred = 0;
    int nb_succ = 0;
    int i = 0;
    int j = 0;

    int last;

    int* longest_simple_path;

    BFT_kmer* succ = NULL;
    BFT_kmer* pred = NULL;
    BFT_kmer* neigh = NULL;

    FILE* file;

    if (get_flag_kmer(kmer, graph) != V_VISITED){ //If the k-mer was not visited before

        set_flag_kmer(V_VISITED, kmer, graph); //Set the k-mer to visited

        succ = get_successors(kmer, graph); //Get the possible successors of the k-mer
        pred = get_predecessors(kmer, graph); //Get the possible predecessors of the k-mer

        for (i = 0; i < 4; i++){ //Compute how many predecessors and successors the k-mer really has
            nb_succ += is_kmer_in_cdbg(&succ[i]);
            nb_pred += is_kmer_in_cdbg(&pred[i]);
        }

        if ((nb_succ < 2) && (nb_pred < 2)){ //If the k-mer has edge in-degree < 2 and edge out-degree < 2

            //We define this k-mer as the extension seed

            //Initiate the ASCII simple path with the extension seed
            l_simple_path = graph->k;

            char* simple_path = malloc(l_simple_path * sizeof(char));
            ASSERT_NULL_PTR(simple_path, "extract_simple_paths()\n");

            memcpy(simple_path, kmer->kmer, l_simple_path * sizeof(char));

            //Test which successor of the extension seed is present
            for (i = 0; i < 4; i++)
                if (is_kmer_in_cdbg(&succ[i])) break;

            // If the extension seed has a single non visited successor
            if ((nb_succ == 1) && (get_flag_kmer(&succ[i], graph) != V_VISITED)){

                neigh = get_successors(&succ[i], graph); //Get possible successors of this (current) successor

                //Count how many are present in the graph
                for (j = 0, nb_succ = 0; j < 4; j++) nb_succ += is_kmer_in_cdbg(&neigh[j]);

                free_BFT_kmer(neigh, 4); //Free them

                //If the current successor has zero or one non visited successor
                while ((nb_succ <= 1) && (get_flag_kmer(&succ[i], graph) != V_VISITED)){

                    set_flag_kmer(V_VISITED, &succ[i], graph); //Set the current successor to visited

                    free_BFT_kmer(pred, 4); //Free the possible predecessors
                    pred = get_predecessors(&succ[i], graph); //Get the possible predecessors of the current successor

                    //Count the number of predecessors the current successor really has
                    for (j=0, nb_pred = 0; j < 4; j++) nb_pred += is_kmer_in_cdbg(&pred[j]);

                    //If the successor has exactly one predecessor, at this point, it means the current successor has
                    //edge in-degree == 1 and edge out-degree < 2
                    if (nb_pred == 1){

                        l_simple_path++; //Increase length of the current simple path

                        //Reallocate the simple path to match new length
                        simple_path = realloc(simple_path, l_simple_path * sizeof(char));
                        ASSERT_NULL_PTR(simple_path, "extract_simple_paths()\n");

                        //Copy last character of the current successor to the end of the simple path
                        simple_path[l_simple_path - 1] = succ[i].kmer[graph->k - 1];

                        //Get the possible successors of the current successor, replace current successor with them
                        neigh = get_successors(&succ[i], graph);
                        free_BFT_kmer(succ, 4);
                        succ = neigh;

                        //Count the number of successors present in the graph and which one is the last one
                        for (i = 0, nb_succ = 0; i < 4; i++){
                            if (is_kmer_in_cdbg(&succ[i])){
                                nb_succ++;
                                last = i;
                            }
                        }

                        i = last;
                    }
                    else break;
                }

                free_BFT_kmer(pred, 4); //Free the current predecessors
                pred = get_predecessors(kmer, graph); //Get the possible predecessors of the extension seed

                //Count how many predecessors of the extension seed are present in the graph
                for (i = 0, nb_pred = 0; i < 4; i++){
                    if (is_kmer_in_cdbg(&pred[i])){
                        nb_pred++;
                        break;
                    }
                }
            }

            // If the extension seed has a single non visited predecessor
            if ((nb_pred == 1) && (get_flag_kmer(&pred[i], graph) != V_VISITED)){

                neigh = get_predecessors(&pred[i], graph); //Get possible predecessors of the (current) predecessor

                //Count how many are present in the graph
                for (j = 0, nb_pred = 0; j < 4; j++) nb_pred += is_kmer_in_cdbg(&neigh[j]);

                free_BFT_kmer(neigh, 4); //Free them

                //If the current predecessor has zero or one non visited predecessor
                while ((nb_pred <= 1) && (get_flag_kmer(&pred[i], graph) != V_VISITED)){

                    set_flag_kmer(V_VISITED, &pred[i], graph); //Set the current predecessor to visited

                    free_BFT_kmer(succ, 4); //Free successors
                    succ = get_successors(&pred[i], graph); //Get the possible successors of the current predecessor

                    //Count the number of successors of the current predecessor present in the graph
                    for (j=0, nb_succ = 0; j < 4; j++) nb_succ += is_kmer_in_cdbg(&succ[j]);

                    //If the current predecessor has exactly one successor, at this point, it means it has
                    //edge out-degree == 1 and edge in-degree < 2
                    if (nb_succ == 1){

                        l_simple_path++; //Increase length of the current simple path

                        //Reallocate the simple path to match new length
                        simple_path = realloc(simple_path, l_simple_path * sizeof(char));
                        ASSERT_NULL_PTR(simple_path, "extract_simple_paths()\n");

                        //Shift the current path of one character on the right
                        memmove(&simple_path[1], simple_path, (l_simple_path - 1) * sizeof(char));

                        //First character of the simple path is the first character of the current predecessor
                        simple_path[0] = pred[i].kmer[0];

                        //Get the possible predecessors of the current predecessor, they become the current predecessors
                        neigh = get_predecessors(&pred[i], graph);
                        free_BFT_kmer(pred, 4);
                        pred = neigh;

                        //Count the number of predecessors present in the graph and which one is the last one
                        for (i=0, nb_pred = 0; i < 4; i++){
                            if (is_kmer_in_cdbg(&pred[i])){
                                nb_pred++;
                                last = i;
                            }
                        }

                        i = last;
                    }
                    else break;
                }
            }

            //Prepare the simple path for file writing
            simple_path = realloc(simple_path, (l_simple_path+1) * sizeof(char));
            ASSERT_NULL_PTR(simple_path, "extract_simple_paths()\n");

            simple_path[l_simple_path] = '\n';

            file = va_arg(args, FILE*); //Extract the file where to write the simple path

            fwrite(simple_path, sizeof(char), l_simple_path+1, file); //Write the simple path in the file

            free(simple_path); //Free the simple path
        }

        //Extract a pointer to the integer storing the max length of extracted simple path
        longest_simple_path = va_arg(args, int*);
        *longest_simple_path = MAX(*longest_simple_path, l_simple_path); //Update this pointer

        //Free all predecessors and successors still allocated
        free_BFT_kmer(succ, 4);
        free_BFT_kmer(pred, 4);
    }

    return 1;
}

/** Function extracting from a colored de Bruijn graph stored as a BFT all simple (non branching) paths.
*   @param graph is a colored de Bruijn graph stored as a BFT.
*   @param filename_output is the name of the file where to write the simple paths.
*/
void extract_simple_paths_to_disk(BFT* graph, char* filename_output){

    ASSERT_NULL_PTR(graph, "extract_simple_paths_to_disk()\n")

    FILE* file;
    int longest_simple_path = 0; //Initiate the number of char. in the longest simple path to 0

    if ((file = fopen(filename_output, "w")) == NULL) //Create a file where to write simple paths
        ERROR("extract_simple_paths_to_disk(): failed to create output file.\n")

    set_neighbors_traversal(graph); //Prepare the graph for traversal
    set_marking(graph); //Prepare the graph to mark the vertices (as visited)

    //Iterate over the k-mers of the graph. extract_simple_paths() will be called on each k-mer
    //with file and  &longest_simple_path as additional arguments
    iterate_over_kmers(graph, extract_simple_paths, file, &longest_simple_path);

    unset_marking(graph); //Delete all marking on vertices, unlock the graph
    unset_neighbors_traversal(graph); //Unset the graph for traversal

    fclose(file); //Close the file

    printf("Longest simple path has %d nuc.\n", longest_simple_path);

    return;
}

/** Function of type BFT_func_ptr extracting from a non-branching core k-mer the simple (non branching)
*   path in which it is. args contains as argument a core_ratio float (between 0 to 1) indicating
*   the ratio of genome ids (compared to the total number of genome ids) that a k-mer annotation must
*   contain to be considered being part of a core simple path.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param args contains all additional parameters given to extract_simple_core_paths_to_disk():
*   A pointer to a file where to write paths and a pointer to the current max. size of a path written in the file.
*/
size_t extract_core_simple_paths(BFT_kmer* kmer, BFT* graph, va_list args){

    int nb_pred = 0;
    int nb_succ = 0;
    int i = 0;
    int j = 0;

    int last;

    BFT_kmer* succ = NULL;
    BFT_kmer* pred = NULL;
    BFT_kmer* neigh = NULL;

    BFT_annotation* curr_annot = NULL;
    BFT_annotation* neigh_annot = NULL;
    BFT_annotation* tmp_annot = NULL;

    FILE* file = va_arg(args, FILE*);
    int* longest_simple_path = va_arg(args, int*);
    double core_ratio = va_arg(args, double);
    int nb_genomes_core = core_ratio * graph->nb_genomes;

    int l_simple_path = 0;

    if (get_flag_kmer(kmer, graph) != V_VISITED){ //If the k-mer was not visited before

        set_flag_kmer(V_VISITED, kmer, graph); //Set the k-mer to visited

        curr_annot = get_annotation(kmer);
        if (get_count_id_genomes(curr_annot, graph) < nb_genomes_core) return 1;

        succ = get_successors(kmer, graph); //Get the successors of the k-mer
        pred = get_predecessors(kmer, graph); //Get the predecessors of the k-mer

        for (i = 0; i < 4; i++){ //Compute how many predecessors and successors for the k-mer
            nb_succ += is_kmer_in_cdbg(&succ[i]);
            nb_pred += is_kmer_in_cdbg(&pred[i]);
        }

        if ((nb_succ < 2) && (nb_pred < 2)){

            l_simple_path = graph->k;

            char* simple_path = malloc(l_simple_path * sizeof(char));
            ASSERT_NULL_PTR(simple_path, "extract_core_simple_paths()\n");

            memcpy(simple_path, kmer->kmer, l_simple_path * sizeof(char));

            for (i = 0; i < 4; i++)
                if (is_kmer_in_cdbg(&succ[i])) break;

            if ((nb_succ == 1) && (get_flag_kmer(&succ[i], graph) != V_VISITED)){

                neigh = get_successors(&succ[i], graph);

                for (j = 0, nb_succ = 0; j < 4; j++) nb_succ += is_kmer_in_cdbg(&neigh[j]);

                free_BFT_kmer(neigh, 4);

                while ((nb_succ <= 1) && (get_flag_kmer(&succ[i], graph) != V_VISITED)){

                    free_BFT_kmer(pred, 4);
                    pred = get_predecessors(&succ[i], graph);

                    for (j=0, nb_pred = 0; j < 4; j++) nb_pred += is_kmer_in_cdbg(&pred[j]);

                    if (nb_pred == 1){

                        neigh_annot = get_annotation(&succ[i]);
                        tmp_annot = intersection_annotations(graph, 2, curr_annot, neigh_annot);

                        if (get_count_id_genomes(tmp_annot, graph) < nb_genomes_core){
                            free_BFT_annotation(neigh_annot);
                            free_BFT_annotation(tmp_annot);
                            break;
                        }

                        set_flag_kmer(V_VISITED, &succ[i], graph);

                        free_BFT_annotation(tmp_annot);
                        free_BFT_annotation(curr_annot);
                        curr_annot = neigh_annot;
                        neigh_annot = NULL;

                        l_simple_path++;

                        simple_path = realloc(simple_path, l_simple_path * sizeof(char));
                        ASSERT_NULL_PTR(simple_path, "extract_core_simple_paths()\n");

                        simple_path[l_simple_path - 1] = succ[i].kmer[graph->k - 1];

                        neigh = get_successors(&succ[i], graph);
                        free_BFT_kmer(succ, 4);
                        succ = neigh;

                        for (i = 0, nb_succ = 0; i < 4; i++){
                            if (is_kmer_in_cdbg(&succ[i])){
                                nb_succ++;
                                last = i;
                            }
                        }

                        i = last;
                    }
                    else{
                        set_flag_kmer(V_VISITED, &succ[i], graph);
                        break;
                    }
                }

                free_BFT_annotation(curr_annot);
                curr_annot = get_annotation(kmer);

                free_BFT_kmer(pred, 4);
                pred = get_predecessors(kmer, graph);

                for (i = 0, nb_pred = 0; i < 4; i++){
                    if (is_kmer_in_cdbg(&pred[i])){
                        nb_pred++;
                        break;
                    }
                }
            }

            if ((nb_pred == 1) && (get_flag_kmer(&pred[i], graph) != V_VISITED)){

                neigh = get_predecessors(&pred[i], graph);

                for (j = 0, nb_pred = 0; j < 4; j++) nb_pred += is_kmer_in_cdbg(&neigh[j]);

                free_BFT_kmer(neigh, 4);

                while ((nb_pred <= 1) && (get_flag_kmer(&pred[i], graph) != V_VISITED)){

                    free_BFT_kmer(succ, 4);
                    succ = get_successors(&pred[i], graph);

                    for (j=0, nb_succ = 0; j < 4; j++) nb_succ += is_kmer_in_cdbg(&succ[j]);

                    if (nb_succ == 1){

                        neigh_annot = get_annotation(&pred[i]);
                        tmp_annot = intersection_annotations(graph, 2, curr_annot, neigh_annot);

                        if (get_count_id_genomes(tmp_annot, graph) < nb_genomes_core){
                            free_BFT_annotation(neigh_annot);
                            free_BFT_annotation(tmp_annot);
                            break;
                        }

                        set_flag_kmer(V_VISITED, &pred[i], graph);

                        free_BFT_annotation(tmp_annot);
                        free_BFT_annotation(curr_annot);
                        curr_annot = neigh_annot;
                        neigh_annot = NULL;

                        l_simple_path++;

                        simple_path = realloc(simple_path, l_simple_path * sizeof(char));
                        ASSERT_NULL_PTR(simple_path, "extract_core_simple_paths()\n");

                        memmove(&simple_path[1], simple_path, (l_simple_path - 1) * sizeof(char));
                        simple_path[0] = pred[i].kmer[0];

                        neigh = get_predecessors(&pred[i], graph);
                        free_BFT_kmer(pred, 4);
                        pred = neigh;

                        for (i=0, nb_pred = 0; i < 4; i++){
                            if (is_kmer_in_cdbg(&pred[i])){
                                nb_pred++;
                                last = i;
                            }
                        }

                        i = last;
                    }
                    else{
                        set_flag_kmer(V_VISITED, &pred[i], graph);
                        break;
                    }
                }
            }

            simple_path = realloc(simple_path, (l_simple_path+1) * sizeof(char));
            ASSERT_NULL_PTR(simple_path, "extract_core_simple_paths()\n");

            simple_path[l_simple_path] = '\n';

            fwrite(simple_path, sizeof(char), l_simple_path+1, file);

            free(simple_path);
        }

        *longest_simple_path = MAX(*longest_simple_path, l_simple_path);

        free_BFT_kmer(succ, 4);
        free_BFT_kmer(pred, 4);

        free_BFT_annotation(curr_annot);
    }

    return 1;
}


/** Function extracting from a colored de Bruijn graph stored as a BFT all simple (non branching) core paths.
*   @param graph is a colored de Bruijn graph stored as a BFT.
*   @param core_ratio is a float (between 0 to 1) indicating the ratio of genome ids (compared to the total
*   number of genome ids inserted in graph) that a k-mer annotation must contain to be considered being part
*   of a core simple path.
*   @param filename_output is the name of the file where to write the simple core paths.
*/
void extract_simple_core_paths_to_disk(BFT* graph, double core_ratio, char* filename_output){

    ASSERT_NULL_PTR(graph, "extract_simple_core_paths_to_disk()\n")

    FILE* file;
    int longest_simple_path = 0; //Initiate the number of char. in the longest simple path to 0

    if ((file = fopen(filename_output, "w")) == NULL) //Create a file where to write simple paths
        ERROR("extract_simple_core_paths_to_disk(): failed to create/open output file.\n")

    set_neighbors_traversal(graph); //Prepare the graph for traversal
    set_marking(graph); //Prepare the graph to mark the vertices (as visited)

    //Iterate over the k-mers of the graph. extract_core_simple_paths() will be called on each k-mer
    //with file, &longest_simple_path and core_ratio as additional arguments
    iterate_over_kmers(graph, extract_core_simple_paths, file, &longest_simple_path, core_ratio);

    unset_marking(graph); //Delete all marking on vertices
    unset_neighbors_traversal(graph); //Unset the graph for traversal

    fclose(file); //Close the file

    printf("Longest simple core path has %d nuc.\n", longest_simple_path);

    return;
}

/** Function of type BFT_func_ptr starting a Breadth-First Search traversal from a k-mer.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param args contains additional parameters given to the calling function (here none).
*   @return true if it is a new connected component, else false.
*/
size_t BFS(BFT_kmer* kmer, BFT* graph, va_list args){

    BFT_kmer* curr_kmer;
    BFT_kmer* kmer_cpy;
    BFT_kmer* neigh;

    List* l;

    if (get_flag_kmer(kmer, graph) == V_NOT_VISITED){ //If the k-mer was not visited before

        set_flag_kmer(V_VISITED, kmer, graph); //Set the k-mer to visited

        l = List_create(); //Queue for the k-mers: mus visit their neighbors
        List_push(l, kmer); //Push k-mer in the queue

        while (List_count(l)){ //While there are k-mers in the queue

            curr_kmer = List_pop(l); //Get k-mer from the queue
            neigh = get_neighbors(curr_kmer, graph); //Get the possible neighbors of this k-mer

            for (int it_neigh = 0; it_neigh < 8; it_neigh++){ //For each possible neighbor

                //If the neighbor exists in the graph but was not visited already
                if (is_kmer_in_cdbg(&neigh[it_neigh]) && (get_flag_kmer(&neigh[it_neigh], graph) == V_NOT_VISITED)){

                    set_flag_kmer(V_VISITED, &neigh[it_neigh], graph); //Set the k-mer to visited

                    //Create an empty k-mer and copy the content of the current neighbor in it
                    kmer_cpy = create_empty_kmer();
                    memcpy(kmer_cpy, &neigh[it_neigh], sizeof(BFT_kmer));

                    List_push(l, kmer_cpy); //Push copied k-mer in the queue
                }
                //Free only the content of the neighbor because the structure neigh is still analyzed
                else free_BFT_kmer_content(&neigh[it_neigh], 1);
            }

            if (curr_kmer != kmer) free_BFT_kmer(curr_kmer, 1); //Free the k-mer (copy) we have got from the queue

            //Free the structure containing the 8 possible neighbors:
            //- Content of existing neighbors is in new k-mers in the queue
            //- Content of non existing neighbors was already freed
            free(neigh);
        }

        List_destroy(l); //Free the queue

        return true; //Return it is a new connected component
    }

    return false; //Return it is not a new connected component
}

/** Function of type BFT_func_ptr starting a Breadth-First Search traversal from a k-mer
*   that is part of a subgraph.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param args contains additional parameters given to the calling function. Here, it contains
*   a number of genome ids, followed by the genome ids, that a k-mer must contain to be considered
*   part of the subgraph.
*   @return true if it is a new connected component, else false.
*/
size_t BFS_subgraph(BFT_kmer* kmer, BFT* graph, va_list args){

    BFT_kmer* curr_kmer;
    BFT_kmer* kmer_cpy;
    BFT_kmer* neigh;

    List* l;

    int nb_id_genomes;

    if (get_flag_kmer(kmer, graph) == V_NOT_VISITED){ //If the k-mer was not visited before

        //Extract the number of genomes that must be present in a k-mer annotation to be
        //considered part of the subgraph
        nb_id_genomes = va_arg(args, int);

        set_flag_kmer(V_VISITED, kmer, graph); //Set the k-mer to visited

        //If the kmer is in the subgraph (it is annotated with the required genome ids)
        if (is_in_subgraph(kmer, graph, nb_id_genomes, args)){

            l = List_create(); //Queue for the k-mers (visit their neighbors)
            List_push(l, kmer); //Push first k-mer

            while (List_count(l)){ //While there are k-mers in the queue

                curr_kmer = List_pop(l); //Get k-mer from the queue
                neigh = get_neighbors(curr_kmer, graph); //Get the possible neighbors of this k-mer

                for (int it_neigh = 0; it_neigh < 8; it_neigh++){ //For each possible neighbor

                    //If the neighbor exists in the graph but was not visited already
                    if (is_kmer_in_cdbg(&neigh[it_neigh]) && (get_flag_kmer(&neigh[it_neigh], graph) == V_NOT_VISITED)){

                        set_flag_kmer(V_VISITED, &neigh[it_neigh], graph); //Set the neighbor to visited

                        if (is_in_subgraph(&neigh[it_neigh], graph, nb_id_genomes, args)){ //If the neighbor is part of the subgraph

                            //Create an empty k-mer and copy the content of the current neighbor in it
                            kmer_cpy = create_empty_kmer();
                            memcpy(kmer_cpy, &neigh[it_neigh], sizeof(BFT_kmer));

                            List_push(l, kmer_cpy); //Push copied neighbor in the queue
                        }
                        //Free only the content of the neighbor because the structure neigh is still analyzed
                        else free_BFT_kmer_content(&neigh[it_neigh], 1);
                    }
                    //Free only the content of the neighbor because the structure neigh is still analyzed
                    else free_BFT_kmer_content(&neigh[it_neigh], 1);
                }

                if (curr_kmer != kmer) free_BFT_kmer(curr_kmer, 1); //Free the copied k-mer we have got from the queue

                //Free the structure containing the 8 possible neighbors:
                //- Content of existing neighbors is in new k-mers in the queue
                //- Content of non existing neighbors was already freed
                free(neigh);
            }

            List_destroy(l); //Free the queue

            return true; //Return it is a new connected component
        }

        return false; //Return it is not a new connected component
    }

    return false; //Return it is not a new connected component
}

/** Function of type BFT_func_ptr starting a Depth-First Search traversal from a k-mer.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param args contains additional parameters given to the calling function (here none).
*   @return true if it is a new connected component, else false.
*/
size_t DFS(BFT_kmer* kmer, BFT* graph, va_list args){

    if (get_flag_kmer(kmer, graph) == V_NOT_VISITED){ //If the k-mer was not visited before

        set_flag_kmer(V_VISITED, kmer, graph); //Set the k-mer to visited

        BFT_kmer* succ = get_neighbors(kmer, graph); //Get the neighbors of the k-mer

        for (int it_succ = 0; it_succ < 8; it_succ++){ //For each possible neighbor
            //If the neighbor is in the graph, recursively DFS-traverse it
            if (is_kmer_in_cdbg(&succ[it_succ])) DFS(&succ[it_succ], graph, args);
        }

        free_BFT_kmer(succ, 8); //Free the allocated possible neighbors

        return true; //Return it is a new connected component
    }

    return false; //Return it is not a new connected component
}

/** Function of type BFT_func_ptr starting a Depth-First Search traversal from a k-mer
*   that is part of a subgraph.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param args contains additional parameters given to the calling function. Here, it contains
*   a number of genome ids, followed by the genome ids, that a k-mer must contain to be considered
*   part of the subgraph.
*   @return true if it is a new connected component, else false.
*/
size_t DFS_subgraph(BFT_kmer* kmer, BFT* graph, va_list args){

    va_list args_cpy;

    int nb_id_genomes;

    if (get_flag_kmer(kmer, graph) == V_NOT_VISITED){ //If the k-mer was not visited before

        //Copy the variable arguments list (because the function is recursive, args is probably
        //already in use in some other DFS_subgraph() function)
        va_copy(args_cpy, args);

        //Extract first integer argument, the number of genome ids that a k-mer must have in
        //its annotation to be considered part of the subgraph
        nb_id_genomes = va_arg(args_cpy, int);

        set_flag_kmer(V_VISITED, kmer, graph); //Set the k-mer to visited

        if (is_in_subgraph(kmer, graph, nb_id_genomes, args_cpy)){ //If the k-mer is part of the subgraph

            BFT_kmer* neigh = get_neighbors(kmer, graph); //Get the neighbors of the k-mer

            for (int it_neigh = 0; it_neigh < 8; it_neigh++){ //For each possible neighbor
                //If neighbor is in the graph, DFS-traverse it
                if (is_kmer_in_cdbg(&neigh[it_neigh])) DFS_subgraph(&neigh[it_neigh], graph, args);
            }

            free_BFT_kmer(neigh, 8); //Free the allocated possible neighbors
            va_end(args_cpy); //End the extraction of the copied variable arguments list

            return true; //Return it is a new connected component
        }

        va_end(args_cpy); //End the extraction of the copied variable arguments list

        return false; //Return it is not a new connected component
    }

    return false; //Return it is not a new connected component
}

/** Function computing if a k-mer is part of a subgraph. A subgraph is determined by k-mers having
*   in their annotation specific genome ids.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param nb_id_genomes is the number of genome ids that a k-mer must contain to be considered
*   part of the subgraph.
*   @param args contains additional parameters given to the calling function. Here, it contains
*   the genome ids that a k-mer must contain to be considered part of the subgraph.
*   @return true if it is a new connected component, else false.
*/
bool is_in_subgraph(BFT_kmer* kmer, BFT* graph, int nb_id_genomes, const va_list args){

    ASSERT_NULL_PTR(kmer, "is_in_subgraph()\n")
    ASSERT_NULL_PTR(graph, "is_in_subgraph()\n")
    ASSERT_NULL_PTR(args, "is_in_subgraph()\n")

    int j = 0;

    BFT_annotation* kmer_annot = get_annotation(kmer); //Get the annotation of the k-mer

    //Get the list (and number) of genome ids stored in this annotation
    uint32_t* list_ids = get_list_id_genomes(kmer_annot, graph);

    uint32_t curr_id;

    va_list args_cpy;

    free_BFT_annotation(kmer_annot); //Free the annotation (no longer necessary)

    //First element of the array returned by get_list_id_genomes() is the number of
    //genome ids stored after it. So the k-mer annotation must contains at least
    //nb_id_genomes genome ids and nb_id_genomes must be > 0
    if ((list_ids[0] >= nb_id_genomes) && (nb_id_genomes > 0)){

        //Copy the variable arguments list (args is probably already in use in the calling function)
        va_copy(args_cpy, args);

        //Extract from it the first genome id that must be present in the k-mer annotation
        curr_id = va_arg(args_cpy, uint32_t);

        for (int i = 1; i <= list_ids[0]; i++){ //Iterate over the array of genome ids

            if (list_ids[i] == curr_id){ //If the current genome ids of the array is the one we extracted

                j++; //Increment the number of genome ids  we are looking for found in the annotation

                if (j == nb_id_genomes) break; //If all of them are found, stop the search
                //Extract next genome id that must be present in the k-mer annotation
                curr_id = va_arg(args_cpy, uint32_t);
            }
            else if (list_ids[i] > curr_id) break; //If current genome id not found in the array, stop search
        }

        va_end(args_cpy); //End the copy of the variable arguments list
        free(list_ids); //Free array of genome ids from the annotation

        if (j == nb_id_genomes) return true; //Return that k-mer is part of the subgraph
    }

    free(list_ids); //Free array of genome ids from the annotation

    return false;
}

/** Function traversing a colored de Bruijn graph stored as a BFT.
*   @param graph is a BFT representing a colored de Bruijn graph.
*   @param f is the traversal function (DFS/BFS, DFS_subgraph/BFS_subgraph) to use.
*   @param ... is the additional arguments to transfer to f() (if there are some).
*/
void cdbg_traversal(BFT* graph, BFT_func_ptr f, ...){

    ASSERT_NULL_PTR(graph, "cdbg_traversal()\n")

    va_list args;

    va_start(args, f); //Start a variable arguments list (parameter ...)

    set_neighbors_traversal(graph); //Prepare the graph for traversal
    set_marking(graph); //Prepare the graph to mark the vertices (as visited)

    /*Iterate over the graph k-mers. f() will be called on each k-mer.
    v_iterate_over_kmers() is used here because cdbg_traversal() cannot
    transfer its variable arguments (parameter ...) to iterate_over_kmers()
    which also has variable arguments (parameter ...)*/
    v_iterate_over_kmers(graph, f, args);

    unset_marking(graph); //Delete all marking on vertices, unlock the graph
    unset_neighbors_traversal(graph); //Unset the graph for traversal

    va_end(args); //End the variable arguments list

    return;
}

/** Function of type BFT_func_ptr calling a traversal method (DFS/BFS, DFS_subgraph/BFS_subgraph)
*   on a k-mer to determine if it is in a new connected component.
*   @param kmer is a k-mer from the BFT graph.
*   @param graph is the BFT from which kmer is from.
*   @param args contains additional parameters given to the calling function. Here, it contains
*   the traversal method to call and its additional arguments.
*/
size_t nb_connected_components(BFT_kmer* kmer, BFT* graph, va_list args){

    //Extract from the variable arguments list args a pointer to the
    //current number of connected components found
    int* nb_connected_comp = va_arg(args, int*);

    //Extract from the variable arguments list args a pointer to the
    //function that iterate over the graph from a vertex (k-mer) and
    //return if it belongs to a new connected component
    BFT_func_ptr f = va_arg(args, BFT_func_ptr);

    //Update the number of connected components
    *nb_connected_comp += f(kmer, graph, args) == true;

    return 1;
}

/** Compute the number of connected components in a colored de-Bruijn graph.
*   @param graph is a BFT representing a colored de Bruijn graph.
*   @param ... is the additional arguments to transfer to nb_connected_components()
*   (traversal method and additional arguments if there are some).
*/
void get_nb_connected_component(BFT* graph, ...){

    ASSERT_NULL_PTR(graph, "get_nb_connected_component()\n")

    va_list args;

    va_start(args, graph); //Start a variable arguments list (parameter ...)

    set_neighbors_traversal(graph); //Prepare the graph for traversal
    set_marking(graph); //Prepare the graph for vertex marking (as visited)

    /*Iterate over the graph k-mers. nb_connected_components() will be called on each k-mer
    v_iterate_over_kmers() is used here because get_nb_connected_component() cannot
    transfer its variable arguments (parameter ...) to iterate_over_kmers() which also
    has variable arguments (parameter ...)*/
    v_iterate_over_kmers(graph, nb_connected_components, args);

    unset_marking(graph); //Delete all marking on vertices
    unset_neighbors_traversal(graph); //Unset the graph for traversal

    va_end(args);

    return;
}
