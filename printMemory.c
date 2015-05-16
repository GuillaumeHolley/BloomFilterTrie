#include "printMemory.h"

/* ---------------------------------------------------------------------------------------------------------------
*  add_memory_Used(mem1, mem2)
*  ---------------------------------------------------------------------------------------------------------------
*  add a structure memory_Used B to a structure memory_Used A
*  ---------------------------------------------------------------------------------------------------------------
*  mem1: pointer to a structure memory_Used
*  mem2: pointer to a structure memory_Used
*  ---------------------------------------------------------------------------------------------------------------
*/
memory_Used* create_memory_Used(){
    memory_Used* mem = calloc(1,sizeof(memory_Used));
    ASSERT_NULL_PTR(mem,"create_memory_Used()")

    return mem;
}

/* ---------------------------------------------------------------------------------------------------------------
*  add_memory_Used(mem1, mem2)
*  ---------------------------------------------------------------------------------------------------------------
*  add a structure memory_Used B to a structure memory_Used A
*  ---------------------------------------------------------------------------------------------------------------
*  mem1: pointer to a structure memory_Used
*  mem2: pointer to a structure memory_Used
*  ---------------------------------------------------------------------------------------------------------------
*/
void add_memory_Used(memory_Used* restrict mem1, memory_Used* restrict mem2){

    ASSERT_NULL_PTR(mem1,"add_memory_Used()")
    ASSERT_NULL_PTR(mem2,"add_memory_Used()")

    mem1->memory += mem2->memory;
    mem1->nb_node_visited += mem2->nb_node_visited;
    mem1->nb_CC_visited += mem2->nb_CC_visited;
    mem1->nb_UCptr_visited += mem2->nb_UCptr_visited;
    mem1->nb_kmers_in_UCptr += mem2->nb_kmers_in_UCptr;
    mem1->nb_pointers_used += mem2->nb_pointers_used;
    mem1->size_biggest_annot = MAX(mem1->size_biggest_annot, mem2->size_biggest_annot);

    mem1->kmers_in_CC63 += mem2->kmers_in_CC63;
    mem1->nb_CC63 += mem2->nb_CC63;
    mem1->kmers_in_CC54 += mem2->kmers_in_CC54;
    mem1->nb_CC54 += mem2->nb_CC54;
    mem1->kmers_in_CC45 += mem2->kmers_in_CC45;
    mem1->nb_CC45 += mem2->nb_CC45;
    mem1->kmers_in_CC36 += mem2->kmers_in_CC36;
    mem1->nb_CC36 += mem2->nb_CC36;
    mem1->kmers_in_CC27 += mem2->kmers_in_CC27;
    mem1->nb_CC27 += mem2->nb_CC27;
    mem1->kmers_in_CC18 += mem2->kmers_in_CC18;
    mem1->nb_CC18 += mem2->nb_CC18;
    mem1->kmers_in_CC9 += mem2->kmers_in_CC9;
    mem1->nb_CC9 += mem2->nb_CC9;
}

/* ---------------------------------------------------------------------------------------------------------------
*  printMemoryUsedFromNode(node, size_kmer, func_on_types)
*  ---------------------------------------------------------------------------------------------------------------
*  Return information about a node and its containers
*  ---------------------------------------------------------------------------------------------------------------
*  node: pointer on a Node structure
*  size_kmer: size of the suffixes stored from this node
*  func_on_types: ptrs_on_func structure, contains information to manipulate CCs field CC->children_type
*  ---------------------------------------------------------------------------------------------------------------
*/
memory_Used* printMemoryUsedFromNode(Node* restrict node, int size_kmer, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(node,"printMemoryUsedFromNode()")
    ASSERT_NULL_PTR(func_on_types,"printMemoryUsedFromNode()")

    memory_Used* mem = create_memory_Used();
    ASSERT_NULL_PTR(mem,"printMemoryUsedFromNode()")

    int i = 0;
    int node_nb_elem = 0;

    mem->nb_node_visited = 1;

    if (node->CC_array != NULL){
        while ((((CC*)node->CC_array)[node_nb_elem].type & 0x1) == 0) node_nb_elem++;
        node_nb_elem++;
    }

    for (i=0; i<node_nb_elem; i++){
        memory_Used* mem_tmp = printMemoryUsedFrom_CC(&(((CC*)node->CC_array)[i]), size_kmer, func_on_types);
        add_memory_Used(mem, mem_tmp);
        free(mem_tmp);
    }

    #if FREE == 1
    free(node->CC_array);
    #endif

    if (node->UC_array.suffixes != NULL){
        memory_Used* mem_tmp = printMemoryUsedFrom_UC(&(node->UC_array));

        if (mem_tmp != NULL){
            add_memory_Used(mem, mem_tmp);
            free(mem_tmp);
        }
        else {
            fprintf (stderr, "printMemoryUsedFrom_Node(): memory usage structure null\n");
            exit (EXIT_FAILURE);
        }

        #if FREE == 1
        free(node->UC_array.substrings);
        #endif
    }

    return mem;
}

/* ---------------------------------------------------------------------------------------------------------------
*  printMemoryUsedFrom_CC(cc, size_kmer, func_on_types)
*  ---------------------------------------------------------------------------------------------------------------
*  Return information about a CC and its children
*  ---------------------------------------------------------------------------------------------------------------
*  cc: pointer on a CC
*  size_kmer: size of the suffixes stored from this CC
*  func_on_types: ptrs_on_func structure, contains information to manipulate CCs field CC->children_type
*  ---------------------------------------------------------------------------------------------------------------
*/
memory_Used* printMemoryUsedFrom_CC(CC* restrict cc, int size_kmer, ptrs_on_func* restrict func_on_types){

    memory_Used* mem = create_memory_Used();
    ASSERT_NULL_PTR(mem,"printMemoryUsedFrom_CC()")

    int nb_skp = CEIL(cc->nb_elem,NB_CHILDREN_PER_SKP);
    int level = (size_kmer/SIZE_SEED)-1;
    int i = 0;
    int pos_Node = 0;

    mem->nb_CC_visited += 1;
    mem->nb_pointers_used = nb_skp;

    UC* uc;

    if (size_kmer == SIZE_SEED){
        //Memory annotations + annotations chunks links
        mem->nb_kmers_in_UCptr = cc->nb_elem;
        mem->kmers_in_CC9 += cc->nb_elem;
        mem->nb_CC9 += 1;
    }
    else {
        if (size_kmer == 63){
            mem->kmers_in_CC63 += cc->nb_elem;
            mem->nb_CC63 += 1;
        }
        else if (size_kmer == 54){
            mem->kmers_in_CC54 += cc->nb_elem;
            mem->nb_CC54 += 1;
        }
        else if (size_kmer == 45){
            mem->kmers_in_CC45 += cc->nb_elem;
            mem->nb_CC45 += 1;
        }
        else if (size_kmer == 36){
            mem->kmers_in_CC36 += cc->nb_elem;
            mem->nb_CC36 += 1;
        }
        else if (size_kmer == 27){
            mem->kmers_in_CC27 += cc->nb_elem;
            mem->nb_CC27 += 1;
        }
        else{
            mem->kmers_in_CC18 += cc->nb_elem;
            mem->nb_CC18 += 1;
        }

        for (i=0; i<nb_skp; i++){
            uc = &(((UC*)cc->children)[i]);
            mem->nb_kmers_in_UCptr += uc->nb_children;
            mem->memory += uc->nb_children * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;
            mem->size_biggest_annot = MAX(mem->size_biggest_annot, uc->size_annot);
        }

        for (i=0; i<cc->nb_elem; i++){
            memory_Used* mem_tmp = NULL;
            if (func_on_types[level].is_child((void*)cc, i) == 0){
                mem_tmp = printMemoryUsedFromNode(&(cc->children_Node_container[pos_Node]), size_kmer-SIZE_SEED, func_on_types);

                ASSERT_NULL_PTR(mem_tmp,"printMemoryUsedFrom_CC()")
                add_memory_Used(mem, mem_tmp);
                free(mem_tmp);

                pos_Node++;
            }
        }
    }

    #if FREE == 1
    freeCC(cc);
    #endif

    return mem;
}

/* ---------------------------------------------------------------------------------------------------------------
*  printMemoryUsedFrom_UC(uc, size_kmer)
*  ---------------------------------------------------------------------------------------------------------------
*  Return information about a UC
*  ---------------------------------------------------------------------------------------------------------------
*  uc: pointer on a UC
*  size_kmer: size of the suffixes stored in this UC
*  ---------------------------------------------------------------------------------------------------------------
*/
memory_Used* printMemoryUsedFrom_UC(UC* uc){
    ASSERT_NULL_PTR(uc,"printMemoryUsedFrom_UC()")

    memory_Used* mem = create_memory_Used();
    ASSERT_NULL_PTR(mem,"printMemoryUsedFrom_UC()")

    mem->nb_UCptr_visited = 1;
    mem->nb_pointers_used = 1;

    mem->memory = (uc->nb_children >> 1) * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;
    mem->nb_kmers_in_UCptr = uc->nb_children >> 1;
    mem->size_biggest_annot = uc->size_annot;

    return mem;
}

void printMemory(memory_Used* mem){

    ASSERT_NULL_PTR(mem,"printMemory()")

    printf("\nSize structure en Bytes = %f\n", mem->memory);
    printf("Nb nodes visited = %f\n", mem->nb_node_visited);
    printf("Nb CC visited = %f\n", mem->nb_CC_visited);
    printf("Nb UC_ptr visited = %f\n", mem->nb_UCptr_visited);
    printf("Nb kmers in UC_ptr = %f\n", mem->nb_kmers_in_UCptr);
    printf("Nb pointers stored = %f\n", mem->nb_pointers_used);
    printf("Size biggest annotation = %f\n", mem->size_biggest_annot);

    printf("Nb container CC63 = %f\n", mem->nb_CC63);
    printf("Nb kmers in CC63 = %f\n", mem->kmers_in_CC63);
    printf("Nb container CC54 = %f\n", mem->nb_CC54);
    printf("Nb kmers in CC54 = %f\n", mem->kmers_in_CC54);
    printf("Nb container CC45 = %f\n", mem->nb_CC45);
    printf("Nb kmers in CC45 = %f\n", mem->kmers_in_CC45);
    printf("Nb container CC36 = %f\n", mem->nb_CC36);
    printf("Nb kmers in CC36 = %f\n", mem->kmers_in_CC36);
    printf("Nb container CC27 = %f\n", mem->nb_CC27);
    printf("Nb kmers in CC27 = %f\n", mem->kmers_in_CC27);
    printf("Nb container CC18 = %f\n", mem->nb_CC18);
    printf("Nb kmers in CC18 = %f\n", mem->kmers_in_CC18);
    printf("Nb container CC9 = %f\n", mem->nb_CC9);
    printf("Nb kmers in CC9 = %f\n\n", mem->kmers_in_CC9);

    return;
}
