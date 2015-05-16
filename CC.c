#include "CC.h"

//Initialize macros used in create_ptrs_on_func()
DEFINE_STRUCT_AND_FUNC(63);
DEFINE_STRUCT_AND_FUNC(54);
DEFINE_STRUCT_AND_FUNC(45);
DEFINE_STRUCT_AND_FUNC(36);
DEFINE_STRUCT_AND_FUNC(27);
DEFINE_STRUCT_AND_FUNC(18);
DEFINE_STRUCT_AND_FUNC(9);

extern CC* createCC(int nb_bits_bf);
extern void initiateCC(CC* cc, int nb_bits_bf);
extern void freeCC(CC* cc);

extern void initializeUC(UC* uc);

extern uint8_t reverse_word_8(uint8_t v);
extern int popcnt_8(uint8_t v);
extern int popcnt_8_par(uint8_t* v, int start, int end);

extern uint16_t hash1(uint8_t* kmer);
extern uint16_t hash2_bis(uint8_t* kmer);

extern uint16_t dbj2(uint8_t c1, uint8_t c2, int size_bf);
extern uint16_t sdbm(uint8_t c1, uint8_t c2, int size_bf);

extern uint8_t* min_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot);
extern int max_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot);
extern int size_annot_sub(uint8_t* annot, int size_substring, int size_annot);

/* ---------------------------------------------------------------------------------------------------------------
*  transform2CC(uc, cc, size_suffix, func_on_types)
*  ---------------------------------------------------------------------------------------------------------------
*  Create a compressed container from an uncompressed container
*  ---------------------------------------------------------------------------------------------------------------
*  uc: pointer to a non-empty uncompressed container
*  cc: pointer to an empty compressed container
*  size_suffix: size of the suffixes in uc
*  func_on_types: ptrs_on_func structure used to manipulate the compressed container field cc->children_type
*  ---------------------------------------------------------------------------------------------------------------
*/
void transform2CC(UC* restrict uc, CC* restrict cc, int size_suffix, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(uc,"transform2CC()")
    ASSERT_NULL_PTR(cc,"transform2CC()")

    UC* uc_tmp;

    int i=0, k=-1, z=0, pos = 0;

    int size_annotation = uc->size_annot*sizeof(uint8_t);
    int nb_substrings_different = 0;

    int nb_cell = func_on_types->size_kmer_in_bytes;
    int size_line_uc = nb_cell + uc->size_annot;

    int nb_cell_children = func_on_types->size_kmer_in_bytes_minus_1;
    int size_line_uc_children;

    int begin_clust = 0;
    int it_children = 0;
    int it_clust = -1;

    int current_posFilter2 = -2;
    int posFilter2 = -1;
    int nb_1 = 0;
    int pos_skp_children = -1;
    int real_pos_i_substring = 0;
    int real_pos_i_uc = 0;

    uint16_t substring_prefix;
    uint16_t hash1_v = 0, hash2_v = 0;
    uint16_t size_bf = cc->type >> 8;

    int* new_order;

    uint8_t* substrings = malloc(func_on_types->nb_kmers_per_cc*NB_CELLS_PER_SEED*sizeof(uint8_t));
    ASSERT_NULL_PTR(substrings,"transform2CC()")

    uint8_t* sub_equal = calloc(func_on_types->nb_kmers_per_cc,sizeof(uint8_t));
    ASSERT_NULL_PTR(sub_equal,"transform2CC()")

    uint8_t* annot_sizes = min_size_per_sub(uc->suffixes, func_on_types->nb_kmers_per_cc, nb_cell, uc->size_annot);
    uint8_t** annot_extend = get_extend_annots(uc, nb_cell, func_on_types->nb_kmers_per_cc, 0, func_on_types->nb_kmers_per_cc-1);
    uint8_t** annot_cplx = get_annots_cplx_nodes(uc, nb_cell, func_on_types->nb_kmers_per_cc, 0, func_on_types->nb_kmers_per_cc-1);
    uint8_t* annot_cplx_sizes = min_size_per_annot_cplx_sub(uc->suffixes, func_on_types->nb_kmers_per_cc, uc->size_annot_cplx_nodes,
                                                            0, func_on_types->nb_kmers_per_cc-1);

    for (i=0; i<func_on_types->nb_kmers_per_cc; i++){

        if ((annot_extend != NULL) && (annot_extend[i] != NULL) && (annot_extend[i][0] != 0)) annot_sizes[i] = uc->size_annot+1;

        pos = i*NB_CELLS_PER_SEED;
        memcpy(&(substrings[pos]), &(uc->suffixes[i*size_line_uc]), NB_CELLS_PER_SEED);

        substrings[pos] = reverse_word_8(substrings[pos]);
        substrings[pos+1] = reverse_word_8(substrings[pos+1]);
        substrings[pos+2] = reverse_word_8(substrings[pos+2]) & 0xc0;

        substring_prefix = substrings[pos] & 0x3f;

        hash1_v = dbj2(substrings[pos+1], substring_prefix, MODULO_HASH);
        hash2_v = sdbm(substrings[pos+1], substring_prefix, MODULO_HASH);

        cc->BF_filter2[hash1_v/SIZE_CELL] |= MASK_POWER_8[hash1_v%SIZE_CELL];
        cc->BF_filter2[hash2_v/SIZE_CELL] |= MASK_POWER_8[hash2_v%SIZE_CELL];

        substring_prefix = (substrings[pos] << SIZE_CELL) | substrings[pos+1];

        substrings[pos] = substring_prefix >> 6;
        substrings[pos+1] = (substring_prefix << 2) | (substrings[pos+2] >> 6);
        substrings[pos+2] = (substring_prefix >> 8) & 0xc0;
    }

    new_order = quicksort_init(substrings, NB_CELLS_PER_SEED, 0, func_on_types->nb_kmers_per_cc-1);

    for (i=0; i<func_on_types->nb_kmers_per_cc; i++){

        if ((i!=0) && (memcmp(&(substrings[i*NB_CELLS_PER_SEED]), &(substrings[(i-1)*NB_CELLS_PER_SEED]), NB_CELLS_PER_SEED) == 0)){
            sub_equal[i] = 1;
            continue;
        }
        else{
            it_clust++;
            nb_substrings_different++;
        }

        if (it_clust == NB_CHILDREN_PER_SKP){

            z = begin_clust;
            while (z < i){
                k = MAX(k, MAX(annot_sizes[new_order[z]], annot_cplx_sizes[new_order[z]]));
                z++;
            }

            cc->children = realloc(cc->children, (it_children+1)*sizeof(UC));
            ASSERT_NULL_PTR(cc->children,"transform2CC()")

            uc_tmp = &(((UC*)cc->children)[it_children]);
            initializeUC(uc_tmp);

            uc_tmp->suffixes = calloc((i-begin_clust)*(nb_cell_children+k), sizeof(uint8_t));
            ASSERT_NULL_PTR(uc_tmp->suffixes,"transform2CC()")

            uc_tmp->size_annot = k;
            uc_tmp->nb_extended_annot = 0;

            k = -1;
            it_children++;
            begin_clust = i;
            it_clust = 0;
        }
    }

    z = begin_clust;
    while (z < i){
        k = MAX(k, MAX(annot_sizes[new_order[z]], annot_cplx_sizes[new_order[z]]));
        z++;
    }

    cc->children = realloc(cc->children, (it_children+1)*sizeof(UC));
    ASSERT_NULL_PTR(cc->children,"transform2CC()")

    uc_tmp = &(((UC*)cc->children)[it_children]);
    initializeUC(uc_tmp);

    uc_tmp->suffixes = calloc((i-begin_clust)*(nb_cell_children+k), sizeof(uint8_t));
    ASSERT_NULL_PTR(uc_tmp->suffixes,"transform2CC()")

    uc_tmp->size_annot = k;
    uc_tmp->nb_extended_annot = 0;

    k = -1;

    //Allocates memory for fields filter3 and extra_filter3 into the CC
    cc->filter3 = malloc(nb_substrings_different*sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->filter3,"transform2CC()")

    if (func_on_types->level_min == 1){
        cc->extra_filter3 = calloc(CEIL(nb_substrings_different, SIZE_CELL),sizeof(uint8_t));
        ASSERT_NULL_PTR(cc->extra_filter3,"transform2CC()")
    }

    //skp_children is not allocated if the CC is a leaf because in that case, there are exactly NB_CHILDREN_PER_SKP
    //annotations (suffix length = 0) in every array of CC->children
    if (size_suffix != SIZE_SEED) (*func_on_types->allocate_children_type)(cc, nb_substrings_different);

    for (i=0; i<func_on_types->nb_kmers_per_cc; i++){
        real_pos_i_substring = i*NB_CELLS_PER_SEED;
        real_pos_i_uc = new_order[i]*size_line_uc;

        //If the prefix of this suffix is the same as the suffix before, no need to update filter2 or filter3
        //We just update to number of suffixes stored in an array of CC->children
        if ((real_pos_i_substring!=0) && (sub_equal[i] == 1)) (*func_on_types->addNewElt)(cc, k, nb_substrings_different);
        else {

            k++;
            pos_skp_children = k/NB_CHILDREN_PER_SKP;
            uc_tmp = &(((UC*)cc->children)[pos_skp_children]);

            //Updates filter2 and filter3
            posFilter2 = (((uint16_t)substrings[real_pos_i_substring]) << 2) | (((uint16_t)substrings[real_pos_i_substring+1]) >> 6);
            cc->BF_filter2[size_bf+posFilter2/SIZE_CELL] |= MASK_POWER_8[posFilter2%SIZE_CELL];
            cc->filter3[k] = (substrings[real_pos_i_substring+1] << 2) | (substrings[real_pos_i_substring+2] >> 6);

            //If the prefix inserted shares the same p_u as the one inserted before, we update extra_filter3
            if (posFilter2 != current_posFilter2) nb_1++;

            //Update SkipFilter3 when a multiple of 248 unique prefix were inserted
            if (k%0xf8 == 0xf7){
                int new_size = (1024/SIZE_CELL)+((k+1)/0xf8);
                cc->BF_filter2 = realloc(cc->BF_filter2, (size_bf+new_size)*sizeof(uint8_t));
                cc->BF_filter2[size_bf+new_size-1] = (uint8_t)nb_1;
                nb_1 = 0;
            }

            if (size_suffix != SIZE_SEED) (*func_on_types->addNewElt)(cc, k, nb_substrings_different);
        }

        if (size_suffix == SIZE_SEED){

            if (annot_sizes[new_order[i]] != 0){

                if ((annot_extend != NULL) && (annot_extend[new_order[i]] != NULL) && (annot_extend[new_order[i]][0] != 0)){
                    memcpy(&(uc_tmp->suffixes[(k%NB_CHILDREN_PER_SKP) * uc_tmp->size_annot]), &(uc->suffixes[real_pos_i_uc + nb_cell]), (annot_sizes[new_order[i]]-1) * sizeof(uint8_t));
                    uc_tmp->suffixes[(k%NB_CHILDREN_PER_SKP) * uc_tmp->size_annot + annot_sizes[new_order[i]] - 1] = annot_extend[new_order[i]][0];
                }
                else memcpy(&(uc_tmp->suffixes[(k%NB_CHILDREN_PER_SKP) * uc_tmp->size_annot]), &(uc->suffixes[real_pos_i_uc + nb_cell]), annot_sizes[new_order[i]] * sizeof(uint8_t));
            }
            else memcpy(&(uc_tmp->suffixes[(k%NB_CHILDREN_PER_SKP) * uc_tmp->size_annot]), annot_cplx[new_order[i]], annot_cplx_sizes[new_order[i]] * sizeof(uint8_t));

            if (posFilter2 != current_posFilter2){
                current_posFilter2 = posFilter2;
                cc->extra_filter3[k/SIZE_CELL] |= MASK_POWER_8[k%SIZE_CELL];
            }
        }
        else{
            //Remove the prefix (which was inserted into filter2 and filter3) from the suffix
            int j=0;

            int nb_cell_to_delete = 2;
            if (size_suffix == 45) nb_cell_to_delete++;

            for (j=0; j < nb_cell - nb_cell_to_delete; j++){
                uc->suffixes[real_pos_i_uc+j] = uc->suffixes[real_pos_i_uc+j+2] >> 2;
                if (j+3 < nb_cell) uc->suffixes[real_pos_i_uc+j] |= uc->suffixes[real_pos_i_uc+j+3] << 6;
            }

            uc->suffixes[real_pos_i_uc+j-1] &= func_on_types->mask_shift_kmer;

            memmove(&(uc->suffixes[real_pos_i_uc+j]), &(uc->suffixes[real_pos_i_uc + j + nb_cell_to_delete]), size_annotation);

            size_line_uc_children = nb_cell_children + uc_tmp->size_annot;

            z = uc_tmp->nb_children * size_line_uc_children;

            memcpy(&(uc_tmp->suffixes[z]), &(uc->suffixes[real_pos_i_uc]), nb_cell_children * sizeof(uint8_t));

            if (annot_sizes[new_order[i]] != 0){

                if ((annot_extend != NULL) && (annot_extend[new_order[i]] != NULL) && (annot_extend[new_order[i]][0] != 0)){

                    memcpy(&(uc_tmp->suffixes[z+nb_cell_children]), &(uc->suffixes[real_pos_i_uc+j]), (annot_sizes[new_order[i]]-1)  * sizeof(uint8_t));
                    uc_tmp->suffixes[z + nb_cell_children + annot_sizes[new_order[i]] - 1] = annot_extend[new_order[i]][0];
                }
                else memcpy(&(uc_tmp->suffixes[z + nb_cell_children]), &(uc->suffixes[real_pos_i_uc + j]), annot_sizes[new_order[i]] * sizeof(uint8_t));
            }
            else{
                memcpy(&(uc_tmp->suffixes[z + nb_cell_children]), &(uc->suffixes[real_pos_i_uc + j]), size_suffix * sizeof(uint8_t));
                memcpy(&(uc_tmp->suffixes[z + nb_cell_children + size_suffix]), annot_cplx[new_order[i]], annot_cplx_sizes[new_order[i]] * sizeof(uint8_t));
            }

            uc_tmp->nb_children++;

            if (posFilter2 != current_posFilter2){
                current_posFilter2 = posFilter2;
                if (func_on_types->level_min == 0) uc_tmp->suffixes[z+nb_cell_children-1] |= 0x80;
                else cc->extra_filter3[k/SIZE_CELL] |= MASK_POWER_8[k%SIZE_CELL];
            }
        }
    }

    cc->nb_elem = nb_substrings_different; //Set the number of unique prefixes inserted into the CC

    free(substrings);
    free(new_order);
    free(sub_equal);
    free(annot_sizes);
    free(annot_cplx);
    free(annot_cplx_sizes);
    if (annot_extend != NULL) free(annot_extend);
}

/* ---------------------------------------------------------------------------------------------------------------
*  transform2CC_from_arraySuffix(array_suffix, cc, size_suffix, size_annot, annot_extend, annot_cplx, size_annot_cplx, func_on_types)
*  ---------------------------------------------------------------------------------------------------------------
*  Create a CC from an array of suffixes and annotations
*  ---------------------------------------------------------------------------------------------------------------
*  array_suffix: array of suffixes and annotations
*  cc: pointer to the CC in which array_suffix will be transformed
*  size_suffix: size of the suffixes in array_suffix
*  size_annot: size of the annotations in array_suffix
*  func_on_types: ptrs_on_func structure used to manipulate cc->children_type
*  ---------------------------------------------------------------------------------------------------------------
*/
void transform2CC_from_arraySuffix(uint8_t* restrict array_suffix, CC* restrict cc, int size_suffix, int size_annot, uint8_t** annot_extend, uint8_t** annot_cplx,
                                int size_annot_cplx, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(array_suffix,"transform2CC_from_arraySuffix()")
    ASSERT_NULL_PTR(cc,"transform2CC_from_arraySuffix()")

    int nb_cell = CEIL(size_suffix*2, SIZE_CELL);
    int size_line_uc = nb_cell + size_annot;

    int nb_cell_children = CEIL((size_suffix-SIZE_SEED)*2, SIZE_CELL);
    int size_line_uc_children;

    uint16_t size_bf = cc->type >> 8;

    uint8_t* substrings = malloc(func_on_types->nb_kmers_per_cc*NB_CELLS_PER_SEED*sizeof(uint8_t));
    ASSERT_NULL_PTR(substrings, "transform2CC_from_arraySuffix()")

    uint8_t* sub_equal = calloc(func_on_types->nb_kmers_per_cc,sizeof(uint8_t));
    ASSERT_NULL_PTR(sub_equal, "transform2CC_from_arraySuffix()")

    uint8_t* annot_sizes = min_size_per_sub(array_suffix, func_on_types->nb_kmers_per_cc, nb_cell, size_annot);
    ASSERT_NULL_PTR(annot_sizes, "transform2CC_from_arraySuffix()")

    uint8_t* annot_sizes_cplx = calloc(func_on_types->nb_kmers_per_cc,sizeof(uint8_t));
    ASSERT_NULL_PTR(annot_sizes_cplx, "transform2CC_from_arraySuffix()")

    int nb_substrings_different = 0;

    int i=0, k=-1, z=0, pos=0;
    uint16_t substring_prefix;
    uint16_t hash1_v = 0, hash2_v = 0;

    for (i=0; i<func_on_types->nb_kmers_per_cc; i++){

        if ((annot_extend != NULL) && (annot_extend[i] != NULL) && (annot_extend[i][0] != 0)) annot_sizes[i] = size_annot+1;
        if (annot_sizes[i] == 0) annot_sizes_cplx[i] = size_annot_sub(annot_cplx[i], 0, size_annot_cplx);

        pos = i*NB_CELLS_PER_SEED;
        memcpy(&(substrings[pos]), &(array_suffix[i*size_line_uc]), NB_CELLS_PER_SEED);

        substrings[pos] = reverse_word_8(substrings[pos]);
        substrings[pos+1] = reverse_word_8(substrings[pos+1]);
        substrings[pos+2] = reverse_word_8(substrings[pos+2]) & 0xc0;

        substring_prefix = substrings[pos] & 0x3f;

        hash1_v = dbj2(substrings[pos+1], substring_prefix, MODULO_HASH);
        hash2_v = sdbm(substrings[pos+1], substring_prefix, MODULO_HASH);

        cc->BF_filter2[hash1_v/SIZE_CELL] |= MASK_POWER_8[hash1_v%SIZE_CELL];
        cc->BF_filter2[hash2_v/SIZE_CELL] |= MASK_POWER_8[hash2_v%SIZE_CELL];

        substring_prefix = (substrings[pos] << SIZE_CELL) | substrings[pos+1];

        substrings[pos] = substring_prefix >> 6;
        substrings[pos+1] = (substring_prefix << 2) | (substrings[pos+2] >> 6);
        substrings[pos+2] = (substring_prefix >> 8) & 0xc0;
    }

    int* new_order = quicksort_init(substrings, NB_CELLS_PER_SEED, 0, func_on_types->nb_kmers_per_cc-1);

    int begin_clust = 0;
    int it_children = 0;
    int it_clust = -1;

    UC* uc;

    for (i=0; i<func_on_types->nb_kmers_per_cc; i++){

        if ((i!=0) && (memcmp(&(substrings[i*NB_CELLS_PER_SEED]), &(substrings[(i-1)*NB_CELLS_PER_SEED]), NB_CELLS_PER_SEED) == 0)){
            sub_equal[i] = 1;
            continue;
        }
        else{
            it_clust++;
            nb_substrings_different++;
        }

        if (it_clust == NB_CHILDREN_PER_SKP){

            z = begin_clust;
            while (z < i){
                k = MAX(k, MAX(annot_sizes[new_order[z]], annot_sizes_cplx[new_order[z]]));
                z++;
            }

            cc->children = realloc(cc->children, (it_children+1)*sizeof(UC));
            ASSERT_NULL_PTR(cc->children,"transform2CC_from_arraySuffix()")

            uc = &(((UC*)cc->children)[it_children]);
            initializeUC(uc);

            uc->suffixes = calloc((i-begin_clust)*(nb_cell_children+k), sizeof(uint8_t));
            ASSERT_NULL_PTR(uc->suffixes,"transform2CC_from_arraySuffix()")

            uc->size_annot = k;
            uc->nb_extended_annot = 0;

            k = -1;
            it_children++;
            begin_clust = i;
            it_clust = 0;
        }
    }

    z = begin_clust;
    while (z < i){
        k = MAX(k, MAX(annot_sizes[new_order[z]], annot_sizes_cplx[new_order[z]]));
        z++;
    }

    cc->children = realloc(cc->children, (it_children+1)*sizeof(UC));
    ASSERT_NULL_PTR(cc->children,"transform2CC_from_arraySuffix()")

    uc = &(((UC*)cc->children)[it_children]);
    initializeUC(uc);

    uc->suffixes = calloc((i-begin_clust)*(nb_cell_children+k), sizeof(uint8_t));
    ASSERT_NULL_PTR(uc->suffixes,"transform2CC_from_arraySuffix()")

    uc->size_annot = k;
    uc->nb_extended_annot = 0;

    k = -1;

    cc->filter3 = malloc(nb_substrings_different*sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->filter3,"transform2CC_from_arraySuffix()")

    if (func_on_types->level_min == 1){
        cc->extra_filter3 = calloc(CEIL(nb_substrings_different, SIZE_CELL),sizeof(uint8_t));
        ASSERT_NULL_PTR(cc->extra_filter3,"transform2CC_from_arraySuffix()")
    }

    if (size_suffix != SIZE_SEED) (*func_on_types->allocate_children_type)(cc, nb_substrings_different);

    int current_posFilter2 = -2;
    int posFilter2 = -1;
    int nb_1 = 0;
    int pos_skp_children = -1;
    int real_pos_i_substring = 0;
    int real_pos_i_uc = 0;

    for (i=0; i<func_on_types->nb_kmers_per_cc; i++){
        real_pos_i_substring = i*NB_CELLS_PER_SEED;
        real_pos_i_uc = new_order[i]*size_line_uc;

        if ((real_pos_i_substring != 0) && (sub_equal[i] == 1)) (*func_on_types->addNewElt)(cc, k, nb_substrings_different);
        else {

            k++;
            pos_skp_children = k/NB_CHILDREN_PER_SKP;
            uc = &(((UC*)cc->children)[pos_skp_children]);

            posFilter2 = (((uint16_t)substrings[real_pos_i_substring]) << 2) | (((uint16_t)substrings[real_pos_i_substring+1]) >> 6);
            cc->BF_filter2[size_bf+posFilter2/8] |= MASK_POWER_8[posFilter2%8];
            cc->filter3[k] = (substrings[real_pos_i_substring+1] << 2) | (substrings[real_pos_i_substring+2] >> 6);

            //If the prefix inserted shares the same p_u as the one inserted before, we update extra_filter3
            if (posFilter2 != current_posFilter2) nb_1++;

            if (k%0xf8 == 0xf7){
                int new_size = (1024/SIZE_CELL)+((k+1)/0xf8);
                cc->BF_filter2 = realloc(cc->BF_filter2, (new_size+size_bf)*sizeof(uint8_t));
                cc->BF_filter2[size_bf+new_size-1] = (uint8_t)nb_1;
                nb_1 = 0;
            }

            if (size_suffix != SIZE_SEED) (*func_on_types->addNewElt)(cc, k, nb_substrings_different);
        }

        if (size_suffix == SIZE_SEED){

            if (annot_sizes[new_order[i]] != 0){

                if ((annot_extend != NULL) && (annot_extend[new_order[i]] != NULL) && (annot_extend[new_order[i]][0] != 0)){
                    memcpy(&(uc->suffixes[(k%NB_CHILDREN_PER_SKP) * uc->size_annot]), &(array_suffix[real_pos_i_uc + nb_cell]), (annot_sizes[new_order[i]]-1) * sizeof(uint8_t));
                    uc->suffixes[(k%NB_CHILDREN_PER_SKP) * uc->size_annot + annot_sizes[new_order[i]] - 1] = annot_extend[new_order[i]][0];
                }
                else memcpy(&(uc->suffixes[(k%NB_CHILDREN_PER_SKP) * uc->size_annot]), &(array_suffix[real_pos_i_uc + nb_cell]), annot_sizes[new_order[i]] * sizeof(uint8_t));
            }
            else memcpy(&(uc->suffixes[(k%NB_CHILDREN_PER_SKP) * uc->size_annot]), annot_cplx[new_order[i]], annot_sizes_cplx[new_order[i]] * sizeof(uint8_t));

            if (posFilter2 != current_posFilter2){
                current_posFilter2 = posFilter2;
                cc->extra_filter3[k/SIZE_CELL] |= MASK_POWER_8[k%SIZE_CELL];
            }
        }
        else{
            int j=0;

            if (size_suffix != 36) array_suffix[real_pos_i_uc+nb_cell-1] &= 0x7f;

            int nb_cell_to_delete = 2;
            if (size_suffix == 45) nb_cell_to_delete++;

            for (j=0; j < nb_cell - nb_cell_to_delete; j++){
                array_suffix[real_pos_i_uc+j] = array_suffix[real_pos_i_uc+j+2] >> 2;
                if (j+3 < nb_cell) array_suffix[real_pos_i_uc+j] |= array_suffix[real_pos_i_uc+j+3] << 6;
            }

            array_suffix[real_pos_i_uc+j-1] &= func_on_types->mask_shift_kmer;

            memmove(&(array_suffix[real_pos_i_uc+j]), &(array_suffix[real_pos_i_uc + j + nb_cell_to_delete]), size_annot);

            size_line_uc_children = nb_cell_children + uc->size_annot;
            z = uc->nb_children * size_line_uc_children;

            //memset(&(uc->suffixes[z]), 0, size_line_uc_children);
            memcpy(&(uc->suffixes[z]), &(array_suffix[real_pos_i_uc]), nb_cell_children * sizeof(uint8_t));

            if (annot_sizes[new_order[i]] != 0){

                if ((annot_extend != NULL) && (annot_extend[new_order[i]] != NULL) && (annot_extend[new_order[i]][0] != 0)){

                    memcpy(&(uc->suffixes[z + nb_cell_children]), &(array_suffix[real_pos_i_uc + j]), (annot_sizes[new_order[i]] - 1) * sizeof(uint8_t));
                    uc->suffixes[z + nb_cell_children + annot_sizes[new_order[i]] - 1] =  annot_extend[new_order[i]][0];
                }
                else memcpy(&(uc->suffixes[z+nb_cell_children]), &(array_suffix[real_pos_i_uc+j]), annot_sizes[new_order[i]] * sizeof(uint8_t));
            }
            else{
                memcpy(&(uc->suffixes[z+nb_cell_children]), &(array_suffix[real_pos_i_uc+j]), size_suffix * sizeof(uint8_t));
                memcpy(&(uc->suffixes[z+nb_cell_children+size_suffix]), annot_cplx[new_order[i]], annot_sizes_cplx[new_order[i]] * sizeof(uint8_t));
            }

            if (posFilter2 != current_posFilter2){
                current_posFilter2 = posFilter2;
                if (func_on_types->level_min == 0) uc->suffixes[z+nb_cell_children-1] |= 0x80;
                else cc->extra_filter3[k/SIZE_CELL] |= MASK_POWER_8[k%SIZE_CELL];
            }

            uc->nb_children++;
        }
    }

    cc->nb_elem = nb_substrings_different;

    free(substrings);
    free(sub_equal);
    free(new_order);
    free(annot_sizes);
    free(annot_sizes_cplx);
}

/* ---------------------------------------------------------------------------------------------------------------
*  insertSP_CC(pres, sp, size_sp, id_genome, func_on_types, ann_inf, annot_sorted)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert a complete suffix prefix into a compressed container.
*  ---------------------------------------------------------------------------------------------------------------
*  pres: ptr to a resultPresence structure, contain information about the suffix prefixes in the container
*  kmer: ptr on the suffix prefix to insert
*  size_sp: size of the suffix prefix to insert, in nb of char.
*  id_genome: genome id of the suffix prefix to insert
*  func_on_types: ptr on a ptrs_on_func structure, used to manipulate CCs field children_type
*  ann_inf: ptr on annotation_inform structure, used to make the transition between reading and modifying an annot
*  ---------------------------------------------------------------------------------------------------------------
*/
void insertSP_CC(resultPresence* restrict pres, uint8_t* restrict sp, int size_sp, int id_genome, ptrs_on_func* restrict func_on_types,
                   annotation_inform* ann_inf, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(pres,"insertSP_CC()")

    CC* cc = pres->container;

    uint16_t size_bf = cc->type >> 8; //Bloom filter size in bytes
    uint8_t suf = (cc->type >> 2) & 0x3f; //Length v of p_v
    uint8_t pref = SIZE_SEED*2-suf; //Length u of p_u

    int pos_extra_filter3 = pres->pos_extra_filter3; //Position where p_v has to be inserted in filter3

    int size_filter2 = size_bf+(MASK_POWER_16[pref]/SIZE_CELL); // Size BF+filter2 in bytes
    int skip_filter2 = MASK_POWER_16[pref]/0xf8; //Size SkipFilter2 in bytes
    int size_filter2_n_skip = size_filter2;
    if (cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM) size_filter2_n_skip += skip_filter2; // Size BF+filter2+SkipFilter2 in bytes

    int i=0;
    int nb_skp = CEIL(cc->nb_elem, NB_CHILDREN_PER_SKP);

    // If the prefix is not present in filter2
    if (pres->presFilter2 == 0){

        //Compute position to set to 1 in filter2
        int posFilter2 = 0;
        if (pref==10) posFilter2 = (((uint16_t)pres->substring[0]) << 2) | ((uint16_t)(pres->substring[1]) >> 6);
        else posFilter2 = (((uint16_t)pres->substring[0]) << 6) | ((uint16_t)(pres->substring[1]) >> 2);

        int skip_posfilter2 = posFilter2/0xf8;

        //Set the position to 1 in filter2, eventually increment one cell in SkipFilter2 if it exists
        cc->BF_filter2[size_bf+posFilter2/SIZE_CELL] |= MASK_POWER_8[posFilter2%SIZE_CELL];
        if ((cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM) && (skip_posfilter2 < skip_filter2)) cc->BF_filter2[size_filter2+skip_posfilter2]++;

        //Compute the Hamming weight between position 0 and the position set to 1 in filter2
        //To know how many p_u are lexicographically inferior or equal to the one we insert
        if (pos_extra_filter3 == INT_MAX){
            int k=size_bf+posFilter2/SIZE_CELL;
            int cnt = 0;
            int hamming_weight = 0;

            //If SkipFilter2 exists, we use it to accelerate the computation of the Hamming weight
            if (cc->nb_elem >= NB_SUBSTRINGS_TRANSFORM){
                while ((cnt<skip_posfilter2) && (cnt < skip_filter2)){
                    hamming_weight += cc->BF_filter2[size_filter2+cnt];
                    cnt++;
                }
            }

            hamming_weight += popcnt_8_par(cc->BF_filter2, size_bf+cnt*31, k);

            uint8_t word_tmp = cc->BF_filter2[k];
            for (k=0; k<=posFilter2%SIZE_CELL; k++)
                if ((word_tmp & MASK_POWER_8[k]) != 0) hamming_weight++;

            pres->posFilter2 = hamming_weight; //pres->posFilter2 now contains this Hamming weight
        }
    }

    //Compute p_v, the v last char. of the prefix we try to insert
    uint8_t suffix = 0;
    if (suf==8) suffix = (pres->substring[1] << 2) | (pres->substring[2] >> 6);
    else if (suf==4) suffix = ((pres->substring[1] & 0x3) << 2) | (pres->substring[2] >> 6);

     //If p_u was also not present in the second filter before the call of this function,
     //we need to determine where to insert it in filter3. For that, we need to determine the position
     //of the first p_v in filter3 which share the same p_u as the prefix we try to insert. This position
     //is the one in extra_filter3 which has the HammingWeightFilter2-th 1. We also call it cluster start position.
    int hamming_weight = 0;
    int hamming_weight_0 = 0;
    int k=0;

    int j=0, m=0;
    int nb_cell_3rdlist = CEIL(cc->nb_elem,SIZE_CELL);
    uint8_t word_tmp;

    int size_line_children = 0;
    int pos_children = 0;
    int cpt_pv = 0;
    int cpt_node = 0;
    int nb_elem_in_pv = 0;

    int end, end_tmp;
    int it = 0;
    int it_tmp = 0;

    uint8_t* children;

    UC* uc;

    if (func_on_types->level_min == 1){

        if (pos_extra_filter3 == INT_MAX){

            //This loop tries to use SkipFilter3 (if it exists) to accelerate the Hamming weight comput.
            //by skipping counting 1s directly in extra_filter3
            int sum = 0;
            while(m<cc->nb_elem/0xf8){
                if ((sum = hamming_weight + cc->BF_filter2[size_filter2_n_skip+m]) < pres->posFilter2){
                    hamming_weight = sum;
                    m++;
                }
                else break;
            }

            for (k=m*31; k<nb_cell_3rdlist; k++){

                if ((sum = hamming_weight+popcnt_8(cc->extra_filter3[k])) >= pres->posFilter2){
                    word_tmp = cc->extra_filter3[k];

                    int size_word = 7;
                    if (k == nb_cell_3rdlist-1) size_word = (cc->nb_elem-1)%8;

                    for (j=0; j<=size_word; j++){
                        if (((word_tmp >> j)&1) == 1){
                            hamming_weight++;
                            if (hamming_weight == pres->posFilter2){
                                pos_extra_filter3 = k*8+j;
                                j++;
                                break;
                            }
                        }
                    }

                    if (k == nb_cell_3rdlist) k--;
                    goto MATCH;
                }
                else hamming_weight = sum;
            }
            //If the substring must be inserted in the last slot of the list
            if (k == nb_cell_3rdlist) k--;
        }

        MATCH: if (pos_extra_filter3 == INT_MAX) pos_extra_filter3 = cc->nb_elem;

        //If the p_v we try to insert in filter3 was not present before the call of this function...
        if (pres->posFilter3 == INT_MAX){
            //...but if p_u was in filter2, it means that the p_v in filter3 corresponding to this p_u
            //can be followed by more p_v having also the same p_u. We need to know how many.
            if ((pres->presFilter2 == 1) && (pos_extra_filter3 != cc->nb_elem)){
                //We know here where is the position of the first p_v in the extra_filter3 where to look for our p_v we want to insert.
                //We need now to know the end position of the last p_v. This is done by computing the number of 0 after the HammingWeightFilter2-th 1
                //in extra_filter3. The position of the last 0 is also called cluster end position. if there is no 0, start position = end position.
                int a=0;
                for (a=k; a<nb_cell_3rdlist; a++){
                    uint8_t word_tmp = cc->extra_filter3[a];
                    while (j<8){
                        if (((word_tmp >> j)&1) == 0){
                            if (a*8+j < cc->nb_elem) hamming_weight_0++;
                            else goto OUT_LOOP;
                        }
                        else goto OUT_LOOP;
                        j++;
                    }
                    j=0;
                }
                //We compare p_vs in the filter3 (in the cluster) to determine where exactly to insert our p_v
        OUT_LOOP:if (pos_extra_filter3 < cc->nb_elem){
                    if (suf==8){
                        for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){
                            if(cc->filter3[k] < suffix) pres->posFilter3 = k+1;
                            else{
                                pres->posFilter3 = k;
                                break;
                            }
                        }
                    }
                    else {
                        uint8_t tmp = 0;
                        for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){
                            if (IS_ODD(k)) tmp = cc->filter3[k/2] >> 4;
                            else tmp = cc->filter3[k/2] & 15;

                            if(tmp < suffix) pres->posFilter3 = k+1;
                            else{
                                pres->posFilter3 = k;
                                break;
                            }
                        }
                    }
                }
            }
            //...but if p_u was not in filter2, no p_vs in the filter3 share the same p_u like the one we have inserted
            else pres->posFilter3 = pos_extra_filter3;
        }
    }
    else{
        if (pos_extra_filter3 == INT_MAX){

            //This loop tries to use SkipFilter3 (if it exists) to accelerate the Hamming weight comput.
            //by skipping counting 1s directly in extra_filter3
            int sum = 0;
            while(m<cc->nb_elem/0xf8){
                if ((sum = hamming_weight + cc->BF_filter2[size_filter2_n_skip+m]) < pres->posFilter2){
                    hamming_weight = sum;
                    m++;
                }
                else break;
            }

            k = m*0xf8;
            pos_children = k/NB_CHILDREN_PER_SKP;
            (*func_on_types->count_Nodes_Children)(pres->container, pos_children*NB_CHILDREN_PER_SKP, k, &cpt_pv, &cpt_node);
            cpt_node += (*func_on_types->count_nodes)(pres->container, 0, pos_children*NB_CHILDREN_PER_SKP);

            nb_elem_in_pv = 0;
            it = k - pos_children * NB_CHILDREN_PER_SKP;
            it_tmp = it;

            int tmp_hamming_weight = hamming_weight;
            UC* uc_tmp;

            pres->count_children = cpt_pv;
            pres->count_nodes = cpt_node;
            pres->pos_children = pos_children;

            while (pres->pos_children < nb_skp){

                uc = &(((UC*)cc->children)[pres->pos_children]);
                size_line_children = func_on_types->size_kmer_in_bytes_minus_1+uc->size_annot;
                children = &(uc->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                if (pres->pos_children == nb_skp - 1) end = cc->nb_elem - pres->pos_children * NB_CHILDREN_PER_SKP;
                else end = NB_CHILDREN_PER_SKP;

                while (it < end){

                    it_tmp = it;

                    for (m=k; m<MIN(k+8, cc->nb_elem); m++){
                        if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, m)) == 0){
                            tmp_hamming_weight += cc->children_Node_container[pres->count_nodes].UC_array.nb_children & 0x1;
                            pres->count_nodes++;
                        }
                        else{
                            tmp_hamming_weight += children[pres->count_children*size_line_children] >> 7;
                            pres->count_children += nb_elem_in_pv;
                        }

                        it++;
                    }

                    if (tmp_hamming_weight >= pres->posFilter2){

                        uc_tmp = &(((UC*)cc->children)[pos_children]);
                        size_line_children = func_on_types->size_kmer_in_bytes_minus_1+uc_tmp->size_annot;
                        children = &(uc_tmp->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                        if (pos_children == nb_skp - 1) end_tmp = cc->nb_elem - pos_children * NB_CHILDREN_PER_SKP;
                        else end_tmp = NB_CHILDREN_PER_SKP;

                        while (it_tmp < end_tmp){

                            if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, k)) == 0){
                                hamming_weight += cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1;
                                cpt_node++;
                            }
                            else{
                                hamming_weight += children[cpt_pv*size_line_children] >> 7;
                                cpt_pv += nb_elem_in_pv;
                            }

                            if (hamming_weight == pres->posFilter2){
                                pos_extra_filter3 = k;
                                goto MATCH2;
                            }

                            k++;
                            it_tmp++;
                        }
                    }
                    else {
                        hamming_weight = tmp_hamming_weight;
                        cpt_pv = pres->count_children;
                        cpt_node = pres->count_nodes;
                        pos_children = pres->pos_children;
                        k = m;
                        //it_tmp = it;
                    }
                }

                it = 0;
                it_tmp = 0;
                pres->count_children = 0;
                cpt_pv = 0;

                pres->pos_children++;
                pos_children = pres->pos_children;
            }
        }
        else{
            pos_children = pres->pos_children;
            cpt_pv = pres->count_children;
            cpt_node = pres->count_nodes;
            k = INT_MAX;

            //it = 0;
            //it_tmp = 0;
        }

        it_tmp = it;

        MATCH2: if (pos_extra_filter3 == INT_MAX) pos_extra_filter3 = cc->nb_elem;

        if (pres->posFilter3 == INT_MAX){
            if ((pres->presFilter2 == 1) && (pos_extra_filter3 != cc->nb_elem)){
                if (pos_children < nb_skp){
                    //if (cpt_pv >= ((UC*)cc->children)[pos_children].nb_children){

                    if (pos_children == nb_skp - 1) end_tmp = cc->nb_elem - pos_children * NB_CHILDREN_PER_SKP;
                    else end_tmp = NB_CHILDREN_PER_SKP;

                    k++;
                    it_tmp++;

                    if (it_tmp >= end_tmp){
                        it = 0;
                        cpt_pv = 0;
                        pos_children++;
                    }
                    else it = it_tmp;

                    while (pos_children < nb_skp){

                        uc = &(((UC*)cc->children)[pos_children]);
                        size_line_children = func_on_types->size_kmer_in_bytes_minus_1 + uc->size_annot;
                        children = &(uc->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                        if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * NB_CHILDREN_PER_SKP;
                        else end = NB_CHILDREN_PER_SKP;

                        while (it < end){

                            if ((nb_elem_in_pv = (*func_on_types->getNbElts)(cc, k)) == 0){
                                if ((cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1) == 0) hamming_weight_0++;
                                else goto OUT_LOOP2;
                                cpt_node++;
                            }
                            else{
                                if (children[cpt_pv*size_line_children] < 0x80) hamming_weight_0++;
                                else goto OUT_LOOP2;
                                cpt_pv += nb_elem_in_pv;
                            }

                            k++;
                            it++;
                        }

                        it = 0;
                        cpt_pv = 0;
                        pos_children++;
                    }
                }

                //We compare p_vs in the filter3 (in the cluster) to determine where exactly to insert our p_v
                OUT_LOOP2:if (pos_extra_filter3 < cc->nb_elem){
                    if (suf==8){
                        for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){
                            if(cc->filter3[k] < suffix) pres->posFilter3 = k+1;
                            else{
                                pres->posFilter3 = k;
                                break;
                            }
                        }
                    }
                    else {
                        uint8_t tmp = 0;
                        for (k=pos_extra_filter3; k<=pos_extra_filter3+hamming_weight_0; k++){
                            if (IS_ODD(k)) tmp = cc->filter3[k/2] >> 4;
                            else tmp = cc->filter3[k/2] & 15;

                            if(tmp < suffix) pres->posFilter3 = k+1;
                            else{
                                pres->posFilter3 = k;
                                break;
                            }
                        }
                    }
                }
            }
            //...but if p_u was not in filter2, no p_vs in the filter3 share the same p_u like the one we have inserted
            else pres->posFilter3 = pos_extra_filter3;
        }
    }

    int to_insert_in_extralist3 = 0;
    int transform_one_2_replace = 0;
    int tmp = 0;

    //Determine how to modify extra_filter3 to insert our p_v in filter3
    if (pres->presFilter2 == 1){
        if (pres->posFilter3 == pos_extra_filter3){
            if (pres->posFilter3 == INT_MAX) pres->posFilter3 = cc->nb_elem;
            to_insert_in_extralist3 = 1;
            transform_one_2_replace = 1;
        }
        else if (pres->posFilter3 == INT_MAX) pres->posFilter3 = cc->nb_elem;
    }
    else{
        to_insert_in_extralist3 = 1;
        tmp = 1;
        if (pres->posFilter3 == INT_MAX) pres->posFilter3 = cc->nb_elem;
    }

    //Realloc the size of the list in the filter3
    if (suf==8){
        cc->filter3 = realloc(cc->filter3, (cc->nb_elem+1)*sizeof(uint8_t));
        ASSERT_NULL_PTR(cc->filter3,"insertSP_CC()")
    }
    else if (cc->nb_elem%2 == 0){
        cc->filter3 = realloc(cc->filter3, ((cc->nb_elem/2)+1)*sizeof(uint8_t));
        ASSERT_NULL_PTR(cc->filter3,"insertSP_CC()")
        cc->filter3[(cc->nb_elem/2)] = 0;
    }

    //Shift the p_vs in filter3 to insert our p_v
    if (suf==8){
        memmove(&(cc->filter3[pres->posFilter3+1]), &(cc->filter3[pres->posFilter3]), cc->nb_elem-pres->posFilter3);
        cc->filter3[pres->posFilter3] = suffix;
    }
    else {
        for (i=cc->nb_elem; i>pres->posFilter3; i--){
            if (IS_ODD(i)) cc->filter3[i/2] <<= 4;
            else cc->filter3[i/2] |= cc->filter3[(i-1)/2] >> 4;
        }
        if (IS_ODD(i)) cc->filter3[i/2] = (cc->filter3[i/2] & 0xf) | (suffix << 4);
        else cc->filter3[i/2] = (cc->filter3[i/2] & 0xf0) | (suffix & 0xf);
    }

    //Realloc extra_filter3 if needed
    if (func_on_types->level_min == 1){
        if (nb_cell_3rdlist*SIZE_CELL == cc->nb_elem){
            nb_cell_3rdlist++;
            cc->extra_filter3 = realloc(cc->extra_filter3, nb_cell_3rdlist*sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->extra_filter3,"insertSP_CC()")
            cc->extra_filter3[nb_cell_3rdlist-1] = 0;
        }

        //Shift and modify the bits after/before the position of insertion in the extra_filter3
        if (nb_cell_3rdlist > 1){
            for (i=nb_cell_3rdlist-1; i>(pres->posFilter3/SIZE_CELL); i--){
                cc->extra_filter3[i] <<= 1;
                cc->extra_filter3[i] |= (cc->extra_filter3[i-1] >> 7);
            }
        }

        i=pres->posFilter3/SIZE_CELL;

        if (pres->posFilter3%SIZE_CELL == 7){
            cc->extra_filter3[i] = (cc->extra_filter3[i] & 127) | (((uint8_t)to_insert_in_extralist3) << 7);

            if ((i+1 < cc->nb_elem) && (transform_one_2_replace == 1)){
                if ((cc->extra_filter3[i+1] & 1) == 1) cc->extra_filter3[i+1] &= 254;
                else cc->extra_filter3[i+1] = (cc->extra_filter3[i+1] & 254) | 1;
            }
        }
        else{
            word_tmp = 0;
            int k=0;
            for (k=7; k>=0; k--){
                uint8_t next = (cc->extra_filter3[i] >> k) & 1;

                if (k == (pres->posFilter3%8)){
                        if (transform_one_2_replace == 1){
                            if (next == 1) word_tmp = word_tmp << 1;
                            else word_tmp = (word_tmp << 1) | 1;
                        }
                        else word_tmp = (word_tmp << 1) | next;
                        word_tmp = (word_tmp << 1) | to_insert_in_extralist3;
                }
                else  word_tmp = (word_tmp << 1) | next;
            }

            cc->extra_filter3[i] = word_tmp;
        }
    }

    compute_best_mode(ann_inf, annot_sorted, NULL, 0, NULL, 0, id_genome); //Compute the size required by the annotation

    int pos_skp = pres->posFilter3/NB_CHILDREN_PER_SKP;
    int pos_sub = 0;

    //We shift the kmer: Its prefix p has been inserted in the filter2, filter3 and extra_filter3 so
    //we only need to keep only the suffix s to insert it into CC->children
    if (size_sp != SIZE_SEED){
        int j=0;
        int nb_cell_to_delete = 2;

        if (size_sp == 45) nb_cell_to_delete++;

        for (j=0; j<func_on_types->size_kmer_in_bytes - nb_cell_to_delete; j++){
            sp[j] = sp[j+2] >> 2;
            if (j+3 < func_on_types->size_kmer_in_bytes) sp[j] |= sp[j+3] << 6;
        }

        sp[j-1] &= func_on_types->mask_shift_kmer;

        int count = 0;

        //pos_sub is the position of insertion in CC->children[pos_skp]
        //count is the number of suffixes stored between the position of insertion in CC->children[pos_skp] and the end of this array
        uc = &(((UC*)cc->children)[pos_skp]);

        if (pres->posFilter3 < cc->nb_elem){
            count = (*func_on_types->count_children)(pres->container, pres->posFilter3, MIN((pos_skp+1)*NB_CHILDREN_PER_SKP, cc->nb_elem));
            pos_sub = uc->nb_children - count;
        }
        else if (pos_skp <= (CEIL(cc->nb_elem,NB_CHILDREN_PER_SKP)-1)) pos_sub = uc->nb_children;

        //Shift suffixes and their annotation in CC->children to insert our suffix
        (*func_on_types->add_skp_children)(pres->container, pres->posFilter3, pos_sub, count, func_on_types->size_kmer_in_bytes_minus_1, ann_inf->min_size);

        uc = &(((UC*)cc->children)[pos_skp]);

        uint8_t* ext_annot = get_extend_annot(uc, func_on_types->size_kmer_in_bytes_minus_1, uc->nb_children, pos_sub);
        if (ext_annot != NULL) memset(ext_annot, 0, sizeof(uint8_t));

        //Insert our suffix at the position computed previously
        size_line_children = pos_sub*(func_on_types->size_kmer_in_bytes_minus_1+uc->size_annot)+func_on_types->size_kmer_in_bytes_minus_1;
        memcpy(&(uc->suffixes[size_line_children-func_on_types->size_kmer_in_bytes_minus_1]), sp, func_on_types->size_kmer_in_bytes_minus_1);
        memset(&(uc->suffixes[size_line_children]), 0, uc->size_annot);

        //Create and insert the annotation next to the suffix
        modify_mode_annotation(ann_inf, &(uc->suffixes[size_line_children]), uc->size_annot, ext_annot, 1, id_genome);

        if ((ext_annot != NULL) && (ext_annot[0] == 0)) delete_extend_annots(uc, func_on_types->size_kmer_in_bytes_minus_1, uc->nb_children, pos_sub, pos_sub, 0, 0, 1);

        //Shift the counters in children_type to insert the new one
        (*func_on_types->realloc_and_int_children_type)(pres->container, cc->nb_elem, pres->posFilter3);

        if (func_on_types->level_min == 0){
            uc->suffixes[size_line_children-1] |= to_insert_in_extralist3 << 7;
            if (transform_one_2_replace == 1){
                if (((pres->posFilter3+1)%NB_CHILDREN_PER_SKP) == 0){
                    ((UC*)cc->children)[pos_skp+1].suffixes[func_on_types->size_kmer_in_bytes_minus_1-1] &= 0x7f;
                }
                else uc->suffixes[size_line_children+func_on_types->size_kmer_in_bytes_minus_1+uc->size_annot-1] &= 0x7f;
            }
        }
    }
    else {
        //prefix p has been inserted and the length of the suffix is 0 so CC->children only contains annotations
        pos_sub = pres->posFilter3 % NB_CHILDREN_PER_SKP;

        int nb_elem = add_skp_annotation(cc, pres->posFilter3, ann_inf->min_size);

        uc = &(((UC*)cc->children)[pos_skp]);

        uint8_t* ext_annot = get_extend_annot(uc, 0, nb_elem, pos_sub);

        if (ext_annot != NULL) ext_annot[0] = 0;
        memset(&(uc->suffixes[pos_sub * uc->size_annot]), 0, uc->size_annot);

        modify_mode_annotation(ann_inf, &(uc->suffixes[pos_sub * uc->size_annot]), uc->size_annot, ext_annot, 1, id_genome);

        if ((ext_annot != NULL) && (ext_annot[0] == 0)) delete_extend_annots(uc, 0, nb_elem, pos_sub, pos_sub, 0, 0, 1);
    }

    //Recompute SkipFilter3
    if (func_on_types->level_min == 1){

        for (i=pres->posFilter3/0xf8; i<cc->nb_elem/0xf8; i++){

            cc->BF_filter2[size_filter2_n_skip+i] += tmp;

            if ((i+1)*31 < nb_cell_3rdlist){

                tmp = cc->extra_filter3[(i+1)*31] & 1;
                cc->BF_filter2[size_filter2_n_skip+i] -= tmp;
            }
        }

        cc->nb_elem += 1;

        //Recompute SkipFilter3
        if ((cc->nb_elem%0xf8 == 0) && (cc->nb_elem >= 0xf8)){

            int new_cell = size_filter2_n_skip+(cc->nb_elem/0xf8);
            cc->BF_filter2 = realloc(cc->BF_filter2, new_cell*sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->BF_filter2,"insertSP_CC()")
            cc->BF_filter2[new_cell-1] = popcnt_8_par(cc->extra_filter3, nb_cell_3rdlist-31, nb_cell_3rdlist);
        }
    }
    else {
        int end = cc->nb_elem/248;
        int last_pos_node = 0;
        uint8_t* children = &(cc->BF_filter2[size_filter2_n_skip]);

        cpt_node = 0;

        for (i=pres->posFilter3/0xf8; i<end; i++){

            children[i] += tmp;
            tmp = (i+1)*0xf8;

            if (i+1 <= end){

                if ((*func_on_types->getNbElts)(cc, tmp) == 0){

                    cpt_node += (*func_on_types->count_nodes)(pres->container, last_pos_node, tmp);
                    last_pos_node = tmp;
                    tmp = cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1;
                }
                else tmp = ((UC*)cc->children)[i+1].suffixes[func_on_types->size_kmer_in_bytes_minus_1-1] >> 7;

                cc->BF_filter2[size_filter2_n_skip+i] -= tmp;
            }
        }

        cc->nb_elem += 1;

        //Recompute SkipFilter3
        if ((cc->nb_elem%0xf8 == 0) && (cc->nb_elem >= 0xf8)){

            int new_cell = size_filter2_n_skip+(cc->nb_elem/0xf8);

            tmp = new_cell-1;
            cc->BF_filter2 = realloc(cc->BF_filter2, new_cell*sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->BF_filter2,"insertSP_CC()")
            cc->BF_filter2[tmp] = 0;

            k = cc->nb_elem-0xf8;
            pos_children = k/NB_CHILDREN_PER_SKP;
            cpt_pv = 0;
            cpt_node = cc->nb_Node_children - (*func_on_types->count_nodes)(pres->container, k, cc->nb_elem);

            uc = &(((UC*)cc->children)[pos_children]);
            size_line_children = func_on_types->size_kmer_in_bytes_minus_1+uc->size_annot;
            children = &(uc->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

            it = 0;

            while (it < NB_CHILDREN_PER_SKP){

                nb_elem_in_pv = (*func_on_types->getNbElts)(cc, k);

                if (nb_elem_in_pv == 0){
                    cc->BF_filter2[tmp] += cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1;
                    cpt_node++;
                }
                else{
                    cc->BF_filter2[tmp] += children[cpt_pv*size_line_children] >> 7;
                    cpt_pv += nb_elem_in_pv;
                }

                k++;
                it++;
            }
        }
    }

    reinit_annotation_inform(ann_inf); //Reinit the annotation_inform structure
}



/* ---------------------------------------------------------------------------------------------------------------
*  transform_Filter2n3(cc, pref_size, suf_size, func_on_types)
*  ---------------------------------------------------------------------------------------------------------------
*  A suffix prefix p in CC is divided in 2 parts: the first pref_size char. (pref) and the remaining puf_size char. (suf).
*  When a CC contains NB_SUBSTRINGS_TRANSFORM suffix prefixes, increasing pref_size and decreasing pref_size decreases the memory used.
*  This is done by increasing the size of filter2 (2^pref_size) and decreasing size of filter3. Both have to be recomputed
*  ---------------------------------------------------------------------------------------------------------------
*  cc: pointer to a CC structure
*  pref_size: Current size of pref in the CC
*  suf_size: Current size of suf in the CC
*  ---------------------------------------------------------------------------------------------------------------
*/
void transform_Filter2n3(CC* cc, int pref_size, int suf_size, ptrs_on_func* restrict func_on_types){

    ASSERT_NULL_PTR(cc,"transform_Filter2n3()")
    ASSERT_NULL_PTR(func_on_types,"transform_Filter2n3()")

    uint16_t it_filter2 = 0;
    int it_filter3_tmp = 0;

    uint8_t current_pref_size = (SIZE_SEED*2)-((cc->type >> 2) & 0x3f);
    uint16_t size_bf = cc->type >> 8;

    int skip_filter_2 = MASK_POWER_16[pref_size]/0xf8;
    int skip_filter_3 = cc->nb_elem/0xf8;
    int size_filter2 = size_bf+(MASK_POWER_16[pref_size]/SIZE_CELL);
    int size_filter2_n_skip = size_filter2+skip_filter_2;

    //Allocation of memory for the new filter 2, we can re-write on the filter 3 and extra filter 3
    //The order of suffixes stays the same in the Filter 3
    uint8_t* filter2_tmp = calloc(size_filter2_n_skip+skip_filter_3, sizeof(uint8_t));
    ASSERT_NULL_PTR(filter2_tmp,"transform_Filter2n3()")

    if (suf_size==4){
        if (func_on_types->level_min == 1){
            //Iterates on every position of the Filter 2
            for (it_filter2=0; it_filter2<MASK_POWER_16[current_pref_size]; it_filter2++){
                //If the position is 0, no cluster of suffix associated
                if ((cc->BF_filter2[size_bf+it_filter2/8] & (MASK_POWER_8[it_filter2%8])) != 0){
                    int first_bit = 1;
                    //Iterates on the cluster of suffix
                    while((it_filter3_tmp < cc->nb_elem) && (((cc->extra_filter3[it_filter3_tmp/8] & MASK_POWER_8[it_filter3_tmp%8]) == 0) || (first_bit == 1))){
                        //Compute the new position of this suffix in the new Filter 2
                        uint16_t new_posFilter2 = (it_filter2 << 4) | ((uint16_t)(cc->filter3[it_filter3_tmp] >> 4));
                        //If the suffix is the smaller one to share the corresponding prefix indicated in Filter 2
                        if ((filter2_tmp[size_bf+new_posFilter2/8] & (MASK_POWER_8[new_posFilter2%8])) == 0){
                            filter2_tmp[size_bf+new_posFilter2/8] |= MASK_POWER_8[new_posFilter2%8];
                            if (new_posFilter2/0xf8 < skip_filter_2) filter2_tmp[size_filter2+new_posFilter2/0xf8]++;
                            if (it_filter3_tmp/0xf8 < skip_filter_3) filter2_tmp[size_filter2_n_skip+it_filter3_tmp/0xf8]++;
                            cc->extra_filter3[it_filter3_tmp/8] |= MASK_POWER_8[it_filter3_tmp%8];
                        }

                        if (IS_ODD(it_filter3_tmp)) cc->filter3[it_filter3_tmp/2] |= cc->filter3[it_filter3_tmp] << 4;
                        else cc->filter3[it_filter3_tmp/2] = cc->filter3[it_filter3_tmp] & 15;

                        it_filter3_tmp++;
                        first_bit=0;
                    }
                }
            }
        }
        else{
            int size_line_children;

            int pos_children = 0;
            int cpt_pv = 0;
            int cpt_node = 0;
            int nb_elem_in_pv;
            int new_posFilter2=0;
            int first_bit = 0;
            int end = 0;
            int it_children = 0;
            int nb_skp = CEIL(cc->nb_elem,NB_CHILDREN_PER_SKP);

            uint8_t* children;

            UC* uc;

            for (it_filter2=0; it_filter2<MASK_POWER_16[current_pref_size]; it_filter2++){
                //If the position is 0, no cluster of suffix associated
                if ((cc->BF_filter2[size_bf+it_filter2/8] & (MASK_POWER_8[it_filter2%8])) != 0){

                    first_bit = 1;

                    while (pos_children < nb_skp){

                        uc = &(((UC*)cc->children)[pos_children]);
                        size_line_children = func_on_types->size_kmer_in_bytes_minus_1+uc->size_annot;
                        children = &(uc->suffixes[func_on_types->size_kmer_in_bytes_minus_1-1]);

                        //it_children = 0;
                        it_children = it_filter3_tmp%NB_CHILDREN_PER_SKP;

                        if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * NB_CHILDREN_PER_SKP;
                        else end = NB_CHILDREN_PER_SKP;

                        while (it_children < end){

                            nb_elem_in_pv = (*func_on_types->getNbElts)(cc, it_filter3_tmp);

                            if (nb_elem_in_pv == 0){

                                if (((cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1) == 0) || (first_bit == 1)){

                                    if (first_bit == 1) first_bit=0;

                                    new_posFilter2 = (int)((((uint16_t)it_filter2) << 4) | ((uint16_t)(cc->filter3[it_filter3_tmp] >> 4)));

                                    if ((filter2_tmp[size_bf+new_posFilter2/8] & (MASK_POWER_8[new_posFilter2%8])) == 0){

                                        filter2_tmp[size_bf+new_posFilter2/8] |= MASK_POWER_8[new_posFilter2%8];
                                        if (new_posFilter2/0xf8 < skip_filter_2) filter2_tmp[size_filter2+new_posFilter2/0xf8]++;
                                        if (it_filter3_tmp/0xf8 < skip_filter_3) filter2_tmp[size_filter2_n_skip+it_filter3_tmp/0xf8]++;
                                        cc->children_Node_container[cpt_node].UC_array.nb_children |= 1;
                                    }

                                    if (IS_ODD(it_filter3_tmp)) cc->filter3[it_filter3_tmp/2] |= cc->filter3[it_filter3_tmp] << 4;
                                    else cc->filter3[it_filter3_tmp/2] = cc->filter3[it_filter3_tmp] & 15;
                                }
                                else goto OUT_LOOP;

                                cpt_node++;
                            }
                            else{
                                if ((children[cpt_pv*size_line_children] < 0x80)  || (first_bit == 1)){

                                    if (first_bit == 1) first_bit=0;

                                    new_posFilter2 = (int)((((uint16_t)it_filter2) << 4) | ((uint16_t)(cc->filter3[it_filter3_tmp] >> 4)));

                                    if ((filter2_tmp[size_bf+new_posFilter2/8] & (MASK_POWER_8[new_posFilter2%8])) == 0){

                                        filter2_tmp[size_bf+new_posFilter2/8] |= MASK_POWER_8[new_posFilter2%8];
                                        if (new_posFilter2/0xf8 < skip_filter_2) filter2_tmp[size_filter2+new_posFilter2/0xf8]++;
                                        if (it_filter3_tmp/0xf8 < skip_filter_3) filter2_tmp[size_filter2_n_skip+it_filter3_tmp/0xf8]++;
                                        children[cpt_pv*size_line_children] |= 0x80;
                                    }

                                    if (IS_ODD(it_filter3_tmp)) cc->filter3[it_filter3_tmp/2] |= cc->filter3[it_filter3_tmp] << 4;
                                    else cc->filter3[it_filter3_tmp/2] = cc->filter3[it_filter3_tmp] & 15;
                                }
                                else goto OUT_LOOP;

                                cpt_pv += nb_elem_in_pv;
                            }

                            it_filter3_tmp++;
                            it_children++;
                        }

                        cpt_pv = 0;
                        pos_children++;
                    }
                }
                OUT_LOOP: continue;
            }
        }

        //Minimizes the memory required by the Filter 3
        if (IS_ODD(cc->nb_elem)) cc->filter3 = realloc(cc->filter3, ((cc->nb_elem/2)+1)*sizeof(uint8_t));
        else cc->filter3 = realloc(cc->filter3, (cc->nb_elem/2)*sizeof(uint8_t));

        ASSERT_NULL_PTR(cc->filter3,"transform_Filter2n3()")

        if ((cc->type & 1) == 0) cc->type = (size_bf << 8) | 0x10; //xxxxxxxx|000100|00 -> s=4
        else cc->type = (size_bf << 8) | 0x11; //xxxxxxxx|000100|01 -> s=4

        memcpy(filter2_tmp, cc->BF_filter2, size_bf);
        free(cc->BF_filter2);
        cc->BF_filter2 = filter2_tmp;
    }
}

/* ---------------------------------------------------------------------------------------------------------------
*  add_skp_annotation(cc, position_type, size_annot)
*  ---------------------------------------------------------------------------------------------------------------
*  Same as add_skp_children() but only for CCs of type CC9. In these CCs, the suffix length is 0 so a suffix prefix can
*  only be linked to a single annotation. Thus, arrays in CC->children contain exactly NB_CHILDREN_PER_SKP annot.
*  ---------------------------------------------------------------------------------------------------------------
*  cc: ptr to a CC
*  position_type: position in cc->filter3 where was inserted the prefix
*  size_annot: annotation size, in bytes
*  ---------------------------------------------------------------------------------------------------------------
*/
int add_skp_annotation(CC* cc, int position_type, int size_annot){

    ASSERT_NULL_PTR(cc, "add_skp_annotation()")

    UC* uc;
    UC* uc_tmp;

    int cpt;
    int cpt_cplx;
    int position_child;

    int z=0;
    int max_size_z = 0;
    int max_size_z_minus1 = 0;
    int max_size_cplx_z_minus1 = 0;
    int start_pos_z_minus1 = 0;
    int nb_children = 0;
    int nb_cell_skp = CEIL(cc->nb_elem+1, NB_CHILDREN_PER_SKP);

    uint8_t* cplx_annot_z_minus1 = NULL;
    uint8_t* extend_annot_z_minus1 = NULL;

    int start = position_type/NB_CHILDREN_PER_SKP;

    if (cc->nb_elem%NB_CHILDREN_PER_SKP == 0){
        cc->children = realloc(cc->children, nb_cell_skp*sizeof(UC));
        ASSERT_NULL_PTR(cc->children, "add_skp_annotation()")

        uc = &(((UC*)cc->children)[nb_cell_skp-1]);
        initializeUC(uc);
        uc->size_annot = 1;
    }

    for (z=nb_cell_skp-1; z>=start; z--){

        uc = &(((UC*)cc->children)[z]);

        if (z == nb_cell_skp-1) nb_children = cc->nb_elem % NB_CHILDREN_PER_SKP;
        else nb_children = NB_CHILDREN_PER_SKP-1;

        if (z == start){

            position_child = position_type % NB_CHILDREN_PER_SKP;

            if (size_annot > uc->size_annot){
                if (nb_children == 0){
                    uc->suffixes = calloc(size_annot, sizeof(uint8_t));
                    uc->size_annot = size_annot;
                }
                else realloc_annotation(uc, 0, nb_children, size_annot, 1, position_child);
            }
            else if (size_annot < uc->size_annot){

                if (uc->nb_extended_annot != 0) max_size_z = uc->size_annot + 1;
                else max_size_z = max_size_per_sub(uc->suffixes, nb_children, 0, uc->size_annot);

                max_size_z = MAX(size_annot, max_size_z);

                if ((max_size_z > 0) && (uc->nb_extended_annot == 0) && (max_size_z < uc->size_annot)){
                    if (nb_children == 0){
                        uc->suffixes = calloc(max_size_z, sizeof(uint8_t));
                        uc->size_annot = max_size_z;
                    }
                    else realloc_annotation(uc, 0, nb_children, max_size_z, 1, position_child);
                }
                else {
                    uc->suffixes = realloc(uc->suffixes, ((nb_children+1) * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                            + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
                    ASSERT_NULL_PTR(uc->suffixes, "add_skp_annotation()")

                    memmove(&(uc->suffixes[(position_child+1) * uc->size_annot]),
                            &(uc->suffixes[position_child * uc->size_annot]),
                            (((nb_children - position_child) * uc->size_annot) + (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT)
                             + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));

                    shift_extended_annot(uc, 0, nb_children+1, position_child);
                }
            }
            else {
                uc->suffixes = realloc(uc->suffixes, ((nb_children + 1) * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                        + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
                ASSERT_NULL_PTR(uc->suffixes, "add_skp_annotation()")

                memmove(&(uc->suffixes[(position_child+1) * uc->size_annot]),
                        &(uc->suffixes[position_child * uc->size_annot]),
                        (((nb_children - position_child) * uc->size_annot) + (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT)
                         + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));

                shift_extended_annot(uc, 0, nb_children+1, position_child);
            }

            shift_annot_cplx_nodes(uc, 0, nb_children+1, position_child);

            return nb_children+1;
        }
        else {
            uc_tmp = &(((UC*)cc->children)[z-1]);

            start_pos_z_minus1 = NB_CHILDREN_PER_SKP - 1;

            cplx_annot_z_minus1 = get_annot_cplx_nodes(uc_tmp, 0, NB_CHILDREN_PER_SKP, start_pos_z_minus1);
            max_size_cplx_z_minus1 = size_annot_sub(cplx_annot_z_minus1, 0, uc_tmp->size_annot_cplx_nodes);

            if (max_size_cplx_z_minus1 != 0) cpt_cplx = 1;
            else cpt_cplx = 0;

            if (max_size_cplx_z_minus1 > uc->size_annot_cplx_nodes) increase_size_annot_cplx_nodes(uc, 0, nb_children, max_size_cplx_z_minus1, 1);
            else if (max_size_cplx_z_minus1 < uc->size_annot_cplx_nodes){
                max_size_z = max_size_annot_cplx_sub(&(uc->suffixes[nb_children * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
                                                        uc->nb_cplx_nodes, uc->size_annot_cplx_nodes, 0, nb_children-1);

                if (max_size_z < max_size_cplx_z_minus1) decrease_size_annot_cplx_nodes(uc, 0, nb_children, max_size_z);
            }

            if ((extend_annot_z_minus1 = get_extend_annot(uc_tmp, 0, NB_CHILDREN_PER_SKP, start_pos_z_minus1)) == NULL){
                max_size_z_minus1 = size_annot_sub(&(uc_tmp->suffixes[start_pos_z_minus1 * uc_tmp->size_annot]), 0, uc_tmp->size_annot);
            }
            else max_size_z_minus1 = uc_tmp->size_annot + 1;

            COMPUTE_ANNOT_EXT:
            if (max_size_z_minus1 > uc->size_annot+1) realloc_annotation(uc, 0, nb_children, max_size_z_minus1, 0, 0);
            else if (max_size_z_minus1 < uc->size_annot){

                if (uc->nb_extended_annot != 0) max_size_z = uc->size_annot + 1;
                else max_size_z = max_size_per_sub(uc->suffixes, nb_children, 0, uc->size_annot);

                max_size_z = MAX(max_size_z_minus1, max_size_z);

                if ((max_size_z > 0) && (uc->nb_extended_annot == 0) && (max_size_z < uc->size_annot))
                    realloc_annotation(uc, 0, nb_children, max_size_z, 0, 0);
            }

            cpt = 0;

            if (max_size_z_minus1 > uc->size_annot){

                if (extend_annot_z_minus1 != NULL) cpt = 1;
                else if (uc_tmp->suffixes[start_pos_z_minus1 * uc_tmp->size_annot + max_size_z_minus1 - 1] != 0) cpt = 1;

                if ((uc->nb_extended_annot + cpt) * SIZE_BYTE_EXT_ANNOT > nb_children + 1){
                    recopy_back_annot_extend(uc, 0, nb_children);
                    goto COMPUTE_ANNOT_EXT;
                }
            }

            nb_children++;

            uc->suffixes = realloc(uc->suffixes, ((nb_children * uc->size_annot) + (uc->nb_extended_annot + cpt) * SIZE_BYTE_EXT_ANNOT
                                    + (uc->nb_cplx_nodes + cpt_cplx) * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
            ASSERT_NULL_PTR(uc->suffixes, "add_skp_annotation()")

            memmove(&(uc->suffixes[uc->size_annot]), uc->suffixes, ((nb_children-1) * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                     + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));

            memset(uc->suffixes, 0, uc->size_annot * sizeof(uint8_t));

            if (max_size_z_minus1 > uc->size_annot){

                if (extend_annot_z_minus1 != NULL){

                    memcpy(uc->suffixes, &(uc_tmp->suffixes[start_pos_z_minus1 * uc_tmp->size_annot]), uc_tmp->size_annot*sizeof(uint8_t));
                    insert_extend_annot(uc, 0, nb_children, 0, extend_annot_z_minus1[0], 1);
                }
                else{

                    memcpy(uc->suffixes, &(uc_tmp->suffixes[start_pos_z_minus1 * uc_tmp->size_annot]), (max_size_z_minus1-1) * sizeof(uint8_t));

                    if (uc_tmp->suffixes[start_pos_z_minus1 * uc_tmp->size_annot + max_size_z_minus1 - 1] == 0) shift_extended_annot(uc, 0, nb_children, 0);
                    else insert_extend_annot(uc, 0, nb_children, 0, uc_tmp->suffixes[start_pos_z_minus1 * uc_tmp->size_annot + max_size_z_minus1 - 1], 1);
                }
            }
            else{
                if (extend_annot_z_minus1 != NULL){

                    memcpy(uc->suffixes, &(uc_tmp->suffixes[start_pos_z_minus1 * uc_tmp->size_annot]), uc_tmp->size_annot * sizeof(uint8_t));
                    uc->suffixes[uc_tmp->size_annot] = extend_annot_z_minus1[0];
                }
                else memcpy(uc->suffixes, &(uc_tmp->suffixes[start_pos_z_minus1 * uc_tmp->size_annot]), max_size_z_minus1 * sizeof(uint8_t));

                shift_extended_annot(uc, 0, nb_children, 0);
            }

            if (cplx_annot_z_minus1 != NULL) insert_annot_cplx_nodes(uc, 0, nb_children, 0, cplx_annot_z_minus1, max_size_cplx_z_minus1, 1);
            else shift_annot_cplx_nodes(uc, 0, nb_children, 0);

            delete_extend_annots(uc_tmp, 0, NB_CHILDREN_PER_SKP, start_pos_z_minus1, start_pos_z_minus1, 0, 1, 0);
            delete_annot_cplx_nodes(uc_tmp, 0, NB_CHILDREN_PER_SKP, start_pos_z_minus1, start_pos_z_minus1, 1, 1, 0);

            extend_annot_z_minus1 = NULL;
            cplx_annot_z_minus1 = NULL;
        }
    }

    return 0;
}

/* ---------------------------------------------------------------------------------------------------------------
*  create_ptrs_on_func(size_min, size_max)
*  ---------------------------------------------------------------------------------------------------------------
*  Creates an array of ptrs_on_func for each possible suffix length between SIZE_SEED and k (the kmers length).
*  Each cell is a set of pointers on functions which manipulate the field children_type of CCs.
*  ---------------------------------------------------------------------------------------------------------------
*  size_min: minimum size of a suffix (in char.)
*  size_max: maximum size of a suffix (in char.)
*  ---------------------------------------------------------------------------------------------------------------
*/
ptrs_on_func* create_ptrs_on_func(int size_min, int size_max){
    ptrs_on_func* ptr = NULL;
    int nb_sizes = 0;
    int i = 0;

    for (i=size_min; i<=size_max; i+=SIZE_SEED){
        ptr = realloc(ptr, (nb_sizes+1)*sizeof(ptrs_on_func));
        ASSERT_NULL_PTR(ptr,"create_ptrs_on_func()")

        ptr[nb_sizes].size_kmer_in_bytes = CEIL(i*2,SIZE_CELL);
        ptr[nb_sizes].size_kmer_in_bytes_minus_1 = MAX(0,CEIL((i-SIZE_SEED)*2,SIZE_CELL));

        if (i == size_max){
            ptr[nb_sizes].level_min = 1;
            ptr[nb_sizes].root = 1;
        }
        else{
            ptr[nb_sizes].level_min = 0;
            ptr[nb_sizes].root = 0;
        }

        if (i == 63){
            ptr[nb_sizes].count_nodes = count_Nodes63;
            ptr[nb_sizes].count_children = count_Children63;
            ptr[nb_sizes].count_Nodes_Children = count_Nodes_Children63;
            ptr[nb_sizes].is_child = isChild63;
            ptr[nb_sizes].add_skp_children = add_skp_children63;
            ptr[nb_sizes].realloc_and_int_children_type = realloc_and_int_children_type63;
            ptr[nb_sizes].addNewElt = addNewElt63;
            ptr[nb_sizes].getNbElts = getNbElts63;
            ptr[nb_sizes].resetChildrenType = resetChildrenType63;
            ptr[nb_sizes].allocate_children_type = allocate_children_type63;
            ptr[nb_sizes].nb_kmers_per_cc = NB_KMERS_PER_CC63;
            ptr[nb_sizes].mask_shift_kmer = 0xf;
            ptr[nb_sizes].size_type_kmer = sizeof(TYPE63);
        }
        else if (i == 54){
            ptr[nb_sizes].count_nodes = count_Nodes54;
            ptr[nb_sizes].count_children = count_Children54;
            ptr[nb_sizes].count_Nodes_Children = count_Nodes_Children54;
            ptr[nb_sizes].is_child = isChild54;
            ptr[nb_sizes].add_skp_children = add_skp_children54;
            ptr[nb_sizes].realloc_and_int_children_type = realloc_and_int_children_type54;
            ptr[nb_sizes].addNewElt = addNewElt54;
            ptr[nb_sizes].getNbElts = getNbElts54;
            ptr[nb_sizes].resetChildrenType = resetChildrenType54;
            ptr[nb_sizes].allocate_children_type = allocate_children_type54;
            ptr[nb_sizes].nb_kmers_per_cc = NB_KMERS_PER_CC54;
            ptr[nb_sizes].mask_shift_kmer = 0x3;
            ptr[nb_sizes].size_type_kmer = sizeof(TYPE54);
        }
        else if (i == 45){
            ptr[nb_sizes].count_nodes = count_Nodes45;
            ptr[nb_sizes].count_children = count_Children45;
            ptr[nb_sizes].count_Nodes_Children = count_Nodes_Children45;
            ptr[nb_sizes].is_child = isChild45;
            ptr[nb_sizes].add_skp_children = add_skp_children45;
            ptr[nb_sizes].realloc_and_int_children_type = realloc_and_int_children_type45;
            ptr[nb_sizes].addNewElt = addNewElt45;
            ptr[nb_sizes].getNbElts = getNbElts45;
            ptr[nb_sizes].resetChildrenType = resetChildrenType45;
            ptr[nb_sizes].allocate_children_type = allocate_children_type45;
            ptr[nb_sizes].nb_kmers_per_cc = NB_KMERS_PER_CC45;
            ptr[nb_sizes].mask_shift_kmer = 0xff;
            ptr[nb_sizes].size_type_kmer = sizeof(TYPE45);
            ptr[nb_sizes].level_min = 1;
        }
        else if (i == 36){
            ptr[nb_sizes].count_nodes = count_Nodes36;
            ptr[nb_sizes].count_children = count_Children36;
            ptr[nb_sizes].count_Nodes_Children = count_Nodes_Children36;
            ptr[nb_sizes].is_child = isChild36;
            ptr[nb_sizes].add_skp_children = add_skp_children36;
            ptr[nb_sizes].realloc_and_int_children_type = realloc_and_int_children_type36;
            ptr[nb_sizes].addNewElt = addNewElt36;
            ptr[nb_sizes].getNbElts = getNbElts36;
            ptr[nb_sizes].resetChildrenType = resetChildrenType36;
            ptr[nb_sizes].allocate_children_type = allocate_children_type36;
            ptr[nb_sizes].nb_kmers_per_cc = NB_KMERS_PER_CC36;
            ptr[nb_sizes].mask_shift_kmer = 0x3f;
            ptr[nb_sizes].size_type_kmer = sizeof(TYPE36);
        }
        else if (i == 27){
            ptr[nb_sizes].count_nodes = count_Nodes27;
            ptr[nb_sizes].count_children = count_Children27;
            ptr[nb_sizes].count_Nodes_Children = count_Nodes_Children27;
            ptr[nb_sizes].is_child = isChild27;
            ptr[nb_sizes].add_skp_children = add_skp_children27;
            ptr[nb_sizes].realloc_and_int_children_type = realloc_and_int_children_type27;
            ptr[nb_sizes].addNewElt = addNewElt27;
            ptr[nb_sizes].getNbElts = getNbElts27;
            ptr[nb_sizes].resetChildrenType = resetChildrenType27;
            ptr[nb_sizes].allocate_children_type = allocate_children_type27;
            ptr[nb_sizes].nb_kmers_per_cc = NB_KMERS_PER_CC27;
            ptr[nb_sizes].mask_shift_kmer = 0xf;
            ptr[nb_sizes].size_type_kmer = sizeof(TYPE27);
        }
        else if (i == 18){
            ptr[nb_sizes].count_nodes = count_Nodes18;
            ptr[nb_sizes].count_children = count_Children18;
            ptr[nb_sizes].count_Nodes_Children = count_Nodes_Children18;
            ptr[nb_sizes].is_child = isChild18;
            ptr[nb_sizes].add_skp_children = add_skp_children18;
            ptr[nb_sizes].realloc_and_int_children_type = realloc_and_int_children_type18;
            ptr[nb_sizes].addNewElt = addNewElt18;
            ptr[nb_sizes].getNbElts = getNbElts18;
            ptr[nb_sizes].resetChildrenType = resetChildrenType18;
            ptr[nb_sizes].allocate_children_type = allocate_children_type18;
            ptr[nb_sizes].nb_kmers_per_cc = NB_KMERS_PER_CC18;
            ptr[nb_sizes].mask_shift_kmer = 0x3;
            ptr[nb_sizes].size_type_kmer = sizeof(TYPE18);
        }
        else{
            ptr[nb_sizes].count_nodes = count_Nodes9;
            ptr[nb_sizes].count_children = count_Children9;
            ptr[nb_sizes].count_Nodes_Children = count_Nodes_Children9;
            ptr[nb_sizes].is_child = isChild9;
            ptr[nb_sizes].add_skp_children = add_skp_children9;
            ptr[nb_sizes].realloc_and_int_children_type = realloc_and_int_children_type9;
            ptr[nb_sizes].addNewElt = addNewElt9;
            ptr[nb_sizes].getNbElts = getNbElts9;
            ptr[nb_sizes].resetChildrenType = resetChildrenType9;
            ptr[nb_sizes].allocate_children_type = allocate_children_type9;
            ptr[nb_sizes].nb_kmers_per_cc = NB_KMERS_PER_CC9;
            ptr[nb_sizes].mask_shift_kmer = 0xff;
            ptr[nb_sizes].size_type_kmer = sizeof(TYPE9);
            ptr[nb_sizes].level_min = 1;
        }

        nb_sizes++;
    }

    return ptr;
}

/* ---------------------------------------------------------------------------------------------------------------
*  build_skip_nodes(node, size_kmer, func_on_types)
*  ---------------------------------------------------------------------------------------------------------------
*  ---------------------------------------------------------------------------------------------------------------
*  size_min: minimum size of a suffix (in char.)
*  size_max: maximum size of a suffix (in char.)
*  ---------------------------------------------------------------------------------------------------------------
*/
uint16_t** build_skip_nodes(Node* node, int size_kmer, ptrs_on_func* func_on_types){

    ASSERT_NULL_PTR(node,"build_skip_nodes()")
    ASSERT_NULL_PTR(node->CC_array,"build_skip_nodes()")

    int node_nb_elem = 0, i = 0, j = 0, nb_skp = 0, level = (size_kmer/SIZE_SEED)-1;
    uint16_t** nodes_skip;
    uint16_t previous_count;

    while ((((CC*)node->CC_array)[node_nb_elem].type & 0x1) == 0) node_nb_elem++;
    node_nb_elem++;

    nodes_skip = malloc(node_nb_elem*sizeof(uint16_t*));
    ASSERT_NULL_PTR(nodes_skip,"build_skip_nodes()")

    for (i=0; i<node_nb_elem; i++){
        previous_count = 0;
        nb_skp = ((CC*)node->CC_array)[i].nb_elem / SIZE_CLUST_SKIP_NODES;

        nodes_skip[i] = malloc(nb_skp*sizeof(uint16_t));
        ASSERT_NULL_PTR(nodes_skip[i],"build_skip_nodes()")

        for (j=0; j<nb_skp; j++){
            nodes_skip[i][j] = previous_count + (*func_on_types[level].count_nodes)(&(((CC*)node->CC_array)[i]), j*SIZE_CLUST_SKIP_NODES, (j+1)*SIZE_CLUST_SKIP_NODES);
            previous_count = nodes_skip[i][j];
        }
    }

    return nodes_skip;
}

void free_skip_nodes(Node* node, uint16_t** skp_nodes){

    ASSERT_NULL_PTR(node,"build_skip_nodes()")
    ASSERT_NULL_PTR(node->CC_array,"build_skip_nodes()")
    ASSERT_NULL_PTR(skp_nodes,"build_skip_nodes()")

    int i = -1;

    do {
        i++;
        free(skp_nodes[i]);
    }
    while ((((CC*)node->CC_array)[i].type & 0x1) == 0);

    free(skp_nodes);

    return;

}
