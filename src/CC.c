#include "./../lib/CC.h"

const uint8_t POPCOUNT_8bit[256] = {
    /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
    /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
    /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
    /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4,
    /* 10 */ 1, /* 11 */ 2, /* 12 */ 2, /* 13 */ 3,
    /* 14 */ 2, /* 15 */ 3, /* 16 */ 3, /* 17 */ 4,
    /* 18 */ 2, /* 19 */ 3, /* 1a */ 3, /* 1b */ 4,
    /* 1c */ 3, /* 1d */ 4, /* 1e */ 4, /* 1f */ 5,
    /* 20 */ 1, /* 21 */ 2, /* 22 */ 2, /* 23 */ 3,
    /* 24 */ 2, /* 25 */ 3, /* 26 */ 3, /* 27 */ 4,
    /* 28 */ 2, /* 29 */ 3, /* 2a */ 3, /* 2b */ 4,
    /* 2c */ 3, /* 2d */ 4, /* 2e */ 4, /* 2f */ 5,
    /* 30 */ 2, /* 31 */ 3, /* 32 */ 3, /* 33 */ 4,
    /* 34 */ 3, /* 35 */ 4, /* 36 */ 4, /* 37 */ 5,
    /* 38 */ 3, /* 39 */ 4, /* 3a */ 4, /* 3b */ 5,
    /* 3c */ 4, /* 3d */ 5, /* 3e */ 5, /* 3f */ 6,
    /* 40 */ 1, /* 41 */ 2, /* 42 */ 2, /* 43 */ 3,
    /* 44 */ 2, /* 45 */ 3, /* 46 */ 3, /* 47 */ 4,
    /* 48 */ 2, /* 49 */ 3, /* 4a */ 3, /* 4b */ 4,
    /* 4c */ 3, /* 4d */ 4, /* 4e */ 4, /* 4f */ 5,
    /* 50 */ 2, /* 51 */ 3, /* 52 */ 3, /* 53 */ 4,
    /* 54 */ 3, /* 55 */ 4, /* 56 */ 4, /* 57 */ 5,
    /* 58 */ 3, /* 59 */ 4, /* 5a */ 4, /* 5b */ 5,
    /* 5c */ 4, /* 5d */ 5, /* 5e */ 5, /* 5f */ 6,
    /* 60 */ 2, /* 61 */ 3, /* 62 */ 3, /* 63 */ 4,
    /* 64 */ 3, /* 65 */ 4, /* 66 */ 4, /* 67 */ 5,
    /* 68 */ 3, /* 69 */ 4, /* 6a */ 4, /* 6b */ 5,
    /* 6c */ 4, /* 6d */ 5, /* 6e */ 5, /* 6f */ 6,
    /* 70 */ 3, /* 71 */ 4, /* 72 */ 4, /* 73 */ 5,
    /* 74 */ 4, /* 75 */ 5, /* 76 */ 5, /* 77 */ 6,
    /* 78 */ 4, /* 79 */ 5, /* 7a */ 5, /* 7b */ 6,
    /* 7c */ 5, /* 7d */ 6, /* 7e */ 6, /* 7f */ 7,
    /* 80 */ 1, /* 81 */ 2, /* 82 */ 2, /* 83 */ 3,
    /* 84 */ 2, /* 85 */ 3, /* 86 */ 3, /* 87 */ 4,
    /* 88 */ 2, /* 89 */ 3, /* 8a */ 3, /* 8b */ 4,
    /* 8c */ 3, /* 8d */ 4, /* 8e */ 4, /* 8f */ 5,
    /* 90 */ 2, /* 91 */ 3, /* 92 */ 3, /* 93 */ 4,
    /* 94 */ 3, /* 95 */ 4, /* 96 */ 4, /* 97 */ 5,
    /* 98 */ 3, /* 99 */ 4, /* 9a */ 4, /* 9b */ 5,
    /* 9c */ 4, /* 9d */ 5, /* 9e */ 5, /* 9f */ 6,
    /* a0 */ 2, /* a1 */ 3, /* a2 */ 3, /* a3 */ 4,
    /* a4 */ 3, /* a5 */ 4, /* a6 */ 4, /* a7 */ 5,
    /* a8 */ 3, /* a9 */ 4, /* aa */ 4, /* ab */ 5,
    /* ac */ 4, /* ad */ 5, /* ae */ 5, /* af */ 6,
    /* b0 */ 3, /* b1 */ 4, /* b2 */ 4, /* b3 */ 5,
    /* b4 */ 4, /* b5 */ 5, /* b6 */ 5, /* b7 */ 6,
    /* b8 */ 4, /* b9 */ 5, /* ba */ 5, /* bb */ 6,
    /* bc */ 5, /* bd */ 6, /* be */ 6, /* bf */ 7,
    /* c0 */ 2, /* c1 */ 3, /* c2 */ 3, /* c3 */ 4,
    /* c4 */ 3, /* c5 */ 4, /* c6 */ 4, /* c7 */ 5,
    /* c8 */ 3, /* c9 */ 4, /* ca */ 4, /* cb */ 5,
    /* cc */ 4, /* cd */ 5, /* ce */ 5, /* cf */ 6,
    /* d0 */ 3, /* d1 */ 4, /* d2 */ 4, /* d3 */ 5,
    /* d4 */ 4, /* d5 */ 5, /* d6 */ 5, /* d7 */ 6,
    /* d8 */ 4, /* d9 */ 5, /* da */ 5, /* db */ 6,
    /* dc */ 5, /* dd */ 6, /* de */ 6, /* df */ 7,
    /* e0 */ 3, /* e1 */ 4, /* e2 */ 4, /* e3 */ 5,
    /* e4 */ 4, /* e5 */ 5, /* e6 */ 5, /* e7 */ 6,
    /* e8 */ 4, /* e9 */ 5, /* ea */ 5, /* eb */ 6,
    /* ec */ 5, /* ed */ 6, /* ee */ 6, /* ef */ 7,
    /* f0 */ 4, /* f1 */ 5, /* f2 */ 5, /* f3 */ 6,
    /* f4 */ 5, /* f5 */ 6, /* f6 */ 6, /* f7 */ 7,
    /* f8 */ 5, /* f9 */ 6, /* fa */ 6, /* fb */ 7,
    /* fc */ 6, /* fd */ 7, /* fe */ 7, /* ff */ 8
};

const uint64_t MASK_POWER_16[17] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};

const uint8_t rev_MSB[16] = {0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15};
const uint8_t rev_LSB[16] = {0, 64, 128, 192, 16, 80, 144, 208, 32, 96, 160, 224, 48, 112, 176, 240};
const uint8_t rev[256] = {0, 64, 128, 192, 16, 80, 144, 208, 32, 96, 160, 224, 48, 112, 176, 240, 4,
                            68, 132, 196, 20, 84, 148, 212, 36, 100, 164, 228, 52, 116, 180, 244, 8,
                            72, 136, 200, 24, 88, 152, 216, 40, 104, 168, 232, 56, 120, 184, 248, 12,
                            76, 140, 204, 28, 92, 156, 220, 44, 108, 172, 236, 60, 124, 188, 252, 1,
                            65, 129, 193, 17, 81, 145, 209, 33, 97, 161, 225, 49, 113, 177, 241, 5,
                            69, 133, 197, 21, 85, 149, 213, 37, 101, 165, 229, 53, 117, 181, 245, 9,
                            73, 137, 201, 25, 89, 153, 217, 41, 105, 169, 233, 57, 121, 185, 249, 13,
                            77, 141, 205, 29, 93, 157, 221, 45, 109, 173, 237, 61, 125, 189, 253, 2,
                            66, 130, 194, 18, 82, 146, 210, 34, 98, 162, 226, 50, 114, 178, 242, 6,
                            70, 134, 198, 22, 86, 150, 214, 38, 102, 166, 230, 54, 118, 182, 246, 10,
                            74, 138, 202, 26, 90, 154, 218, 42, 106, 170, 234, 58, 122, 186, 250, 14,
                            78, 142, 206, 30, 94, 158, 222, 46, 110, 174, 238, 62, 126, 190, 254, 3,
                            67, 131, 195, 19, 83, 147, 211, 35, 99, 163, 227, 51, 115, 179, 243, 7,
                            71, 135, 199, 23, 87, 151, 215, 39, 103, 167, 231, 55, 119, 183, 247, 11,
                            75, 139, 203, 27, 91, 155, 219, 43, 107, 171, 235, 59, 123, 187, 251, 15,
                            79, 143, 207, 31, 95, 159, 223, 47, 111, 175, 239, 63, 127, 191, 255};

extern CC* createCC(int nb_bits_bf);
extern void initiateCC(CC* cc, int nb_bits_bf);
extern void freeCC(CC* cc, int lvl_cc, info_per_level* restrict info_per_level);
extern void freeNode(Node* restrict node, int lvl_node, info_per_level* restrict info_per_lvl);

extern int count_nodes(CC* cc, int start, int end, uint8_t type);
extern int count_children(CC* cc, int start, int end, uint8_t type);

extern void initializeUC(UC* uc);

extern uint8_t reverse_word_8(uint8_t v);
extern int popcnt_8(uint8_t v);
extern int popcnt_8_par(const uint8_t* v, int start, int end);

extern uint16_t hash1(uint8_t* kmer);
extern uint16_t hash2_bis(uint8_t* kmer);

extern uint16_t dbj2(uint8_t c1, uint8_t c2, int size_bf);
extern uint16_t sdbm(uint8_t c1, uint8_t c2, int size_bf);

extern UC_SIZE_ANNOT_T *min_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot);
extern int max_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot);
extern int size_annot_sub(uint8_t* annot, int size_substring, int size_annot);

/* ---------------------------------------------------------------------------------------------------------------
*  transform2CC(uc, cc, size_suffix, info_per_lvl)
*  ---------------------------------------------------------------------------------------------------------------
*  Create a compressed container from an uncompressed container
*  ---------------------------------------------------------------------------------------------------------------
*  uc: pointer to a non-empty uncompressed container
*  cc: pointer to an empty compressed container
*  size_suffix: size of the suffixes in uc
*  info_per_lvl: info_per_level structure used to manipulate the compressed container field cc->children_type
*  ---------------------------------------------------------------------------------------------------------------
*/
void transform2CC(UC* restrict uc, CC* restrict cc, Root* root, int lvl_cc, int size_suffix){

    ASSERT_NULL_PTR(uc,"transform2CC()")
    ASSERT_NULL_PTR(cc,"transform2CC()")

    UC* uc_tmp;

    int i=0, k=-1, z=0, pos = 0;

    int size_annotation = uc->size_annot*sizeof(uint8_t);
    int nb_substrings_different = 0;
    int pos_start_cluster = 0;

    int nb_cell = root->info_per_lvl[lvl_cc].size_kmer_in_bytes;
    int size_line_uc = nb_cell + uc->size_annot;

    int nb_cell_children = root->info_per_lvl[lvl_cc].size_kmer_in_bytes_minus_1;
    int size_line_uc_children;

    int begin_clust = 0;
    int it_children = 0;
    int it_clust = -1;

    int current_posFilter2 = -2;
    int posFilter2 = -1;
    int nb_1 = 0;
    int cmp = 0;
    int pos_skp_children = -1;
    int real_pos_i_substring = 0;
    int real_pos_i_uc = 0;

    uint8_t bit_new_suf;
    uint8_t mask;
    uint8_t type = (cc->type >> 6) & 0x1;

    uint32_t substring_prefix;

    uint16_t hash1_v = 0, hash2_v = 0;
    uint16_t size_bf = cc->type >> 7;

    int* new_order;

    uint8_t* substrings = malloc(root->info_per_lvl[lvl_cc].nb_kmers_uc * SIZE_BYTES_SUF_PREF * sizeof(uint8_t));
    ASSERT_NULL_PTR(substrings,"transform2CC()")

    uint8_t* sub_equal = calloc(root->info_per_lvl[lvl_cc].nb_kmers_uc, sizeof(uint8_t));
    ASSERT_NULL_PTR(sub_equal,"transform2CC()")

    UC_SIZE_ANNOT_T *annot_sizes = min_size_per_sub(uc->suffixes, root->info_per_lvl[lvl_cc].nb_kmers_uc, nb_cell, uc->size_annot);
    uint8_t** annot_extend = get_extend_annots(uc, nb_cell, root->info_per_lvl[lvl_cc].nb_kmers_uc, 0, root->info_per_lvl[lvl_cc].nb_kmers_uc-1);
    uint8_t** annot_cplx = get_annots_cplx_nodes(uc, nb_cell, root->info_per_lvl[lvl_cc].nb_kmers_uc, 0, root->info_per_lvl[lvl_cc].nb_kmers_uc-1);
    UC_SIZE_ANNOT_CPLX_T *annot_cplx_sizes = min_size_per_annot_cplx_sub(uc, root->info_per_lvl[lvl_cc].nb_kmers_uc, nb_cell, 0, root->info_per_lvl[lvl_cc].nb_kmers_uc-1);

    if (root->info_per_lvl[lvl_cc].level_min == 0) mask = 0x7f;
    else mask = 0xff;

    for (i=0; i<root->info_per_lvl[lvl_cc].nb_kmers_uc; i++){

        if ((annot_extend != NULL) && (annot_extend[i] != NULL) && (annot_extend[i][0] != 0)) annot_sizes[i] = uc->size_annot+1;

        pos = i*SIZE_BYTES_SUF_PREF;
        memcpy(&(substrings[pos]), &(uc->suffixes[i*size_line_uc]), SIZE_BYTES_SUF_PREF);

        substrings[pos] = reverse_word_8(substrings[pos]);
        substrings[pos+1] = reverse_word_8(substrings[pos+1]);
        substrings[pos+2] = reverse_word_8(substrings[pos+2]) & 0xc0;

        if (root->compressed > 0){

            substring_prefix = (substrings[pos] << 10) | (substrings[pos+1] << 2) | (substrings[pos+2] >> 6);

            //hash1_v = KnutDivision(substring_prefix, root->info_per_lvl[lvl_cc].modulo_hash);
            //hash2_v = MutiplicationMethod(substring_prefix, root->info_per_lvl[lvl_cc].modulo_hash);

            hash1_v = root->hash_v[substring_prefix * 2] % root->info_per_lvl[lvl_cc].modulo_hash;
            hash2_v = root->hash_v[substring_prefix * 2 + 1] % root->info_per_lvl[lvl_cc].modulo_hash;

            cc->BF_filter2[hash1_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash1_v%SIZE_BITS_UINT_8T];
            cc->BF_filter2[hash2_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash2_v%SIZE_BITS_UINT_8T];
        }
        else {

            substring_prefix = ((substrings[pos] & 0x3f) << SIZE_BITS_UINT_8T) | substrings[pos+1];

            hash1_v = root->hash_v[substring_prefix * 2] % root->info_per_lvl[lvl_cc].modulo_hash;
            hash2_v = root->hash_v[substring_prefix * 2 + 1] % root->info_per_lvl[lvl_cc].modulo_hash;

            /*substring_prefix = substrings[pos] & 0x3f;

            hash1_v = dbj2(substrings[pos+1], substring_prefix, root->info_per_lvl[lvl_cc].modulo_hash);
            hash2_v = sdbm(substrings[pos+1], substring_prefix, root->info_per_lvl[lvl_cc].modulo_hash);*/

            cc->BF_filter2[hash1_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash1_v%SIZE_BITS_UINT_8T];
            cc->BF_filter2[hash2_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash2_v%SIZE_BITS_UINT_8T];

            substring_prefix = (substrings[pos] << SIZE_BITS_UINT_8T) | substrings[pos+1];

            substrings[pos] = substring_prefix >> 6;
            substrings[pos+1] = (substring_prefix << 2) | (substrings[pos+2] >> 6);
            substrings[pos+2] = (substring_prefix >> 8) & 0xc0;
        }
    }

    new_order = quicksort_init(substrings, SIZE_BYTES_SUF_PREF, 0, root->info_per_lvl[lvl_cc].nb_kmers_uc-1);

    for (i=0; i<root->info_per_lvl[lvl_cc].nb_kmers_uc; i++){

        if ((i!=0) && (memcmp(&(substrings[i*SIZE_BYTES_SUF_PREF]), &(substrings[(i-1)*SIZE_BYTES_SUF_PREF]), SIZE_BYTES_SUF_PREF) == 0)){
            sub_equal[i] = 1;
            continue;
        }
        else{
            it_clust++;
            nb_substrings_different++;
        }

        if (it_clust == root->info_per_lvl[lvl_cc].nb_ucs_skp){

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

    if (root->info_per_lvl[lvl_cc].level_min == 1){
        cc->extra_filter3 = calloc(CEIL(nb_substrings_different, SIZE_BITS_UINT_8T),sizeof(uint8_t));
        ASSERT_NULL_PTR(cc->extra_filter3,"transform2CC()")
    }

    //skp_children is not allocated if the CC is a leaf because in that case, there are exactly NB_UC_PER_SKP
    //annotations (suffix length = 0) in every array of CC->children
    if (size_suffix != NB_CHAR_SUF_PREF) allocate_children_type(cc, nb_substrings_different);

    for (i=0; i<root->info_per_lvl[lvl_cc].nb_kmers_uc; i++){
        real_pos_i_substring = i*SIZE_BYTES_SUF_PREF;
        real_pos_i_uc = new_order[i]*size_line_uc;

        //If the prefix of this suffix is the same as the suffix before, no need to update filter2 or filter3
        //We just update to number of suffixes stored in an array of CC->children
        if ((real_pos_i_substring!=0) && (sub_equal[i] == 1)) type = addNewElt(cc, k, nb_substrings_different, type);
        else {

            k++;
            pos_skp_children = k / root->info_per_lvl[lvl_cc].nb_ucs_skp;
            uc_tmp = &(((UC*)cc->children)[pos_skp_children]);

            pos_start_cluster = uc_tmp->nb_children;

            //Updates filter2 and filter3
            posFilter2 = (((uint16_t)substrings[real_pos_i_substring]) << 2) | (((uint16_t)substrings[real_pos_i_substring+1]) >> 6);
            cc->BF_filter2[size_bf+posFilter2/SIZE_BITS_UINT_8T] |= MASK_POWER_8[posFilter2%SIZE_BITS_UINT_8T];
            cc->filter3[k] = (substrings[real_pos_i_substring+1] << 2) | (substrings[real_pos_i_substring+2] >> 6);

            //If the prefix inserted shares the same p_u as the one inserted before, we update extra_filter3
            if (posFilter2 != current_posFilter2) nb_1++;

            //Update SkipFilter3 when a multiple of 248 unique prefix were inserted
            if (k % root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter3 == root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter3 - 1){
                int new_size = (SIZE_FILTER2_DEFAULT/SIZE_BITS_UINT_8T)+((k+1)/root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter3);
                cc->BF_filter2 = realloc(cc->BF_filter2, (size_bf+new_size)*sizeof(uint8_t));
                cc->BF_filter2[size_bf+new_size-1] = (uint8_t)nb_1;
                nb_1 = 0;
            }

            if (size_suffix != NB_CHAR_SUF_PREF) type = addNewElt(cc, k, nb_substrings_different, type);
        }

        if (size_suffix == NB_CHAR_SUF_PREF){

            if (annot_sizes[new_order[i]] != 0){

                if ((annot_extend != NULL) && (annot_extend[new_order[i]] != NULL) && (annot_extend[new_order[i]][0] != 0)){

                    memcpy(&(uc_tmp->suffixes[(k % root->info_per_lvl[lvl_cc].nb_ucs_skp) * uc_tmp->size_annot]),
                           &(uc->suffixes[real_pos_i_uc + nb_cell]),
                           (annot_sizes[new_order[i]]-1) * sizeof(uint8_t));

                    uc_tmp->suffixes[(k % root->info_per_lvl[lvl_cc].nb_ucs_skp) * uc_tmp->size_annot + annot_sizes[new_order[i]] - 1] = annot_extend[new_order[i]][0];
                }
                else memcpy(&(uc_tmp->suffixes[(k % root->info_per_lvl[lvl_cc].nb_ucs_skp) * uc_tmp->size_annot]),
                            &(uc->suffixes[real_pos_i_uc + nb_cell]),
                            annot_sizes[new_order[i]] * sizeof(uint8_t));
            }
            else memcpy(&(uc_tmp->suffixes[(k % root->info_per_lvl[lvl_cc].nb_ucs_skp) * uc_tmp->size_annot]),
                        annot_cplx[new_order[i]],
                        annot_cplx_sizes[new_order[i]] * sizeof(uint8_t));

            if (posFilter2 != current_posFilter2){
                current_posFilter2 = posFilter2;
                cc->extra_filter3[k/SIZE_BITS_UINT_8T] |= MASK_POWER_8[k%SIZE_BITS_UINT_8T];
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

            uc->suffixes[real_pos_i_uc+j-1] &= root->info_per_lvl[lvl_cc].mask_shift_kmer;

            memmove(&(uc->suffixes[real_pos_i_uc+j]), &(uc->suffixes[real_pos_i_uc + j + nb_cell_to_delete]), size_annotation);

            size_line_uc_children = nb_cell_children + uc_tmp->size_annot;

            //z = uc_tmp->nb_children * size_line_uc_children;

            bit_new_suf = 0;

            if (pos_start_cluster > uc_tmp->nb_children-1) z = pos_start_cluster;
            else{

                z = binary_search_UC_array(uc_tmp->suffixes, uc_tmp->size_annot, pos_start_cluster, uc_tmp->nb_children-1,
                                           &(uc->suffixes[real_pos_i_uc]), nb_cell_children, mask);

                if (root->info_per_lvl[lvl_cc].level_min == 1){
                    if (memcmp(&(uc->suffixes[real_pos_i_uc]),
                                &(uc_tmp->suffixes[z * size_line_uc_children]),
                                nb_cell_children * sizeof(uint8_t)) > 0) z++;
                }
                else{

                    cmp = memcmp(&(uc->suffixes[real_pos_i_uc]), &(uc_tmp->suffixes[z * size_line_uc_children]), (nb_cell_children-1) * sizeof(uint8_t));

                    if (cmp > 0) z++;
                    else if ((cmp == 0) && (uc->suffixes[real_pos_i_uc + nb_cell_children-1]
                                            > (uc_tmp->suffixes[z * size_line_uc_children + nb_cell_children-1] & 0x7f))) z++;

                    if ((z == pos_start_cluster) && ((uc_tmp->suffixes[z * size_line_uc_children + nb_cell_children - 1] >> 7) == 1)){
                        bit_new_suf = 0x80;
                        uc_tmp->suffixes[z * size_line_uc_children + nb_cell_children - 1] &= 0x7f;
                    }
                }
            }

            memmove(&(uc_tmp->suffixes[(z+1) * size_line_uc_children]),
                    &(uc_tmp->suffixes[z * size_line_uc_children]),
                    (uc_tmp->nb_children - z) * size_line_uc_children * sizeof(uint8_t));

            z *= size_line_uc_children;

            memset(&(uc_tmp->suffixes[z]), 0, size_line_uc_children * sizeof(uint8_t));
            memcpy(&(uc_tmp->suffixes[z]), &(uc->suffixes[real_pos_i_uc]), nb_cell_children * sizeof(uint8_t));

            uc_tmp->suffixes[z + nb_cell_children - 1] |= bit_new_suf;

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
                if (root->info_per_lvl[lvl_cc].level_min == 0) uc_tmp->suffixes[z+nb_cell_children-1] |= 0x80;
                else cc->extra_filter3[k/SIZE_BITS_UINT_8T] |= MASK_POWER_8[k%SIZE_BITS_UINT_8T];
            }
        }
    }

    cc->nb_elem = nb_substrings_different; //Set the number of unique prefixes inserted into the CC
    cc->type = (cc->type & 0xffbf) | (type << 6);

    free(substrings);
    free(new_order);
    free(sub_equal);
    free(annot_sizes);
    free(annot_cplx);
    free(annot_cplx_sizes);
    if (annot_extend != NULL) free(annot_extend);
}

/* ---------------------------------------------------------------------------------------------------------------
*  transform2CC_from_arraySuffix(array_suffix, cc, size_suffix, size_annot, annot_extend, annot_cplx, size_annot_cplx, info_per_lvl)
*  ---------------------------------------------------------------------------------------------------------------
*  Create a CC from an array of suffixes and annotations
*  ---------------------------------------------------------------------------------------------------------------
*  array_suffix: array of suffixes and annotations
*  cc: pointer to the CC in which array_suffix will be transformed
*  size_suffix: size of the suffixes in array_suffix
*  size_annot: size of the annotations in array_suffix
*  info_per_lvl: info_per_level structure used to manipulate cc->children_type
*  ---------------------------------------------------------------------------------------------------------------
*/
void transform2CC_from_arraySuffix(uint8_t* restrict array_suffix, CC* restrict cc, Root* root, int lvl_cc, int size_suffix,
                                   int size_annot, uint8_t** annot_extend, uint8_t** annot_cplx, int size_annot_cplx){

    ASSERT_NULL_PTR(array_suffix,"transform2CC_from_arraySuffix()")
    ASSERT_NULL_PTR(cc,"transform2CC_from_arraySuffix()")

    int nb_cell = CEIL(size_suffix*2, SIZE_BITS_UINT_8T);
    int size_line_uc = nb_cell + size_annot;

    int nb_cell_children = CEIL((size_suffix-NB_CHAR_SUF_PREF)*2, SIZE_BITS_UINT_8T);
    int size_line_uc_children;

    uint16_t size_bf = cc->type >> 7;

    uint8_t mask;
    uint8_t bit_new_suf;
    uint8_t type = (cc->type >> 6) & 0x1;

    uint8_t* substrings = malloc(root->info_per_lvl[lvl_cc].nb_kmers_uc*SIZE_BYTES_SUF_PREF*sizeof(uint8_t));
    ASSERT_NULL_PTR(substrings, "transform2CC_from_arraySuffix()")

    uint8_t* sub_equal = calloc(root->info_per_lvl[lvl_cc].nb_kmers_uc,sizeof(uint8_t));
    ASSERT_NULL_PTR(sub_equal, "transform2CC_from_arraySuffix()")

    UC_SIZE_ANNOT_T *annot_sizes = min_size_per_sub(array_suffix, root->info_per_lvl[lvl_cc].nb_kmers_uc, nb_cell, size_annot);
    ASSERT_NULL_PTR(annot_sizes, "transform2CC_from_arraySuffix()")

    UC_SIZE_ANNOT_CPLX_T *annot_sizes_cplx = calloc(root->info_per_lvl[lvl_cc].nb_kmers_uc,sizeof( UC_SIZE_ANNOT_CPLX_T ));
    ASSERT_NULL_PTR(annot_sizes_cplx, "transform2CC_from_arraySuffix()")

    int cmp = 0;
    int nb_substrings_different = 0;
    int pos_start_cluster = 0;

    int i=0, k=-1, z=0, pos=0;
    uint32_t substring_prefix;
    uint16_t hash1_v = 0, hash2_v = 0;

    if (root->info_per_lvl[lvl_cc].level_min == 0) mask = 0x7f;
    else mask = 0xff;

    for (i=0; i<root->info_per_lvl[lvl_cc].nb_kmers_uc; i++){

        if ((annot_extend != NULL) && (annot_extend[i] != NULL) && (annot_extend[i][0] != 0)) annot_sizes[i] = size_annot+1;
        if (annot_sizes[i] == 0) annot_sizes_cplx[i] = size_annot_sub(annot_cplx[i], 0, size_annot_cplx);

        pos = i*SIZE_BYTES_SUF_PREF;
        memcpy(&(substrings[pos]), &(array_suffix[i*size_line_uc]), SIZE_BYTES_SUF_PREF);

        substrings[pos] = reverse_word_8(substrings[pos]);
        substrings[pos+1] = reverse_word_8(substrings[pos+1]);
        substrings[pos+2] = reverse_word_8(substrings[pos+2]) & 0xc0;

        if (root->compressed > 0){

            substring_prefix = (substrings[pos] << 10) | (substrings[pos+1] << 2) | (substrings[pos+2] >> 6);

            //hash1_v = KnutDivision(substring_prefix, root->info_per_lvl[lvl_cc].modulo_hash);
            //hash2_v = MutiplicationMethod(substring_prefix, root->info_per_lvl[lvl_cc].modulo_hash);

            hash1_v = root->hash_v[substring_prefix * 2] % root->info_per_lvl[lvl_cc].modulo_hash;
            hash2_v = root->hash_v[substring_prefix * 2 + 1] % root->info_per_lvl[lvl_cc].modulo_hash;

            cc->BF_filter2[hash1_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash1_v%SIZE_BITS_UINT_8T];
            cc->BF_filter2[hash2_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash2_v%SIZE_BITS_UINT_8T];
        }
        else{

            /*substring_prefix = substrings[pos] & 0x3f;

            hash1_v = dbj2(substrings[pos+1], substring_prefix, root->info_per_lvl[lvl_cc].modulo_hash);
            hash2_v = sdbm(substrings[pos+1], substring_prefix, root->info_per_lvl[lvl_cc].modulo_hash);*/

            substring_prefix = ((substrings[pos] & 0x3f) << SIZE_BITS_UINT_8T) | substrings[pos+1];

            hash1_v = root->hash_v[substring_prefix * 2] % root->info_per_lvl[lvl_cc].modulo_hash;
            hash2_v = root->hash_v[substring_prefix * 2 + 1] % root->info_per_lvl[lvl_cc].modulo_hash;

            cc->BF_filter2[hash1_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash1_v%SIZE_BITS_UINT_8T];
            cc->BF_filter2[hash2_v/SIZE_BITS_UINT_8T] |= MASK_POWER_8[hash2_v%SIZE_BITS_UINT_8T];

            substring_prefix = (substrings[pos] << SIZE_BITS_UINT_8T) | substrings[pos+1];

            substrings[pos] = substring_prefix >> 6;
            substrings[pos+1] = (substring_prefix << 2) | (substrings[pos+2] >> 6);
            substrings[pos+2] = (substring_prefix >> 8) & 0xc0;
        }
    }

    int* new_order = quicksort_init(substrings, SIZE_BYTES_SUF_PREF, 0, root->info_per_lvl[lvl_cc].nb_kmers_uc-1);

    int begin_clust = 0;
    int it_children = 0;
    int it_clust = -1;

    UC* uc;

    for (i=0; i<root->info_per_lvl[lvl_cc].nb_kmers_uc; i++){

        if ((i!=0) && (memcmp(&(substrings[i*SIZE_BYTES_SUF_PREF]), &(substrings[(i-1)*SIZE_BYTES_SUF_PREF]), SIZE_BYTES_SUF_PREF) == 0)){
            sub_equal[i] = 1;
            continue;
        }
        else{
            it_clust++;
            nb_substrings_different++;
        }

        if (it_clust == root->info_per_lvl[lvl_cc].nb_ucs_skp){

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

    if (root->info_per_lvl[lvl_cc].level_min == 1){
        cc->extra_filter3 = calloc(CEIL(nb_substrings_different, SIZE_BITS_UINT_8T),sizeof(uint8_t));
        ASSERT_NULL_PTR(cc->extra_filter3,"transform2CC_from_arraySuffix()")
    }

    if (size_suffix != NB_CHAR_SUF_PREF) allocate_children_type(cc, nb_substrings_different);

    int current_posFilter2 = -2;
    int posFilter2 = -1;
    int nb_1 = 0;
    int pos_skp_children = -1;
    int real_pos_i_substring = 0;
    int real_pos_i_uc = 0;

    for (i=0; i<root->info_per_lvl[lvl_cc].nb_kmers_uc; i++){
        real_pos_i_substring = i*SIZE_BYTES_SUF_PREF;
        real_pos_i_uc = new_order[i]*size_line_uc;

        if ((real_pos_i_substring != 0) && (sub_equal[i] == 1)) type = addNewElt(cc, k, nb_substrings_different, type);
        else {

            k++;
            pos_skp_children = k / root->info_per_lvl[lvl_cc].nb_ucs_skp;
            uc = &(((UC*)cc->children)[pos_skp_children]);

            pos_start_cluster = uc->nb_children;

            posFilter2 = (((uint16_t)substrings[real_pos_i_substring]) << 2) | (((uint16_t)substrings[real_pos_i_substring+1]) >> 6);
            cc->BF_filter2[size_bf+posFilter2/8] |= MASK_POWER_8[posFilter2%8];
            cc->filter3[k] = (substrings[real_pos_i_substring+1] << 2) | (substrings[real_pos_i_substring+2] >> 6);

            //If the prefix inserted shares the same p_u as the one inserted before, we update extra_filter3
            if (posFilter2 != current_posFilter2) nb_1++;

            if (k % root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter3 == root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter3 - 1){
                int new_size = (SIZE_FILTER2_DEFAULT/SIZE_BITS_UINT_8T)+((k+1)/root->info_per_lvl[lvl_cc].nb_bits_per_cell_skip_filter3);
                cc->BF_filter2 = realloc(cc->BF_filter2, (new_size+size_bf)*sizeof(uint8_t));
                cc->BF_filter2[size_bf+new_size-1] = (uint8_t)nb_1;
                nb_1 = 0;
            }

            if (size_suffix != NB_CHAR_SUF_PREF) type = addNewElt(cc, k, nb_substrings_different, type);
        }

        if (size_suffix == NB_CHAR_SUF_PREF){

            if (annot_sizes[new_order[i]] != 0){

                if ((annot_extend != NULL) && (annot_extend[new_order[i]] != NULL) && (annot_extend[new_order[i]][0] != 0)){

                    memcpy(&(uc->suffixes[(k % root->info_per_lvl[lvl_cc].nb_ucs_skp) * uc->size_annot]),
                           &(array_suffix[real_pos_i_uc + nb_cell]),
                           (annot_sizes[new_order[i]]-1) * sizeof(uint8_t));

                    uc->suffixes[(k % root->info_per_lvl[lvl_cc].nb_ucs_skp) * uc->size_annot + annot_sizes[new_order[i]] - 1] = annot_extend[new_order[i]][0];
                }
                else memcpy(&(uc->suffixes[(k % root->info_per_lvl[lvl_cc].nb_ucs_skp) * uc->size_annot]),
                            &(array_suffix[real_pos_i_uc + nb_cell]),
                            annot_sizes[new_order[i]] * sizeof(uint8_t));
            }
            else memcpy(&(uc->suffixes[(k % root->info_per_lvl[lvl_cc].nb_ucs_skp) * uc->size_annot]),
                        annot_cplx[new_order[i]],
                        annot_sizes_cplx[new_order[i]] * sizeof(uint8_t));

            if (posFilter2 != current_posFilter2){
                current_posFilter2 = posFilter2;
                cc->extra_filter3[k/SIZE_BITS_UINT_8T] |= MASK_POWER_8[k%SIZE_BITS_UINT_8T];
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

            array_suffix[real_pos_i_uc+j-1] &= root->info_per_lvl[lvl_cc].mask_shift_kmer;

            memmove(&(array_suffix[real_pos_i_uc+j]), &(array_suffix[real_pos_i_uc + j + nb_cell_to_delete]), size_annot);

            size_line_uc_children = nb_cell_children + uc->size_annot;

            //z = uc->nb_children * size_line_uc_children;

            bit_new_suf = 0;

            if (pos_start_cluster > uc->nb_children-1) z = pos_start_cluster;
            else{

                z = binary_search_UC_array(uc->suffixes, uc->size_annot, pos_start_cluster, uc->nb_children-1,
                                           &(array_suffix[real_pos_i_uc]), nb_cell_children, mask);

                if (root->info_per_lvl[lvl_cc].level_min == 1){
                    if (memcmp(&(array_suffix[real_pos_i_uc]),
                                &(uc->suffixes[z * size_line_uc_children]),
                                nb_cell_children * sizeof(uint8_t)) > 0) z++;
                }
                else{

                    cmp = memcmp(&(array_suffix[real_pos_i_uc]), &(uc->suffixes[z * size_line_uc_children]), (nb_cell_children-1) * sizeof(uint8_t));

                    if (cmp > 0) z++;
                    else if ((cmp == 0) && (array_suffix[real_pos_i_uc + nb_cell_children-1]
                                            > (uc->suffixes[z * size_line_uc_children + nb_cell_children-1] & 0x7f))) z++;

                    if ((z == pos_start_cluster) && ((uc->suffixes[z * size_line_uc_children + nb_cell_children - 1] >> 7) == 1)){
                        bit_new_suf = 0x80;
                        uc->suffixes[z * size_line_uc_children + nb_cell_children - 1] &= 0x7f;
                    }
                }
            }

            memmove(&(uc->suffixes[(z+1) * size_line_uc_children]),
                    &(uc->suffixes[z * size_line_uc_children]),
                    (uc->nb_children - z) * size_line_uc_children * sizeof(uint8_t));

            z *= size_line_uc_children;

            memset(&(uc->suffixes[z]), 0, size_line_uc_children * sizeof(uint8_t));
            memcpy(&(uc->suffixes[z]), &(array_suffix[real_pos_i_uc]), nb_cell_children * sizeof(uint8_t));

            uc->suffixes[z + nb_cell_children - 1] |= bit_new_suf;

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
                if (root->info_per_lvl[lvl_cc].level_min == 0) uc->suffixes[z+nb_cell_children-1] |= 0x80;
                else cc->extra_filter3[k/SIZE_BITS_UINT_8T] |= MASK_POWER_8[k%SIZE_BITS_UINT_8T];
            }

            uc->nb_children++;
        }
    }

    cc->nb_elem = nb_substrings_different;
    cc->type = (cc->type & 0xffbf) | (type << 6);

    free(substrings);
    free(sub_equal);
    free(new_order);
    free(annot_sizes);
    free(annot_sizes_cplx);
}

/* ---------------------------------------------------------------------------------------------------------------
*  insertSP_CC(pres, sp, size_sp, id_genome, info_per_lvl, ann_inf, annot_sorted)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert a complete suffix prefix into a compressed container.
*  ---------------------------------------------------------------------------------------------------------------
*  pres: ptr to a resultPresence structure, contain information about the suffix prefixes in the container
*  kmer: ptr on the suffix prefix to insert
*  size_sp: size of the suffix prefix to insert, in nb of char.
*  id_genome: genome id of the suffix prefix to insert
*  info_per_lvl: ptr on a info_per_level structure, used to manipulate CCs field children_type
*  ann_inf: ptr on annotation_inform structure, used to make the transition between reading and modifying an annot
*  ---------------------------------------------------------------------------------------------------------------
*/
void insertSP_CC(resultPresence* restrict pres, uint8_t* restrict sp, int size_sp, uint32_t id_genome, int size_id_genome,
                 info_per_level* restrict info_per_lvl, annotation_inform* ann_inf, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(pres,"insertSP_CC()")

    CC* cc = pres->container;

    uint16_t size_bf = cc->type >> 7; //Bloom filter size in bytes

    uint8_t suf = (cc->type >> 1) & 0x1f; //Length v of p_v
    uint8_t pref = NB_CHAR_SUF_PREF*2-suf; //Length u of p_u

    uint8_t type = (cc->type >> 6) & 0x1;

    int pos_extra_filter3 = pres->pos_extra_filter3; //Position where p_v has to be inserted in filter3

    int size_filter2 = size_bf+(MASK_POWER_16[pref]/SIZE_BITS_UINT_8T); // Size BF+filter2 in bytes
    int skip_filter2 = MASK_POWER_16[pref]/info_per_lvl->nb_bits_per_cell_skip_filter2; //Size SkipFilter2 in bytes
    int size_filter2_n_skip = size_filter2;
    if (cc->nb_elem >= info_per_lvl->tresh_suf_pref) size_filter2_n_skip += skip_filter2; // Size BF+filter2+SkipFilter2 in bytes

    int i=0;
    int nb_skp = CEIL(cc->nb_elem, info_per_lvl->nb_ucs_skp);

    int imid, imin, imax;

    int hamming_weight = 0;
    int hamming_weight_0 = 0;

    int j=0, k=0, m=0;

    int size_line_children = 0;
    int pos_children = 0;
    int cpt_pv = 0;
    int cpt_node = 0;
    int nb_elem_in_pv = 0;
    int it = 0;
    int sum = 0;

    int end;

    int nb_cell_3rdlist = CEIL(cc->nb_elem,SIZE_BITS_UINT_8T);

    uint8_t* children;

    uint8_t word_tmp;
    uint8_t suffix = 0;

    UC* uc;

    // If the prefix is not present in filter2
    if (pres->presFilter2 == 0){

        //Compute position to set to 1 in filter2
        int posFilter2 = 0;
        if (pref==10) posFilter2 = (((uint16_t)pres->substring[0]) << 2) | ((uint16_t)(pres->substring[1]) >> 6);
        else posFilter2 = (((uint16_t)pres->substring[0]) << 6) | ((uint16_t)(pres->substring[1]) >> 2);

        int skip_posfilter2 = posFilter2/info_per_lvl->nb_bits_per_cell_skip_filter2;

        //Set the position to 1 in filter2, eventually increment one cell in SkipFilter2 if it exists
        cc->BF_filter2[size_bf+posFilter2/SIZE_BITS_UINT_8T] |= MASK_POWER_8[posFilter2%SIZE_BITS_UINT_8T];

        if ((cc->nb_elem >= info_per_lvl->tresh_suf_pref) && (skip_posfilter2 < skip_filter2))
            cc->BF_filter2[size_filter2+skip_posfilter2]++;

        //Compute the Hamming weight between position 0 and the position set to 1 in filter2
        //To know how many p_u are lexicographically inferior or equal to the one we insert
        if (pos_extra_filter3 == INT_MAX){
            int k = size_bf+posFilter2/SIZE_BITS_UINT_8T;
            int cnt = 0;
            int hamming_weight = 0;

            //If SkipFilter2 exists, we use it to accelerate the computation of the Hamming weight
            if (cc->nb_elem >= info_per_lvl->tresh_suf_pref){

                /*while ((cnt < skip_posfilter2) && (cnt < skip_filter2)){
                    hamming_weight += cc->BF_filter2[size_filter2+cnt];
                    cnt++;
                }*/

                for (uint8_t* it_ptr = &(cc->BF_filter2[size_filter2]); cnt < MIN(skip_posfilter2, skip_filter2); cnt++, it_ptr++)
                    hamming_weight += *it_ptr;
            }

            hamming_weight += popcnt_8_par(cc->BF_filter2, size_bf + cnt * info_per_lvl->nb_bytes_per_cell_skip_filter2, k);

            word_tmp = cc->BF_filter2[k];
            for (k=0; k<=posFilter2%SIZE_BITS_UINT_8T; k++, word_tmp >>= 1) hamming_weight += word_tmp & 1;

            pres->posFilter2 = hamming_weight; //pres->posFilter2 now contains this Hamming weight
        }
    }

    //Compute p_v, the v last char. of the prefix we try to insert
    if (suf==8) suffix = (pres->substring[1] << 2) | (pres->substring[2] >> 6);
    else if (suf==4) suffix = ((pres->substring[1] & 0x3) << 2) | (pres->substring[2] >> 6);

    if (info_per_lvl->level_min == 1){

        if (pos_extra_filter3 == INT_MAX){

            //This loop tries to use SkipFilter3 (if it exists) to accelerate the Hamming weight comput.
            //by skipping counting 1s directly in extra_filter3

            /*while(m < cc->nb_elem / info_per_lvl->nb_bits_per_cell_skip_filter3){
                if ((sum = hamming_weight + cc->BF_filter2[size_filter2_n_skip+m]) < pres->posFilter2){
                    hamming_weight = sum;
                    m++;
                }
                else break;
            }*/

            while ((m < cc->nb_elem / info_per_lvl->nb_bits_per_cell_skip_filter3)
                   && ((hamming_weight += cc->BF_filter2[size_filter2_n_skip+m]) < pres->posFilter2)) m++;

            if (hamming_weight >= pres->posFilter2) hamming_weight -= cc->BF_filter2[size_filter2_n_skip+m];

            for (k = m * info_per_lvl->nb_bytes_per_cell_skip_filter3; k < nb_cell_3rdlist; k++){

                if ((sum = hamming_weight+popcnt_8(cc->extra_filter3[k])) >= pres->posFilter2){
                    word_tmp = cc->extra_filter3[k];

                    int size_word = 7;
                    if (k == nb_cell_3rdlist-1) size_word = (cc->nb_elem-1)%8;

                    for (j=0; j<=size_word; j++, word_tmp >>= 1){
                        if ((hamming_weight += word_tmp & 1) == pres->posFilter2){
                            pos_extra_filter3 = k*8+j;
                            j++;
                            break;
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

                /*int a=0;
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
                }*/

                for (sum = k * SIZE_BITS_UINT_8T, word_tmp = cc->extra_filter3[k] >> j; sum < cc->nb_elem; sum++, word_tmp >>= 1){
                    if (sum%SIZE_BITS_UINT_8T == 0) word_tmp = cc->extra_filter3[sum/SIZE_BITS_UINT_8T];
                    if (IS_ODD(word_tmp)) break;
                    else hamming_weight_0++;
                }

                /*a = MIN(SIZE_BITS_UINT_8T, cc->nb_elem - k * SIZE_BITS_UINT_8T);
                word_tmp = cc->extra_filter3[k] >> j;

                while ((j < a) && ((word_tmp & 1) == 0)){
                    hamming_weight_0++;
                    j++;
                    word_tmp >>= 1;
                }

                if (j == SIZE_BITS_UINT_8T){

                    for (a = k+1; a < nb_cell_3rdlist-1; a++){
                        if (cc->extra_filter3[a] == 0) hamming_weight_0 += SIZE_BITS_UINT_8T;
                        else break;
                    }

                    word_tmp = cc->extra_filter3[a];
                    j = a * SIZE_BITS_UINT_8T;

                    while ((j < cc->nb_elem) && ((word_tmp & 1) == 0)){
                        hamming_weight_0++;
                        j++;
                        word_tmp >>= 1;
                    }
                }*/

                //We compare p_vs in the filter3 (in the cluster) to determine where exactly to insert our p_v
                if (pos_extra_filter3 < cc->nb_elem){

                    imin = pos_extra_filter3;
                    imax = pos_extra_filter3 + hamming_weight_0;

                    if (suf==8){

                        while (imin < imax){
                            imid = (imin + imax)/2;

                            if (cc->filter3[imid] < suffix) imin = imid + 1;
                            else imax = imid;
                        }

                        if (cc->filter3[imin] < suffix) pres->posFilter3 = imin+1;
                        else pres->posFilter3 = imin;
                    }
                    else {

                        uint8_t tmp = 0;

                        while (imin < imax){
                            imid = (imin + imax)/2;

                            if (IS_ODD(imid)) tmp = cc->filter3[imid/2] >> 4;
                            else tmp = cc->filter3[imid/2] & 0xf;

                            if (tmp < suffix) imin = imid + 1;
                            else imax = imid;
                        }

                        if (IS_ODD(imin)) tmp = cc->filter3[imin/2] >> 4;
                        else tmp = cc->filter3[imin/2] & 0xf;

                        if (tmp < suffix) pres->posFilter3 = imin+1;
                        else pres->posFilter3 = imin;
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

            /*int sum = 0;
            while(m < cc->nb_elem/info_per_lvl->nb_bits_per_cell_skip_filter3){
                if ((sum = hamming_weight + cc->BF_filter2[size_filter2_n_skip+m]) < pres->posFilter2){
                    hamming_weight = sum;
                    m++;
                }
                else break;
            }*/
            while ((m < cc->nb_elem / info_per_lvl->nb_bits_per_cell_skip_filter3)
                   && ((hamming_weight += cc->BF_filter2[size_filter2_n_skip+m]) < pres->posFilter2)) m++;

            if (hamming_weight >= pres->posFilter2) hamming_weight -= cc->BF_filter2[size_filter2_n_skip+m];

            k = m * info_per_lvl->nb_bits_per_cell_skip_filter3;
            pos_children = k/info_per_lvl->nb_ucs_skp;

            //count_Nodes_Children(pres->container, pos_children * info_per_lvl->nb_ucs_skp, k, &cpt_pv, &cpt_node);
            //cpt_node += count_nodes(pres->container, 0, pos_children * info_per_lvl->nb_ucs_skp);

            nb_elem_in_pv = pos_children * info_per_lvl->nb_ucs_skp;

            if (k <= cc->nb_elem - k){
                count_Nodes_Children(cc, nb_elem_in_pv, k, &cpt_pv, &cpt_node, type);
                cpt_node += count_nodes(cc, 0, nb_elem_in_pv, type);
            }
            else{
                end = MIN(cc->nb_elem - nb_elem_in_pv, info_per_lvl->nb_ucs_skp);

                if (k - nb_elem_in_pv <= nb_elem_in_pv + end - k) cpt_pv = count_children(cc, nb_elem_in_pv, k, type);
                else cpt_pv = (((UC*)cc->children)[pos_children]).nb_children - count_children(cc, k, nb_elem_in_pv + end, type);

                cpt_node = cc->nb_Node_children - count_nodes(cc, k, cc->nb_elem, type);
            }

            nb_elem_in_pv = 0;
            it = k - pos_children * info_per_lvl->nb_ucs_skp;

            pres->count_children = cpt_pv;
            pres->count_nodes = cpt_node;
            pres->pos_children = pos_children;

            while (pres->pos_children < nb_skp){

                uc = &(((UC*)cc->children)[pres->pos_children]);
                size_line_children = info_per_lvl->size_kmer_in_bytes_minus_1+uc->size_annot;
                children = &(uc->suffixes[info_per_lvl->size_kmer_in_bytes_minus_1-1]);

                if (pres->pos_children == nb_skp - 1) end = cc->nb_elem - pres->pos_children * info_per_lvl->nb_ucs_skp;
                else end = info_per_lvl->nb_ucs_skp;

                while (it < end){

                    if ((nb_elem_in_pv = getNbElts(cc, k, type)) == 0){
                        hamming_weight += cc->children_Node_container[pres->count_nodes].UC_array.nb_children & 0x1;
                        pres->count_nodes++;
                    }
                    else{
                        hamming_weight += children[pres->count_children*size_line_children] >> 7;
                        pres->count_children += nb_elem_in_pv;
                    }

                    if (hamming_weight == pres->posFilter2){
                        pos_extra_filter3 = k;
                        goto MATCH2;
                    }

                    it++;
                    k++;
                }

                it = 0;
                pres->count_children = 0;
                pres->pos_children++;
            }
        }
        else{
            pos_children = pres->pos_children;
            cpt_pv = pres->count_children;
            cpt_node = pres->count_nodes;
            k = INT_MAX;
        }

        MATCH2: if (pos_extra_filter3 == INT_MAX) pos_extra_filter3 = cc->nb_elem;

        if (pres->posFilter3 == INT_MAX){

            if ((pres->presFilter2 == 1) && (pos_extra_filter3 != cc->nb_elem)){

                if (pos_children < nb_skp){

                    pos_children = pres->pos_children;
                    cpt_pv = pres->count_children;
                    cpt_node = pres->count_nodes;

                    if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * info_per_lvl->nb_ucs_skp;
                    else end = info_per_lvl->nb_ucs_skp;

                    k++;
                    it++;

                    if (it >= end){
                        it = 0;
                        cpt_pv = 0;
                        pos_children++;
                    }

                    while (pos_children < nb_skp){

                        uc = &(((UC*)cc->children)[pos_children]);
                        size_line_children = info_per_lvl->size_kmer_in_bytes_minus_1 + uc->size_annot;
                        children = &(uc->suffixes[info_per_lvl->size_kmer_in_bytes_minus_1-1]);

                        if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * info_per_lvl->nb_ucs_skp;
                        else end = info_per_lvl->nb_ucs_skp;

                        while (it < end){

                            if ((nb_elem_in_pv = getNbElts(cc, k, type)) == 0){
                                if (IS_ODD(cc->children_Node_container[cpt_node].UC_array.nb_children)) goto OUT_LOOP2;
                                else hamming_weight_0++;;
                                cpt_node++;
                            }
                            else{
                                if (IS_ODD(children[cpt_pv*size_line_children] >> 7)) goto OUT_LOOP2;
                                else hamming_weight_0++;
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

                    imin = pos_extra_filter3;
                    imax = pos_extra_filter3 + hamming_weight_0;

                    if (suf==8){

                        while (imin < imax){
                            imid = (imin + imax)/2;

                            if (cc->filter3[imid] < suffix) imin = imid + 1;
                            else imax = imid;
                        }

                        if (cc->filter3[imin] < suffix) pres->posFilter3 = imin+1;
                        else pres->posFilter3 = imin;
                    }
                    else {

                        uint8_t tmp = 0;

                        while (imin < imax){
                            imid = imin + (imax-imin)/2;

                            if (IS_ODD(imid)) tmp = cc->filter3[imid/2] >> 4;
                            else tmp = cc->filter3[imid/2] & 0xf;

                            if (tmp < suffix) imin = imid + 1;
                            else imax = imid;
                        }

                        if (IS_ODD(imin)) tmp = cc->filter3[imin/2] >> 4;
                        else tmp = cc->filter3[imin/2] & 0xf;

                        if (tmp < suffix) pres->posFilter3 = imin+1;
                        else pres->posFilter3 = imin;
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
            if (pres->posFilter3 < cc->nb_elem) transform_one_2_replace = 1;
            to_insert_in_extralist3 = 1;
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
        /*for (i=cc->nb_elem; i>pres->posFilter3; i--){
            if (IS_ODD(i)) cc->filter3[i/2] <<= 4;
            else cc->filter3[i/2] |= cc->filter3[(i-1)/2] >> 4;
        }

        if (IS_ODD(i)) cc->filter3[i/2] = (cc->filter3[i/2] & 0xf) | (suffix << 4);
        else cc->filter3[i/2] = (cc->filter3[i/2] & 0xf0) | (suffix & 0xf);*/

        i = cc->nb_elem/2;

        if (IS_EVEN(cc->nb_elem) && (pres->posFilter3 != cc->nb_elem)) {
            cc->filter3[i] = cc->filter3[i-1] >> 4;
            i--;
        }

        for (; i>pres->posFilter3/2; i--) cc->filter3[i] = (cc->filter3[i] << 4) | (cc->filter3[i-1] >> 4);

        if (IS_ODD(pres->posFilter3)) cc->filter3[i] = (suffix << 4) | (cc->filter3[i] & 0xf);
        else cc->filter3[i] = (cc->filter3[i] << 4) | (suffix & 0xf);

    }

    //Realloc extra_filter3 if needed
    if (info_per_lvl->level_min == 1){
        if (nb_cell_3rdlist*SIZE_BITS_UINT_8T == cc->nb_elem){
            nb_cell_3rdlist++;
            cc->extra_filter3 = realloc(cc->extra_filter3, nb_cell_3rdlist*sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->extra_filter3,"insertSP_CC()")
            cc->extra_filter3[nb_cell_3rdlist-1] = 0;
        }

        //Shift and modify the bits after/before the position of insertion in the extra_filter3
        if (nb_cell_3rdlist > 1){
            for (i=nb_cell_3rdlist-1; i>(pres->posFilter3/SIZE_BITS_UINT_8T); i--)
                cc->extra_filter3[i] = (cc->extra_filter3[i] << 1) | (cc->extra_filter3[i-1] >> 7);
        }

        i=pres->posFilter3/SIZE_BITS_UINT_8T;

        if (pres->posFilter3%SIZE_BITS_UINT_8T == 7){
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

    compute_best_mode(ann_inf, annot_sorted, NULL, 0, NULL, 0, id_genome, size_id_genome); //Compute the size required by the annotation

    int pos_skp = pres->posFilter3/info_per_lvl->nb_ucs_skp;
    int pos_sub = 0;

    //We shift the kmer: Its prefix p has been inserted in the filter2, filter3 and extra_filter3 so
    //we only need to keep only the suffix s to insert it into CC->children
    if (size_sp != NB_CHAR_SUF_PREF){
        int j=0;
        int nb_cell_to_delete = 2;

        if (size_sp == 45) nb_cell_to_delete++;

        for (j=0; j<info_per_lvl->size_kmer_in_bytes - nb_cell_to_delete; j++){
            sp[j] = sp[j+2] >> 2;
            if (j+3 < info_per_lvl->size_kmer_in_bytes) sp[j] |= sp[j+3] << 6;
        }

        sp[j-1] &= info_per_lvl->mask_shift_kmer;

        int count = 0;

        //pos_sub is the position of insertion in CC->children[pos_skp]
        //count is the number of suffixes stored between the position of insertion in CC->children[pos_skp] and the end of this array
        uc = &(((UC*)cc->children)[pos_skp]);

        if (pres->posFilter3 < cc->nb_elem){
            count = count_children(cc, pres->posFilter3, MIN((pos_skp+1) * info_per_lvl->nb_ucs_skp, cc->nb_elem), type);
            pos_sub = uc->nb_children - count;
        }
        else if (pos_skp <= (CEIL(cc->nb_elem, info_per_lvl->nb_ucs_skp)-1)) pos_sub = uc->nb_children;

        //Shift suffixes and their annotation in CC->children to insert our suffix
        add_skp_children(cc, pres->posFilter3, pos_sub, count, info_per_lvl->size_kmer_in_bytes_minus_1,
                         ann_inf->min_size, info_per_lvl);

        uc = &(((UC*)cc->children)[pos_skp]);

        uint8_t* ext_annot = get_extend_annot(uc, info_per_lvl->size_kmer_in_bytes_minus_1, uc->nb_children, pos_sub);
        if (ext_annot != NULL) memset(ext_annot, 0, sizeof(uint8_t));

        //Insert our suffix at the position computed previously
        size_line_children = pos_sub * (info_per_lvl->size_kmer_in_bytes_minus_1 + uc->size_annot) + info_per_lvl->size_kmer_in_bytes_minus_1;
        memcpy(&(uc->suffixes[size_line_children - info_per_lvl->size_kmer_in_bytes_minus_1]), sp, info_per_lvl->size_kmer_in_bytes_minus_1);
        memset(&(uc->suffixes[size_line_children]), 0, uc->size_annot);

        //Create and insert the annotation next to the suffix
        modify_mode_annotation(ann_inf, &(uc->suffixes[size_line_children]), uc->size_annot, ext_annot, 1, id_genome, size_id_genome);

        if ((ext_annot != NULL) && (ext_annot[0] == 0))
            delete_extend_annots(uc, info_per_lvl->size_kmer_in_bytes_minus_1, uc->nb_children, pos_sub, pos_sub, 0, 0, 1);

        //Shift the counters in children_type to insert the new one
        realloc_and_int_children_type(cc, cc->nb_elem, pres->posFilter3, type);

        if (info_per_lvl->level_min == 0){

            uc->suffixes[size_line_children-1] |= to_insert_in_extralist3 << 7;

            if (transform_one_2_replace == 1){

                if (getNbElts(cc, pres->posFilter3+1, type) == 0){
                    //cc->children_Node_container[count_nodes(cc, 0, pres->posFilter3+1)].UC_array.nb_children &= 0xfffe;

                    if (pres->posFilter3+1 < cc->nb_elem - pres->posFilter3 - 1){
                        cc->children_Node_container[count_nodes(cc, 0, pres->posFilter3+1, type)].UC_array.nb_children &= 0xfffe;
                    }
                    else cc->children_Node_container[cc->nb_Node_children - count_nodes(cc, pres->posFilter3+1, cc->nb_elem, type)].UC_array.nb_children &= 0xfffe;
                }
                else if (((pres->posFilter3+1)%info_per_lvl->nb_ucs_skp) == 0){
                    ((UC*)cc->children)[pos_skp+1].suffixes[info_per_lvl->size_kmer_in_bytes_minus_1-1] &= 0x7f;
                }
                else uc->suffixes[size_line_children+info_per_lvl->size_kmer_in_bytes_minus_1+uc->size_annot-1] &= 0x7f;
            }
        }
    }
    else {
        //prefix p has been inserted and the length of the suffix is 0 so CC->children only contains annotations
        pos_sub = pres->posFilter3 % info_per_lvl->nb_ucs_skp;

        int nb_elem = add_skp_annotation(cc, pres->posFilter3, ann_inf->min_size, info_per_lvl);

        uc = &(((UC*)cc->children)[pos_skp]);

        uint8_t* ext_annot = get_extend_annot(uc, 0, nb_elem, pos_sub);

        if (ext_annot != NULL) ext_annot[0] = 0;
        memset(&(uc->suffixes[pos_sub * uc->size_annot]), 0, uc->size_annot);

        modify_mode_annotation(ann_inf, &(uc->suffixes[pos_sub * uc->size_annot]), uc->size_annot, ext_annot, 1, id_genome, size_id_genome);

        if ((ext_annot != NULL) && (ext_annot[0] == 0)) delete_extend_annots(uc, 0, nb_elem, pos_sub, pos_sub, 0, 0, 1);
    }

    //Recompute SkipFilter3
    if (info_per_lvl->level_min == 1){

        for (i = size_filter2_n_skip + pres->posFilter3/info_per_lvl->nb_bits_per_cell_skip_filter3,
             j=(pres->posFilter3/info_per_lvl->nb_bits_per_cell_skip_filter3 + 1) * info_per_lvl->nb_bytes_per_cell_skip_filter3;
                i < size_filter2_n_skip + cc->nb_elem/info_per_lvl->nb_bits_per_cell_skip_filter3;
                i++, j += info_per_lvl->nb_bytes_per_cell_skip_filter3){

            cc->BF_filter2[i] += tmp;
            if (j < nb_cell_3rdlist) cc->BF_filter2[i] -= (tmp = cc->extra_filter3[j] & 0x1);
        }

        cc->nb_elem += 1;

        //Recompute SkipFilter3
        if ((cc->nb_elem%info_per_lvl->nb_bits_per_cell_skip_filter3 == 0) && (cc->nb_elem >= info_per_lvl->nb_bits_per_cell_skip_filter3)){

            int new_cell = size_filter2_n_skip+(cc->nb_elem/info_per_lvl->nb_bits_per_cell_skip_filter3);
            cc->BF_filter2 = realloc(cc->BF_filter2, new_cell*sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->BF_filter2,"insertSP_CC()")
            cc->BF_filter2[new_cell-1] = popcnt_8_par(cc->extra_filter3, nb_cell_3rdlist - info_per_lvl->nb_bytes_per_cell_skip_filter3, nb_cell_3rdlist);
        }
    }
    else {
        int end = cc->nb_elem/info_per_lvl->nb_bits_per_cell_skip_filter3;
        int last_pos_node = 0;
        int last_it_uc = -1;

        int it_uc;
        int last_tmp_uc;

        uint8_t* children = &(cc->BF_filter2[size_filter2_n_skip]);

        cpt_node = 0;

        for (i=pres->posFilter3/info_per_lvl->nb_bits_per_cell_skip_filter3; i<end; i++){

            children[i] += tmp;
            tmp = (i+1) * info_per_lvl->nb_bits_per_cell_skip_filter3;

            if (i+1 <= end){

                if (getNbElts(cc, tmp, type) == 0){
                    cpt_node += count_nodes(cc, last_pos_node, tmp, type);
                    last_pos_node = tmp;
                    tmp = cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1;
                }
                else{

                    if (info_per_lvl->nb_ucs_skp != info_per_lvl->nb_bits_per_cell_skip_filter3){
                        it_uc = tmp/info_per_lvl->nb_ucs_skp;
                        last_tmp_uc = tmp;

                        if (it_uc != last_it_uc){
                            last_it_uc = it_uc;
                            uc = &(((UC*)cc->children)[it_uc]);
                            size_line_children = info_per_lvl->size_kmer_in_bytes_minus_1 + uc->size_annot;
                            cpt_pv = count_children(cc, it_uc * info_per_lvl->nb_ucs_skp, tmp, type);
                        }
                        else cpt_pv += count_children(cc, last_tmp_uc, tmp, type);

                        tmp = uc->suffixes[cpt_pv * size_line_children + info_per_lvl->size_kmer_in_bytes_minus_1-1] >> 7;
                    }
                    else tmp = ((UC*)cc->children)[i+1].suffixes[info_per_lvl->size_kmer_in_bytes_minus_1-1] >> 7;
                }

                cc->BF_filter2[size_filter2_n_skip+i] -= tmp;
            }
        }

        cc->nb_elem += 1;

        //Recompute SkipFilter3
        if ((cc->nb_elem%info_per_lvl->nb_bits_per_cell_skip_filter3 == 0) && (cc->nb_elem >= info_per_lvl->nb_bits_per_cell_skip_filter3)){

            int new_cell = size_filter2_n_skip + cc->nb_elem / info_per_lvl->nb_bits_per_cell_skip_filter3;

            cc->BF_filter2 = realloc(cc->BF_filter2, new_cell*sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->BF_filter2,"insertSP_CC()")

            tmp = new_cell-1;
            cc->BF_filter2[tmp] = 0;

            k = cc->nb_elem - info_per_lvl->nb_bits_per_cell_skip_filter3;
            cpt_pv = -1;
            cpt_node = -1;

            last_tmp_uc = -1;

            it = 0;

            while (it < info_per_lvl->nb_bits_per_cell_skip_filter3){

                pos_children = k/info_per_lvl->nb_ucs_skp;

                if (last_tmp_uc != pos_children){

                    last_tmp_uc = pos_children;
                    uc = &(((UC*)cc->children)[k/info_per_lvl->nb_ucs_skp]);

                    size_line_children = info_per_lvl->size_kmer_in_bytes_minus_1 + uc->size_annot;

                    if (info_per_lvl->nb_ucs_skp == info_per_lvl->nb_bits_per_cell_skip_filter3) cpt_pv = 0;
                    else cpt_pv = count_children(cc, pos_children * info_per_lvl->nb_ucs_skp, k, type);

                    children = &(uc->suffixes[cpt_pv * size_line_children + info_per_lvl->size_kmer_in_bytes_minus_1-1]);
                }

                if ((nb_elem_in_pv = getNbElts(cc, k, type)) == 0){

                    if (cpt_node == -1) cpt_node = cc->nb_Node_children - count_nodes(cc, k, cc->nb_elem, type);

                    cc->BF_filter2[tmp] += cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1;
                    cpt_node++;
                }
                else{
                    cc->BF_filter2[tmp] += children[cpt_pv * size_line_children] >> 7;
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
*  transform_Filter2n3(cc, pref_size, suf_size, info_per_lvl)
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
void transform_Filter2n3(CC* cc, int pref_size, int suf_size, info_per_level* restrict info_per_lvl){

    ASSERT_NULL_PTR(cc,"transform_Filter2n3()")
    ASSERT_NULL_PTR(info_per_lvl,"transform_Filter2n3()")

    uint16_t it_filter2 = 0;
    uint16_t size_bf = cc->type >> 7;

    uint8_t current_pref_size = (NB_CHAR_SUF_PREF*2)-((cc->type >> 1) & 0x1f);
    uint8_t type = (cc->type >> 6) & 0x1;

    int it_filter3_tmp = 0;
    int skip_filter_2 = MASK_POWER_16[pref_size]/info_per_lvl->nb_bits_per_cell_skip_filter2;
    int skip_filter_3 = cc->nb_elem/info_per_lvl->nb_bits_per_cell_skip_filter3;
    int size_filter2 = size_bf+(MASK_POWER_16[pref_size]/SIZE_BITS_UINT_8T);
    int size_filter2_n_skip = size_filter2+skip_filter_2;

    //Allocation of memory for the new filter 2, we can re-write on the filter 3 and extra filter 3
    //The order of suffixes stays the same in the Filter 3
    uint8_t* filter2_tmp = calloc(size_filter2_n_skip+skip_filter_3, sizeof(uint8_t));
    ASSERT_NULL_PTR(filter2_tmp,"transform_Filter2n3()")

    if (suf_size==4){
        if (info_per_lvl->level_min == 1){
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

                            if (new_posFilter2/info_per_lvl->nb_bits_per_cell_skip_filter2 < skip_filter_2)
                                filter2_tmp[size_filter2+new_posFilter2/info_per_lvl->nb_bits_per_cell_skip_filter2]++;

                            if (it_filter3_tmp/info_per_lvl->nb_bits_per_cell_skip_filter3 < skip_filter_3)
                                filter2_tmp[size_filter2_n_skip+it_filter3_tmp/info_per_lvl->nb_bits_per_cell_skip_filter3]++;

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
            int nb_skp = CEIL(cc->nb_elem, info_per_lvl->nb_ucs_skp);

            uint8_t* children;

            UC* uc;

            for (it_filter2=0; it_filter2<MASK_POWER_16[current_pref_size]; it_filter2++){
                //If the position is 0, no cluster of suffix associated
                if ((cc->BF_filter2[size_bf+it_filter2/8] & (MASK_POWER_8[it_filter2%8])) != 0){

                    first_bit = 1;

                    while (pos_children < nb_skp){

                        uc = &(((UC*)cc->children)[pos_children]);
                        size_line_children = info_per_lvl->size_kmer_in_bytes_minus_1+uc->size_annot;
                        children = &(uc->suffixes[info_per_lvl->size_kmer_in_bytes_minus_1-1]);

                        //it_children = 0;
                        it_children = it_filter3_tmp % info_per_lvl->nb_ucs_skp;

                        if (pos_children == nb_skp - 1) end = cc->nb_elem - pos_children * info_per_lvl->nb_ucs_skp;
                        else end = info_per_lvl->nb_ucs_skp;

                        while (it_children < end){

                            if ((nb_elem_in_pv = getNbElts(cc, it_filter3_tmp, type)) == 0){

                                if (((cc->children_Node_container[cpt_node].UC_array.nb_children & 0x1) == 0) || (first_bit == 1)){

                                    if (first_bit == 1) first_bit=0;

                                    new_posFilter2 = (int)((((uint16_t)it_filter2) << 4) | ((uint16_t)(cc->filter3[it_filter3_tmp] >> 4)));

                                    if ((filter2_tmp[size_bf+new_posFilter2/8] & (MASK_POWER_8[new_posFilter2%8])) == 0){

                                        filter2_tmp[size_bf+new_posFilter2/8] |= MASK_POWER_8[new_posFilter2%8];

                                        if (new_posFilter2/info_per_lvl->nb_bits_per_cell_skip_filter2 < skip_filter_2)
                                            filter2_tmp[size_filter2+new_posFilter2/info_per_lvl->nb_bits_per_cell_skip_filter2]++;

                                        if (it_filter3_tmp/info_per_lvl->nb_bits_per_cell_skip_filter3 < skip_filter_3)
                                            filter2_tmp[size_filter2_n_skip+it_filter3_tmp/info_per_lvl->nb_bits_per_cell_skip_filter3]++;

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

                                        if (new_posFilter2/info_per_lvl->nb_bits_per_cell_skip_filter2 < skip_filter_2)
                                            filter2_tmp[size_filter2+new_posFilter2/info_per_lvl->nb_bits_per_cell_skip_filter2]++;

                                        if (it_filter3_tmp/info_per_lvl->nb_bits_per_cell_skip_filter3 < skip_filter_3)
                                            filter2_tmp[size_filter2_n_skip+it_filter3_tmp/info_per_lvl->nb_bits_per_cell_skip_filter3]++;

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

        if (IS_ODD(cc->type)) cc->type = (size_bf << 7) | (type << 6) | 0x9; //xxxxxxxx|000100|01 -> s=4
        else cc->type = (size_bf << 7) | (type << 6) | 0x8; //xxxxxxxx|000100|00 -> s=4

        memcpy(filter2_tmp, cc->BF_filter2, size_bf);
        free(cc->BF_filter2);
        cc->BF_filter2 = filter2_tmp;
    }
}

/* ---------------------------------------------------------------------------------------------------------------
*  add_skp_annotation(cc, position_type, size_annot)
*  ---------------------------------------------------------------------------------------------------------------
*  Same as add_skp_children() but only for CCs of type CC9. In these CCs, the suffix length is 0 so a suffix prefix can
*  only be linked to a single annotation. Thus, arrays in CC->children contain exactly NB_UC_PER_SKP annot.
*  ---------------------------------------------------------------------------------------------------------------
*  cc: ptr to a CC
*  position_type: position in cc->filter3 where was inserted the prefix
*  size_annot: annotation size, in bytes
*  ---------------------------------------------------------------------------------------------------------------
*/
int add_skp_annotation(CC* cc, int position_type, int size_annot, info_per_level* restrict info_per_level){

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
    int nb_cell_skp = CEIL(cc->nb_elem+1, info_per_level->nb_ucs_skp);

    uint8_t* cplx_annot_z_minus1 = NULL;
    uint8_t* extend_annot_z_minus1 = NULL;

    int start = position_type/info_per_level->nb_ucs_skp;

    if (cc->nb_elem%info_per_level->nb_ucs_skp == 0){
        cc->children = realloc(cc->children, nb_cell_skp*sizeof(UC));
        ASSERT_NULL_PTR(cc->children, "add_skp_annotation()")

        uc = &(((UC*)cc->children)[nb_cell_skp-1]);
        initializeUC(uc);
        uc->size_annot = 1;
    }

    for (z=nb_cell_skp-1; z>=start; z--){

        uc = &(((UC*)cc->children)[z]);

        if (z == nb_cell_skp-1) nb_children = cc->nb_elem % info_per_level->nb_ucs_skp;
        else nb_children = info_per_level->nb_ucs_skp - 1;

        if (z == start){

            position_child = position_type % info_per_level->nb_ucs_skp;

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

            start_pos_z_minus1 = info_per_level->nb_ucs_skp - 1;

            cplx_annot_z_minus1 = get_annot_cplx_nodes(uc_tmp, 0, info_per_level->nb_ucs_skp, start_pos_z_minus1);
            max_size_cplx_z_minus1 = size_annot_sub(cplx_annot_z_minus1, 0, uc_tmp->size_annot_cplx_nodes);

            if (max_size_cplx_z_minus1 != 0) cpt_cplx = 1;
            else cpt_cplx = 0;

            if (max_size_cplx_z_minus1 > uc->size_annot_cplx_nodes) increase_size_annot_cplx_nodes(uc, 0, nb_children, max_size_cplx_z_minus1, 1);
            else if (max_size_cplx_z_minus1 < uc->size_annot_cplx_nodes){
                max_size_z = max_size_annot_cplx_sub(&(uc->suffixes[nb_children * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
                                                        uc->nb_cplx_nodes, uc->size_annot_cplx_nodes, 0, nb_children-1);

                if (max_size_z < max_size_cplx_z_minus1) decrease_size_annot_cplx_nodes(uc, 0, nb_children, max_size_z);
            }

            if ((extend_annot_z_minus1 = get_extend_annot(uc_tmp, 0, info_per_level->nb_ucs_skp, start_pos_z_minus1)) == NULL){
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

            delete_extend_annots(uc_tmp, 0, info_per_level->nb_ucs_skp, start_pos_z_minus1, start_pos_z_minus1, 0, 1, 0);
            delete_annot_cplx_nodes(uc_tmp, 0, info_per_level->nb_ucs_skp, start_pos_z_minus1, start_pos_z_minus1, 1, 1, 0);

            extend_annot_z_minus1 = NULL;
            cplx_annot_z_minus1 = NULL;
        }
    }

    return 0;
}

/* ---------------------------------------------------------------------------------------------------------------
*  create_info_per_level(size_min, size_max)
*  ---------------------------------------------------------------------------------------------------------------
*  Creates an array of info_per_level for each possible suffix length between NB_CHAR_SUF_PREF and k (the kmers length).
*  Each cell is a set of pointers on functions which manipulate the field children_type of CCs.
*  ---------------------------------------------------------------------------------------------------------------
*  size_min: minimum size of a suffix (in char.)
*  size_max: maximum size of a suffix (in char.)
*  ---------------------------------------------------------------------------------------------------------------
*/
info_per_level* create_info_per_level(int size_max){

    info_per_level* ptr = NULL;

    int nb_sizes, i;

    for (i = NB_CHAR_SUF_PREF, nb_sizes = 0; i <= size_max; i += NB_CHAR_SUF_PREF, nb_sizes++){

        ptr = realloc(ptr, (nb_sizes+1)*sizeof(info_per_level));
        ASSERT_NULL_PTR(ptr,"create_info_per_level()")

        ptr[nb_sizes].nb_ucs_skp = NB_UC_PER_SKP;
        ptr[nb_sizes].modulo_hash = MODULO_HASH;
        ptr[nb_sizes].tresh_suf_pref = TRESH_SUF_PREF;
        ptr[nb_sizes].nb_bits_per_cell_skip_filter2 = NB_BITS_IN_CELL_SKIP_FILTER2;
        ptr[nb_sizes].nb_bits_per_cell_skip_filter3 = NB_BITS_IN_CELL_SKIP_FILTER3;
        ptr[nb_sizes].nb_bytes_per_cell_skip_filter2 = NB_BYTES_IN_CELL_SKIP_FILTER2;
        ptr[nb_sizes].nb_bytes_per_cell_skip_filter3 = NB_BYTES_IN_CELL_SKIP_FILTER3;
        ptr[nb_sizes].size_kmer_in_bytes = CEIL(i*2,SIZE_BITS_UINT_8T);
        ptr[nb_sizes].size_kmer_in_bytes_minus_1 = MAX(0,CEIL((i-NB_CHAR_SUF_PREF)*2,SIZE_BITS_UINT_8T));

        if (i == size_max){
            ptr[nb_sizes].level_min = 1;
            ptr[nb_sizes].root = 1;
        }
        else{
            ptr[nb_sizes].level_min = 0;
            ptr[nb_sizes].root = 0;
        }

        switch(i){
            case 63: {
                ptr[nb_sizes].nb_kmers_uc = NB_KMERS_PER_UC63;
                ptr[nb_sizes].mask_shift_kmer = 0xf;
                break;
            }
            case 54: {
                ptr[nb_sizes].nb_kmers_uc = NB_KMERS_PER_UC54;
                ptr[nb_sizes].mask_shift_kmer = 0x3;
                break;
            }
            case 45: {
                ptr[nb_sizes].nb_kmers_uc = NB_KMERS_PER_UC45;
                ptr[nb_sizes].mask_shift_kmer = 0xff;
                ptr[nb_sizes].level_min = 1;
                break;
            }
            case 36: {
                ptr[nb_sizes].nb_kmers_uc = NB_KMERS_PER_UC36;
                ptr[nb_sizes].mask_shift_kmer = 0x3f;
                break;
            }
            case 27: {
                ptr[nb_sizes].nb_kmers_uc = NB_KMERS_PER_UC27;
                ptr[nb_sizes].mask_shift_kmer = 0xf;
                break;
            }
            case 18: {
                ptr[nb_sizes].nb_kmers_uc = NB_KMERS_PER_UC18;
                ptr[nb_sizes].mask_shift_kmer = 0x3;
                break;
            }
            case 9: {
                ptr[nb_sizes].nb_kmers_uc = NB_KMERS_PER_UC9;
                ptr[nb_sizes].mask_shift_kmer = 0xff;
                ptr[nb_sizes].level_min = 1;
                break;
            }
            default: ERROR("create_info_per_level()")
        }
    }

    return ptr;
}

/* ---------------------------------------------------------------------------------------------------------------
*  add_skp_children(obj, position_type, position_child, count_before_child, size_substrings, size_annot)
*  ---------------------------------------------------------------------------------------------------------------
*  Prepare insertion of a new suffix in the field CC->children corresponding to a newly inserted prefix in the CC.
*  CC->children and its content is eventually reallocated and suffixes already inserted a eventually shifted.
*  ---------------------------------------------------------------------------------------------------------------
*  obj: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  position_type: position in CC->filter3 where was inserted the prefix
*  position_child: position of the suffix to insert in the array of insertion in CC->children
*  count_before_child: number of suffixes between the end the array of insertion in CC->children and
*                      position_child (position_child included)
*  size_substrings: size of the suffixe to insert
*  size_annot: size of the annotation attached to the suffix to insert
*  ---------------------------------------------------------------------------------------------------------------
*/

void add_skp_children(CC* cc, int position_type, int position_child, int count_before_child, int size_substrings,
                      int size_annot, info_per_level* restrict info_per_lvl){
    ASSERT_NULL_PTR(cc, "add_skp_children()")
    ASSERT_NULL_PTR(info_per_lvl, "add_skp_children()")


    UC* uc;
    UC* uc_tmp;

    uint8_t** extend_annot_z_minus1 = NULL;
    uint8_t** cplx_annot_z_minus1 = NULL;
    uint8_t* tab;

    uint8_t type = (cc->type >> 6) & 0x1;

    int i=0, z=0;
    int count2push = 0;
    int count2delete = 0;
    int max_size_z = 0;
    int max_size_z_minus1 = 0;
    int max_size_cplx_z_minus1 = 0;
    int size_line_children = 0;
    int size_line_children_z_minus1 = 0;
    int cpt_annot_extend_z = 0;
    int cpt_annot_cplx_z = 0;
    int start_pos_z_minus1 = 0;

    int nb_cell_skp = CEIL(cc->nb_elem+1, info_per_lvl->nb_ucs_skp);
    int start = position_type/info_per_lvl->nb_ucs_skp;

    if (cc->nb_elem%info_per_lvl->nb_ucs_skp == 0){
        cc->children = realloc(cc->children, nb_cell_skp*sizeof(UC));
        ASSERT_NULL_PTR(cc->children, "add_skp_children()")

        uc = &(((UC*)cc->children)[nb_cell_skp-1]);
        initializeUC(uc);
        uc->size_annot = 1;
    }

    for (z=nb_cell_skp-1; z>=start; z--){

        uc = &(((UC*)cc->children)[z]);
        uc->nb_children -= count2delete;

        if (z == start){

            if (size_annot > uc->size_annot){
                if (uc->nb_children == 0){
                    uc->suffixes = calloc(size_substrings + size_annot, sizeof(uint8_t));
                    uc->size_annot = size_annot;
                }
                else realloc_annotation(uc, size_substrings, uc->nb_children, size_annot, 1, position_child);
            }
            else if (size_annot < uc->size_annot){

                if (uc->nb_extended_annot != 0) max_size_z = uc->size_annot + 1;
                else max_size_z = max_size_per_sub(uc->suffixes, uc->nb_children, size_substrings, uc->size_annot);

                max_size_z = MAX(size_annot, max_size_z);

                if ((max_size_z > 0) && (uc->nb_extended_annot == 0) && (max_size_z < uc->size_annot)){
                    if (uc->nb_children == 0){
                        uc->suffixes = calloc(size_substrings + max_size_z, sizeof(uint8_t));
                        uc->size_annot = max_size_z;
                    }
                    else realloc_annotation(uc, size_substrings, uc->nb_children, max_size_z, 1, position_child);
                }
                else {
                    size_line_children = size_substrings + uc->size_annot;
                    uc->suffixes = realloc(uc->suffixes, ((uc->nb_children+1) * size_line_children + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                                          + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
                    ASSERT_NULL_PTR(uc->suffixes, "add_skp_children()")

                    count_before_child -= count2delete;

                    memmove(&(uc->suffixes[(position_child+1)*size_line_children]),
                            &(uc->suffixes[position_child*size_line_children]),
                            ((count_before_child * size_line_children) + (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT)
                             + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));

                    shift_extended_annot(uc, size_substrings, uc->nb_children+1, position_child);
                }
            }
            else {
                size_line_children = size_substrings + uc->size_annot;
                uc->suffixes = realloc(uc->suffixes, ((uc->nb_children + 1) * size_line_children + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                                      + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
                ASSERT_NULL_PTR(uc->suffixes, "add_skp_children()")

                count_before_child -= count2delete;
                memmove(&(uc->suffixes[(position_child+1)*size_line_children]),
                        &(uc->suffixes[position_child*size_line_children]),
                        ((count_before_child * size_line_children) + (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT)
                         + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));

                shift_extended_annot(uc, size_substrings, uc->nb_children+1, position_child);
            }

            uc->nb_children++;
            shift_annot_cplx_nodes(uc, size_substrings, uc->nb_children, position_child);
        }
        else {
            uc_tmp = &(((UC*)cc->children)[z-1]);

            count2push = getNbElts(cc, z*info_per_lvl->nb_ucs_skp-1, type);

            if (count2push != 0){

                start_pos_z_minus1 = uc_tmp->nb_children - count2push;

                size_line_children_z_minus1 = size_substrings + uc_tmp->size_annot;
                size_line_children = size_substrings + uc->size_annot;

                cplx_annot_z_minus1 = get_annots_cplx_nodes(uc_tmp, size_substrings, uc_tmp->nb_children, start_pos_z_minus1, uc_tmp->nb_children-1);
                max_size_cplx_z_minus1 = max_size_annot_cplx_sub(&(uc_tmp->suffixes[uc_tmp->nb_children * size_line_children_z_minus1 +
                                                                    uc_tmp->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
                                                                    uc_tmp->nb_cplx_nodes, uc_tmp->size_annot_cplx_nodes, start_pos_z_minus1, uc_tmp->nb_children-1);

                if (max_size_cplx_z_minus1 > uc->size_annot_cplx_nodes){
                    increase_size_annot_cplx_nodes(uc, size_substrings, uc->nb_children, max_size_cplx_z_minus1, 1);
                }
                else if (max_size_cplx_z_minus1 < uc->size_annot_cplx_nodes){
                    max_size_z = max_size_annot_cplx_sub(&(uc->suffixes[uc->nb_children * size_line_children + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
                                                            uc->nb_cplx_nodes, uc->size_annot_cplx_nodes, 0, uc->nb_children-1);

                    if (max_size_z <= max_size_cplx_z_minus1) decrease_size_annot_cplx_nodes(uc, size_substrings, uc->nb_children, max_size_z);
                }

                cpt_annot_cplx_z = 0;

                if (cplx_annot_z_minus1 != NULL){
                    for (i = 0; i < count2push; i++) cpt_annot_cplx_z += cplx_annot_z_minus1[i] != NULL;
                }

                if ((extend_annot_z_minus1 = get_extend_annots(uc_tmp, size_substrings, uc_tmp->nb_children, start_pos_z_minus1, uc_tmp->nb_children-1)) != NULL){
                    max_size_z_minus1 = uc_tmp->size_annot+1;
                }
                else max_size_z_minus1 = max_size_per_sub(&(uc_tmp->suffixes[start_pos_z_minus1 * size_line_children_z_minus1]),
                                                            count2push,
                                                            size_substrings,
                                                            uc_tmp->size_annot);

                COMPUTE_ANNOT_EXT:if (max_size_z_minus1 > uc->size_annot+1) realloc_annotation(uc, size_substrings, uc->nb_children, max_size_z_minus1, 0, 0);
                else if (max_size_z_minus1 < uc->size_annot){

                    if (uc->nb_extended_annot != 0) max_size_z = uc->size_annot + 1;
                    else max_size_z = max_size_per_sub(uc->suffixes, uc->nb_children, size_substrings, uc->size_annot);

                    max_size_z = MAX(max_size_z_minus1, max_size_z);

                    if ((max_size_z > 0) && (uc->nb_extended_annot == 0) && (max_size_z < uc->size_annot))
                        realloc_annotation(uc, size_substrings, uc->nb_children, max_size_z, 0, 0);
                }

                size_line_children_z_minus1 = size_substrings + uc_tmp->size_annot;
                size_line_children = size_substrings + uc->size_annot;
                cpt_annot_extend_z = 0;

                if (max_size_z_minus1 > uc->size_annot){
                    if (extend_annot_z_minus1 != NULL){
                        for (i=0; i<count2push; i++) cpt_annot_extend_z += extend_annot_z_minus1[i] != NULL;
                    }
                    else {
                        tab = uc_tmp->suffixes + start_pos_z_minus1 * size_line_children_z_minus1;
                        for (i = size_substrings + max_size_z_minus1 - 1; i < count2push*size_line_children_z_minus1; i += size_line_children_z_minus1)
                            cpt_annot_extend_z += tab[i] != 0;
                    }

                    if ((uc->nb_extended_annot + cpt_annot_extend_z) * SIZE_BYTE_EXT_ANNOT > uc->nb_children + count2push){
                        recopy_back_annot_extend(uc, size_substrings, uc->nb_children);
                        goto COMPUTE_ANNOT_EXT;
                    }
                }

                uc->nb_children += count2push;

                uc->suffixes = realloc(uc->suffixes,
                                    ((uc->nb_children * size_line_children) + (uc->nb_extended_annot + cpt_annot_extend_z) * SIZE_BYTE_EXT_ANNOT
                                     + (uc->nb_cplx_nodes + cpt_annot_cplx_z) * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
                ASSERT_NULL_PTR(uc->suffixes, "add_skp_children()")

                memmove(&(uc->suffixes[count2push*size_line_children]),
                        uc->suffixes,
                        ((uc->nb_children - count2push) * size_line_children + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                         + (uc->nb_cplx_nodes + cpt_annot_cplx_z) * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));
                memset(uc->suffixes, 0, count2push * size_line_children * sizeof(uint8_t));

                if (max_size_z_minus1 > uc->size_annot){
                    if (extend_annot_z_minus1 != NULL){

                        for (i=0; i<count2push; i++){
                             memcpy(&(uc->suffixes[i*size_line_children]),
                                    &(uc_tmp->suffixes[(start_pos_z_minus1+i) * size_line_children_z_minus1]),
                                    size_line_children_z_minus1*sizeof(uint8_t));
                        }

                        for (i=0; i<count2push; i++){
                            if (extend_annot_z_minus1[i] == NULL) shift_extended_annot(uc, size_substrings, uc->nb_children, i);
                            else insert_extend_annot(uc, size_substrings, uc->nb_children, i, extend_annot_z_minus1[i][0], 1);
                        }

                        for (i=0; i<count2push; i++){
                            if ((cplx_annot_z_minus1 == NULL) || (cplx_annot_z_minus1[i] == NULL)) shift_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i);
                            else insert_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i, cplx_annot_z_minus1[i], max_size_cplx_z_minus1, 1);
                        }
                    }
                    else{
                        tab = &(uc_tmp->suffixes[start_pos_z_minus1 * size_line_children_z_minus1]);

                        for (i=0; i<count2push; i++){
                             memcpy(&(uc->suffixes[i*size_line_children]),
                                    &(uc_tmp->suffixes[(start_pos_z_minus1+i)*size_line_children_z_minus1]),
                                    (size_substrings + max_size_z_minus1 - 1) * sizeof(uint8_t));
                        }

                        for (i=0; i<count2push; i++){
                            if (tab[i*size_line_children_z_minus1 + size_substrings + max_size_z_minus1 - 1] != 0){
                                insert_extend_annot(uc, size_substrings, uc->nb_children, i,
                                                    tab[i*size_line_children_z_minus1 + size_substrings + max_size_z_minus1 - 1], 1);
                            }
                            else shift_extended_annot(uc, size_substrings, uc->nb_children, i);
                        }

                        for (i=0; i<count2push; i++){
                            if ((cplx_annot_z_minus1 == NULL) || (cplx_annot_z_minus1[i] == NULL)) shift_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i);
                            else insert_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i, cplx_annot_z_minus1[i], max_size_cplx_z_minus1, 1);
                        }

                    }
                }
                else{

                    for (i=0; i<count2push; i++){
                        if (extend_annot_z_minus1 != NULL){
                             memcpy(&(uc->suffixes[i*size_line_children]),
                                    &(uc_tmp->suffixes[(start_pos_z_minus1+i)*size_line_children_z_minus1]),
                                    size_line_children_z_minus1*sizeof(uint8_t));

                            if (extend_annot_z_minus1[i] != NULL)
                                uc->suffixes[i*size_line_children+size_line_children_z_minus1] = extend_annot_z_minus1[i][0];
                        }
                        else{
                             memcpy(&(uc->suffixes[i*size_line_children]),
                                    &(uc_tmp->suffixes[(start_pos_z_minus1 + i) * size_line_children_z_minus1]),
                                    (size_substrings + max_size_z_minus1) * sizeof(uint8_t));
                        }
                    }

                    for (i=0; i<count2push; i++) shift_extended_annot(uc, size_substrings, uc->nb_children, i);

                    for (i=0; i<count2push; i++){
                        if ((cplx_annot_z_minus1 == NULL) || (cplx_annot_z_minus1[i] == NULL)) shift_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i);
                        else insert_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i, cplx_annot_z_minus1[i], max_size_cplx_z_minus1, 1);
                    }
                }

                delete_extend_annots(uc_tmp, size_substrings, uc_tmp->nb_children , start_pos_z_minus1,  uc_tmp->nb_children-1, 0, 1, 0);
                delete_annot_cplx_nodes(uc_tmp, size_substrings, uc_tmp->nb_children, start_pos_z_minus1, uc_tmp->nb_children-1, 1, 1, 0);

            }

            count2delete = count2push;

            if (extend_annot_z_minus1 != NULL){
                free(extend_annot_z_minus1);
                extend_annot_z_minus1 = NULL;
            }

            if (cplx_annot_z_minus1 != NULL){
                free(cplx_annot_z_minus1);
                cplx_annot_z_minus1 = NULL;
            }
        }
    }

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  build_skip_nodes(node, size_kmer, info_per_lvl)
*  ---------------------------------------------------------------------------------------------------------------
*  ---------------------------------------------------------------------------------------------------------------
*  size_min: minimum size of a suffix (in char.)
*  size_max: maximum size of a suffix (in char.)
*  ---------------------------------------------------------------------------------------------------------------
*/
uint16_t** build_skip_nodes(Node* node){

    ASSERT_NULL_PTR(node,"build_skip_nodes()")
    ASSERT_NULL_PTR(node->CC_array,"build_skip_nodes()")

    CC* cc;

    int node_nb_elem = 0, i = 0, j = 0, nb_skp = 0;

    uint16_t** nodes_skip;

    uint16_t previous_count;

    while ((((CC*)node->CC_array)[node_nb_elem].type & 0x1) == 0) node_nb_elem++;
    node_nb_elem++;

    nodes_skip = malloc(node_nb_elem*sizeof(uint16_t*));
    ASSERT_NULL_PTR(nodes_skip,"build_skip_nodes()")

    for (i=0; i<node_nb_elem; i++){
        previous_count = 0;
        cc = &(((CC*)node->CC_array)[i]);
        nb_skp = cc->nb_elem / SIZE_CLUST_SKIP_NODES;

        nodes_skip[i] = malloc(nb_skp*sizeof(uint16_t));
        ASSERT_NULL_PTR(nodes_skip[i],"build_skip_nodes()")

        for (j=0; j<nb_skp; j++){
            nodes_skip[i][j] = previous_count + count_nodes(cc, j*SIZE_CLUST_SKIP_NODES, (j+1)*SIZE_CLUST_SKIP_NODES, (cc->type >> 6) & 0x1);
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
