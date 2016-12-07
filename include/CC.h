#pragma once

/* ===================================================================================================================================
*  INCLUDES AND DEFINES
*  ===================================================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "default_param.h"
#include "useful_macros.h"
#include "Node.h"

#include "hash.h"
#include "quicksort.h"
#include "popcnt.h"
#include "annotation.h"

extern const uint64_t MASK_POWER_16[19];

/* ===================================================================================================================================
*  STRUCTURES DECLARATION
*  ===================================================================================================================================
*/

// A Compress container (CC) stores prefixes (fixed length) ps of k-mers suffixes contained in a UC.
// A prefix p stored in a CC is itself divided into a prefix p_u and a suffix p_v
typedef struct{
    uint16_t type; //8bits -> size BF in bytes, 5bits -> v(length of p_v), 1bit: full or not, 1bit: 1 if the CC is the last one of its node
    uint16_t nb_elem; //Number of prefix stored in the CC
    uint16_t nb_Node_children; //Number of nodes in children_Node_container
    // The array BF_filter2 contains several arrays. In order:
    // 1 - The Bloom Filter (fixed size) is a bit array used to record inexactly the presence or absence of a set of ps
    // 2 - The Second Filter (fixed size) is a bit array used to record exactly the presence of p_u
    // 3 - The Skip Filter 2 (fixed size) is an array in which each cell records the Hamming weight of 248 bits in the Second Filter
    //      so for example, SkipFilter2[0] = HammingWeight(Filter2[i]) for 0 <= i <= 31, because 31 cells of 8 bits = 248 bits
    //      This array only exists if the CC stores NB_SUBSTRINGS_TRANSFORM prefixes or more.
    // 4 - The Skip Filter 3 is an array in which each cell records the Hamming weight of 248 bits in the Third Filter
    uint8_t*  BF_filter2;
    // The array filter3 stores explicitly the p_v, in lexicographic ascending order of p_u
    uint8_t*  filter3;
    // The bit array extra_filter3 is a "list of cluster position": One or more consecutive p_v in filter3 belong to the same cluster
    // if they share the same p_u in the Second Filter. A 1 at a position in extra_filter3 signal a new cluster for the prefix at the
    // same position in filter3.
    uint8_t*  extra_filter3;
    // The array children_type contains for each prefix stored the number of suffixes linked to it, on 4 bits or 8 bits.
    // It must only be accessed using functions stored in the structure info_per_level because the type of the array is dynamic,
    //depending on the level of the tree. A 0 at a position indicates that the suffixes linked to the prefix at this position
    // in filter3 are stored in a new node. First cell of this array indicates if it is 4 bits mode or 8 bits mode.

    uint8_t* children_type;

    // The suffixes linked to the CC prefixes are not really stored in new UCs in new children nodes.
    // Instead they are compactly stored in the CC. children is a list of arrays. Each of these arrays contains the suffixes for
    // NB_UC_PER_SKP prefixes, so children[0] contains suffixes and annotations for the first NB_UC_PER_SKP prefixes
    // stored in filter3. A suffix is always immediately followed by its annotation. The length of the annotation is stored in
    // the first cell of each array.
    void* children;
    // children_Node_container is an array of the current CC children nodes
    Node*  children_Node_container;
} __attribute__ ((__packed__)) CC;

/* ===================================================================================================================================
*  INLINE FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

inline CC* createCC(int nb_bits_bf);
inline void initiateCC(CC* cc, int nb_bits_bf);
inline void freeCC(CC* cc, int lvl_cc, info_per_level*  info_per_level);

inline void freeNode(Node* node, int lvl_node, info_per_level*  info_per_lvl);

inline BFT_Root* createBFT_Root(int k, int treshold_compression, uint8_t compressed);
inline void initialize_BFT_Root(BFT_Root* root, int k, int treshold_compression, uint8_t compressed);
inline void freeBFT_Root(BFT_Root* root);
inline void free_content_BFT_Root(BFT_Root* root);
inline BFT_Root* copy_BFT_Root(BFT_Root* root_src);

inline void add_genomes_BFT_Root(int nb_files, char** filenames, BFT_Root* root);

inline void allocate_children_type (CC* cc, int nb_elt);
inline int is_child(CC* cc, int position, uint8_t type);
inline int getNbElts(CC* cc, int position, uint8_t type);
inline uint8_t addNewElt(CC* cc, int position, int current_nb_elem, uint8_t type);

inline void realloc_and_int_children_type(CC* cc, int current_nb_elem, int position_insert, uint8_t type);
inline void resetChildrenType (CC* cc, int position, uint8_t type);
inline int count_nodes(CC* cc, int start, int end, uint8_t type);
inline void count_Nodes_Children(CC* cc, int start, int end, int* count_children, int* count_nodes, uint8_t type);
inline int count_children(CC* cc, int start, int end, uint8_t type);

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

void transform2CC(UC*  uc, CC*  cc, BFT_Root* root, int lvl_cc, int size_suffix);
void transform2CC_from_arraySuffix(uint8_t*  array_suffix, CC*  cc, BFT_Root* root, int lvl_cc, int size_suffix,
                                   int size_annot, uint8_t** annot_extend, uint8_t** annot_cplx, int size_annot_cplx);

void insertSP_CC(resultPresence*  pres, BFT_Root* root, int lvl_node_insert,
                 uint8_t*  kmer, int size_sp, uint32_t id_genome, int size_id_genome);

void transform_Filter2n3(CC* cc, int pref_size, int suf_size, info_per_level*  info_per_lvl);

int add_skp_annotation(CC*, int position_type, int size_annot, info_per_level*  info_per_lvl);
void add_skp_children(CC* cc, int position_type, int position_child, int count_before_child, int size_substrings,
                      int size_annot, info_per_level*  info_per_lvl);

info_per_level* create_info_per_level(int size_max);
uint16_t** build_skip_nodes(Node* node);
void free_skip_nodes(Node* node, uint16_t** skp_nodes);

inline CC* createCC(int nb_bits_bf){
    CC* cc = malloc(sizeof(CC));
    ASSERT_NULL_PTR(cc,"createCC()")

    cc->BF_filter2 = calloc((nb_bits_bf/SIZE_BITS_UINT_8T)+(SIZE_FILTER2_DEFAULT/SIZE_BITS_UINT_8T), sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->BF_filter2,"createCC()")

    cc->filter3 = NULL;
    cc->extra_filter3 = NULL;

    cc->children_type = NULL;
    cc->children = NULL;
    cc->children_Node_container = NULL;

    cc->type = (((uint16_t)nb_bits_bf/SIZE_BITS_UINT_8T) << 7) | 0x11; //xxxxxxxx|001000|01 -> s=8
    cc->nb_elem = 0;
    cc->nb_Node_children = 0;

    return cc;
}

inline void initiateCC(CC* cc, int nb_bits_bf){
    ASSERT_NULL_PTR(cc,"createCC()")

    cc->BF_filter2 = calloc((nb_bits_bf/SIZE_BITS_UINT_8T)+(SIZE_FILTER2_DEFAULT/SIZE_BITS_UINT_8T), sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->BF_filter2,"createCC()")

    cc->filter3 = NULL;
    cc->extra_filter3 = NULL;

    cc->children_type = NULL;
    cc->children = NULL;
    cc->children_Node_container = NULL;

    cc->type = (((uint16_t)nb_bits_bf/SIZE_BITS_UINT_8T) << 7) | 0x11; //xxxxxxxx|001000|01 -> s=8
    cc->nb_elem = 0;
    cc->nb_Node_children = 0;

    return;
}

inline void freeCC(CC* cc, int lvl_cc, info_per_level*  info_per_level){

    if (cc == NULL) return;

    free(cc->BF_filter2);
    free(cc->filter3);

    if (cc->extra_filter3 != NULL) free(cc->extra_filter3);
    if (cc->children_type != NULL) free(cc->children_type);

    int i = 0;
    for (i=0; i<CEIL(cc->nb_elem, info_per_level[lvl_cc].nb_ucs_skp); i++)
        if (((UC*)cc->children)[i].suffixes != NULL) free(((UC*)cc->children)[i].suffixes);

    free(cc->children);

    if (cc->children_Node_container != NULL){
        for (i=0; i<cc->nb_Node_children; i++) freeNode(&(cc->children_Node_container[i]), lvl_cc-1, info_per_level);
        free(cc->children_Node_container);
    }

    return;
}

inline void freeNode(Node* node, int lvl_node, info_per_level*  info_per_lvl){

    if (node == NULL) return;

    if (node->CC_array != NULL){

        int i = -1;
        int it;

        do {
            i++;
            it = ((CC*)node->CC_array)[i].type & 0x1;
            freeCC(&(((CC*)node->CC_array)[i]), lvl_node, info_per_lvl);
        }
        while (it == 0);

        free(node->CC_array);
        node->CC_array = NULL;
    }

    if (node->UC_array.suffixes != NULL){
        free(node->UC_array.suffixes);
        node->UC_array.suffixes = NULL;
    }

    return;
}

inline BFT_Root* createBFT_Root(int k, int treshold_compression, uint8_t compressed){

    BFT_Root* root = malloc(sizeof(BFT_Root));
    ASSERT_NULL_PTR(root,"createBFT_Root()")

    initialize_BFT_Root(root, k, treshold_compression, compressed);

    return root;
}

inline void initialize_BFT_Root(BFT_Root* root, int k, int treshold_compression, uint8_t compressed){

    ASSERT_NULL_PTR(root, "initialize_BFT_Root()")

    root->compressed = compressed;
    root->treshold_compression = treshold_compression;
    root->k = k;
    root->nb_genomes = 0;
    root->length_comp_set_colors = 0;
    root->marked = 0;
    root->filenames = NULL;
    root->comp_set_colors = NULL;
    root->skip_sp = NULL;
    root->ann_inf = NULL;

    if (k == 0){
        root->r1 = 0;
        root->r2 = 0;
        root->hash_v = NULL;
        root->info_per_lvl = NULL;
        root->res = NULL;
    }
    else{
        root->r1 = rand();
        while (root->r1 == (root->r2 = rand())){};

        root->hash_v = create_hash_v_array(root->r1, root->r2);
        root->info_per_lvl = create_info_per_level(root->k);
        root->res = create_resultPresence();
    }

    initiateNode(&(root->node));

    return;
}

inline void freeBFT_Root(BFT_Root* root){

    if (root == NULL) return;

    free_content_BFT_Root(root);
    free(root);

    return;
}

inline void free_content_BFT_Root(BFT_Root* root){

    if (root == NULL) return;

    if (root->filenames != NULL){
        for (int i = 0; i < root->nb_genomes; i++) free(root->filenames[i]);
        free(root->filenames);
    }

    if (root->comp_set_colors != NULL) free_annotation_array_elem(&root->comp_set_colors, &root->length_comp_set_colors);

    free_annotation_inform(root->ann_inf);

    freeNode(&(root->node), (root->k / NB_CHAR_SUF_PREF) - 1, root->info_per_lvl);

    free(root->res);
    free(root->info_per_lvl);
    free(root->hash_v);

    return;
}

inline BFT_Root* copy_BFT_Root(BFT_Root* root_src){

    ASSERT_NULL_PTR(root_src,"copy_BFT_Root()")

    BFT_Root* root_dest = malloc(sizeof(BFT_Root));
    ASSERT_NULL_PTR(root_dest,"copy_BFT_Root()")

    memcpy(root_dest, root_src, sizeof(BFT_Root));

    if (root_dest->ann_inf != NULL) root_dest->ann_inf = create_annotation_inform(root_dest->nb_genomes, root_dest->compressed == 1);
    if (root_dest->res != NULL) root_dest->res = create_resultPresence();

    return root_dest;
}

inline void add_genomes_BFT_Root(int nb_files, char** filenames, BFT_Root* root){

    if (nb_files < 0) ERROR("add_genomes_BFT_Root(): the number of genomes to insert cannot be less than 0.\n")

    if (nb_files > 0){

        ASSERT_NULL_PTR(filenames, "add_genomes_BFT_Root()\n")

        if (root->filenames == NULL) root->filenames = malloc(nb_files * sizeof(char*));
        else root->filenames = realloc(root->filenames, (root->nb_genomes + nb_files) * sizeof(char*));

        ASSERT_NULL_PTR(root->filenames, "add_genomes_BFT_Root()\n")

        for (int i = 0; i < nb_files; i++){

            ASSERT_NULL_PTR(filenames[i], "add_genomes_BFT_Root()\n")

            root->filenames[i + root->nb_genomes] = calloc(strlen(filenames[i]) + 1, sizeof(char));
            ASSERT_NULL_PTR(root->filenames[i + root->nb_genomes] , "add_genomes_BFT_Root()\n")

            strcpy(root->filenames[i + root->nb_genomes], filenames[i]);
        }

        root->nb_genomes += nb_files;

        if (root->ann_inf != NULL) free_annotation_inform(root->ann_inf);
        root->ann_inf = create_annotation_inform(root->nb_genomes, root->compressed == 1);
    }

    return;
}

inline void allocate_children_type (CC* cc, int nb_elt){

    ASSERT_NULL_PTR(cc, "allocate_children_type()")

    cc->children_type = calloc(CEIL(nb_elt,2), sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->children_type, "allocate_children_type()")

    return;
}

inline int is_child(CC* cc, int position, uint8_t type){
    ASSERT_NULL_PTR(cc, "isChild()")

    if (type) return cc->children_type[position] != 0;
    else if (IS_ODD(position)) return cc->children_type[position/2] > 0xf;

    return (cc->children_type[position/2] & 0xf) != 0;
}

inline int getNbElts(CC* cc, int position, uint8_t type){

    ASSERT_NULL_PTR(cc, "getNbElts()")

    if (type) return cc->children_type[position];
    else if (IS_ODD(position)) return cc->children_type[position/2] >> 4;

    return cc->children_type[position/2] & 0xf;
}

inline uint8_t addNewElt(CC* cc, int position, int current_nb_elem, uint8_t type){

    ASSERT_NULL_PTR(cc, "addNewElt()")

    if (type == 0){

        if (IS_ODD(position)){

            if ((cc->children_type[position/2] >> 4) + 1 == 0x10){

                uint8_t* tmp_children_type = calloc(current_nb_elem, sizeof(uint8_t));
                ASSERT_NULL_PTR(tmp_children_type, "addNewElt()")

                int i, j;
                for (i = 0, j = 0; j < current_nb_elem/2; i += 2, j++){
                    tmp_children_type[i] = cc->children_type[j] & 0xf;
                    tmp_children_type[i+1] = cc->children_type[j] >> 4;
                }

                if (IS_ODD(current_nb_elem)) tmp_children_type[i] = cc->children_type[j] & 0xf;

                free(cc->children_type);

                cc->children_type = tmp_children_type;
                cc->children_type[position]++;

                return 1;
            }
            else cc->children_type[position/2] += 0x10;
        }
        else{
            if ((cc->children_type[position/2] & 0xf) + 1 == 0x10){

                uint8_t* tmp_children_type = calloc(current_nb_elem, sizeof(uint8_t));
                ASSERT_NULL_PTR(tmp_children_type, "addNewElt()")

                int i, j;
                for (i = 0, j = 0; j < current_nb_elem/2; i += 2, j++){
                    tmp_children_type[i] = cc->children_type[j] & 0xf;
                    tmp_children_type[i+1] = cc->children_type[j] >> 4;
                }

                if (IS_ODD(current_nb_elem)) tmp_children_type[i] = cc->children_type[j] & 0xf;

                free(cc->children_type);

                cc->children_type = tmp_children_type;
                cc->children_type[position]++;

                return 1;
            }
            else cc->children_type[position/2]++;
        }
    }
    else cc->children_type[position]++;

    return type;
}

inline void realloc_and_int_children_type(CC* cc, int current_nb_elem, int position_insert, uint8_t type){

    ASSERT_NULL_PTR(cc, "realloc_and_int_children_type ()")

    if (type == 1){
        cc->children_type = realloc(cc->children_type, (current_nb_elem+1)*sizeof(uint8_t));
        ASSERT_NULL_PTR(cc->children_type, "realloc_and_int_children_type ()")

        memmove(&(cc->children_type[position_insert+1]),
                &(cc->children_type[position_insert]),
                (current_nb_elem-position_insert)*sizeof(uint8_t));

        cc->children_type[position_insert] = 1;
    }
    else{
        if (IS_EVEN(current_nb_elem)){

            cc->children_type = realloc(cc->children_type, (current_nb_elem / 2 + 1) * sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->children_type, "realloc_and_int_children_type ()")

            cc->children_type[current_nb_elem / 2] = 0;
        }

        int i = current_nb_elem;
        for (; i > position_insert; i--){
            if (IS_ODD(i)) cc->children_type[i/2] <<= 4;
            else cc->children_type[i/2] |= cc->children_type[i/2-1] >> 4;
        }

        if (IS_ODD(i)) cc->children_type[i/2] = (cc->children_type[i/2] & 0xf) + 16;
        else cc->children_type[i/2] = (cc->children_type[i/2] & 0xf0) + 1;
    }
}

inline void resetChildrenType (CC* cc, int position, uint8_t type){
    ASSERT_NULL_PTR(cc, "resetChildrenType()")

    if (type) cc->children_type[position] = 0;
    else if (IS_ODD(position)) cc->children_type[position/2] &= 0xf;
    else cc->children_type[position/2] &= 0xf0;

    return;
}

inline int count_nodes(CC* cc, int start, int end, uint8_t type){

    ASSERT_NULL_PTR(cc, "count_Nodes()")

    int count = 0;
    uint8_t* z;

    if (type){
        for (z = cc->children_type + start; z < cc->children_type + end; z++) count += *z == 0;
    }
    else{
        z = cc->children_type + start/2;
        if (IS_ODD(start)) count -= (*z & 0xf) == 0;
        for (; z < cc->children_type + end/2; z++) count += (*z < 0x10) + ((*z & 0xf) == 0);
        if (IS_ODD(end)) count += (*z & 0xf) == 0;
    }

    return count;
}

inline void count_Nodes_Children(CC* cc, int start, int end, int* count_children, int* count_nodes, uint8_t type){
    ASSERT_NULL_PTR(cc, "count_Nodes_Children()")
    ASSERT_NULL_PTR(count_children, "count_Nodes_Children()")
    ASSERT_NULL_PTR(count_nodes, "count_Nodes_Children()")

    uint8_t* z;
    int count_nodes_tmp = 0;
    int count_children_tmp = 0;

    if (type){

        for (z = cc->children_type + start; z < cc->children_type + end; z++){
            count_children_tmp += *z;
            count_nodes_tmp += *z == 0;
        }
    }
    else {

        z = cc->children_type + start/2;

        if (IS_ODD(start)){
            count_nodes_tmp -= (*z & 0xf) == 0;
            count_children_tmp -= *z & 0xf;
        }

        for (; z < cc->children_type + end/2; z++){
            count_nodes_tmp += (*z < 0x10) + ((*z & 0xf) == 0);
            count_children_tmp += (*z >> 4) + (*z & 0xf);
        }

        if (IS_ODD(end)){
            count_nodes_tmp += (*z & 0xf) == 0;
            count_children_tmp += *z & 0xf;
        }
    }

    *count_children = count_children_tmp;
    *count_nodes = count_nodes_tmp;

    return;
}

inline int count_children(CC* cc, int start, int end, uint8_t type){
    ASSERT_NULL_PTR(cc, "count_Children()")

    int count = 0;
    uint8_t* z;

    if (type){
        for (z = cc->children_type + start; z < cc->children_type + end; z++) count += *z;
    }
    else {
        z = cc->children_type + start/2;
        if (IS_ODD(start)) count -= *z & 0xf;
        for (; z < cc->children_type + end/2; z++) count += (*z >> 4) + (*z & 0xf);
        if (IS_ODD(end)) count += *z & 0xf;
    }

    return count;
}
