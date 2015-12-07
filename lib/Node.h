#ifndef DEF_NODE
#define DEF_NODE

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

#include <Judy.h>

#include "./../lib/default_param.h"
#include "./../lib/useful_macros.h"
#include "./../lib/UC.h"
#include "./../lib/xxhash.h"

/* ===================================================================================================================================
*  STRUCTURES DECLARATION
*  ===================================================================================================================================
*/

// info_per_level is a structure which contains pointers on macro used to manipulate the field children_type of a CC.
// The structure is used in an array, initialized in create_info_per_level()
typedef struct{
    int nb_ucs_skp;
    int nb_kmers_uc; //Number of prefixes a CC, at a specific level of the tree, can contain.
    int mask_shift_kmer; //Suffixes are encoded in arrays of 8bits cells: mask_shift_kmer covers only the bits used on the last cell
    int size_kmer_in_bytes; //Size of the suffixes represented at a given level of the tree, in bytes
    int size_kmer_in_bytes_minus_1; //Size of the suffixes represented at a given level of the tree, minus the size of the prefixes ps, in bytes
    int nb_bits_per_cell_skip_filter2;
    int nb_bits_per_cell_skip_filter3;
    int nb_bytes_per_cell_skip_filter2;
    int nb_bytes_per_cell_skip_filter3;
    int modulo_hash;
    int tresh_suf_pref;
    int level_min;
    int root;
} info_per_level;

// A node is a list of containers of two types: Compressed Containers (CC) and Uncompressed Containers (UC).
// It can contain 0 or more CCs in CC_array, plus always one UC which can be empty (UC_array.substrings == NULL) or not.
typedef struct {
    void* CC_array;
    UC UC_array;
} __attribute__ ((__packed__)) Node;

// The BFT root is a special node: it contains a node, the id of the last genome inserted and the size in bytes of these ideas
// TODO: Include the file names associated to each genome ID
typedef struct {
    Node node;

    char** filenames;

    annotation_array_elem* comp_set_colors;

    info_per_level* info_per_lvl;

    uint64_t* hash_v;

    uint8_t compressed;

    int k;
    int r1;
    int r2;
    int nb_genomes;
    int length_comp_set_colors;
    int treshold_compression;
} Root;

//resultPresence is a structure produced by presenceKmer(). It contains information about the presence of a prefix p into a given node.
typedef struct{
    void* restrict node;
    void* restrict container; //Ptr to the container (UC or CC) which contain the prefix p or cc->children that contain the substring we are looking for
    void* restrict link_child; //Ptr to the container (Node or uint8_t*) having potentially the suffix linked to the prefix p

    int level_node;
    int pos_container;

    int bucket; //Position of the array containing the suffixes linked to prefix p in children
    int pos_sub_bucket; // Position (in term of suffix+annotation) of the first suffix linked to prefix p in children[bucket]

    int children_type_leaf; //Boolean indicating that container is a leaf of the tree
    int container_is_UC; //Boolean indicating if container is a UC

    int posFilter2; //position of p in filter2 (where it is or where it should be) or size of suffix
    int posFilter3; //position of p in filter3 (where it is or where it should be) or position of match in link_child
    int pos_extra_filter3; //position of p in extra_filter3 (where it is or where it should be) or size of suffix

    int pos_children;
    int count_children;
    int count_nodes;

    uint8_t substring[SIZE_BYTES_SUF_PREF]; // the prefix p

    uint8_t presBF; //Boolean indicating if p is said present in the Bloom Filter
    uint8_t presFilter2; //Boolean indicating if p_u is present in the Second Filter
    uint8_t presFilter3; //Boolean indicating if p_v is present in the Third Filter

} resultPresence;

typedef struct{

    int start_pos_filter2;
    int end_pos_filter2;
    int start_pos_filter3;
    int end_pos_filter3;

    uint8_t* pres_in_filter2;
    uint8_t* pres_in_filter3;

} resultPresenceSeedCC;

typedef struct{
    Node* n;
    resultPresenceSeedCC* res_pres_ccs;
    uint8_t* pres_uc;

    int nb_res_pres_ccs;
} resultPresenceSeed;

/* ===================================================================================================================================
*  INLINE FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

inline uint64_t* create_hash_v_array(int rand_seed1, int rand_seed2){

    int j, nb_bits;
    uint32_t nb_hash_v = pow(4, NB_CHAR_SUF_PREF);

    uint32_t i;

    uint64_t* hash_v = malloc(nb_hash_v * 2 * sizeof(uint64_t));
    ASSERT_NULL_PTR(hash_v, "create_hash_v_array()")

    uint8_t gen_sub[SIZE_BYTES_SUF_PREF];

    for (i = 0; i < nb_hash_v; i++){

        nb_bits = NB_CHAR_SUF_PREF * 2;

        for (j = 0; j < SIZE_BYTES_SUF_PREF; j++){

            nb_bits -= SIZE_BITS_UINT_8T;

            if (nb_bits >= 0) gen_sub[j] = (i >> nb_bits) & 0xff;
            else gen_sub[j] = (i << (-nb_bits)) & 0xff;
        }

        hash_v[i * 2] = XXH64(gen_sub, SIZE_BYTES_SUF_PREF, rand_seed1);
        hash_v[i * 2 + 1] = XXH64(gen_sub, SIZE_BYTES_SUF_PREF, rand_seed2);
    }

    return hash_v;
}

inline void init_resultPresenceSeed(resultPresenceSeed* res_seed){

    ASSERT_NULL_PTR(res_seed,"init_resultPresenceSeed()")

    res_seed->n = NULL;
    res_seed->res_pres_ccs = NULL;
    res_seed->pres_uc = NULL;

    res_seed->nb_res_pres_ccs = 0;
}

inline void init_resultPresenceSeedCC(resultPresenceSeedCC* res_pres_CC){

    ASSERT_NULL_PTR(res_pres_CC,"init_resultPresenceSeedCC()")

    res_pres_CC->start_pos_filter2 = -1;
    res_pres_CC->end_pos_filter2 = -1;
    res_pres_CC->start_pos_filter3 = -1;
    res_pres_CC->end_pos_filter3 = -1;

    res_pres_CC->pres_in_filter2 = NULL;
    res_pres_CC->pres_in_filter3 = NULL;
}

inline void free_resultPresenceSeedCC(resultPresenceSeedCC* res_pres_CC, int free_res_pres_CC_or_not){

    ASSERT_NULL_PTR(res_pres_CC,"free_resultPresenceSeedCC()")

    if (res_pres_CC->pres_in_filter2 != NULL) free(res_pres_CC->pres_in_filter2);
    if (res_pres_CC->pres_in_filter3 != NULL) free(res_pres_CC->pres_in_filter3);

    if (free_res_pres_CC_or_not != 0) free(res_pres_CC);
}

inline void free_resultPresenceSeed(resultPresenceSeed* res_seed, int free_res_seed_or_not){

    ASSERT_NULL_PTR(res_seed,"free_resultPresenceSeed()")

    int i = 0;

    if (res_seed->res_pres_ccs != NULL){
        for (i = 0; i < res_seed->nb_res_pres_ccs; i++) free_resultPresenceSeedCC(&(res_seed->res_pres_ccs[i]), 0);
        free(res_seed->res_pres_ccs);
    }

    if (res_seed->pres_uc != NULL) free(res_seed->pres_uc);
    if (free_res_seed_or_not != 0) free(res_seed);
}

/* ---------------------------------------------------------------------------------------------------------------
*  createNode()
*  ---------------------------------------------------------------------------------------------------------------
*  Allocate and initialize a node
*  ---------------------------------------------------------------------------------------------------------------
*  ---------------------------------------------------------------------------------------------------------------
*/
inline Node* createNode(void){
    Node* node = malloc(sizeof(Node));
    ASSERT_NULL_PTR(node,"createNode()")

    node->CC_array = NULL;
    initializeUC(&(node->UC_array));

    return node;
}

/* ---------------------------------------------------------------------------------------------------------------
*  initiateNode(node)
*  ---------------------------------------------------------------------------------------------------------------
*  Initialize a node
*  ---------------------------------------------------------------------------------------------------------------
*  node: a pointer on a Node structure
*  ---------------------------------------------------------------------------------------------------------------
*/
inline void initiateNode(Node* node){

    ASSERT_NULL_PTR(node,"initiateNode()")

    node->CC_array = NULL;
    initializeUC(&(node->UC_array));

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  createRoot()
*  ---------------------------------------------------------------------------------------------------------------
*  Create the root of a BFT
*  ---------------------------------------------------------------------------------------------------------------
*  ---------------------------------------------------------------------------------------------------------------
*/
inline Root* createRoot(char** filenames, int nb_files, int k, int treshold_compression, uint8_t compressed,
                        int r1, int r2, info_per_level* info_per_lvl){

    Root* root = malloc(sizeof(Root));
    ASSERT_NULL_PTR(root,"createRoot()")

    root->filenames = filenames;
    root->nb_genomes = nb_files;
    root->compressed = compressed;
    root->comp_set_colors = NULL;
    root->length_comp_set_colors = 0;
    root->k = k;
    root->r1 = r1;
    root->r2 = r2;
    root->treshold_compression = treshold_compression;
    root->info_per_lvl = info_per_lvl;

    root->hash_v = create_hash_v_array(root->r1, root->r2);

    initiateNode(&(root->node));

    return root;
}

inline resultPresence* create_resultPresence(){

    resultPresence* res = calloc(1,sizeof(resultPresence));
    ASSERT_NULL_PTR(res,"create_resultPresence()")

    res->node = NULL;
    res->container = NULL;
    res->link_child = NULL;
    res->posFilter2 = INT_MAX;
    res->posFilter3 = INT_MAX;
    res->pos_extra_filter3 = INT_MAX;

    return res;
}

inline void initialize_resultPresence(resultPresence* res){

    ASSERT_NULL_PTR(res,"create_resultPresence()")

    res->node = NULL;
    res->container = NULL;
    res->link_child = NULL;
    res->pos_container = 0;
    res->level_node = 0;
    res->bucket = 0;
    res->pos_sub_bucket = 0;
    res->presBF = 0;
    res->presFilter2 = 0;
    res->presFilter3 = 0;
    res->count_children = 0;
    res->count_nodes = 0;
    res->pos_children = 0;
    res->children_type_leaf = 0;
    res->container_is_UC = 0;
    res->posFilter2 = INT_MAX;
    res->posFilter3 = INT_MAX;
    res->pos_extra_filter3 = INT_MAX;

    return;
}

#endif
