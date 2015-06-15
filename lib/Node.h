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

#include "./../lib/useful_macros.h"
#include "./../lib/UC.h"

#define SIZE_SEED 9 //Length of prefixes (in characters) stored in CCs
#define NB_CELLS_PER_SEED CEIL(SIZE_SEED*2,SIZE_CELL) //Length of prefixes (in 8bits cells) stored in CCs

#define SIZE_CLUST_SKIP_NODES NB_CHILDREN_PER_SKP //Value max

/* ===================================================================================================================================
*  STRUCTURES DECLARATION
*  ===================================================================================================================================
*/

// Power of 2, up to 2^16
static const uint64_t MASK_POWER_16[17] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};

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
    annotation_array_elem* comp_set_colors;
    char** filenames;
    int nb_genomes;
    int length_comp_set_colors;
    int k;
    int treshold_compression;
} Root;

//resultPresence is a structure produced by presenceKmer(). It contains information about the presence of a prefix p into a given node.
typedef struct{
    void* restrict container; //Ptr to the container (UC or CC) which contain the prefix p or cc->children that contain the substring we are looking for
    void* restrict link_child; //Ptr to the container (Node or uint8_t*) having potentially the suffix linked to the prefix p

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

    uint8_t substring[NB_CELLS_PER_SEED]; // the prefix p

    uint8_t presBF; //Boolean indicating if p is said present in the Bloom Filter
    uint8_t presFilter2; //Boolean indicating if p_u is present in the Second Filter
    uint8_t presFilter3; //Boolean indicating if p_v is present in the Third Filter

} resultPresence;

typedef struct{
    double double1;
    double double2;
} Duo;

/* ===================================================================================================================================
*  INLINE FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

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
inline Root* createRoot(char** filenames, int nb_files, int k, int treshold_compression){

    Root* root = malloc(sizeof(Root));
    ASSERT_NULL_PTR(root,"createRoot()")

    initiateNode(&(root->node));

    root->filenames = filenames;
    root->nb_genomes = nb_files;
    root->comp_set_colors = NULL;
    root->length_comp_set_colors = 0;
    root->k = k;
    root->treshold_compression = treshold_compression;

    return root;
}

inline resultPresence* create_resultPresence(){

    resultPresence* res = calloc(1,sizeof(resultPresence));
    ASSERT_NULL_PTR(res,"create_resultPresence()")

    res->container = NULL;
    res->link_child = NULL;
    res->posFilter2 = INT_MAX;
    res->posFilter3 = INT_MAX;
    res->pos_extra_filter3 = INT_MAX;

    return res;
}

inline void initialize_resultPresence(resultPresence* res){

    ASSERT_NULL_PTR(res,"create_resultPresence()")

    res->container = NULL;
    res->link_child = NULL;
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
