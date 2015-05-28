#ifndef DEF_CC
#define DEF_CC

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

#include "./../lib/useful_macros.h"
#include "./../lib/Node.h"

#include "./../lib/hash.h"
#include "./../lib/quicksort.h"
#include "./../lib/popcnt.h"
#include "./../lib/annotation.h"

#define TRESHOLD_COMPRESSION 0 //if a CC has less than TRESHOLD_COMPRESSION prefixes stored, children_type is compressed
#define NB_SUBSTRINGS_TRANSFORM 3584 //if a CC has exactly NB_SUBSTRINGS_TRANSFORM prefixes, it is recomputed with transform_Filter2n3()
#define MODULO_HASH 1504 //Size bloom filter in bits
#define SIZE_FILTER2_DEFAULT 1024 //Nb cells in filter2 by default (when s=8)
#define NB_CHILDREN_PER_SKP 248 //An array of CC->children will stores suffixes corresponding to NB_CHILDREN_PER_SKP prefixes

#define NB_KMERS_PER_CC63 255 //Maximum capacity of UCs containing suffixes of length 63 is NB_KMERS_PER_CC63
#define NB_KMERS_PER_CC54 255 //Maximum capacity of UCs containing suffixes of length 54 is NB_KMERS_PER_CC54
#define NB_KMERS_PER_CC45 255 //Maximum capacity of UCs containing suffixes of length 45 is NB_KMERS_PER_CC45
#define NB_KMERS_PER_CC36 255 //Maximum capacity of UCs containing suffixes of length 36 is NB_KMERS_PER_CC36
#define NB_KMERS_PER_CC27 255 //Maximum capacity of UCs containing suffixes of length 27 is NB_KMERS_PER_CC27
#define NB_KMERS_PER_CC18 255 //Maximum capacity of UCs containing suffixes of length 18 is NB_KMERS_PER_CC18
#define NB_KMERS_PER_CC9 255 //Maximum capacity of UCs containing suffixes of length 9 is NB_KMERS_PER_CC9
//Note that the maximum capacity of a UC is also a bound for the number of different prefixes in a CC.

//Define the type of the field children_type in a CC
#if NB_KMERS_PER_CC63 > 255
    #define TYPE63 uint16_t
#else
    #define TYPE63 uint8_t
#endif

#if NB_KMERS_PER_CC54 > 255
    #define TYPE54 uint16_t
#else
    #define TYPE54 uint8_t
#endif

#if NB_KMERS_PER_CC45 > 255
    #define TYPE45 uint16_t
#else
    #define TYPE45 uint8_t
#endif

#if NB_KMERS_PER_CC36 > 255
    #define TYPE36 uint16_t
#else
    #define TYPE36 uint8_t
#endif

#if NB_KMERS_PER_CC27 > 255
    #define TYPE27 uint16_t
#else
    #define TYPE27 uint8_t
#endif

#if NB_KMERS_PER_CC18 > 255
    #define TYPE18 uint16_t
#else
    #define TYPE18 uint8_t
#endif

#if NB_KMERS_PER_CC9 > 255
    #define TYPE9 uint16_t
#else
    #define TYPE9 uint8_t
#endif

//Define how many uint16_t (16bits) are used by a pointer, 2 for 32bits system and 4 for 64bits system
#if (defined(__x86_64__) || defined(_M_X64) || defined(_WIN64) \
  || defined(__powerpc64__) || defined(__ppc64__) || defined(__PPC64__) \
  || defined(__64BIT__) || defined(_LP64) || defined(__LP64__) \
  || defined(__ia64) || defined(__itanium__) || defined(_M_IA64) )
#  define NB_BYTE_PTR 4
#else
#  define NB_BYTE_PTR 2
#endif

/* ===================================================================================================================================
*  STRUCTURES DECLARATION
*  ===================================================================================================================================
*/

// A Compress container (CC) stores prefixes (fixed length) ps of k-mers suffixes contained in a UC.
// A prefix p stored in a CC is itself divided into a prefix p_u and a suffix p_v
typedef struct{
    // The array BF_filter2 contains several arrays. In order:
    // 1 - The Bloom Filter (fixed size) is a bit array used to record inexactly the presence or absence of a set of ps
    // 2 - The Second Filter (fixed size) is a bit array used to record exactly the presence of p_u
    // 3 - The Skip Filter 2 (fixed size) is an array in which each cell records the Hamming weight of 248 bits in the Second Filter
    //      so for example, SkipFilter2[0] = HammingWeight(Filter2[i]) for 0 <= i <= 31, because 31 cells of 8 bits = 248 bits
    //      This array only exists if the CC stores NB_SUBSTRINGS_TRANSFORM prefixes or more.
    // 4 - The Skip Filter 3 is an array in which each cell records the Hamming weight of 248 bits in the Third Filter
    uint8_t* restrict BF_filter2;
    // The array filter3 stores explicitly the p_v, in lexicographic ascending order of p_u
    uint8_t* restrict filter3;
    // The bit array extra_filter3 is a "list of cluster position": One or more consecutive p_v in filter3 belong to the same cluster
    // if they share the same p_u in the Second Filter. A 1 at a position in extra_filter3 signal a new cluster for the prefix at the
    // same position in filter3.
    uint8_t* restrict extra_filter3;
    // The array children_type contains for each prefix stored the number of suffixes linked to it, on 4 bits or 8 bits.
    // It must only be accessed using functions stored in the structure ptrs_on_func because the type of the array is dynamic,
    //depending on the level of the tree. A 0 at a position indicates that the suffixes linked to the prefix at this position
    // in filter3 are stored in a new node. First cell of this array indicates if it is 4 bits mode or 8 bits mode.
    void* restrict children_type;
    // The suffixes linked to the CC prefixes are not really stored in new UCs in new children nodes.
    // Instead they are compactly stored in the CC. children is a list of arrays. Each of these arrays contains the suffixes for
    // NB_CHILDREN_PER_SKP prefixes, so children[0] contains suffixes and annotations for the first NB_CHILDREN_PER_SKP prefixes
    // stored in filter3. A suffix is always immediately followed by its annotation. The length of the annotation is stored in
    // the first cell of each array.
    void* children;
    // children_Node_container is an array of the current CC children nodes
    Node* restrict children_Node_container;

    uint16_t type; //8bits -> size BF in bytes, 5bits -> v(length of p_v), 1bit: full or not, 1bit: 1 if the CC is the last one of its node
    uint16_t nb_elem; //Number of prefix stored in the CC
    uint16_t nb_Node_children; //Number of nodes in children_Node_container
} __attribute__ ((__packed__)) CC;

#define DECLARE_CC( SIZE_KMER )                     \
    typedef struct{                                 \
                                                    \
        uint8_t* restrict BF_filter2;               \
        uint8_t* restrict filter3;                  \
        uint8_t* restrict extra_filter3;            \
                                                    \
        TYPE##SIZE_KMER * restrict children_type;   \
        void* children;                            \
        Node* restrict children_Node_container;     \
                                                    \
        uint16_t type;                              \
        uint16_t nb_elem;                           \
        uint16_t nb_Node_children;                  \
    } __attribute__ ((__packed__)) CC##SIZE_KMER;

// ptrs_on_func is a structure which contains pointers on macro used to manipulate the field children_type of a CC.
// The structure is used in an array, initialized in create_ptrs_on_func()
typedef struct{
    int (*count_nodes)(void*,int, int); //Count the number of nodes in the CC
    int (*count_children)(void*,int, int); //Count the number of suffixes linked to prefixes in the CC
    void (*count_Nodes_Children) (void*, int, int, int*, int*);
    int (*is_child)(void*,int); //Determine if a position corresponds to a suffix (stored in children) or a Node
    void (*add_skp_children)(void*, int, int, int, int, int); //Realloc children arrays and shift the kmers and annotations
    void (*realloc_and_int_children_type)(void*, int, int); //Realloc children_type, shift the counters and initialize one counter
    void (*addNewElt)(void*, int, int); //Increment the counter of one cell in children_type of 1
    int (*getNbElts)(void*, int); //Output a counter of children_type for a specific position
    void (*resetChildrenType)(void*, int); //Set a position in children_type to 0 (transformation to a node)
    void (*allocate_children_type)(void*, int); //Allocate memory for children_type
    int nb_kmers_per_cc; //Number of prefixes a CC, at a specific level of the tree, can contain.
    int mask_shift_kmer; //Suffixes are encoded in arrays of 8bits cells: mask_shift_kmer covers only the bits used on the last cell
    int size_kmer_in_bytes; //Size of the suffixes represented at a given level of the tree, in bytes
    int size_kmer_in_bytes_minus_1; //Size of the suffixes represented at a given level of the tree, minus the size of the prefixes ps, in bytes
    size_t size_type_kmer; //Size in bytes of the type used for children_type, at a specific level of the tree
    int level_min;
    int root;
} ptrs_on_func;

/* ===================================================================================================================================
*  INLINE FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

inline CC* createCC(int nb_bits_bf);
inline void initiateCC(CC* cc, int nb_bits_bf);

inline void freeNode(Node* node);
inline void freeRoot(Root* root);
inline void freeCC(CC* cc);

/* ---------------------------------------------------------------------------------------------------------------
*  createCC(nb_bits_bf)
*  ---------------------------------------------------------------------------------------------------------------
*  Create a CC
*  ---------------------------------------------------------------------------------------------------------------
*  nb_bits_bf: size of the Bloom Filter, in bits
*  ---------------------------------------------------------------------------------------------------------------
*/
inline CC* createCC(int nb_bits_bf){
    CC* cc = malloc(sizeof(CC));
    ASSERT_NULL_PTR(cc,"createCC()")

    cc->BF_filter2 = calloc((nb_bits_bf/SIZE_CELL)+(SIZE_FILTER2_DEFAULT/SIZE_CELL), sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->BF_filter2,"createCC()")

    cc->filter3 = NULL;
    cc->extra_filter3 = NULL;

    cc->children_type = NULL;
    cc->children = NULL;
    cc->children_Node_container = NULL;

    cc->type = (((uint16_t)nb_bits_bf/SIZE_CELL) << SIZE_CELL) | 0x21; //xxxxxxxx|001000|01 -> s=8
    cc->nb_elem = 0;
    cc->nb_Node_children = 0;

    return cc;
}

/* ---------------------------------------------------------------------------------------------------------------
*  initiateCC(cc, nb_bits_bf)
*  ---------------------------------------------------------------------------------------------------------------
*  Initialize a CC
*  ---------------------------------------------------------------------------------------------------------------
*  cc: pointer to a CC
*  nb_bits_bf: size of the Bloom Filter, in bits
*  ---------------------------------------------------------------------------------------------------------------
*/
inline void initiateCC(CC* cc, int nb_bits_bf){
    ASSERT_NULL_PTR(cc,"createCC()")

    cc->BF_filter2 = calloc((nb_bits_bf/SIZE_CELL)+(SIZE_FILTER2_DEFAULT/SIZE_CELL), sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->BF_filter2,"createCC()")

    cc->filter3 = NULL;
    cc->extra_filter3 = NULL;

    cc->children_type = NULL;
    cc->children = NULL;
    cc->children_Node_container = NULL;

    cc->type = (((uint16_t)nb_bits_bf/SIZE_CELL) << 8) | 0x21; //xxxxxxxx|001000|01 -> s=8
    cc->nb_elem = 0;
    cc->nb_Node_children = 0;

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  freeCC(cc)
*  ---------------------------------------------------------------------------------------------------------------
*  free the memory used by a CC
*  ---------------------------------------------------------------------------------------------------------------
*  cc: pointer to a CC
*  ---------------------------------------------------------------------------------------------------------------
*/
inline void freeCC(CC* cc){

    ASSERT_NULL_PTR(cc,"freeCC()")

    free(cc->BF_filter2);
    free(cc->filter3);

    if (cc->extra_filter3 != NULL) free(cc->extra_filter3);
    if (cc->children_type != NULL) free(cc->children_type);

    int i = 0;
    for (i=0; i<CEIL(cc->nb_elem, NB_CHILDREN_PER_SKP); i++)
        if (((UC*)cc->children)[i].suffixes != NULL) free(((UC*)cc->children)[i].suffixes);

    free(cc->children);

    for (i=0; i<cc->nb_Node_children; i++) freeNode(&(cc->children_Node_container[i]));
    if (cc->children_Node_container != NULL) free(cc->children_Node_container);

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  freeNode(node)
*  ---------------------------------------------------------------------------------------------------------------
*  free the memory used by a node
*  ---------------------------------------------------------------------------------------------------------------
*  node: pointer to a Node
*  ---------------------------------------------------------------------------------------------------------------
*/
inline void freeNode(Node* restrict node){

    ASSERT_NULL_PTR(node,"freeNode()")

    if (node->CC_array != NULL){

        int i = -1;
        int it;

        do {
            i++;
            it = ((CC*)node->CC_array)[i].type & 0x1;
            freeCC(&(((CC*)node->CC_array)[i]));
        }
        while (it == 0);

        free(node->CC_array);
    }

    if (node->UC_array.suffixes != NULL) free(node->UC_array.suffixes);

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  freeRoot()
*  ---------------------------------------------------------------------------------------------------------------
*  Free the root of a BFT
*  ---------------------------------------------------------------------------------------------------------------
*  ---------------------------------------------------------------------------------------------------------------
*/
inline void freeRoot(Root* root){

    ASSERT_NULL_PTR(root,"freeRoot()")

    if (root->filenames != NULL){
        int i = 0;
        for (i=0; i<root->nb_genomes; i++) free(root->filenames[i]);
        free(root->filenames);
    }

    free_annotation_array_elem(root->comp_set_colors, root->length_comp_set_colors);

    freeNode(&(root->node));
    free(root);

    return;
}

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

void transform2CC(UC* restrict uc, CC* restrict cc, int size_suffix, ptrs_on_func* restrict func_on_types);
void transform2CC_from_arraySuffix(uint8_t* restrict array_suffix, CC* restrict cc, int size_suffix, int size_annot, uint8_t** annot_extend, uint8_t** annot_cplx,
                                int size_annot_cplx, ptrs_on_func* restrict func_on_types);

void insertSP_CC(resultPresence* restrict pres, uint8_t* restrict kmer, int size_sp, int id_genome, ptrs_on_func* restrict func_on_types,
                   annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

void transform_Filter2n3(CC* cc, int pref_size, int suf_size, ptrs_on_func* restrict func_on_types);
int add_skp_annotation(CC*, int position_type, int size_annot);

ptrs_on_func* create_ptrs_on_func(int size_min, int size_max);
uint16_t** build_skip_nodes(Node* node, int size_kmer, ptrs_on_func* func_on_types);
void free_skip_nodes(Node* node, uint16_t** skp_nodes);

/* ===================================================================================================================================
*  MACROS DECLARATION AND DEFINITIONS
*  ===================================================================================================================================
*/


/* ---------------------------------------------------------------------------------------------------------------
*  allocate_children_type(cc, nb_elt)
*  ---------------------------------------------------------------------------------------------------------------
*  Allocate and initialize CC->children_type for nb_elt counters
*  ---------------------------------------------------------------------------------------------------------------
*  cc: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  nb_elt: number of counters to allocate
*  ---------------------------------------------------------------------------------------------------------------
*/
#define ALLOCATE_CHILDREN_TYPE( SIZE_KMER )                                                      \
void allocate_children_type##SIZE_KMER (void* cc, int nb_elt){                                   \
    ((CC##SIZE_KMER *)cc)->children_type = calloc(CEIL(nb_elt,2)+1, sizeof(TYPE##SIZE_KMER ));   \
    ASSERT_NULL_PTR(((CC##SIZE_KMER *)cc)->children_type, "allocate_children_type##SIZE_KMER ()")\
    ((CC##SIZE_KMER *)cc)->children_type[0] = 4;                                                 \
    return;                                                                                      \
}

/* ---------------------------------------------------------------------------------------------------------------
*  isChild (cc, position)
*  ---------------------------------------------------------------------------------------------------------------
*  Determine if a position in CC->children_type corresponds to suffixes in CC->children (a child) or to a Node
*  ---------------------------------------------------------------------------------------------------------------
*  cc: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  position: counter position
*  ---------------------------------------------------------------------------------------------------------------
*/
#define IS_CHILD( SIZE_KMER )                                                                   \
int isChild##SIZE_KMER (void* cc, int position){                                                \
    ASSERT_NULL_PTR(cc, "isChild##SIZE_KMER ()")                                                \
    if (((CC##SIZE_KMER *)cc)->children_type[0] == 4){                                          \
        if (IS_ODD(position)){                                                                  \
           if ((((CC##SIZE_KMER *)cc)->children_type[(position/2)+1] >> 4) != 0) return 1;      \
        }                                                                                       \
        else if ((((CC##SIZE_KMER *)cc)->children_type[(position/2)+1] & 0xf) != 0) return 1;   \
    }                                                                                           \
    else if (((CC##SIZE_KMER *)cc)->children_type[position+1] != 0) return 1;                   \
    return 0;                                                                                   \
}

/* ---------------------------------------------------------------------------------------------------------------
*  getNbElts (cc, position)
*  ---------------------------------------------------------------------------------------------------------------
*  Get a counter
*  ---------------------------------------------------------------------------------------------------------------
*  cc: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  position: counter position
*  ---------------------------------------------------------------------------------------------------------------
*/
#define GET_NB_ELTS( SIZE_KMER )                                                                    \
int getNbElts##SIZE_KMER (void* cc, int position){                                                  \
    ASSERT_NULL_PTR(cc, "getNbElts##SIZE_KMER ()")                                                  \
    if (((CC##SIZE_KMER *)cc)->children_type[0] == 4){                                              \
        if (IS_ODD(position)) return ((CC##SIZE_KMER *)cc)->children_type[(position/2)+1] >> 4;     \
        else return ((CC##SIZE_KMER *)cc)->children_type[(position/2)+1] & 0xf;                     \
    }                                                                                               \
    else return ((CC##SIZE_KMER *)cc)->children_type[position+1];                                   \
}

/* ---------------------------------------------------------------------------------------------------------------
*  addNewElt(obj, position, current_nb_elem)
*  ---------------------------------------------------------------------------------------------------------------
*  Increment a counter of 1
*  ---------------------------------------------------------------------------------------------------------------
*  obj: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  position: counter position
*  current_nb_elem: number of counters in CC->children_type
*  ---------------------------------------------------------------------------------------------------------------
*/
#define ADD_NEW_ELT( SIZE_KMER )                                                                            \
void addNewElt##SIZE_KMER (void* obj, int position, int current_nb_elem){                                   \
    ASSERT_NULL_PTR(obj, "addNewElt##SIZE_KMER ()")                                                         \
    CC##SIZE_KMER * cc = (CC##SIZE_KMER *)obj;                                                              \
    if (cc->children_type[0] == 4){                                                                         \
        if (IS_ODD(position)){                                                                              \
            if ((((uint16_t)(cc->children_type[(position/2)+1] >> 4)) + 1) == MASK_POWER_8[4]){             \
                TYPE##SIZE_KMER *tmp_children_type = calloc((current_nb_elem+1),sizeof(TYPE##SIZE_KMER ));  \
                ASSERT_NULL_PTR(tmp_children_type, "addNewElt##SIZE_KMER ()")                              \
                int i = 0;                                                                                  \
                for (i=0; i<current_nb_elem; i++){                                                          \
                    if (IS_ODD(i)) tmp_children_type[i+1] = cc->children_type[(i/2)+1] >> 4;                \
                    else tmp_children_type[i+1] = cc->children_type[(i/2)+1] & 0xf;                         \
                }                                                                                           \
                free(cc->children_type);                                                                    \
                cc->children_type = tmp_children_type;                                                      \
                cc->children_type[position+1]++;                                                            \
                cc->children_type[0] = 8;                                                                   \
            }                                                                                               \
            else cc->children_type[(position/2)+1] += 16;                                                   \
        }                                                                                                   \
        else{                                                                                               \
            if ((((uint16_t)(cc->children_type[(position/2)+1] & 0xf)) + 1) == MASK_POWER_8[4]){            \
                TYPE##SIZE_KMER *tmp_children_type = calloc((current_nb_elem+1),sizeof(TYPE##SIZE_KMER ));  \
                ASSERT_NULL_PTR(tmp_children_type, "addNewElt##SIZE_KMER ()")                              \
                int i = 0;                                                                                  \
                for (i=0; i<current_nb_elem; i++){                                                          \
                    if (IS_ODD(i)) tmp_children_type[i+1] = cc->children_type[(i/2)+1] >> 4;                \
                    else tmp_children_type[i+1] = cc->children_type[(i/2)+1] & 0xf;                         \
                }                                                                                           \
                free(cc->children_type);                                                                    \
                cc->children_type = tmp_children_type;                                                      \
                cc->children_type[position+1]++;                                                            \
                cc->children_type[0] = 8;                                                                   \
            }                                                                                               \
            else cc->children_type[(position/2)+1]++;                                                       \
        }                                                                                                   \
    }                                                                                                       \
    else cc->children_type[position+1]++;                                                                   \
}

/* ---------------------------------------------------------------------------------------------------------------
*  realloc_and_int_children_type(obj, current_nb_elem, position_insert)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert a new counter initialized to 1 at position_insert in CC->children_type
*  ---------------------------------------------------------------------------------------------------------------
*  obj: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  current_nb_elem: number of counters in CC->children_type
*  position_insert: position to insert the new counter
*  ---------------------------------------------------------------------------------------------------------------
*/
#define REALLOC_INIT_CHILDREN_TYPE( SIZE_KMER )                                                                                                                     \
void realloc_and_int_children_type##SIZE_KMER (void* obj, int current_nb_elem, int position_insert){                                                                \
    ASSERT_NULL_PTR(obj, "realloc_and_int_children_type##SIZE_KMER ()")                                                                                             \
    CC##SIZE_KMER * cc = (CC##SIZE_KMER *)obj;                                                                                                                      \
                                                                                                                                                                    \
    if (cc->children_type[0] == 8){                                                                                                                                 \
        cc->children_type = realloc(cc->children_type, (current_nb_elem+2)*sizeof( TYPE##SIZE_KMER ));                                                              \
        ASSERT_NULL_PTR(cc->children_type, "realloc_and_int_children_type##SIZE_KMER ()")                                                                           \
        memmove(&(cc->children_type[position_insert+2]), &(cc->children_type[position_insert+1]), (current_nb_elem-position_insert)*sizeof( TYPE##SIZE_KMER ));     \
        cc->children_type[position_insert+1] = 1;                                                                                                                   \
    }                                                                                                                                                               \
    else{                                                                                                                                                           \
        if (CEIL(current_nb_elem+1,2) > CEIL(current_nb_elem,2)){                                                                                                   \
            cc->children_type = realloc(cc->children_type, (CEIL(current_nb_elem+1,2)+1)*sizeof( TYPE##SIZE_KMER ));                                                \
            ASSERT_NULL_PTR(cc->children_type, "realloc_and_int_children_type##SIZE_KMER ()")                                                                       \
            cc->children_type[CEIL(current_nb_elem+1,2)] = 0;                                                                                                       \
        }                                                                                                                                                           \
                                                                                                                                                                    \
        int i = 0;                                                                                                                                                  \
        for (i=current_nb_elem+2; i>position_insert+2; i--){                                                                                                        \
            if (IS_ODD(i)) cc->children_type[i/2] <<= 4;                                                                                                            \
            else cc->children_type[i/2] |= cc->children_type[(i-1)/2] >> 4;                                                                                         \
        }                                                                                                                                                           \
                                                                                                                                                                    \
        if (IS_ODD(i)) cc->children_type[i/2] = (cc->children_type[i/2] & 0xf) + 16;                                                                                \
        else cc->children_type[i/2] = (cc->children_type[i/2] & 0xf0) + 1;                                                                                          \
    }                                                                                                                                                               \
}

/* ---------------------------------------------------------------------------------------------------------------
*  resetChildrenType (cc, position)
*  ---------------------------------------------------------------------------------------------------------------
*  Set a position to 0 in CC->children_type
*  ---------------------------------------------------------------------------------------------------------------
*  cc: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  position: position to set to 0
*  ---------------------------------------------------------------------------------------------------------------
*/
#define RESET_CHILDREN_TYPE( SIZE_KMER )                                                                    \
void resetChildrenType##SIZE_KMER (void* cc, int position){                                                 \
    ASSERT_NULL_PTR(cc, "resetChildrenType##SIZE_KMER ()")                                                  \
    if (((CC##SIZE_KMER *)cc)->children_type[0] == 8) ((CC##SIZE_KMER *)cc)->children_type[position+1] = 0; \
    else if (IS_ODD(position)) ((CC##SIZE_KMER *)cc)->children_type[(position/2)+1] &= 0xf;                 \
    else ((CC##SIZE_KMER *)cc)->children_type[(position/2)+1] &= 0xf0;                                      \
}

/* ---------------------------------------------------------------------------------------------------------------
*  count_Nodes (obj, start, end)
*  ---------------------------------------------------------------------------------------------------------------
*  Count the number of Nodes between two positions in CC->children_Node_container
*  ---------------------------------------------------------------------------------------------------------------
*  obj: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  start: position of start (including itself)
*  end: position of end (not including itself)
*  ---------------------------------------------------------------------------------------------------------------
*/
#define DECLARE_COUNT_NODES( SIZE_KMER )                                            \
int count_Nodes##SIZE_KMER (void* obj, int start, int end){                         \
    ASSERT_NULL_PTR(obj, "count_Nodes##SIZE_KMER ()")                               \
    int count = 0;                                                                  \
    int z=0;                                                                        \
    CC##SIZE_KMER * cc = (CC##SIZE_KMER *)obj;                                      \
    if (cc->children_type[0] == 8){                                                 \
        for (z=start+1; z<end+1; z++){                                              \
            if (cc->children_type[z] == 0) count++;                                 \
        }                                                                           \
    }                                                                               \
    else {                                                                          \
        for (z=start+2; z<end+2; z++){                                              \
            if (IS_ODD(z)){                                                         \
                if (cc->children_type[z/2] < 0x10) count++;                         \
            }                                                                       \
            else if ((cc->children_type[z/2] & 0xf) == 0) count++;                  \
        }                                                                           \
    }                                                                               \
                                                                                    \
    return count;                                                                   \
}

/* ---------------------------------------------------------------------------------------------------------------
*  count_Nodes_Children (obj, start, end, count_children, count_nodes)
*  ---------------------------------------------------------------------------------------------------------------
*  Count the number of nodes and children between two positions in CC->children_Node_container
*  ---------------------------------------------------------------------------------------------------------------
*  obj: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  start: position of start (including itself)
*  end: position of end (not including itself)
*  count_children: pointer to an integer that will contain the count of children
*  count_nodes: pointer to an integer that will contain the count of nodes
*  ---------------------------------------------------------------------------------------------------------------
*/
#define DECLARE_COUNT_NODES_CHILDREN( SIZE_KMER )                                                            \
void count_Nodes_Children##SIZE_KMER (void* obj, int start, int end, int* count_children, int* count_nodes){ \
    ASSERT_NULL_PTR(obj, "count_Nodes_Children##SIZE_KMER ()")                                               \
    ASSERT_NULL_PTR(count_children, "count_Nodes_Children##SIZE_KMER ()")                                    \
    ASSERT_NULL_PTR(count_nodes, "count_Nodes_Children##SIZE_KMER ()")                                       \
                                                                                                             \
    int z=0;                                                                                                 \
    CC##SIZE_KMER * cc = (CC##SIZE_KMER *)obj;                                                               \
                                                                                                             \
    if (cc->children_type[0] == 8){                                                                          \
        for (z=start+1; z<end+1; z++){                                                                       \
            if (cc->children_type[z] == 0) *count_nodes += 1;                                                \
            else *count_children += cc->children_type[z];                                                    \
        }                                                                                                    \
    }                                                                                                        \
    else {                                                                                                   \
        for (z=start+2; z<end+2; z++){                                                                       \
            if (IS_ODD(z)){                                                                                  \
                if ((cc->children_type[z/2] < 0x10) == 0) *count_nodes += 1;                                 \
                else *count_children += cc->children_type[z/2] >> 4;                                         \
            }                                                                                                \
            else {                                                                                           \
                if ((cc->children_type[z/2] & 0xf) == 0) *count_nodes += 1;                                  \
                else *count_children += cc->children_type[z/2] & 0xf;                                        \
            }                                                                                                \
        }                                                                                                    \
    }                                                                                                        \
                                                                                                             \
    return;                                                                                                  \
}
/* ---------------------------------------------------------------------------------------------------------------
*  count_Children (obj, start, end)
*  ---------------------------------------------------------------------------------------------------------------
*  Count the number of suffixes between two positions in CC->children
*  ---------------------------------------------------------------------------------------------------------------
*  obj: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  start: position of start (including itself)
*  end: position of end (not including itself)
*  ---------------------------------------------------------------------------------------------------------------
*/
#define DECLARE_COUNT_CHILDREN( SIZE_KMER )                                         \
int count_Children##SIZE_KMER (void* obj, int start, int end){                      \
    ASSERT_NULL_PTR(obj, "count_Children##SIZE_KMER ()")                            \
                                                                                    \
    int count = 0;                                                                  \
    int z=0;                                                                        \
                                                                                    \
    CC##SIZE_KMER * cc = (CC##SIZE_KMER *)obj;                                      \
    if (cc->children_type[0] == 8){                                                 \
        for (z=start+1; z<end+1; z++){                                              \
            count += cc->children_type[z];                                          \
        }                                                                           \
    }                                                                               \
    else {                                                                          \
        for (z=start+2; z<end+2; z++){                                              \
            if (IS_ODD(z)) count += cc->children_type[z/2] >> 4;                    \
            else count += cc->children_type[z/2] & 0xf;                             \
        }                                                                           \
    }                                                                               \
                                                                                    \
    return count;                                                                   \
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

#define DECLARE_ADD_SKP_CHILDREN( SIZE_KMER )                                                                                                                                               \
void add_skp_children##SIZE_KMER (void* obj, int position_type, int position_child, int count_before_child, int size_substrings, int size_annot){                                       \
    ASSERT_NULL_PTR(obj, "add_skp_children##SIZE_KMER ()") \
\
    CC##SIZE_KMER * cc = (CC##SIZE_KMER *)obj; \
\
    UC* uc; \
    UC* uc_tmp; \
\
    uint8_t** extend_annot_z_minus1 = NULL;\
    uint8_t** cplx_annot_z_minus1 = NULL;\
\
    int i=0, z=0;\
    int count2push = 0;\
    int count2delete = 0;\
    int max_size_z = 0;\
    int max_size_z_minus1 = 0;\
    int max_size_cplx_z_minus1 = 0;\
    int size_line_children = 0;\
    int size_line_children_z_minus1 = 0;\
    int cpt_annot_extend_z = 0;\
    int cpt_annot_cplx_z = 0;\
    int start_pos_z_minus1 = 0;\
\
    int nb_cell_skp = CEIL(cc->nb_elem+1, NB_CHILDREN_PER_SKP);\
    int start = position_type/NB_CHILDREN_PER_SKP;\
\
    if (cc->nb_elem%NB_CHILDREN_PER_SKP == 0){\
        cc->children = realloc(cc->children, nb_cell_skp*sizeof(UC));\
        ASSERT_NULL_PTR(cc->children, "add_skp_children##SIZE_KMER ()") \
\
        uc = &(((UC*)cc->children)[nb_cell_skp-1]); \
        initializeUC(uc); \
        uc->size_annot = 1;\
    }\
\
    for (z=nb_cell_skp-1; z>=start; z--){\
\
        uc = &(((UC*)cc->children)[z]); \
        uc->nb_children -= count2delete;\
        \
        if (z == start){\
\
            if (size_annot > uc->size_annot){ \
                if (uc->nb_children == 0){ \
                    uc->suffixes = calloc(size_substrings + size_annot, sizeof(uint8_t)); \
                    uc->size_annot = size_annot; \
                } \
                else realloc_annotation(uc, size_substrings, uc->nb_children, size_annot, 1, position_child);\
            } \
            else if (size_annot < uc->size_annot){\
\
                if (uc->nb_extended_annot != 0) max_size_z = uc->size_annot + 1;\
                else max_size_z = max_size_per_sub(uc->suffixes, uc->nb_children, size_substrings, uc->size_annot);\
\
                max_size_z = MAX(size_annot, max_size_z);\
\
                if ((max_size_z > 0) && (uc->nb_extended_annot == 0) && (max_size_z < uc->size_annot)){\
                    if (uc->nb_children == 0){ \
                        uc->suffixes = calloc(size_substrings + max_size_z, sizeof(uint8_t)); \
                        uc->size_annot = max_size_z; \
                    } \
                    else realloc_annotation(uc, size_substrings, uc->nb_children, max_size_z, 1, position_child);\
                }\
                else {\
                    size_line_children = size_substrings + uc->size_annot;\
                    uc->suffixes = realloc(uc->suffixes, ((uc->nb_children+1) * size_line_children + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT \
                                                          + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));\
                    ASSERT_NULL_PTR(uc->suffixes, "add_skp_children##SIZE_KMER ()") \
\
                    count_before_child -= count2delete;\
\
                    memmove(&(uc->suffixes[(position_child+1)*size_line_children]),\
                            &(uc->suffixes[position_child*size_line_children]),\
                            ((count_before_child * size_line_children) + (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT) \
                             + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));\
\
                    shift_extended_annot(uc, size_substrings, uc->nb_children+1, position_child);\
                }\
            }\
            else {\
                size_line_children = size_substrings + uc->size_annot;\
                uc->suffixes = realloc(uc->suffixes, ((uc->nb_children + 1) * size_line_children + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT \
                                                      + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));\
                ASSERT_NULL_PTR(uc->suffixes, "add_skp_children##SIZE_KMER ()") \
\
                count_before_child -= count2delete;\
                memmove(&(uc->suffixes[(position_child+1)*size_line_children]),\
                        &(uc->suffixes[position_child*size_line_children]),\
                        ((count_before_child * size_line_children) + (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT) \
                         + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));\
\
                shift_extended_annot(uc, size_substrings, uc->nb_children+1, position_child);\
            }\
\
            uc->nb_children++;\
            shift_annot_cplx_nodes(uc, size_substrings, uc->nb_children, position_child); \
        }\
        else {\
            uc_tmp = &(((UC*)cc->children)[z-1]); \
\
            count2push = getNbElts##SIZE_KMER (obj, z*NB_CHILDREN_PER_SKP-1);\
\
            if (count2push != 0){ \
\
                start_pos_z_minus1 = uc_tmp->nb_children - count2push;\
\
                size_line_children_z_minus1 = size_substrings + uc_tmp->size_annot;\
                size_line_children = size_substrings + uc->size_annot;\
\
                cplx_annot_z_minus1 = get_annots_cplx_nodes(uc_tmp, size_substrings, uc_tmp->nb_children, start_pos_z_minus1, uc_tmp->nb_children-1); \
                max_size_cplx_z_minus1 = max_size_annot_cplx_sub(&(uc_tmp->suffixes[uc_tmp->nb_children * size_line_children_z_minus1 + \
                                                                    uc_tmp->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),\
                                                                    uc_tmp->nb_cplx_nodes, uc_tmp->size_annot_cplx_nodes, start_pos_z_minus1, uc_tmp->nb_children-1); \
\
                if (max_size_cplx_z_minus1 > uc->size_annot_cplx_nodes){ \
                    increase_size_annot_cplx_nodes(uc, size_substrings, uc->nb_children, max_size_cplx_z_minus1, 1); \
                } \
                else if (max_size_cplx_z_minus1 < uc->size_annot_cplx_nodes){ \
                    max_size_z = max_size_annot_cplx_sub(&(uc->suffixes[uc->nb_children * size_line_children + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),\
                                                            uc->nb_cplx_nodes, uc->size_annot_cplx_nodes, 0, uc->nb_children-1); \
\
                    if (max_size_z <= max_size_cplx_z_minus1) decrease_size_annot_cplx_nodes(uc, size_substrings, uc->nb_children, max_size_z); \
                } \
\
                cpt_annot_cplx_z = 0; \
\
                if (cplx_annot_z_minus1 != NULL){ \
                    for (i = 0; i < count2push; i++)\
                        if (cplx_annot_z_minus1[i] != NULL) cpt_annot_cplx_z++; \
                } \
\
                if ((extend_annot_z_minus1 = get_extend_annots(uc_tmp, size_substrings, uc_tmp->nb_children, start_pos_z_minus1, uc_tmp->nb_children-1)) != NULL){\
                    max_size_z_minus1 = uc_tmp->size_annot+1;\
                }\
                else max_size_z_minus1 = max_size_per_sub(&(uc_tmp->suffixes[start_pos_z_minus1 * size_line_children_z_minus1]),\
                                                            count2push,\
                                                            size_substrings,\
                                                            uc_tmp->size_annot);\
    \
                COMPUTE_ANNOT_EXT:if (max_size_z_minus1 > uc->size_annot+1) realloc_annotation(uc, size_substrings, uc->nb_children, max_size_z_minus1, 0, 0);\
                else if (max_size_z_minus1 < uc->size_annot){\
    \
                    if (uc->nb_extended_annot != 0) max_size_z = uc->size_annot + 1;\
                    else max_size_z = max_size_per_sub(uc->suffixes, uc->nb_children, size_substrings, uc->size_annot);\
    \
                    max_size_z = MAX(max_size_z_minus1, max_size_z);\
    \
                    if ((max_size_z > 0) && (uc->nb_extended_annot == 0) && (max_size_z < uc->size_annot)){\
                        realloc_annotation(uc, size_substrings, uc->nb_children, max_size_z, 0, 0);\
                    }\
                }\
    \
                size_line_children_z_minus1 = size_substrings + uc_tmp->size_annot;\
                size_line_children = size_substrings + uc->size_annot;\
                cpt_annot_extend_z = 0;\
    \
                if (max_size_z_minus1 > uc->size_annot){\
                    if (extend_annot_z_minus1 != NULL){\
                        for (i=0; i<count2push; i++)\
                            if (extend_annot_z_minus1[i] != NULL) cpt_annot_extend_z++;\
                    }\
                    else {\
                        uint8_t* tab = &(uc_tmp->suffixes[start_pos_z_minus1*size_line_children_z_minus1]);\
                        for (i=0; i<count2push; i++){\
                            if (tab[i*size_line_children_z_minus1 + size_substrings + max_size_z_minus1 - 1] != 0) cpt_annot_extend_z++;\
                        }\
                    }\
    \
                    if ((uc->nb_extended_annot + cpt_annot_extend_z) * SIZE_BYTE_EXT_ANNOT > uc->nb_children + count2push){\
                        recopy_back_annot_extend(uc, size_substrings, uc->nb_children);\
                        goto COMPUTE_ANNOT_EXT;\
                    }\
                }\
    \
                uc->nb_children += count2push;\
    \
                uc->suffixes = realloc(uc->suffixes, \
                                    ((uc->nb_children * size_line_children) + (uc->nb_extended_annot + cpt_annot_extend_z) * SIZE_BYTE_EXT_ANNOT \
                                     + (uc->nb_cplx_nodes + cpt_annot_cplx_z) * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));\
                ASSERT_NULL_PTR(uc->suffixes, "add_skp_children##SIZE_KMER ()") \
    \
                memmove(&(uc->suffixes[count2push*size_line_children]),\
                        uc->suffixes,\
                        ((uc->nb_children - count2push) * size_line_children + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT \
                         + (uc->nb_cplx_nodes + cpt_annot_cplx_z) * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes)) * sizeof(uint8_t));\
                memset(uc->suffixes, 0, count2push * size_line_children * sizeof(uint8_t));\
    \
                if (max_size_z_minus1 > uc->size_annot){\
                    if (extend_annot_z_minus1 != NULL){\
\
                        for (i=0; i<count2push; i++){\
                             memcpy(&(uc->suffixes[i*size_line_children]),\
                                    &(uc_tmp->suffixes[(start_pos_z_minus1+i) * size_line_children_z_minus1]),\
                                    size_line_children_z_minus1*sizeof(uint8_t));\
                        }\
\
                        for (i=0; i<count2push; i++){\
                            if (extend_annot_z_minus1[i] == NULL) shift_extended_annot(uc, size_substrings, uc->nb_children, i); \
                            else insert_extend_annot(uc, size_substrings, uc->nb_children, i, extend_annot_z_minus1[i][0], 1); \
                        }\
\
                        for (i=0; i<count2push; i++){\
                            if ((cplx_annot_z_minus1 == NULL) || (cplx_annot_z_minus1[i] == NULL)) shift_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i);\
                            else insert_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i, cplx_annot_z_minus1[i], max_size_cplx_z_minus1, 1);\
                        }\
                    }\
                    else{\
\
                        uint8_t* tab = &(uc_tmp->suffixes[start_pos_z_minus1 * size_line_children_z_minus1]);\
\
                        for (i=0; i<count2push; i++){\
                             memcpy(&(uc->suffixes[i*size_line_children]),\
                                    &(uc_tmp->suffixes[(start_pos_z_minus1+i)*size_line_children_z_minus1]),\
                                    (size_substrings + max_size_z_minus1 - 1) * sizeof(uint8_t));\
                        }\
\
                        for (i=0; i<count2push; i++){\
                            if (tab[i*size_line_children_z_minus1 + size_substrings + max_size_z_minus1 - 1] != 0){\
                                insert_extend_annot(uc, size_substrings, uc->nb_children, i, \
                                                    tab[i*size_line_children_z_minus1 + size_substrings + max_size_z_minus1 - 1], 1);\
                            }\
                            else shift_extended_annot(uc, size_substrings, uc->nb_children, i);\
                        }\
\
                        for (i=0; i<count2push; i++){\
                            if ((cplx_annot_z_minus1 == NULL) || (cplx_annot_z_minus1[i] == NULL)) shift_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i);\
                            else insert_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i, cplx_annot_z_minus1[i], max_size_cplx_z_minus1, 1);\
                        }\
\
                    }\
                }\
                else{\
    \
                    for (i=0; i<count2push; i++){\
                        if (extend_annot_z_minus1 != NULL){\
                             memcpy(&(uc->suffixes[i*size_line_children]),\
                                    &(uc_tmp->suffixes[(start_pos_z_minus1+i)*size_line_children_z_minus1]),\
                                    size_line_children_z_minus1*sizeof(uint8_t));\
    \
                            if (extend_annot_z_minus1[i] != NULL) \
                                uc->suffixes[i*size_line_children+size_line_children_z_minus1] = extend_annot_z_minus1[i][0];\
                        }\
                        else{\
                             memcpy(&(uc->suffixes[i*size_line_children]),\
                                    &(uc_tmp->suffixes[(start_pos_z_minus1 + i) * size_line_children_z_minus1]),\
                                    (size_substrings + max_size_z_minus1) * sizeof(uint8_t));\
                        }\
                    }\
\
                    for (i=0; i<count2push; i++) shift_extended_annot(uc, size_substrings, uc->nb_children, i);\
\
                    for (i=0; i<count2push; i++){ \
                        if ((cplx_annot_z_minus1 == NULL) || (cplx_annot_z_minus1[i] == NULL)) shift_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i);\
                        else insert_annot_cplx_nodes(uc, size_substrings, uc->nb_children, i, cplx_annot_z_minus1[i], max_size_cplx_z_minus1, 1);\
                    } \
                }\
    \
                delete_extend_annots(uc_tmp, size_substrings, uc_tmp->nb_children , start_pos_z_minus1,  uc_tmp->nb_children-1, 0, 1, 0);\
                delete_annot_cplx_nodes(uc_tmp, size_substrings, uc_tmp->nb_children, start_pos_z_minus1, uc_tmp->nb_children-1, 1, 1, 0);\
                \
            } \
\
            count2delete = count2push;\
\
            if (extend_annot_z_minus1 != NULL){\
                free(extend_annot_z_minus1);\
                extend_annot_z_minus1 = NULL;\
            }\
\
            if (cplx_annot_z_minus1 != NULL){\
                free(cplx_annot_z_minus1);\
                cplx_annot_z_minus1 = NULL;\
            }\
        }\
    }\
\
    return;\
}

/* --------------------------------------------------------------------------------------
*  DEFINE_STRUCT_AND_FUNC( SIZE_KMER )
*  --------------------------------------------------------------------------------------
*  Initialize CCs and functions manipulating the children_type field.
*  The CC initialized will be used to transform UC containing suffix of length SIZE_KMER.
*  --------------------------------------------------------------------------------------
*  SIZE_KMER : size of the UC suffixes (in char) the CC will have to transform
*  --------------------------------------------------------------------------------------
*/

#define DEFINE_STRUCT_AND_FUNC( SIZE_KMER )   \
    DECLARE_CC( SIZE_KMER )                   \
    GET_NB_ELTS( SIZE_KMER )                  \
    DECLARE_ADD_SKP_CHILDREN( SIZE_KMER )     \
    IS_CHILD( SIZE_KMER )                     \
    DECLARE_COUNT_NODES( SIZE_KMER )          \
    DECLARE_COUNT_CHILDREN( SIZE_KMER )       \
    DECLARE_COUNT_NODES_CHILDREN( SIZE_KMER ) \
    REALLOC_INIT_CHILDREN_TYPE( SIZE_KMER )   \
    ADD_NEW_ELT( SIZE_KMER )                  \
    ALLOCATE_CHILDREN_TYPE( SIZE_KMER )       \
    RESET_CHILDREN_TYPE( SIZE_KMER )

#endif
