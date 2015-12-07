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

#include "./../lib/default_param.h"
#include "./../lib/useful_macros.h"
#include "./../lib/Node.h"

#include "./../lib/hash.h"
#include "./../lib/quicksort.h"
#include "./../lib/popcnt.h"
#include "./../lib/annotation.h"

const uint64_t MASK_POWER_16[17];

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
    uint8_t* restrict BF_filter2;
    // The array filter3 stores explicitly the p_v, in lexicographic ascending order of p_u
    uint8_t* restrict filter3;
    // The bit array extra_filter3 is a "list of cluster position": One or more consecutive p_v in filter3 belong to the same cluster
    // if they share the same p_u in the Second Filter. A 1 at a position in extra_filter3 signal a new cluster for the prefix at the
    // same position in filter3.
    uint8_t* restrict extra_filter3;
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
    Node* restrict children_Node_container;
} __attribute__ ((__packed__)) CC;

/* ===================================================================================================================================
*  INLINE FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

inline CC* createCC(int nb_bits_bf);
inline void initiateCC(CC* cc, int nb_bits_bf);

inline void freeNode(Node* node, int lvl_node, info_per_level* restrict info_per_lvl);
inline void freeRoot(Root* root);
inline void freeCC(CC* cc, int lvl_cc, info_per_level* restrict info_per_lvl);

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

/* ---------------------------------------------------------------------------------------------------------------
*  freeCC(cc)
*  ---------------------------------------------------------------------------------------------------------------
*  free the memory used by a CC
*  ---------------------------------------------------------------------------------------------------------------
*  cc: pointer to a CC
*  ---------------------------------------------------------------------------------------------------------------
*/
inline void freeCC(CC* cc, int lvl_cc, info_per_level* restrict info_per_level){

    ASSERT_NULL_PTR(cc,"freeCC()")

    free(cc->BF_filter2);
    free(cc->filter3);

    if (cc->extra_filter3 != NULL) free(cc->extra_filter3);
    if (cc->children_type != NULL) free(cc->children_type);

    int i = 0;
    for (i=0; i<CEIL(cc->nb_elem, info_per_level[lvl_cc].nb_ucs_skp); i++)
        if (((UC*)cc->children)[i].suffixes != NULL) free(((UC*)cc->children)[i].suffixes);

    free(cc->children);

    for (i=0; i<cc->nb_Node_children; i++) freeNode(&(cc->children_Node_container[i]), lvl_cc-1, info_per_level);
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
inline void freeNode(Node* restrict node, int lvl_node, info_per_level* restrict info_per_lvl){

    ASSERT_NULL_PTR(node,"freeNode()")

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

    freeNode(&(root->node), (root->k / NB_CHAR_SUF_PREF) - 1, root->info_per_lvl);

    free(root->info_per_lvl);
    free(root->hash_v);
    free(root);

    return;
}

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

void transform2CC(UC* restrict uc, CC* restrict cc, Root* root, int lvl_cc, int size_suffix);
void transform2CC_from_arraySuffix(uint8_t* restrict array_suffix, CC* restrict cc, Root* root, int lvl_cc, int size_suffix,
                                   int size_annot, uint8_t** annot_extend, uint8_t** annot_cplx, int size_annot_cplx);

void insertSP_CC(resultPresence* restrict pres, uint8_t* restrict kmer, int size_sp, uint32_t id_genome, int size_id_genome,
                 info_per_level* restrict info_per_lvl, annotation_inform* ann_inf, annotation_array_elem* annot_sorted);

void transform_Filter2n3(CC* cc, int pref_size, int suf_size, info_per_level* restrict info_per_lvl);

int add_skp_annotation(CC*, int position_type, int size_annot, info_per_level* restrict info_per_lvl);
void add_skp_children(CC* cc, int position_type, int position_child, int count_before_child, int size_substrings,
                      int size_annot, info_per_level* restrict info_per_lvl);

info_per_level* create_info_per_level(int size_max);
uint16_t** build_skip_nodes(Node* node);
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
/*inline void allocate_children_type (CC* cc, int nb_elt){

    ASSERT_NULL_PTR(cc, "allocate_children_type()")

    cc->children_type = calloc(CEIL(nb_elt,2)+1, sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->children_type, "allocate_children_type()")

    cc->children_type[0] = 4;
    return;
}*/

inline void allocate_children_type (CC* cc, int nb_elt){

    ASSERT_NULL_PTR(cc, "allocate_children_type()")

    cc->children_type = calloc(CEIL(nb_elt,2), sizeof(uint8_t));
    ASSERT_NULL_PTR(cc->children_type, "allocate_children_type()")

    return;
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
/*inline int is_child(CC* cc, int position){
    ASSERT_NULL_PTR(cc, "isChild()")

    if (cc->children_type[0] == 4){
        if (IS_ODD(position)){
           if ((cc->children_type[(position/2)+1] >> 4) != 0) return 1;
        }
        else if ((cc->children_type[(position/2)+1] & 0xf) != 0) return 1;
    }
    else if (cc->children_type[position+1] != 0) return 1;

    return 0;
}*/

inline int is_child(CC* cc, int position, uint8_t type){
    ASSERT_NULL_PTR(cc, "isChild()")

    if (type == 0){
        if (IS_ODD(position)) return cc->children_type[position/2] > 0xf;
        else return (cc->children_type[position/2] & 0xf) != 0;
    }
    else return cc->children_type[position] != 0;
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
/*inline int getNbElts(CC* cc, int position){

    ASSERT_NULL_PTR(cc, "getNbElts()")

    if (cc->children_type[0] == 4){
        if (IS_ODD(position)) return cc->children_type[(position/2)+1] >> 4;
        else return cc->children_type[(position/2)+1] & 0xf;
    }
    else return cc->children_type[position+1];
}*/

inline int getNbElts(CC* cc, int position, uint8_t type){

    ASSERT_NULL_PTR(cc, "getNbElts()")

    if (type == 0){
        if (IS_ODD(position)) return cc->children_type[position/2] >> 4;
        else return cc->children_type[position/2] & 0xf;
    }
    else return cc->children_type[position];
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
/*inline void addNewElt(CC* cc, int position, int current_nb_elem){
    ASSERT_NULL_PTR(cc, "addNewElt()")

    if (cc->children_type[0] == 4){

        if (IS_ODD(position)){

            if ((cc->children_type[(position/2)+1] >> 4) + 1 == 0x10){

                uint8_t* tmp_children_type = calloc((current_nb_elem+1),sizeof(uint8_t ));
                ASSERT_NULL_PTR(tmp_children_type, "addNewElt()")

                int i = 1, j = 1;
                for (; j < (current_nb_elem/2)+1; i += 2, j++){
                    tmp_children_type[i] = cc->children_type[j] & 0xf;
                    tmp_children_type[i+1] = cc->children_type[j] >> 4;
                }

                if (IS_ODD(current_nb_elem)) tmp_children_type[i] = cc->children_type[j] & 0xf;

                free(cc->children_type);

                cc->children_type = tmp_children_type;
                cc->children_type[position+1]++;
                cc->children_type[0] = 8;
            }
            else cc->children_type[(position/2)+1] += 16;
        }
        else{
            if ((cc->children_type[(position/2)+1] & 0xf) + 1 == 0x10){

                uint8_t* tmp_children_type = calloc((current_nb_elem+1),sizeof(uint8_t ));
                ASSERT_NULL_PTR(tmp_children_type, "addNewElt()")

                int i = 1, j = 1;
                for (; j < (current_nb_elem/2)+1; i += 2, j++){
                    tmp_children_type[i] = cc->children_type[j] & 0xf;
                    tmp_children_type[i+1] = cc->children_type[j] >> 4;
                }

                if (IS_ODD(current_nb_elem)) tmp_children_type[i] = cc->children_type[j] & 0xf;

                free(cc->children_type);

                cc->children_type = tmp_children_type;
                cc->children_type[position+1]++;
                cc->children_type[0] = 8;
            }
            else cc->children_type[(position/2)+1]++;
        }
    }
    else cc->children_type[position+1]++;
}*/

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
/*inline void realloc_and_int_children_type(CC* cc, int current_nb_elem, int position_insert){
    ASSERT_NULL_PTR(cc, "realloc_and_int_children_type ()")

    if (cc->children_type[0] == 8){
        cc->children_type = realloc(cc->children_type, (current_nb_elem+2)*sizeof(uint8_t));
        ASSERT_NULL_PTR(cc->children_type, "realloc_and_int_children_type ()")

        memmove(&(cc->children_type[position_insert+2]),
                &(cc->children_type[position_insert+1]),
                (current_nb_elem-position_insert)*sizeof(uint8_t));

        cc->children_type[position_insert+1] = 1;
    }
    else{
        if (CEIL(current_nb_elem+1,2) > CEIL(current_nb_elem,2)){

            cc->children_type = realloc(cc->children_type, (CEIL(current_nb_elem+1,2)+1)*sizeof(uint8_t));
            ASSERT_NULL_PTR(cc->children_type, "realloc_and_int_children_type ()")

            cc->children_type[CEIL(current_nb_elem+1,2)] = 0;
        }

        int i = 0;
        for (i=current_nb_elem+2; i>position_insert+2; i--){
            if (IS_ODD(i)) cc->children_type[i/2] <<= 4;
            else cc->children_type[i/2] |= cc->children_type[(i-1)/2] >> 4;
        }

        if (IS_ODD(i)) cc->children_type[i/2] = (cc->children_type[i/2] & 0xf) + 16;
        else cc->children_type[i/2] = (cc->children_type[i/2] & 0xf0) + 1;
    }
}*/

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

/* ---------------------------------------------------------------------------------------------------------------
*  resetChildrenType (cc, position)
*  ---------------------------------------------------------------------------------------------------------------
*  Set a position to 0 in CC->children_type
*  ---------------------------------------------------------------------------------------------------------------
*  cc: a CC for which we do not know the type (CC63, CC27, etc. ?)
*  position: position to set to 0
*  ---------------------------------------------------------------------------------------------------------------
*/
/*inline void resetChildrenType (CC* cc, int position){
    ASSERT_NULL_PTR(cc, "resetChildrenType()")

    if (cc->children_type[0] == 8) cc->children_type[position+1] = 0;
    else if (IS_ODD(position)) cc->children_type[(position/2)+1] &= 0xf;
    else cc->children_type[(position/2)+1] &= 0xf0;
}*/

inline void resetChildrenType (CC* cc, int position, uint8_t type){
    ASSERT_NULL_PTR(cc, "resetChildrenType()")

    if (type == 1) cc->children_type[position] = 0;
    else if (IS_ODD(position)) cc->children_type[position/2] &= 0xf;
    else cc->children_type[position/2] &= 0xf0;

    return;
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
/*inline int count_nodes(CC* cc, int start, int end){

    ASSERT_NULL_PTR(cc, "count_Nodes()")

    int count = 0;
    uint8_t* z;

    if (cc->children_type[0] == 8){
        for (z = cc->children_type + start + 1; z < cc->children_type + end + 1; z++) count += *z == 0;
    }
    else{
        for (z = cc->children_type + (start+2)/2; z <= cc->children_type + (end+1)/2; z++)
            count += (*z < 0x10) + ((*z & 0xf) == 0);

        if (IS_ODD(start)) count -= (cc->children_type[(start+2)/2] & 0xf) == 0;
        if (IS_EVEN(end-1)) count -= cc->children_type[(end+1)/2] < 0x10;
    }

    return count;
}*/

inline int count_nodes(CC* cc, int start, int end, uint8_t type){

    ASSERT_NULL_PTR(cc, "count_Nodes()")

    int count = 0;
    uint8_t* z;

    if (type == 1){
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
/*inline void count_Nodes_Children(CC* cc, int start, int end, int* count_children, int* count_nodes){
    ASSERT_NULL_PTR(cc, "count_Nodes_Children()")
    ASSERT_NULL_PTR(count_children, "count_Nodes_Children()")
    ASSERT_NULL_PTR(count_nodes, "count_Nodes_Children()")

    uint8_t* z;

    *count_nodes = 0;
    *count_children = 0;

    if (cc->children_type[0] == 8){
        for (z = cc->children_type + start + 1; z < cc->children_type + end + 1; z++){
            *count_children += *z;
            *count_nodes += *z == 0;
        }
    }
    else {

        int start_r = (start+2)/2;
        int end_r = (end+1)/2;

        for (z = cc->children_type + start_r; z <= cc->children_type + end_r; z++){
            *count_nodes += (*z < 0x10) + ((*z & 0xf) == 0);
            *count_children += (*z >> 4) + (*z & 0xf);
        }

        if (IS_ODD(start)){
            *count_nodes -= (cc->children_type[start_r] & 0xf) == 0;
            *count_children -= cc->children_type[start_r] & 0xf;
        }
        if (IS_EVEN(end-1)){
            *count_nodes -= cc->children_type[end_r] < 0x10;
            *count_children -= cc->children_type[end_r] >> 4;
        }
    }

    return;
}*/

inline void count_Nodes_Children(CC* cc, int start, int end, int* count_children, int* count_nodes, uint8_t type){
    ASSERT_NULL_PTR(cc, "count_Nodes_Children()")
    ASSERT_NULL_PTR(count_children, "count_Nodes_Children()")
    ASSERT_NULL_PTR(count_nodes, "count_Nodes_Children()")

    uint8_t* z;
    int count_nodes_tmp = 0;
    int count_children_tmp = 0;

    if (type == 1){

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
/*inline int count_children(CC* cc, int start, int end){
    ASSERT_NULL_PTR(cc, "count_Children()")

    int count = 0;
    uint8_t* z;

    if (cc->children_type[0] == 8){
        for (z= cc->children_type + start + 1; z < cc->children_type + end + 1; z++) count += *z;
    }
    else {
        for (z= cc->children_type + (start+2)/2; z <= cc->children_type + (end+1)/2; z++)
            count += (*z >> 4) + (*z & 0xf);

        if (IS_ODD(start)) count -= cc->children_type[(start+2)/2] & 0xf;
        if (!IS_ODD(end-1)) count -= cc->children_type[(end+1)/2] >> 4;
    }

    return count;
}*/

inline int count_children(CC* cc, int start, int end, uint8_t type){
    ASSERT_NULL_PTR(cc, "count_Children()")

    int count = 0;
    uint8_t* z;

    if (type == 1){
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

#endif
