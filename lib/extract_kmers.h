#ifndef DEF_EXTRACT_KMERS
#define DEF_EXTRACT_KMERS

/* ===================================================================================================================================
*  INCLUDES
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
#include "./../lib/CC.h"
#include "./../lib/retrieveAnnotation.h"

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

int extract_kmers_from_node(Node* n, int lvl_node, uint8_t* kmer, int size_kmer, int bucket, int pos_in_bucket,
                            int size_kmer_root, int compression, info_per_level* restrict info_per_lvl, FILE* file_output);

#endif
