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

int extract_kmers_from_node(Node* n, uint8_t* kmer, int size_kmer, int bucket, int pos_in_bucket, int size_kmer_root,
                             ptrs_on_func* restrict func_on_types, FILE* file_output);

#endif
