#pragma once

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

#include <stdarg.h>

#include "./../lib/useful_macros.h"
#include "./../lib/fasta.h"
#include "./../lib/CC.h"
#include "./../lib/retrieveAnnotation.h"
#include "./../lib/write_to_disk.h"

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

void iterate_over_kmers_from_node(Node* n, BFT_Root* root, int lvl_node, uint8_t* kmer, BFT_kmer* bft_kmer, int size_kmer,
                                  int bucket, int pos_in_bucket, size_t (*f)(BFT_kmer*, BFT_Root*, va_list), va_list args);
