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

#include "useful_macros.h"
#include "fasta.h"
#include "CC.h"
#include "retrieveAnnotation.h"
#include "write_to_disk.h"

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

//int extract_kmers_from_node(Node* n, BFT_Root* root, int lvl_node, uint8_t* kmer, int size_kmer, int bucket, int pos_in_bucket, FILE* file_output);

size_t iterate_over_kmers_from_node(Node* n, BFT_Root* root, int lvl_node, uint8_t* kmer, BFT_kmer* bft_kmer, int size_kmer,
                                    int bucket, int pos_in_bucket, size_t (*f)(BFT_kmer*, BFT_Root*, va_list), va_list args);

size_t iterate_over_prefixes_from_node(Node* n, BFT_Root* root, int lvl_node, uint8_t* kmer, BFT_kmer* bft_kmer, int size_kmer,
                                       uint8_t* prefix_comp, int length_prefix, int bucket, int pos_in_bucket,
                                       size_t (*f)(BFT_kmer*, BFT_Root*, va_list), va_list args);
