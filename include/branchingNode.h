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

#include "useful_macros.h"
#include "CC.h"
#include "presenceNode.h"
#include "insertNode.h"

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

int isBranchingLeft(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer);
int isBranchingRight(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer);

resultPresence* getRightNeighbors(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer);
resultPresence* getLeftNeighbors(Node*  node, BFT_Root* root, int lvl_node, uint8_t*  kmer, int size_kmer);

void extractSuffix(uint8_t* kmer_tmp, int size_kmer, int size_kmer_array, int shifting_prefix,
                   int it_bucket, uint8_t* suffix_start, int size_suffix_start_bytes, int shifting_suffix);
