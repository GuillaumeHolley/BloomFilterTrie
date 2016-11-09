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
#include "presenceNode.h"
#include "CC.h"
#include "retrieveAnnotation.h"

/* ===================================================================================================================================
*  FUNCTIONS DECLARATION
*  ===================================================================================================================================
*/

void insertKmers(BFT_Root*  root, uint8_t*  array_kmers, int nb_kmers, uint32_t id_genome,
                 int size_id_genome);

void insertKmer_Node(Node*  node, BFT_Root*  root, int lvl_node, uint8_t*  suffix, int size_suffix,
                     uint8_t* kmer, uint32_t id_genome, int size_id_genome, int pos_CC_start_search);

Node* insertKmer_Node_special(BFT_Root*  root, int lvl_cont, uint8_t*  suffix, int size_suffix,
                              uint8_t*  kmer, uint32_t id_genome, int size_id_genome);
