#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <stdarg.h>
#include <stdbool.h>

#include <omp.h>

#include "bft.h"

void merging_BFT(char* prefix_bft1, char* prefix_bft2, char* output_prefix, int cut_lvl,
                 bool packed_in_subtries);
size_t l_insert_kmer(BFT_kmer* kmer, BFT* graph, va_list args);
void l_insert_kmer_bis(BFT* bft, int lvl_root, uint8_t* kmer_comp, uint8_t* kmer_comp_cpy,
                        int nb_bytes, int old_nb_genome_ids_bft, bool overlap_genome_ids,
                        BFT_annotation* bft_annot, int pos_start_search);
