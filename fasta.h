#ifndef DEF_FASTA
#define DEF_FASTA

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "inttypes.h"
#include "annotation.h"
#include "useful_macros.h"

static const uint8_t MASK_INSERT[3][4] = {{1, 4, 16, 64}, {2, 8, 32, 128}, {3, 12, 48, 192}};

int parseKmerCount(char* line, int size_kmer, uint8_t* tab, int pos_tab);
void parseSequenceBuffer(char* buf, uint8_t* tab, int* nb_kmers, int size_kmers, int nb_cell);

#endif
