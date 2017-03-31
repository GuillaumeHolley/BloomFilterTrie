#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include "annotation.h"
#include "useful_macros.h"

static const char COMP_TO_ASCII[4] = {'A', 'C', 'G', 'T'};

static const char COMP_IUPAC_TO_ASCII[15] = {'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'};

static const uint8_t MASK_INSERT[3][4] = {{1, 4, 16, 64}, {2, 8, 32, 128}, {3, 12, 48, 192}};


static const uint8_t MASK_INSERT_IUPAC[15][2] = {{1, 16}, {2, 32}, {3, 48}, {4, 64}, {5, 80},
                                                {6, 96}, {7, 112}, {8, 128}, {9, 144}, {10, 160},
                                                {11, 176}, {12, 192}, {13, 208}, {14, 224}, {15, 240}};

static const uint8_t MASK_INSERT_IUPAC_REV[15][2] = {{16, 1}, {32, 2}, {48, 3}, {64, 4}, {80, 5},
                                                {96, 6}, {112, 7}, {128, 8}, {144, 9}, {160, 10},
                                                {176, 11}, {192, 12}, {208, 13}, {224, 14}, {240, 15}};

int parseKmerCount(const char* line, int size_kmer, uint8_t* tab, int pos_tab);
void kmer_comp_to_ascii(const uint8_t* kmer_comp, int k, char* kmer);
void kmer_iupac_comp_to_ascii(const uint8_t* kmer_comp, int k, char* kmer);
void parseSequenceBuffer(char* buf, uint8_t* tab, int* nb_kmers, int size_kmers, int nb_cell);

int parseKmerCount_IUPAC(char* line, int size_kmer, uint8_t* tab, int pos_tab);
int parseKmerCount_IUPAC_rev(char* line, int size_kmer, uint8_t* tab, int pos_tab);
void parseSequenceBuffer_IUPAC(char* buf, uint8_t* tab, int* nb_kmers, int size_kmers, int nb_cell);

int parseKmerCount_xIUPAC(char* line, int length_to_parse, uint8_t* kmers_arr_no_iupac, uint8_t* kmers_arr_iupac,
                          int pos_kmers_arr_no_iupac, int pos_kmers_arr_iupac);

int is_substring_IUPAC(const char* line);
int is_substring_nonACGT(const char* str, int len, bool check_last_char_only);

void reverse_complement(const char* s1, char* s2, int length);
char reverse_complement_char(char c);
