#pragma once

#define SIZE_BITS_UINT_8T 8
#define KMER_LENGTH_MAX 126

#define SIZE_BYTE_EXT_ANNOT 3
#define SIZE_BYTE_CPLX_N 2

#define NB_MAX_ID_GENOMES 100000000
#define SIZE_MAX_BYTE_ANNOT CEIL(NB_MAX_ID_GENOMES+2,SIZE_BITS_UINT_8T)

#define NB_CHAR_SUF_PREF 9 //Length of suffix prefixes (in characters) stored in CCs
#define SIZE_BYTES_SUF_PREF CEIL(NB_CHAR_SUF_PREF*2,SIZE_BITS_UINT_8T) //Length of suffix prefixes (in bytes) stored in CCs

#define SIZE_FILTER2_DEFAULT 1024//Nb cells in filter2 by default (when s=8), SIZE_FILTER2_DEFAULT % 8 = 0

#define NB_KMERS_PER_UC126 255 //Maximum capacity of UCs containing suffixes of length 45 is NB_KMERS_PER_CC45
#define NB_KMERS_PER_UC117 255 //Maximum capacity of UCs containing suffixes of length 36 is NB_KMERS_PER_CC36
#define NB_KMERS_PER_UC108 255 //Maximum capacity of UCs containing suffixes of length 27 is NB_KMERS_PER_CC27
#define NB_KMERS_PER_UC99 255 //Maximum capacity of UCs containing suffixes of length 18 is NB_KMERS_PER_CC18
#define NB_KMERS_PER_UC90 255 //Maximum capacity of UCs containing suffixes of length 9 is NB_KMERS_PER_CC9
#define NB_KMERS_PER_UC81 255 //Maximum capacity of UCs containing suffixes of length 9 is NB_KMERS_PER_CC
#define NB_KMERS_PER_UC72 255 //Maximum capacity of UCs containing suffixes of length 63 is NB_KMERS_PER_CC63
#define NB_KMERS_PER_UC63 255 //Maximum capacity of UCs containing suffixes of length 63 is NB_KMERS_PER_CC63
#define NB_KMERS_PER_UC54 255 //Maximum capacity of UCs containing suffixes of length 54 is NB_KMERS_PER_CC54
#define NB_KMERS_PER_UC45 255 //Maximum capacity of UCs containing suffixes of length 45 is NB_KMERS_PER_CC45
#define NB_KMERS_PER_UC36 255 //Maximum capacity of UCs containing suffixes of length 36 is NB_KMERS_PER_CC36
#define NB_KMERS_PER_UC27 255 //Maximum capacity of UCs containing suffixes of length 27 is NB_KMERS_PER_CC27
#define NB_KMERS_PER_UC18 255 //Maximum capacity of UCs containing suffixes of length 18 is NB_KMERS_PER_CC18
#define NB_KMERS_PER_UC9 255 //Maximum capacity of UCs containing suffixes of length 9 is NB_KMERS_PER_CC9
#define NB_KMERS_PER_UC 255 //Maximum capacity of UCs containing suffixes of length 9 is NB_KMERS_PER_CC

#define NB_UC_PER_SKP 128 //SHOULD BE <= 248 && NB_UC_PER_SKP%8 = 0
#define NB_BITS_IN_CELL_SKIP_FILTER2 NB_UC_PER_SKP //SHOULD BE <= 248 & NB_BITS_IN_CELL_SKIP_FILTER2%8 = 0
#define NB_BITS_IN_CELL_SKIP_FILTER3 NB_UC_PER_SKP //SHOULD BE <= 248 & NB_BITS_IN_CELL_SKIP_FILTER3%8 = 0

#define NB_BYTES_IN_CELL_SKIP_FILTER2 CEIL(NB_BITS_IN_CELL_SKIP_FILTER2, SIZE_BITS_UINT_8T)
#define NB_BYTES_IN_CELL_SKIP_FILTER3 CEIL(NB_BITS_IN_CELL_SKIP_FILTER3, SIZE_BITS_UINT_8T)

#define SIZE_CLUST_SKIP_NODES NB_UC_PER_SKP //Value max

#define MODULO_HASH 1504 //Size bloom filter in bits, MODULO_HASH%8 = 0
#define TRESH_SUF_PREF 3584 //if a CC has exactly NB_SUBSTRINGS_TRANSFORM prefixes, it is recomputed with transform_Filter2n3()

#define SIZE_BUFFER 4096
#define PRINT_EVERY_X_KMERS 1000000
