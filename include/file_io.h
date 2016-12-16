#pragma once

#define _XOPEN_SOURCE 500
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>

#include <libgen.h>
#include <omp.h>

#include "bft.h"
#include "getRSS.h"
#include "insertNode.h"
#include "printMemory.h"
#include "fasta.h"
#include "replaceAnnotation.h"

void insert_Genomes_from_KmerFiles(BFT_Root* root, int nb_files, char** filenames, int binary_files, char* filename_bft);
void insert_Genomes_from_FASTxFiles(BFT_Root* root, int nb_files, char** filenames);

int queryBFT_kmerPresences_from_KmerFiles(BFT_Root* root, char* query_filename, int binary_file, char* output_filename);
int queryBFT_kmerBranching_from_KmerFiles(BFT_Root* root, char* query_filename, int binary_file);

void compress_annotations_disk(BFT_Root* bft, char* filename_bft);

void query_sequences_outputCSV(BFT_Root* root, char* query_filename, char* output_filename, double threshold, bool canonical_search);

/*void par_insert_Genomes_from_KmerFiles(int nb_files, char** filenames, int binary_files, int size_kmer,
                                       int treshold_compression, char* prefix_output, int cut_lvl, int memory_limit);*/
