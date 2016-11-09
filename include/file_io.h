#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <libgen.h>
#include <omp.h>

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

/*void par_insert_Genomes_from_KmerFiles(int nb_files, char** filenames, int binary_files, int size_kmer,
                                       int treshold_compression, char* prefix_output, int cut_lvl, int memory_limit);*/
