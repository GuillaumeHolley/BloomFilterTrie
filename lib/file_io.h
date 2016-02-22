#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <libgen.h>
#include <omp.h>

#include "./../lib/insertNode.h"
#include "./../lib/printMemory.h"
#include "./../lib/fasta.h"
#include "./../lib/replaceAnnotation.h"

void insert_Genomes_from_KmerFiles(BFT_Root* root, int nb_files, char** filenames, int binary_files);
void insert_Genomes_from_FASTxFiles(BFT_Root* root, int nb_files, char** filenames);

int queryBFT_kmerPresences_from_KmerFiles(BFT_Root* root, char* query_filename, int binary_file, char* output_filename);
int queryBFT_kmerBranching_from_KmerFiles(BFT_Root* root, char* query_filename, int binary_file);
