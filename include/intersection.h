#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <immintrin.h>

#include "useful_macros.h"

void prepare_shuffling_dictionary();
int getBit(int value, int position);

int comp_uint32(const void *a, const void *b);
int comp_uint64(const void *a, const void *b);

uint32_t is_intersecting_uint32(uint32_t* list_a, uint32_t* list_b, uint32_t nb_inter_max);

uint32_t* intersection_uint32(uint32_t* list_a, uint32_t* list_b);
uint32_t* intersection_uint32_SIMD(uint32_t* list_a, uint32_t* list_b);

uint64_t is_intersecting_uint64(uint64_t* list_a, uint64_t* list_b, uint32_t nb_inter_max);
uint64_t is_intersecting_uint64_SIMD(uint64_t* list_a, uint64_t* list_b, uint64_t nb_inter_max);

uint64_t is_intersecting_uint8(uint8_t* list_a, uint64_t nb_elem_a, uint8_t* list_b, uint64_t nb_elem_b, int size_line, uint32_t nb_inter_max);

uint64_t is_intersecting2x2(uint64_t* list_a, uint64_t* list_b, uint64_t* list_c, uint64_t* list_d, uint32_t nb_inter_max);

uint32_t* union_lists_uint32(uint32_t* list_a, uint32_t* list_b);
