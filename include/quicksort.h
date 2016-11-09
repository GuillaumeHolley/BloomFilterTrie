#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "useful_macros.h"

inline void swap_mem_uint8_t(uint8_t* a, uint8_t* b, int size_line_sub);
inline void swap_mem_int(int* a, int* b);
inline void swap_mem_uint64(uint64_t* a, uint64_t* b);

inline void swap_mem_uint8_t(uint8_t* a, uint8_t* b, int size_line_sub){

    uint8_t temp[size_line_sub];
    memcpy(&temp, a, size_line_sub);
    memmove(a, b, size_line_sub);
    memcpy(b, &temp, size_line_sub);
    return;
}

inline void swap_mem_int(int* a, int* b){

    int tmp = *a;
    *a = *b;
    *b = tmp;
    return;
}

inline void swap_mem_uint64(uint64_t* a, uint64_t* b){

    uint64_t tmp = *a;
    *a = *b;
    *b = tmp;
    return;
}

uint8_t* median3(uint8_t *a, int size_line, int left, int right, int* tab_ind);
void manualSort(uint8_t *a, int size_line, int left, int right, int* tab_ind);
void quicksort(uint8_t* substrings, int size_line_sub, int p, int r, int* tab_ind);
int* quicksort_init(uint8_t* substrings, int size_line_sub, int p, int r);

uint64_t median3_uint64(uint64_t *a, uint64_t* b, uint64_t left, uint64_t right);
void manualSort_uint64(uint64_t *a, uint64_t* b, uint64_t left, uint64_t right);
void quicksort_uint64(uint64_t* a, uint64_t* b, uint64_t p, uint64_t r);
