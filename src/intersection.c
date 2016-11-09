
#include "intersection.h"

static __m128i shuffle_mask[16];

// a simple implementation, we don't care about performance here
void prepare_shuffling_dictionary() {
    for(int i = 0; i < 16; i++) {
        int counter = 0;
        char permutation[16];
        memset(permutation, 0xFF, sizeof(permutation));
        for(char b = 0; b < 4; b++) {
            if(getBit(i, b)) {
                permutation[counter++] = 4*b;
                permutation[counter++] = 4*b + 1;
                permutation[counter++] = 4*b + 2;
                permutation[counter++] = 4*b + 3;
            }
        }
        __m128i mask = _mm_loadu_si128((const __m128i*)permutation);
        shuffle_mask[i] = mask;
    }
}

int getBit(int value, int position) {
    return ( ( value & (1 << position) ) >> position);
}

uint32_t is_intersecting_uint32(uint32_t* list_a, uint32_t* list_b, uint32_t nb_inter_max){

    ASSERT_NULL_PTR(list_a, "is_intersecting_uint32()\n")
    ASSERT_NULL_PTR(list_b, "is_intersecting_uint32()\n")

    uint32_t i = 1, j = 1, nb = 0;
    uint32_t size_a = list_a[0]+1, size_b = list_b[0]+1;

    while (i < size_a && j < size_b){
        if (list_a[i] > list_b[j]) j++;
        else if (list_b[j] > list_a[i]) i++;
        else {
            nb++; i++; j++;
            if (nb > nb_inter_max) return nb;
        }
    }

    return nb;
}

uint32_t* intersection_uint32(uint32_t* list_a, uint32_t* list_b){

    ASSERT_NULL_PTR(list_a, "is_intersecting_uint32()\n")
    ASSERT_NULL_PTR(list_b, "is_intersecting_uint32()\n")

    uint32_t i = 1, j = 1;
    uint32_t size_a = list_a[0]+1, size_b = list_b[0]+1;

    uint32_t* list_c = malloc(MIN(size_a, size_b) * sizeof(uint32_t));
    ASSERT_NULL_PTR(list_c, "union_lists_uint32() 1\n")

    list_c[0] = 0;

    while (i < size_a && j < size_b){
        if (list_a[i] > list_b[j]) j++;
        else if (list_b[j] > list_a[i]) i++;
        else {
            list_c[0]++;
            list_c[list_c[0]] = list_a[i];
            i++; j++;
        }
    }

    return list_c;
}

int comp_uint32(const void *a, const void *b)  {
    if(*((uint32_t *)a) > *((uint32_t *)b))  return(+1);
    if(*((uint32_t *)a) < *((uint32_t *)b))  return(-1);
    return(0);
}

int comp_uint64(const void *a, const void *b)  {
    if(*((uint64_t *)a) > *((uint64_t *)b))  return(+1);
    if(*((uint64_t *)a) < *((uint64_t *)b))  return(-1);
    return(0);
}

/**
 * Taken almost verbatim from http://highlyscalable.wordpress.com/2012/06/05/fast-intersection-sorted-lists-sse/
 *
 * It is not safe for out to be either A or B.
 */
uint32_t* intersection_uint32_SIMD(uint32_t* list_a, uint32_t* list_b) {

    uint32_t* list_c = calloc(CEIL(MIN(list_a[0], list_b[0]) + 1, 4) * 4, sizeof(uint32_t));
    ASSERT_NULL_PTR(list_c, "intersection_uint32_SIMD()");

    uint32_t* A = &(list_a[1]);
    uint32_t* B = &(list_b[1]);
    uint32_t* C = &(list_c[1]);

    uint32_t count = 0;

    size_t i_a = 0, i_b = 0;

    // trim lengths to be a multiple of 4
    size_t st_a = (list_a[0] / 4) * 4;
    size_t st_b = (list_b[0] / 4) * 4;

    while (i_a < st_a && i_b < st_b) {

        //[ load segments of four 32-bit elements
        __m128i v_a = _mm_loadu_si128((__m128i *) &A[i_a]);
        __m128i v_b = _mm_loadu_si128((__m128i *) &B[i_b]);
        //]

        //[ move pointers
        const uint32_t a_max = A[i_a + 3];
        const uint32_t b_max = B[i_b + 3];
        i_a += (a_max <= b_max) * 4;
        i_b += (a_max >= b_max) * 4;
        //]

        //[ compute mask of common elements
        //const uint32_t cyclic_shift = _MM_SHUFFLE(0, 3, 2, 1);
        __m128i cmp_mask1 = _mm_cmpeq_epi32(v_a, v_b); // pairwise comparison
        v_b = _mm_shuffle_epi32(v_b, _MM_SHUFFLE(0, 3, 2, 1)); // shuffling
        __m128i cmp_mask2 = _mm_cmpeq_epi32(v_a, v_b); // again...
        v_b = _mm_shuffle_epi32(v_b, _MM_SHUFFLE(0, 3, 2, 1));
        __m128i cmp_mask3 = _mm_cmpeq_epi32(v_a, v_b); // and again...
        v_b = _mm_shuffle_epi32(v_b, _MM_SHUFFLE(0, 3, 2, 1));
        __m128i cmp_mask4 = _mm_cmpeq_epi32(v_a, v_b); // and again.
        __m128i cmp_mask = _mm_or_si128(_mm_or_si128(cmp_mask1, cmp_mask2), _mm_or_si128(cmp_mask3, cmp_mask4)); // OR-ing of comparison masks
        // convert the 128-bit mask to the 4-bit mask
        const int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp_mask));
        //]

        //[ copy out common elements
        __m128i p = _mm_shuffle_epi8(v_a, shuffle_mask[mask]);
        _mm_storeu_si128((__m128i*)&C[count], p);
        count += _mm_popcnt_u32(mask);
        //]
    }

    //qsort(C, count, sizeof(uint32_t), comp_uint32);

    // intersect the tail using scalar intersection
    while (i_a < list_a[0] && i_b < list_b[0]) {
        if (A[i_a] < B[i_b]) i_a++;
        else if (B[i_b] < A[i_a]) i_b++;
        else {
            C[count] = B[i_b];
            count++;
            i_a++;
            i_b++;
        }
    }

    list_c[0] = count;

    list_c = realloc(list_c, (count + 1) * sizeof(uint32_t));
    ASSERT_NULL_PTR(list_c, "intersection_uint32_SIMD()");

    return list_c;
}

uint64_t is_intersecting_uint8(uint8_t* list_a, uint64_t nb_elem_a, uint8_t* list_b, uint64_t nb_elem_b, int size_line, uint32_t nb_inter_max){

    ASSERT_NULL_PTR(list_a, "is_intersecting_uint8()\n")
    ASSERT_NULL_PTR(list_b, "is_intersecting_uint8()\n")

    uint64_t i = 0, j = 0, nb = 0;

    int cmp;

    nb_elem_a *= size_line;
    nb_elem_b *= size_line;

    while (i < nb_elem_a && j < nb_elem_b){

        cmp = memcmp(&(list_a[i]), &(list_b[j]), size_line * sizeof(uint8_t));

        if (cmp > 0) j += size_line;
        else if (cmp < 0) i += size_line;
        else {
            nb++;
            i += size_line;
            j += size_line;

            if (nb > nb_inter_max) return nb;
        }
    }

    return nb;
}

uint64_t is_intersecting_uint64(uint64_t* list_a, uint64_t* list_b, uint32_t nb_inter_max){

    ASSERT_NULL_PTR(list_a, "is_intersecting_uint64()\n")
    ASSERT_NULL_PTR(list_b, "is_intersecting_uint64()\n")

    uint64_t i = 1, j = 1, nb = 0;
    uint64_t size_a = list_a[0]+1, size_b = list_b[0]+1;

    while (i < size_a && j < size_b){
        if (list_a[i] > list_b[j]) j++;
        else if (list_b[j] > list_a[i]) i++;
        else {
            nb++; i++; j++;
            if (nb > nb_inter_max) return nb;
        }
    }

    return nb;
}

uint64_t is_intersecting_uint64_SIMD(uint64_t* list_a, uint64_t* list_b, uint64_t nb_inter_max) {

    uint64_t* A = &(list_a[1]);
    uint64_t* B = &(list_b[1]);

    size_t i_a = 0, i_b = 0, nb = 0;

    // trim lengths to be a multiple of 4
    size_t st_a = (list_a[0] / 2) * 2;
    size_t st_b = (list_b[0] / 2) * 2;

    while (i_a < st_a && i_b < st_b) {

        // Load segments of two 64-bit elements
        __m128i v_a = _mm_loadu_si128((__m128i *) &A[i_a]);
        __m128i v_b = _mm_loadu_si128((__m128i *) &B[i_b]);

        // Move pointers
        const uint64_t a_max = A[i_a + 1];
        const uint64_t b_max = B[i_b + 1];
        i_a += (a_max <= b_max) * 2;
        i_b += (a_max >= b_max) * 2;

        // Compute mask of common elements
        //const uint64_t cyclic_shift = _MM_SHUFFLE(1, 0, 3, 2);
        __m128i cmp_mask1 = _mm_cmpeq_epi64(v_a, v_b); // pairwise comparison
        v_b = _mm_shuffle_epi32(v_b, _MM_SHUFFLE(1, 0, 3, 2)); // shuffling
        __m128i cmp_mask2 = _mm_cmpeq_epi64(v_a, v_b); // again...
        __m128i cmp_mask = _mm_or_si128(cmp_mask1, cmp_mask2); // OR-ing of comparison masks

        // Convert the 128-bit mask to the 4-bit mask
        const int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp_mask));

        nb += _mm_popcnt_u64(mask)/* / sizeof(uint64_t)*/;
        if (nb > nb_inter_max) return nb;
    }

    // intersect the tail using scalar intersection
    while (i_a < list_a[0] && i_b < list_b[0]) {
        if (A[i_a] < B[i_b]) i_a++;
        else if (B[i_b] < A[i_a]) i_b++;
        else {
            nb++; i_a++; i_b++;
            if (nb > nb_inter_max) return nb;
        }
    }

    return nb;
}

uint64_t is_intersecting2x2(uint64_t* list_a, uint64_t* list_b, uint64_t* list_c, uint64_t* list_d, uint32_t nb_inter_max){

    ASSERT_NULL_PTR(list_a, "is_intersecting()\n")
    ASSERT_NULL_PTR(list_b, "is_intersecting()\n")
    ASSERT_NULL_PTR(list_c, "is_intersecting()\n")
    ASSERT_NULL_PTR(list_d, "is_intersecting()\n")

    if (list_a[0] != list_b[0]) ERROR("is_intersecting2x2(): size of list a is not the same as the one of list b.\n")
    if (list_c[0] != list_d[0]) ERROR("is_intersecting2x2(): size of list c is not the same as the one of list d.\n")

    uint64_t i = 1, j = 1, nb = 0;
    uint64_t size_a = list_a[0]+1, size_c = list_c[0]+1;

    while (i < size_a && j < size_c){
        if (list_a[i] > list_c[j]) j++;
        else if (list_c[j] > list_a[i]) i++;
        else if (list_b[i] > list_d[j]) j++;
        else if (list_d[j] > list_b[i]) i++;
        else{
            nb++; i++; j++;
            if (nb > nb_inter_max) return nb;
        }
    }

    return nb;
}

uint64_t* union_lists_uint64(uint64_t* list_a, uint64_t* list_b){

    ASSERT_NULL_PTR(list_a, "union_lists_uint64()\n")
    ASSERT_NULL_PTR(list_b, "union_lists_uint64()\n")

    int i = 1, j = 1;

    uint64_t size_a = list_a[0]+1, size_b = list_b[0]+1;

    uint64_t* list_c = malloc((size_a + size_b) * sizeof(uint64_t));
    ASSERT_NULL_PTR(list_c, "union_lists_uint64() 1\n")

    list_c[0] = 0;

    while (i < size_a && j < size_b){

        list_c[0]++;

        if (list_a[i] > list_b[j]){
            list_c[list_c[0]] = list_b[j];
            j++;
        }
        else if (list_b[j] > list_a[i]){
            list_c[list_c[0]] = list_a[i];
            i++;
        }
        else {
            list_c[list_c[0]] = list_a[i];
            i++; j++;
        }
    }

    while (i < size_a){
        list_c[0]++;
        list_c[list_c[0]] = list_a[i];
        i++;
    }

    while (j < size_b){
        list_c[0]++;
        list_c[list_c[0]] = list_b[j];
        j++;
    }

    list_c = realloc(list_c, (list_c[0] + 1) * sizeof(uint64_t));
    ASSERT_NULL_PTR(list_c, "union_lists_uint64() 2\n")

    return list_c;
}

uint32_t* union_lists_uint32(uint32_t* list_a, uint32_t* list_b){

    ASSERT_NULL_PTR(list_a, "union_lists_uint32()\n")
    ASSERT_NULL_PTR(list_b, "union_lists_uint32()\n")

    int i = 1, j = 1;

    uint32_t size_a = list_a[0], size_b = list_b[0];

    uint32_t* list_c = malloc((size_a + size_b) * sizeof(uint32_t));
    ASSERT_NULL_PTR(list_c, "union_lists_uint32() 1\n")

    list_c[0] = 0;

    while (i <= size_a && j <= size_b){

        list_c[0]++;

        if (list_a[i] > list_b[j]){
            list_c[list_c[0]] = list_b[j];
            j++;
        }
        else if (list_b[j] > list_a[i]){
            list_c[list_c[0]] = list_a[i];
            i++;
        }
        else {
            list_c[list_c[0]] = list_a[i];
            i++; j++;
        }
    }

    while (i <= size_a){
        list_c[0]++;
        list_c[list_c[0]] = list_a[i];
        i++;
    }

    while (j <= size_b){
        list_c[0]++;
        list_c[list_c[0]] = list_b[j];
        j++;
    }

    list_c = realloc(list_c, (list_c[0] + 1) * sizeof(uint32_t));
    ASSERT_NULL_PTR(list_c, "union_lists_uint32() 2\n")

    return list_c;
}
