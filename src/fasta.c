#include "./../lib/fasta.h"

int parseKmerCount(const char* line, int size_kmer, uint8_t* tab, int pos_tab){
    ASSERT_NULL_PTR(line,"parseKmerCount()")
    ASSERT_NULL_PTR(tab,"parseKmerCount()")

    uint8_t* tab_tmp = &(tab[pos_tab]);

    int pos_car = 0;
    for (pos_car=0; pos_car < size_kmer; pos_car++){

        switch(line[pos_car]){
            case 'a':
            case 'A': break;
            case 'c':
            case 'C': tab_tmp[pos_car/4] |= MASK_INSERT[0][pos_car%4]; break;
            case 'g':
            case 'G': tab_tmp[pos_car/4] |= MASK_INSERT[1][pos_car%4]; break;
            case 'u':
            case 'U':
            case 't':
            case 'T': tab_tmp[pos_car/4] |= MASK_INSERT[2][pos_car%4]; break;
            case 'r':
            case 'R':
            case 'y':
            case 'Y':
            case 's':
            case 'S':
            case 'w':
            case 'W':
            case 'k':
            case 'K':
            case 'm':
            case 'M':
            case 'b':
            case 'B':
            case 'd':
            case 'D':
            case 'h':
            case 'H':
            case 'v':
            case 'V':
            case 'n':
            case 'N': memset(tab_tmp, 0, ((pos_car+1)/4)*sizeof(uint8_t)); goto END_LOOP;
            case '\n': break;
            default: memset(tab_tmp, 0, ((pos_car+1)/4)*sizeof(uint8_t)); goto END_LOOP;
        }
    }

    END_LOOP: if (pos_car == size_kmer) return 1;
    return 0;
}

void kmer_comp_to_ascii(const uint8_t* kmer_comp, int k, char* kmer){
    ASSERT_NULL_PTR(kmer_comp,"kmer_comp_to_ascii()")
    ASSERT_NULL_PTR(kmer,"kmer_comp_to_ascii()")

    uint8_t tmp;

    int i = 0, j = 0;

    for (i = 0; i < (k * 2) / SIZE_BITS_UINT_8T; i++){
        tmp = kmer_comp[i];
        *kmer = COMP_TO_ASCII[tmp & 0x3];
        kmer++; tmp >>= 2;
        *kmer = COMP_TO_ASCII[tmp & 0x3];
        kmer++; tmp >>= 2;
        *kmer = COMP_TO_ASCII[tmp & 0x3];
        kmer++; tmp >>= 2;
        *kmer = COMP_TO_ASCII[tmp & 0x3];
        kmer++;
    }

    if (k % 4){
        for (j = 0, tmp = kmer_comp[i]; j < k % 4; j++, kmer++, tmp >>= 2)
            *kmer = COMP_TO_ASCII[tmp & 0x3];
    }

    *kmer = '\0';

    return;
}

void parseSequenceBuffer(char* buf, uint8_t* tab, int* nb_kmers, int size_kmers, int nb_cell){
    int nb_ACGT_kmers = 0;

    ASSERT_NULL_PTR(buf,"parseSequenceBuffer()")
    ASSERT_NULL_PTR(tab,"parseSequenceBuffer()")

    int i = 0, k = 0;
    while (i < *nb_kmers){
        int pos_car = 0;
        for (pos_car=0; pos_car < size_kmers; pos_car++){
            switch(buf[i+pos_car]){
                case 'a':
                case 'A': break;
                case 'c':
                case 'C': tab[k + pos_car/4] |= MASK_INSERT[0][pos_car%4]; break;
                case 'g':
                case 'G': tab[k + pos_car/4] |= MASK_INSERT[1][pos_car%4]; break;
                case 'u':
                case 'U':
                case 't':
                case 'T': tab[k + pos_car/4] |= MASK_INSERT[2][pos_car%4]; break;
                case 'r':
                case 'R':
                case 'y':
                case 'Y':
                case 's':
                case 'S':
                case 'w':
                case 'W':
                case 'k':
                case 'K':
                case 'm':
                case 'M':
                case 'b':
                case 'B':
                case 'd':
                case 'D':
                case 'h':
                case 'H':
                case 'v':
                case 'V':
                case 'n':
                case 'N': memset(&(tab[k]), 0, (pos_car+1)/4); goto END_LOOP;
                case '\n': break;
                default: {  memset(&(tab[k]), 0, (pos_car+1)/4);
                            goto END_LOOP;
                        }
            }
        }

        k += nb_cell;
        nb_ACGT_kmers++;

        END_LOOP:i++;
    }

    *nb_kmers = nb_ACGT_kmers;

    return;
}

int parseKmerCount_IUPAC(char* line, int size_kmer, uint8_t* tab, int pos_tab){
    ASSERT_NULL_PTR(line,"parseKmerCount_IUPAC()")
    ASSERT_NULL_PTR(tab,"parseKmerCount_IUPAC()")

    int pos_car = 0;
    for (pos_car=0; pos_car < size_kmer; pos_car++){
        switch(line[pos_car]){
            case 'a':
            case 'A': break;
            case 'c':
            case 'C': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[0][pos_car%2]; break;
            case 'g':
            case 'G': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[1][pos_car%2]; break;
            case 'u':
            case 'U':
            case 't':
            case 'T': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[2][pos_car%2]; break;
            case 'r':
            case 'R': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[3][pos_car%2]; break;
            case 'y':
            case 'Y': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[4][pos_car%2]; break;
            case 's':
            case 'S': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[5][pos_car%2]; break;
            case 'w':
            case 'W': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[6][pos_car%2]; break;
            case 'k':
            case 'K': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[7][pos_car%2]; break;
            case 'm':
            case 'M': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[8][pos_car%2]; break;
            case 'b':
            case 'B': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[9][pos_car%2]; break;
            case 'd':
            case 'D': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[10][pos_car%2]; break;
            case 'h':
            case 'H': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[11][pos_car%2]; break;
            case 'v':
            case 'V': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[12][pos_car%2]; break;
            case 'n':
            case 'N': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[13][pos_car%2]; break;
            case '.':
            case '-': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC[14][pos_car%2]; break;
            case '\n': break;
            default: {  memset(&(tab[pos_tab]), 0, ((pos_car+1)/2)*sizeof(uint8_t));
                        goto END_LOOP;
                    }
        }
    }

    END_LOOP: if (pos_car == size_kmer) return 1;
    return 0;
}

void parseSequenceBuffer_IUPAC(char* buf, uint8_t* tab, int* nb_kmers, int size_kmers, int nb_cell){
    int nb_IUPAC_kmers = 0;

    ASSERT_NULL_PTR(buf,"parseSequenceBuffer_IUPAC()")
    ASSERT_NULL_PTR(tab,"parseSequenceBuffer_IUPAC()")

    int i = 0, k = 0;
    while (i < *nb_kmers){
        int pos_car = 0;
        for (pos_car=0; pos_car < size_kmers; pos_car++){
            switch(buf[i+pos_car]){
                case 'a':
                case 'A': break;
                case 'c':
                case 'C': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[0][pos_car%2]; break;
                case 'g':
                case 'G': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[1][pos_car%2]; break;
                case 'u':
                case 'U':
                case 't':
                case 'T': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[2][pos_car%2]; break;
                case 'r':
                case 'R': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[3][pos_car%2]; break;
                case 'y':
                case 'Y': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[4][pos_car%2]; break;
                case 's':
                case 'S': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[5][pos_car%2]; break;
                case 'w':
                case 'W': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[6][pos_car%2]; break;
                case 'k':
                case 'K': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[7][pos_car%2]; break;
                case 'm':
                case 'M': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[8][pos_car%2]; break;
                case 'b':
                case 'B': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[9][pos_car%2]; break;
                case 'd':
                case 'D': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[10][pos_car%2]; break;
                case 'h':
                case 'H': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[11][pos_car%2]; break;
                case 'v':
                case 'V': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[12][pos_car%2]; break;
                case 'n':
                case 'N': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[13][pos_car%2]; break;
                case '.':
                case '-': tab[k + pos_car/2] |= MASK_INSERT_IUPAC[14][pos_car%2]; break;
                case '\n': break;
                default: {  memset(&(tab[k]), 0, ((pos_car+1)/2)*sizeof(uint8_t));
                            goto END_LOOP;
                        }
            }
        }

        k += nb_cell;
        nb_IUPAC_kmers++;

        END_LOOP:i++;
    }

    *nb_kmers = nb_IUPAC_kmers;

    return;
}

int parseKmerCount_IUPAC_rev(char* line, int size_kmer, uint8_t* tab, int pos_tab){
    ASSERT_NULL_PTR(line,"parseKmerCount_IUPAC_rev()")
    ASSERT_NULL_PTR(tab,"parseKmerCount_IUPAC_rev()")

    int pos_car = 0;
    for (pos_car=0; pos_car < size_kmer; pos_car++){
        switch(line[pos_car]){
            case 'a':
            case 'A': break;
            case 'c':
            case 'C': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[0][pos_car%2]; break;
            case 'g':
            case 'G': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[1][pos_car%2]; break;
            case 'u':
            case 'U':
            case 't':
            case 'T': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[2][pos_car%2]; break;
            case 'r':
            case 'R': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[3][pos_car%2]; break;
            case 'y':
            case 'Y': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[4][pos_car%2]; break;
            case 's':
            case 'S': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[5][pos_car%2]; break;
            case 'w':
            case 'W': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[6][pos_car%2]; break;
            case 'k':
            case 'K': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[7][pos_car%2]; break;
            case 'm':
            case 'M': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[8][pos_car%2]; break;
            case 'b':
            case 'B': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[9][pos_car%2]; break;
            case 'd':
            case 'D': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[10][pos_car%2]; break;
            case 'h':
            case 'H': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[11][pos_car%2]; break;
            case 'v':
            case 'V': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[12][pos_car%2]; break;
            case 'n':
            case 'N': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[13][pos_car%2]; break;
            case '.':
            case '-': tab[pos_tab + pos_car/2] |= MASK_INSERT_IUPAC_REV[14][pos_car%2]; break;
            case '\n': break;
            default: {  memset(&(tab[pos_tab]), 0, ((pos_car+1)/2)*sizeof(uint8_t));
                        goto END_LOOP;
                    }
        }
    }

    END_LOOP: if (pos_car == size_kmer) return 1;
    return 0;
}

int parseKmerCount_xIUPAC(char* line, int size_kmer, uint8_t* kmers_arr_no_iupac, uint8_t* kmers_arr_iupac,
                          int pos_kmers_arr_no_iupac, int pos_kmers_arr_iupac){

    ASSERT_NULL_PTR(line,"parseKmerCount_xIUPAC()")
    ASSERT_NULL_PTR(kmers_arr_no_iupac,"parseKmerCount_xIUPAC()")
    ASSERT_NULL_PTR(kmers_arr_iupac,"parseKmerCount_xIUPAC()")

    if (strpbrk(line, "rRyYsSwWkKmMbBdDhHvVnN.-") != NULL){
        if (parseKmerCount_IUPAC(line, size_kmer, kmers_arr_iupac, pos_kmers_arr_iupac) < 1)
            ERROR("parseKmerCount_xIUPAC(): invalid character encountered in the input data")
        return 1;
    }
    else{
        if (parseKmerCount(line, size_kmer, kmers_arr_no_iupac, pos_kmers_arr_no_iupac) < 1)
            ERROR("parseKmerCount_xIUPAC(): invalid character encountered in the input data")
        return 0;
    }

    return 0;
}

int is_substring_IUPAC(char* line){

    ASSERT_NULL_PTR(line,"parseKmerCount_xIUPAC()")

    if (strpbrk(line, "rRyYsSwWkKmMbBdDhHvVnN.-") != NULL) return 1;
    return 0;
}

void reverse_complement(const char* s1, char* s2, int length){

    ASSERT_NULL_PTR(s1,"reverse_complement()")
    ASSERT_NULL_PTR(s2,"reverse_complement()")

    int pos, i;

    for (pos = length-1, i = 0; pos >= 0; pos--, i++){
        switch(s1[pos]){
            case 'a': s2[i] = 't'; break;
            case 'A': s2[i] = 'T'; break;
            case 'c': s2[i] = 'g'; break;
            case 'C': s2[i] = 'G'; break;
            case 'g': s2[i] = 'c'; break;
            case 'G': s2[i] = 'C'; break;
            case 'u': s2[i] = 'a'; break;
            case 'U': s2[i] = 'A'; break;
            case 't': s2[i] = 'a'; break;
            case 'T': s2[i] = 'A'; break;
            case 'n': s2[i] = 'n'; break;
            case 'N': s2[i] = 'N'; break;
            default: {
                printf("%s\n", s1);
                ERROR("reverse_complement(): encountered an non-IUPAC character");
            }
        }
    }

    return;
}
