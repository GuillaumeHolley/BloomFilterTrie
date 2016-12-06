#include "fasta.h"

int parseKmerCount(const char* line, int size_kmer, uint8_t* tab, int pos_tab){
    ASSERT_NULL_PTR(line,"parseKmerCount()")
    ASSERT_NULL_PTR(tab,"parseKmerCount()")

    uint8_t* tab_tmp = &tab[pos_tab];

    int pos_car = 0;
    for (; pos_car < size_kmer; pos_car++){

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

void kmer_iupac_comp_to_ascii(const uint8_t* kmer_comp, int k, char* kmer){
    ASSERT_NULL_PTR(kmer_comp,"parseKmerCount()")
    ASSERT_NULL_PTR(kmer,"parseKmerCount()")

    int i = 0;

    for (i = 0; i < (k * 4) / SIZE_BITS_UINT_8T; i++){
        *kmer = COMP_IUPAC_TO_ASCII[kmer_comp[i] & 0xf];
        kmer++;
        *kmer = COMP_IUPAC_TO_ASCII[kmer_comp[i] >> 4];
        kmer++;
    }

    if (k % 2){
        *kmer = COMP_IUPAC_TO_ASCII[kmer_comp[i] & 0xf];
        kmer++;
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
    for (pos_car = 0; pos_car < size_kmer; pos_car++){

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

    ASSERT_NULL_PTR(buf,"parseSequenceBuffer_IUPAC()")
    ASSERT_NULL_PTR(tab,"parseSequenceBuffer_IUPAC()")

    int nb_IUPAC_kmers = 0;
    int i = 0, k = 0;

    while (i < *nb_kmers){

        for (int pos_car = 0; pos_car < size_kmers; pos_car++){

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

int parseKmerCount_xIUPAC(char* line, int length_to_parse, uint8_t* kmers_arr_no_iupac, uint8_t* kmers_arr_iupac,
                          int pos_kmers_arr_no_iupac, int pos_kmers_arr_iupac){

    ASSERT_NULL_PTR(line,"parseKmerCount_xIUPAC()")
    ASSERT_NULL_PTR(kmers_arr_no_iupac,"parseKmerCount_xIUPAC()")
    ASSERT_NULL_PTR(kmers_arr_iupac,"parseKmerCount_xIUPAC()")

    const char* pos = strpbrk(line, "rRyYsSwWkKmMbBdDhHvVnN.-");

    if ((pos != NULL) && (pos < line + length_to_parse)){
        if (parseKmerCount_IUPAC(line, length_to_parse, kmers_arr_iupac, pos_kmers_arr_iupac) < 1)
            ERROR("parseKmerCount_xIUPAC(): invalid character encountered in the input data")
        return 1;
    }
    else{
        if (parseKmerCount(line, length_to_parse, kmers_arr_no_iupac, pos_kmers_arr_no_iupac) < 1)
            ERROR("parseKmerCount_xIUPAC(): invalid character encountered in the input data")
        return 0;
    }

    return 0;
}

int is_substring_IUPAC(char* line){

    ASSERT_NULL_PTR(line,"is_substring_IUPAC()\n")

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
            case 'm': s2[i] = 'k'; break;
            case 'M': s2[i] = 'K'; break;
            case 'r': s2[i] = 'y'; break;
            case 'R': s2[i] = 'Y'; break;
            case 'w': s2[i] = 'w'; break;
            case 'W': s2[i] = 'W'; break;
            case 's': s2[i] = 's'; break;
            case 'S': s2[i] = 'S'; break;
            case 'y': s2[i] = 'r'; break;
            case 'Y': s2[i] = 'R'; break;
            case 'k': s2[i] = 'm'; break;
            case 'K': s2[i] = 'M'; break;
            case 'v': s2[i] = 'b'; break;
            case 'V': s2[i] = 'B'; break;
            case 'h': s2[i] = 'd'; break;
            case 'H': s2[i] = 'D'; break;
            case 'd': s2[i] = 'h'; break;
            case 'D': s2[i] = 'H'; break;
            case 'b': s2[i] = 'v'; break;
            case 'B': s2[i] = 'V'; break;
            case 'n': s2[i] = 'n'; break;
            case 'N': s2[i] = 'N'; break;
            default: {
                printf("%s\n", s1);
                ERROR("reverse_complement(): encountered an non-IUPAC character");
            }
        }
    }

    s2[length] = '\0';

    return;
}

char reverse_complement_char(char c){

    switch(c){

        case 'a': return 't';
        case 'A': return 'T';
        case 'c': return 'g';
        case 'C': return 'G';
        case 'g': return 'c';
        case 'G': return 'C';
        case 'u': return 'a';
        case 'U': return 'A';
        case 't': return 'a';
        case 'T': return 'A';
        case 'm': return 'k';
        case 'M': return 'K';
        case 'r': return 'y';
        case 'R': return 'Y';
        case 'w': return 'w';
        case 'W': return 'W';
        case 's': return 's';
        case 'S': return 'S';
        case 'y': return 'r';
        case 'Y': return 'R';
        case 'k': return 'm';
        case 'K': return 'M';
        case 'v': return 'b';
        case 'V': return 'B';
        case 'h': return 'd';
        case 'H': return 'D';
        case 'd': return 'h';
        case 'D': return 'H';
        case 'b': return 'v';
        case 'B': return 'V';
        case 'n': return 'n';
        case 'N': return 'N';
        default: {
            printf("%c\n", c);
            ERROR("reverse_complement(): encountered an non-IUPAC character");
        }
    }

    return c;
}
