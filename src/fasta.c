#include "./../lib/fasta.h"

int parseKmerCount(char* line, int size_kmer, uint8_t* tab, int pos_tab){
    ASSERT_NULL_PTR(line,"parseKmerCount()")
    ASSERT_NULL_PTR(tab,"parseKmerCount()")

    int pos_car = 0;
    for (pos_car=0; pos_car < size_kmer; pos_car++){
        switch(line[pos_car]){
            case 'a':
            case 'A': break;
            case 'c':
            case 'C': tab[pos_tab + pos_car/4] |= MASK_INSERT[0][pos_car%4]; break;
            case 'g':
            case 'G': tab[pos_tab + pos_car/4] |= MASK_INSERT[1][pos_car%4]; break;
            case 'u':
            case 'U':
            case 't':
            case 'T': tab[pos_tab + pos_car/4] |= MASK_INSERT[2][pos_car%4]; break;
            case 'm':
            case 'M':
            case 'r':
            case 'R':
            case 'w':
            case 'W':
            case 's':
            case 'S':
            case 'y':
            case 'Y':
            case 'k':
            case 'K':
            case 'v':
            case 'V':
            case 'h':
            case 'H':
            case 'd':
            case 'D':
            case 'b':
            case 'B':
            case 'n':
            case 'N': memset(&(tab[pos_tab]), 0, ((pos_car+1)/4)*sizeof(uint8_t)); goto END_LOOP;
            case '\n': break;
            default: {  memset(&(tab[pos_tab]), 0, ((pos_car+1)/4)*sizeof(uint8_t));
                        goto END_LOOP;
                    }
        }
    }

    END_LOOP: if (pos_car == size_kmer) return 1;
    return 0;
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
                case 'm':
                case 'M':
                case 'r':
                case 'R':
                case 'w':
                case 'W':
                case 's':
                case 'S':
                case 'y':
                case 'Y':
                case 'k':
                case 'K':
                case 'v':
                case 'V':
                case 'h':
                case 'H':
                case 'd':
                case 'D':
                case 'b':
                case 'B':
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
