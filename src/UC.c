#include "./../lib/UC.h"

extern void reinit_annotation_inform(annotation_inform* ann_inf);
extern UC* createUC();
extern void initiateUC(UC* uc);
extern void freeUC(UC* uc);

void insertKmer_UC(UC*  uc, uint8_t*  kmer, uint32_t id_genome, int size_id_genome, int size_kmer, int pos_insertion,
                   annotation_inform* ann_inf, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(uc, "insertKmer_UC()")

    if (id_genome > NB_MAX_ID_GENOMES)
        ERROR("insertKmer_UC(): insertion of a genome with unauthorized ID (<1 or >NB_MAX_ID_GENOMES)")

    int nb_cell = CEIL(size_kmer*2, SIZE_BITS_UINT_8T);
    int size_line = nb_cell;
    int nb_elem = uc->nb_children >> 1;
    int beginning_cluster = uc->nb_children & 0x1;

    compute_best_mode(ann_inf, annot_sorted, NULL, 0, NULL, 0, id_genome, size_id_genome);

    if (nb_elem == 0){

        uc->suffixes = calloc(nb_cell + ann_inf->min_size, sizeof(uint8_t));
        ASSERT_NULL_PTR(uc->suffixes, "insertKmer_UC()")

        uc->nb_children = 0x2 | beginning_cluster;

        uc->size_annot = ann_inf->min_size;
        memcpy(uc->suffixes, kmer, nb_cell * sizeof(uint8_t));
        modify_mode_annotation(ann_inf, &(uc->suffixes[nb_cell]),
                               uc->size_annot, NULL, 0, id_genome, size_id_genome);
    }
    else{

        uint8_t* extend_annot = NULL;

        if (ann_inf->min_size > uc->size_annot){
            extend_annot = realloc_annotation(uc, nb_cell, nb_elem, ann_inf->min_size, 1, pos_insertion);
            size_line += uc->size_annot;
        }
        else{
            size_line += uc->size_annot;
            beginning_cluster = uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                + uc->nb_cplx_nodes * (uc->size_annot_cplx_nodes + SIZE_BYTE_CPLX_N);

            uc->suffixes = realloc(uc->suffixes,
                                   ((nb_elem+1) * size_line + beginning_cluster) * sizeof(uint8_t));

            memmove(&(uc->suffixes[(pos_insertion+1) * size_line]), &(uc->suffixes[pos_insertion * size_line]),
                    (beginning_cluster + (nb_elem - pos_insertion) * size_line) * sizeof(uint8_t));

            shift_extended_annot(uc, nb_cell, nb_elem+1, pos_insertion);
        }

        if (uc->suffixes != NULL){
            size_line = pos_insertion * size_line;

            memcpy(&(uc->suffixes[size_line]), kmer, nb_cell * sizeof(uint8_t));
            memset(&(uc->suffixes[size_line + nb_cell]), 0, uc->size_annot * sizeof(uint8_t));

            modify_mode_annotation(ann_inf, &(uc->suffixes[size_line + nb_cell]),
                                   uc->size_annot, extend_annot, 1, id_genome, size_id_genome);

            uc->nb_children += 2;
        }
        else ERROR("insertKmer_UC(): out of memory\n")
    }

    reinit_annotation_inform(ann_inf);

    return;
}

int binary_search_UC(const UC* uc, int pos_start, int pos_end, const uint8_t* suf, int size_suf_byte, uint8_t mask_for_last_byte){

    ASSERT_NULL_PTR(uc, "binary_search_UC()")
    ASSERT_NULL_PTR(suf, "binary_search_UC()")

    if (pos_start > pos_end) ERROR("binary_search_UC()")
    if (uc->suffixes == NULL) return 0;

    int imid;

    int imin = pos_start;
    int imax = pos_end;
    int size_line = size_suf_byte + uc->size_annot;

    if (mask_for_last_byte == 0xff){

        while (imin < imax){

            imid = imin + (imax-imin)/2;

            if (memcmp(&(uc->suffixes[imid * size_line]), suf, size_suf_byte * sizeof(uint8_t)) < 0) imin = imid+1;
            else imax = imid;
        }
    }
    else{

        int cmp = 0;

        while (imin < imax){

            imid = imin + (imax-imin)/2;

            cmp = memcmp(&(uc->suffixes[imid * size_line]), suf, (size_suf_byte-1) * sizeof(uint8_t));

            if (cmp < 0) imin = imid+1;
            else if ((cmp == 0) && ((uc->suffixes[imid * size_line + size_suf_byte - 1] & mask_for_last_byte) < suf[size_suf_byte - 1])){
                imin = imid+1;
            }
            else imax = imid;
        }
    }

    return imin;
}

int binary_search_UC_array(const uint8_t* uc_array, int size_annot, int pos_start, int pos_end, const uint8_t* suf,
                           int size_suf_byte, uint8_t mask_for_last_byte){

    ASSERT_NULL_PTR(suf, "binary_search_UC()")

    if (pos_start > pos_end) ERROR("binary_search_UC()")
    if (uc_array == NULL) return 0;

    int imid;

    int imin = pos_start;
    int imax = pos_end;
    int size_line = size_suf_byte + size_annot;

    if (mask_for_last_byte == 0xff){

        while (imin < imax){

            imid = imin + (imax-imin)/2;

            if (memcmp(&(uc_array[imid * size_line]), suf, size_suf_byte * sizeof(uint8_t)) < 0) imin = imid+1;
            else imax = imid;
        }
    }
    else{

        int cmp = 0;

        while (imin < imax){

            imid = imin + (imax-imin)/2;

            cmp = memcmp(&(uc_array[imid * size_line]), suf, (size_suf_byte-1) * sizeof(uint8_t));

            if (cmp < 0) imin = imid+1;
            else if ((cmp == 0) && ((uc_array[imid * size_line + size_suf_byte - 1] & mask_for_last_byte) < suf[size_suf_byte - 1])){
                imin = imid+1;
            }
            else imax = imid;
        }
    }

    return imin;
}
