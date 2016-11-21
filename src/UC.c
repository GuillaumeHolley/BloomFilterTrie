#include "UC.h"

extern inline UC* createUC();
extern inline void initiateUC(UC* uc);
extern inline void initializeUC(UC* uc);
extern inline void resetUC(UC* uc);
extern inline void freeUC(UC* uc);

extern inline UC_SIZE_ANNOT_CPLX_T * min_size_per_annot_cplx(UC* uc, int nb_substrings, int size_substring);
extern inline UC_SIZE_ANNOT_CPLX_T * min_size_per_annot_cplx_sub(UC* uc, int nb_substrings, int size_substring, int pos_start, int pos_end);
extern inline int max_size_annot_cplx_sub(uint8_t* annot_cplx, int nb_cplx, int size_cplx, int pos_start, int pos_end);

void insertKmer_UC(UC* uc, uint8_t*  kmer, uint32_t id_genome, int size_id_genome, int size_kmer, int pos_insertion,
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

int get_annot(UC* uc, uint8_t** annot, uint8_t** annot_ext, uint8_t** annot_cplx,
                   int* size_annot, int* size_annot_cplx, int size_substring, int nb_substring, int position){

    if (position >= nb_substring) ERROR("get_annot(): position >= nb_substring")

    ASSERT_NULL_PTR(uc,"get_annot()")
    ASSERT_NULL_PTR(uc->suffixes,"get_annot()")

    uint8_t* ext_annot = NULL;
    uint8_t* cplx_annot = NULL;

    if (uc->size_annot == 0){ // No standard annotation

        *annot = NULL;
        *size_annot = 0;

        if ((ext_annot = get_extend_annot(uc, size_substring, nb_substring, position)) == NULL){ // No extended annotation

            *annot_ext = NULL;

            if ((cplx_annot = get_annot_cplx_nodes(uc, size_substring, nb_substring, position)) == NULL){ // No complex annotation

                *annot_cplx = NULL;
                *size_annot_cplx = 0;

                return 0;
            }
            else {
                *annot_cplx = cplx_annot;
                *size_annot_cplx = size_annot_sub(*annot_cplx, 0, uc->size_annot_cplx_nodes);
            }
        }
        else{
            *annot_ext = ext_annot;
            *size_annot_cplx = 0;
        }
    }
    else{

        *annot = &(uc->suffixes[position * (size_substring + uc->size_annot) + size_substring]);

        if ((ext_annot = get_extend_annot(uc, size_substring, nb_substring, position)) == NULL){// No extended annotation

            *annot_ext = NULL;

            if ((cplx_annot = get_annot_cplx_nodes(uc, size_substring, nb_substring, position)) == NULL){ // No complex annotation

                *annot_cplx = NULL;
                *size_annot_cplx = 0;
                *size_annot = size_annot_sub(*annot, 0, uc->size_annot);

                if (*size_annot == 0) return 0;
            }
            else {
                *annot_cplx = cplx_annot;
                *size_annot_cplx = size_annot_sub(*annot_cplx, 0, uc->size_annot_cplx_nodes);
                *size_annot = 0;
            }
        }
        else{
            *annot_ext = ext_annot;
            *annot_cplx = NULL;
            *size_annot = uc->size_annot;
            *size_annot_cplx = 0;
        }
    }

    return 1;
}

void get_annots(UC* uc, uint8_t*** annots, uint8_t*** annots_ext, uint8_t*** annots_cplx,
                   int** size_annots, int** size_annots_cplx, int size_substring, int nb_substring,
                   int position_start, int position_end){

    if (nb_substring <= 0) return;

    if ((position_start >= nb_substring) || (position_start < 0)){
        ERROR("get_annots(): position_start >= nb_substring or < 0")
    }

    if ((position_end >= nb_substring) || (position_end < 0)){
        ERROR("get_annots(): position_end >= nb_substring or < 0")
    }

    ASSERT_NULL_PTR(uc,"get_annot()")
    ASSERT_NULL_PTR(uc->suffixes,"get_annot()")

    int size_line;

    int i = 0;
    int nb_suffixes = position_end - position_start + 1;

    *annots = malloc(nb_suffixes * sizeof(uint8_t*));
    ASSERT_NULL_PTR(*annots, "get_annots()")

    *size_annots_cplx = calloc(nb_suffixes, sizeof(int));
    ASSERT_NULL_PTR(*size_annots_cplx, "get_annots()")

    *size_annots = calloc(nb_suffixes, sizeof(int));
    ASSERT_NULL_PTR(*size_annots, "get_annots()")

    *annots_ext = NULL;
    *annots_cplx = NULL;

    if (uc->size_annot == 0){

        for (i = 0; i < nb_suffixes; i++) (*annots)[i] = NULL;

        if ((*annots_ext = get_extend_annots(uc, size_substring, nb_substring,
                                             position_start, position_end)) == NULL){

            if ((*annots_cplx = get_annots_cplx_nodes(uc, size_substring, nb_substring,
                                                     position_start, position_end)) != NULL){
                free(*size_annots_cplx);
                *size_annots_cplx = min_size_per_annot_cplx_sub(uc, size_substring, nb_substring,
                                                                position_start, position_end);
            }
        }
    }
    else{

        size_line = size_substring + uc->size_annot;

        for (i = position_start; i <= position_end; i++)
            (*annots)[i - position_start] = &(uc->suffixes[i * size_line + size_substring]);

        if ((*annots_ext = get_extend_annots(uc, size_substring, nb_substring,
                                             position_start, position_end)) == NULL){

            for (i = 0; i < nb_suffixes; i++)
                (*size_annots)[i] = size_annot_sub((*annots)[i], 0, uc->size_annot);

            if ((*annots_cplx = get_annots_cplx_nodes(uc, size_substring, nb_substring,
                                                     position_start, position_end)) != NULL){

                for (i = 0; i < nb_suffixes; i++)
                    (*size_annots_cplx)[i] = size_annot_sub((*annots_cplx)[i], 0, uc->size_annot_cplx_nodes);
            }
        }
        else{
            for (i = 0; i < nb_suffixes; i++){
                if ((*annots_ext)[i] != NULL) (*size_annots)[i] = uc->size_annot;
                else (*size_annots)[i] = size_annot_sub((*annots)[i], 0, uc->size_annot);
            }
        }
    }

    return;
}

void shift_extended_annot(UC* uc, int size_substring, int nb_substring, int pos_insert){

    if (uc->nb_extended_annot != 0){

        uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot)]);

        int end_ext = uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;
        int pos = 0;
        int i = 0;

        pos += ((((uint16_t)extend_annot[0]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[1]));

        while (pos < pos_insert){
            i += SIZE_BYTE_EXT_ANNOT;
            if (i == end_ext) break;
            pos += ((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1]));
        }

        if (i < end_ext){
            if (extend_annot[i+1] != 0xff) extend_annot[i+1]++;
            else {
                extend_annot[i]++;
                extend_annot[i+1] = 0x0;
            }
        }
    }
}

void insert_extend_annot(UC* uc, int size_substring, int nb_substring, int pos_insert, uint8_t annot, int shift_or_not){

    ASSERT_NULL_PTR(uc,"insert_extend_annot()")

    uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot)]);
    int tot = uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

    if (uc->nb_extended_annot != 0){

        int pos = 0;
        int i = 0;
        int old_pos = 0;
        int delta;
        int sum;

        pos += ((((uint16_t)extend_annot[0]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[1]));

        while (pos < pos_insert){
            i++;
            if (i == uc->nb_extended_annot) break;
            old_pos = pos;
            pos += ((((uint16_t)extend_annot[i*SIZE_BYTE_EXT_ANNOT]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i*SIZE_BYTE_EXT_ANNOT+1]));
        }

        delta = pos_insert - pos;

        if (i == 0){
            memmove(&(extend_annot[SIZE_BYTE_EXT_ANNOT]), &(extend_annot[0]), (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + tot) * sizeof(uint8_t));

            extend_annot[0] = pos_insert >> SIZE_BITS_UINT_8T;
            extend_annot[1] = pos_insert & 0xff;
            extend_annot[2] = annot;

            delta = pos - pos_insert;
        }
        else if (i == uc->nb_extended_annot) extend_annot[(i+1)*SIZE_BYTE_EXT_ANNOT-1] = annot;
        else {
            sum = i*SIZE_BYTE_EXT_ANNOT;
            memmove(&(extend_annot[sum+SIZE_BYTE_EXT_ANNOT]), &(extend_annot[sum]), ((uc->nb_extended_annot - i) * SIZE_BYTE_EXT_ANNOT + tot) * sizeof(uint8_t));

            delta = pos_insert - old_pos;

            extend_annot[sum] = delta >> SIZE_BITS_UINT_8T;
            extend_annot[sum+1] = delta & 0xff;
            extend_annot[sum+2] = annot;

            delta = pos - pos_insert;
        }

        if (i != uc->nb_extended_annot){
            if (shift_or_not == 1) delta++;
        }
        else i--;

        extend_annot[(i+1)*SIZE_BYTE_EXT_ANNOT] = delta >> SIZE_BITS_UINT_8T;
        extend_annot[(i+1)*SIZE_BYTE_EXT_ANNOT+1] = delta & 0xff;

        uc->nb_extended_annot++;
    }
    else{
        extend_annot[0] = pos_insert >> SIZE_BITS_UINT_8T;
        extend_annot[1] = pos_insert & 0xff;
        extend_annot[2] = annot;
        uc->nb_extended_annot = 1;
    }

    return;
}

void delete_extend_annots(UC* uc, int size_substring, int nb_substring, int pos_sub_start, int pos_sub_end, int delete_sub,
                          int delete_ext_sub_array, int realloc_table){

    ASSERT_NULL_PTR(uc,"delete_extend_annots()")

    int size_line = size_substring + uc->size_annot;
    int tot = uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

    if (uc->nb_extended_annot != 0){

        uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot)]);
        int pos = 0;
        int old_pos = pos;
        int size_delta_list = uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;
        int i = 0;
        int it_pos = pos_sub_start;
        int tmp_pos_sub_end = pos_sub_end;
        int shift;
        int deleted = 0;

        while(1){
            pos += ((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1]));

            if (pos > tmp_pos_sub_end){

                if (delete_ext_sub_array == 1){
                    pos = ((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1])) - (pos_sub_end - pos_sub_start + 1 - deleted);
                    extend_annot[i] = pos >> SIZE_BITS_UINT_8T;
                    extend_annot[i+1] = pos & 0xff;
                }

                break;
            }

            if (it_pos <= pos){

                shift = pos - it_pos + 1;

                if (i < size_delta_list-SIZE_BYTE_EXT_ANNOT){

                    memmove(&(extend_annot[i]), &(extend_annot[i+SIZE_BYTE_EXT_ANNOT]), (size_delta_list-i-SIZE_BYTE_EXT_ANNOT + tot)*sizeof(uint8_t));
                    pos += (((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1])));
                    if (delete_ext_sub_array == 1) pos -= shift;
                    extend_annot[i] = (pos - old_pos) >> SIZE_BITS_UINT_8T;
                    extend_annot[i+1] = (pos - old_pos) & 0xff;
                    it_pos = pos;
                    pos = old_pos;
                }

                deleted += shift;
                if (delete_ext_sub_array == 1) tmp_pos_sub_end -= shift;
                uc->nb_extended_annot--;
                size_delta_list -= SIZE_BYTE_EXT_ANNOT;

                if (i >= size_delta_list) break;
            }
            else{
                i += SIZE_BYTE_EXT_ANNOT;
                if (i >= size_delta_list) break;
                old_pos = pos;
            }
        }
    }

    if (delete_sub == 1){
        memmove(&(uc->suffixes[pos_sub_start * size_line]),
                &(uc->suffixes[(pos_sub_end+1) * size_line]),
                ((nb_substring - pos_sub_end  - 1) * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + tot) * sizeof(uint8_t));

        if (realloc_table == 1){
            uc->suffixes = realloc(uc->suffixes, ((nb_substring - (pos_sub_end - pos_sub_start + 1)) * size_line +
                                    uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + tot) * sizeof(uint8_t));
            ASSERT_NULL_PTR(uc->suffixes,"delete_extend_annots()")
        }
    }
    else if (realloc_table == 1){
        uc->suffixes = realloc(uc->suffixes, (nb_substring * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + tot) * sizeof(uint8_t));
        ASSERT_NULL_PTR(uc->suffixes,"delete_extend_annots()")
    }

    return;
}

uint8_t* get_extend_annot(UC* uc, int size_substring, int nb_substring, int pos_substring){

    ASSERT_NULL_PTR(uc,"get_extend_annot()")

    if (uc->nb_extended_annot != 0){

        uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot)]);

        int pos = (extend_annot[0] << SIZE_BITS_UINT_8T) | extend_annot[1];
        int i = SIZE_BYTE_EXT_ANNOT;

        while ((pos < pos_substring) && (i < uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT)){
            pos += (extend_annot[i] << SIZE_BITS_UINT_8T) | extend_annot[i+1];
            i += SIZE_BYTE_EXT_ANNOT;
        }

        if (pos == pos_substring) return &(extend_annot[i-1]);
    }

    return NULL;
}

uint8_t** get_extend_annots(UC* uc, int size_substring, int nb_substring, int pos_substring_begin, int pos_substring_end){

    ASSERT_NULL_PTR(uc,"get_extend_annots()")

    int i;
    int count_sub = pos_substring_end - pos_substring_begin + 1;
    uint8_t** ptr_extend_annot = NULL;

    if (uc->nb_extended_annot != 0){
        ptr_extend_annot = malloc(count_sub*sizeof(uint8_t*));
        ASSERT_NULL_PTR(ptr_extend_annot,"get_extend_annots()")

        for (i = 0; i < count_sub; i++) ptr_extend_annot[i] = NULL;

        uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot)]);
        int pos = 0;
        int it_pos = pos_substring_begin;
        int bool_pres_extend_annot = 0;

        for (i = 0; i < uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT; i += SIZE_BYTE_EXT_ANNOT){

            pos += ((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1]));

            if (pos > pos_substring_end) break;
            if (it_pos <= pos){
                it_pos = pos;
                ptr_extend_annot[it_pos-pos_substring_begin] = &(extend_annot[i+2]);
                bool_pres_extend_annot = 1;
            }
        }

        if (bool_pres_extend_annot == 1) return ptr_extend_annot;
        else free(ptr_extend_annot);
    }

    return NULL;
}

void recopy_back_annot_extend(UC* uc, int size_substring, int nb_substring){

    ASSERT_NULL_PTR(uc,"recopy_back_annot_extend()")

    if (nb_substring == 0){
        uc->size_annot++;
        return;
    }

    int i, j;
    int old_size_line = size_substring + uc->size_annot;
    int new_size_line = old_size_line + 1;
    int tot_size_line = nb_substring * old_size_line;
    int pos = old_size_line;

    uint8_t* new_suffixes;

    new_suffixes = calloc(nb_substring * new_size_line + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes), sizeof(uint8_t));
    ASSERT_NULL_PTR(new_suffixes,"recopy_back_annot_extend()")

    for (i = 0, j = 0; i < nb_substring * new_size_line; i += new_size_line, j += old_size_line)
        memcpy(&(new_suffixes[i]), &(uc->suffixes[j]), old_size_line*sizeof(uint8_t));

    if (uc->nb_extended_annot != 0){

        for (i = tot_size_line; i < tot_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT; i += SIZE_BYTE_EXT_ANNOT){

            pos += ((((uint16_t)uc->suffixes[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)uc->suffixes[i + 1])) * new_size_line;
            new_suffixes[pos] = uc->suffixes[i + SIZE_BYTE_EXT_ANNOT - 1];
        }

        uc->nb_extended_annot = 0;
    }

    memcpy(&(new_suffixes[nb_substring * new_size_line]), &(uc->suffixes[tot_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
           uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes) * sizeof(uint8_t));

    uc->size_annot++;

    free(uc->suffixes);
    uc->suffixes = new_suffixes;

    return;
}

void create_annot_extended(UC* uc, int size_substring, int nb_substring){

    if (nb_substring == 0) return;

    ASSERT_NULL_PTR(uc,"create_annot_extended()")
    ASSERT_NULL_PTR(uc->suffixes,"create_annot_extended()")

    if (uc->size_annot <= 1) return;
    if (uc->nb_extended_annot > 0) return;

    int i = 0;
    int old_pos = 0;
    int it_annot_extend = 0;
    int nb_possible_annot_extend = 0;
    int size_line = size_substring + uc->size_annot;

    int delta;
    int new_size_line;
    int tot_new_size_line;

    uint8_t* new_tab_suffixes;
    uint8_t* extend_annot;

    for (i = size_line - 1; i < nb_substring * size_line; i += size_line)
        nb_possible_annot_extend += uc->suffixes[i] != 0;

    if ((nb_possible_annot_extend * SIZE_BYTE_EXT_ANNOT) >= nb_substring) return;

    new_size_line = size_substring + uc->size_annot - 1;
    tot_new_size_line = nb_substring * new_size_line;

    new_tab_suffixes = calloc(tot_new_size_line + nb_possible_annot_extend * SIZE_BYTE_EXT_ANNOT
                              + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes), sizeof(uint8_t));
    ASSERT_NULL_PTR(new_tab_suffixes,"create_annot_extended()")

    extend_annot = &(new_tab_suffixes[tot_new_size_line]);

    for (i = 0; i < nb_substring; i++){

        memcpy(&(new_tab_suffixes[i*new_size_line]), &(uc->suffixes[i*size_line]), new_size_line*sizeof(uint8_t));

        if (uc->suffixes[(i+1) * size_line - 1] != 0){

            delta = i - old_pos;

            extend_annot[it_annot_extend] = delta >> SIZE_BITS_UINT_8T;
            extend_annot[it_annot_extend+1] = delta & 0xff;
            extend_annot[it_annot_extend+2] = uc->suffixes[(i+1) * size_line - 1];

            old_pos = i;
            it_annot_extend += SIZE_BYTE_EXT_ANNOT;
        }
    }

    memcpy(&(new_tab_suffixes[tot_new_size_line + it_annot_extend]),
           &(uc->suffixes[nb_substring * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
           uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes) * sizeof(uint8_t));

    uc->size_annot--;
    uc->nb_extended_annot = nb_possible_annot_extend;

    free(uc->suffixes);
    uc->suffixes = new_tab_suffixes;

    return;
}

uint8_t* realloc_annotation(UC* uc, int size_substring, int nb_substring, int new_size_annotation, int new_insertion, int pos_insert_extend){

    ASSERT_NULL_PTR(uc,"realloc_annotation()")

    if (nb_substring == 0){
        uc->size_annot = new_size_annotation;
        uc->nb_extended_annot = 0;
        return NULL;
    }

    ASSERT_NULL_PTR(uc->suffixes,"realloc_annotation()")

    int i = 0, j;
    int old_size_line = size_substring + uc->size_annot;
    int new_size_line = size_substring + new_size_annotation;
    int tot_size_line_cplx = uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

    int tot_size_line;

    uint8_t* new_tab_suffixes;
    uint8_t* to_return = NULL;

    if (new_size_annotation > uc->size_annot+1){

        if (new_insertion == 1) new_size_line--;

        new_tab_suffixes = calloc(((nb_substring + new_insertion) * new_size_line + new_insertion * SIZE_BYTE_EXT_ANNOT + tot_size_line_cplx), sizeof(uint8_t));
        ASSERT_NULL_PTR(new_tab_suffixes,"realloc_annotation()")

        tot_size_line = nb_substring * old_size_line;

        for (i = 0, j = 0; i < nb_substring * new_size_line; i += new_size_line, j += old_size_line)
            memcpy(&(new_tab_suffixes[i]), &(uc->suffixes[j]), old_size_line*sizeof(uint8_t));

        if (uc->nb_extended_annot != 0){

            int pos = old_size_line;

            for (i = tot_size_line; i < tot_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT; i += SIZE_BYTE_EXT_ANNOT){
                pos += ((((uint16_t)uc->suffixes[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)uc->suffixes[i + 1])) * new_size_line;
                new_tab_suffixes[pos] = uc->suffixes[i + SIZE_BYTE_EXT_ANNOT - 1];
            }
        }

        if (new_insertion == 1){

            memmove(&(new_tab_suffixes[(pos_insert_extend+1) * new_size_line]),
                    &(new_tab_suffixes[pos_insert_extend * new_size_line]),
                    (nb_substring - pos_insert_extend) * new_size_line * sizeof(uint8_t));

            nb_substring = (nb_substring+1) * new_size_line;
            new_tab_suffixes[nb_substring] = pos_insert_extend >> SIZE_BITS_UINT_8T;
            new_tab_suffixes[nb_substring + 1] = pos_insert_extend & 0xff;
            to_return = &(new_tab_suffixes[nb_substring + SIZE_BYTE_EXT_ANNOT - 1]);

            memcpy(&(to_return[1]), &(uc->suffixes[tot_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]), tot_size_line_cplx * sizeof(uint8_t));
        }
        else memcpy(&(new_tab_suffixes[nb_substring * new_size_line]),
                    &(uc->suffixes[tot_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
                    tot_size_line_cplx * sizeof(uint8_t));

        free(uc->suffixes);
        uc->suffixes = new_tab_suffixes;
        uc->size_annot = new_size_line - size_substring;
        uc->nb_extended_annot = new_insertion;
    }
    else if (new_size_annotation > uc->size_annot){

        if (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT > nb_substring){
            recopy_back_annot_extend(uc, size_substring, nb_substring);
            return realloc_annotation(uc, size_substring, nb_substring, new_size_annotation, new_insertion, pos_insert_extend);
        }

        uc->suffixes = realloc(uc->suffixes, ((nb_substring + new_insertion) * old_size_line + (uc->nb_extended_annot+1) * SIZE_BYTE_EXT_ANNOT
                                + tot_size_line_cplx) * sizeof(uint8_t));
        ASSERT_NULL_PTR(uc->suffixes,"realloc_annotation()")

        tot_size_line = nb_substring*old_size_line;

        if (new_insertion == 1){

            memmove(&(uc->suffixes[(pos_insert_extend+1) * old_size_line]),
                    &(uc->suffixes[pos_insert_extend * old_size_line]),
                    (((nb_substring - pos_insert_extend) * old_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + tot_size_line_cplx) * sizeof(uint8_t)));

            memset(&(uc->suffixes[pos_insert_extend * old_size_line]), 0, old_size_line * sizeof(uint8_t));

            tot_size_line += old_size_line;
        }

        if (uc->nb_extended_annot != 0){

            uint8_t* extended_annot = &(uc->suffixes[tot_size_line]);

            int pos = 0;
            int old_pos = 0;
            int delta;
            int sum;

            i = 0;

            pos += ((((uint16_t)extended_annot[0]) << SIZE_BITS_UINT_8T) | ((uint16_t)extended_annot[1]));

            while (pos < pos_insert_extend){
                i++;
                if (i == uc->nb_extended_annot) break;
                old_pos = pos;
                pos += ((((uint16_t)extended_annot[i*SIZE_BYTE_EXT_ANNOT]) << SIZE_BITS_UINT_8T) | ((uint16_t)extended_annot[i*SIZE_BYTE_EXT_ANNOT+1]));
            }

            sum = i*SIZE_BYTE_EXT_ANNOT;

            if (i==0){

                memmove(&(extended_annot[SIZE_BYTE_EXT_ANNOT]), extended_annot, (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + tot_size_line_cplx) * sizeof(uint8_t));

                extended_annot[0] = pos_insert_extend >> SIZE_BITS_UINT_8T;
                extended_annot[1] = pos_insert_extend & 0xff;

                delta = pos - pos_insert_extend + new_insertion;
            }
            else if (i != uc->nb_extended_annot){

                memmove(&(extended_annot[sum + SIZE_BYTE_EXT_ANNOT]), &(extended_annot[sum]),
                        ((uc->nb_extended_annot - i) * SIZE_BYTE_EXT_ANNOT + tot_size_line_cplx) * sizeof(uint8_t));

                delta = pos_insert_extend - old_pos;

                extended_annot[sum] = delta >> SIZE_BITS_UINT_8T;
                extended_annot[sum+1] = delta & 0xff;

                delta = pos - pos_insert_extend + new_insertion;
            }
            else delta = pos_insert_extend - pos;

            if (i == uc->nb_extended_annot){

                extended_annot[sum] = delta >> SIZE_BITS_UINT_8T;
                extended_annot[sum+1] = delta & 0xff;
            }
            else{

                extended_annot[sum+SIZE_BYTE_EXT_ANNOT] = delta >> SIZE_BITS_UINT_8T;
                extended_annot[sum+SIZE_BYTE_EXT_ANNOT+1] = delta & 0xff;
            }

            to_return = &(extended_annot[sum+SIZE_BYTE_EXT_ANNOT-1]);
        }
        else {

            uc->suffixes[tot_size_line] = pos_insert_extend >> SIZE_BITS_UINT_8T;
            uc->suffixes[tot_size_line+1] = pos_insert_extend & 0xff;
            to_return = &(uc->suffixes[tot_size_line+2]);
        }

        uc->nb_extended_annot++;
    }
    else if (new_size_annotation == uc->size_annot){

        if (new_insertion == 1){

            uc->suffixes = realloc(uc->suffixes, ((nb_substring + 1) * old_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                                    + tot_size_line_cplx) * sizeof(uint8_t));
            ASSERT_NULL_PTR(uc->suffixes,"realloc_annotation()")

            memmove(&(uc->suffixes[(pos_insert_extend+1) * old_size_line]),
                    &(uc->suffixes[pos_insert_extend * old_size_line]),
                    (((nb_substring - pos_insert_extend) * old_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + tot_size_line_cplx) * sizeof(uint8_t)));

            memset(&(uc->suffixes[pos_insert_extend * old_size_line]), 0, old_size_line * sizeof(uint8_t));
        }
    }
    else{

        new_tab_suffixes = calloc(((nb_substring + new_insertion) * new_size_line + tot_size_line_cplx), sizeof(uint8_t));
        ASSERT_NULL_PTR(new_tab_suffixes,"realloc_annotation()")

        if (new_insertion == 1){

            for (i = 0, j = 0; j < pos_insert_extend * old_size_line; i += new_size_line, j += old_size_line)
                memcpy(&(new_tab_suffixes[i]), &(uc->suffixes[j]), new_size_line * sizeof(uint8_t));

            i += new_size_line;

            for (; j < nb_substring * old_size_line; i += new_size_line,j += old_size_line)
                memcpy(&(new_tab_suffixes[i]), &(uc->suffixes[j]), new_size_line * sizeof(uint8_t));


        }
        else{
            for (i=0; i<nb_substring; i++)
                memcpy(&(new_tab_suffixes[i*new_size_line]), &(uc->suffixes[i*old_size_line]), new_size_line*sizeof(uint8_t));
        }

        /*if (new_insertion == 1){
            memmove(&(new_tab_suffixes[(pos_insert_extend+1) * new_size_line]),
                    &(new_tab_suffixes[pos_insert_extend * new_size_line]),
                    ((nb_substring - pos_insert_extend) * new_size_line)* sizeof(uint8_t));

            memset(&(new_tab_suffixes[pos_insert_extend * new_size_line]), 0, new_size_line*sizeof(uint8_t));
        }*/

        memcpy(&(new_tab_suffixes[(nb_substring + new_insertion) * new_size_line]),
               &(uc->suffixes[nb_substring * old_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
               tot_size_line_cplx * sizeof(uint8_t));

        free(uc->suffixes);
        uc->suffixes = new_tab_suffixes;
        uc->size_annot = new_size_annotation;
        uc->nb_extended_annot = 0;
    }

    if (to_return != NULL) memset(to_return, 0, sizeof(uint8_t));

    return to_return;
}

void shift_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_insert){

    ASSERT_NULL_PTR(uc,"shift_annot_cplx_nodes()")

    if (uc->nb_cplx_nodes != 0){

        uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]);
        int pos = 0;
        int i = 0;
        int size_line = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;

        while (1){
            pos += ((((uint16_t)extend_annot[i * size_line]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i * size_line + 1]));
            if (pos < pos_insert){
                i++;
                if (i == uc->nb_cplx_nodes) break;
            }
            else break;
        }

        if (i < uc->nb_cplx_nodes){
            i *= size_line;
            uint16_t delta = ((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1]));
            delta++;
            extend_annot[i] = delta >> SIZE_BITS_UINT_8T;
            extend_annot[i+1] = delta & 0xff;
        }
    }
}

void insert_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_insert, uint8_t* annot, int size_annot, int shift_or_not){

    ASSERT_NULL_PTR(uc,"insert_annot_cplx_nodes()")
    ASSERT_NULL_PTR(annot,"insert_annot_cplx_nodes()")

    uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]);

    if (uc->nb_cplx_nodes != 0){

        int pos = 0;
        int i = 0;
        int old_pos;
        int delta;
        int sum;
        int size_line = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;

        while (1){
            old_pos = pos;
            pos += ((((uint16_t)extend_annot[i * size_line]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i * size_line + 1]));
            if (pos < pos_insert){
                i++;
                if (i == uc->nb_cplx_nodes) break;
            }
            else break;
        }

        delta = pos_insert - pos;

        if (i == 0){
            memmove(&(extend_annot[size_line]), extend_annot, uc->nb_cplx_nodes * size_line * sizeof(uint8_t));

            extend_annot[0] = pos_insert >> SIZE_BITS_UINT_8T;
            extend_annot[1] = pos_insert & 0xff;
            memcpy(&(extend_annot[SIZE_BYTE_CPLX_N]), annot, size_annot * sizeof(uint8_t));

            delta = pos - pos_insert;
        }
        else if (i == uc->nb_cplx_nodes) memcpy(&(extend_annot[i * size_line + SIZE_BYTE_CPLX_N]), annot, size_annot * sizeof(uint8_t));
        else {
            sum = i * size_line;
            memmove(&(extend_annot[sum + size_line]), &(extend_annot[sum]), (uc->nb_cplx_nodes - i) * size_line * sizeof(uint8_t));

            delta = pos_insert - old_pos;

            extend_annot[sum] = delta >> SIZE_BITS_UINT_8T;
            extend_annot[sum+1] = delta & 0xff;
            memcpy(&(extend_annot[sum + SIZE_BYTE_CPLX_N]), annot, size_annot * sizeof(uint8_t));

            delta = pos - pos_insert;
        }

        if (i != uc->nb_cplx_nodes){
            if (shift_or_not == 1) delta++;
        }
        else i--;

        extend_annot[(i+1) * size_line] = delta >> SIZE_BITS_UINT_8T;
        extend_annot[(i+1) * size_line + 1] = delta & 0xff;

        uc->nb_cplx_nodes++;
    }
    else{
        extend_annot[0] = pos_insert >> SIZE_BITS_UINT_8T;
        extend_annot[1] = pos_insert & 0xff;
        memcpy(&(extend_annot[SIZE_BYTE_CPLX_N]), annot, size_annot * sizeof(uint8_t));
        uc->nb_cplx_nodes = 1;
    }

    return;
}

void delete_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_sub_start, int pos_sub_end, int delete_sub, int delete_ext_sub_array, int realloc_table){

    ASSERT_NULL_PTR(uc,"delete_annot_cplx_nodes()")

    int size_line = size_substring + uc->size_annot;
    int size_line_cplx_n = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;

    if (uc->nb_cplx_nodes != 0){

        uint8_t* extend_annot = &(uc->suffixes[nb_substring * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]);

        int pos = 0;
        int old_pos = pos;
        int size_delta_list = uc->nb_cplx_nodes * size_line_cplx_n;
        int i = 0;
        int it_pos = pos_sub_start;
        int tmp_pos_sub_end = pos_sub_end;
        int shift;
        int deleted = 0;

        while(1){
            pos += ((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1]));

            if (pos > tmp_pos_sub_end){
                if (delete_ext_sub_array == 1){
                    pos = ((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1])) - (pos_sub_end - pos_sub_start + 1 - deleted);
                    extend_annot[i] = pos >> SIZE_BITS_UINT_8T;
                    extend_annot[i+1] = pos & 0xff;
                }

                break;
            }

            if (it_pos <= pos){

                shift = pos - it_pos + 1;

                if (i < size_delta_list - size_line_cplx_n){
                    memmove(&(extend_annot[i]), &(extend_annot[i + size_line_cplx_n]), (size_delta_list - i - size_line_cplx_n) * sizeof(uint8_t));
                    pos += (((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1])));
                    if (delete_ext_sub_array == 1) pos -= shift;
                    extend_annot[i] = (pos - old_pos) >> SIZE_BITS_UINT_8T;
                    extend_annot[i+1] = (pos - old_pos) & 0xff;
                    it_pos = pos;
                    pos = old_pos;
                }

                deleted += shift;
                if (delete_ext_sub_array == 1) tmp_pos_sub_end -= shift;
                uc->nb_cplx_nodes--;
                size_delta_list -= size_line_cplx_n;

                if (i >= size_delta_list) break;
            }
            else{
                i += size_line_cplx_n;
                if (i >= size_delta_list) break;
                old_pos = pos;
            }
        }
    }

    if (delete_sub == 1){
        memmove(&(uc->suffixes[pos_sub_start * size_line]),
                &(uc->suffixes[(pos_sub_end+1) * size_line]),
                ((nb_substring - pos_sub_end  -1) * size_line +
                 uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + uc->nb_cplx_nodes * size_line_cplx_n) * sizeof(uint8_t));

        if (realloc_table == 1){
            uc->suffixes = realloc(uc->suffixes,
                                   ((nb_substring - (pos_sub_end - pos_sub_start + 1)) * size_line +
                                    uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + uc->nb_cplx_nodes * size_line_cplx_n) * sizeof(uint8_t));

            ASSERT_NULL_PTR(uc->suffixes,"delete_annot_cplx_nodes()")
        }
    }
    else if (realloc_table == 1){
        uc->suffixes = realloc(uc->suffixes, (nb_substring * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT +
                                              uc->nb_cplx_nodes * size_line_cplx_n) * sizeof(uint8_t));

        ASSERT_NULL_PTR(uc->suffixes,"delete_annot_cplx_nodes()")
    }

    return;
}

uint8_t* get_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_substring){

    ASSERT_NULL_PTR(uc,"get_annot_cplx_nodes()")

    if (uc->nb_cplx_nodes != 0){
        uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]);
        int size_line = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;
        int pos = 0;
        int i = 0;

        for (i=0; i<uc->nb_cplx_nodes * size_line; i += size_line){
            pos += ((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1]));
            if (pos == pos_substring) return &(extend_annot[i+2]);
            else if (pos > pos_substring) return NULL;
        }
    }

    return NULL;
}

uint8_t** get_annots_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_substring_begin, int pos_substring_end){

    ASSERT_NULL_PTR(uc,"get_annots_cplx_nodes()")

    int i;
    int size_line = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;
    int count_sub = pos_substring_end - pos_substring_begin + 1;
    uint8_t** ptr_extend_annot = NULL;

    if (uc->nb_cplx_nodes != 0){
        ptr_extend_annot = malloc(count_sub*sizeof(uint8_t*));
        ASSERT_NULL_PTR(ptr_extend_annot,"get_annots_cplx_nodes()")

        for (i=0; i<count_sub; i++) ptr_extend_annot[i] = NULL;

        uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]);
        int pos = 0;
        int it_pos = pos_substring_begin;
        int bool_pres_extend_annot = 0;

        for (i=0; i < uc->nb_cplx_nodes * size_line; i += size_line){
            pos += ((((uint16_t)extend_annot[i]) << SIZE_BITS_UINT_8T) | ((uint16_t)extend_annot[i+1]));

            if (pos > pos_substring_end) break;
            if (it_pos <= pos){
                it_pos = pos;
                ptr_extend_annot[it_pos-pos_substring_begin] = &(extend_annot[i+2]);
                bool_pres_extend_annot = 1;
            }
        }

        if (bool_pres_extend_annot == 1) return ptr_extend_annot;
        else free(ptr_extend_annot);
    }

    return NULL;
}

void recopy_back_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring){

    ASSERT_NULL_PTR(uc,"recopy_back_annot_cplx_nodes()")

    if (nb_substring == 0){
        uc->size_annot += uc->size_annot_cplx_nodes;
        uc->nb_cplx_nodes = 0;
        return;
    }

    int i = 0;
    int old_size_line = size_substring + uc->size_annot;
    int new_size_line = old_size_line + uc->size_annot_cplx_nodes;
    int tot_size_line = nb_substring * old_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;

    uint8_t* new_suffixes = calloc(nb_substring * new_size_line,sizeof(uint8_t));
    ASSERT_NULL_PTR(new_suffixes,"recopy_back_annot_cplx_nodes()")

    for (i=0; i<nb_substring; i++)
        memcpy(&(new_suffixes[i*new_size_line]), &(uc->suffixes[i*old_size_line]), old_size_line*sizeof(uint8_t));

    if (uc->nb_cplx_nodes != 0){

        int pos = 0;
        int size_line_cplx_nodes = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;

        for (i=0; i < uc->nb_cplx_nodes * size_line_cplx_nodes; i += size_line_cplx_nodes){
            pos += ((((uint16_t)uc->suffixes[tot_size_line + i]) << SIZE_BITS_UINT_8T) | ((uint16_t)uc->suffixes[tot_size_line + i + 1]));

            memcpy(&(new_suffixes[pos * new_size_line + old_size_line]),
                   &(uc->suffixes[tot_size_line + i + SIZE_BYTE_CPLX_N]),
                   uc->size_annot_cplx_nodes * sizeof(uint8_t*));
        }

        uc->nb_cplx_nodes = 0;
        uc->size_annot_cplx_nodes = 0;
    }

    uc->size_annot += uc->size_annot_cplx_nodes;

    free(uc->suffixes);
    uc->suffixes = new_suffixes;
}

void increase_size_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int new_size, int to_realloc){

    ASSERT_NULL_PTR(uc,"increase_size_annot_cplx_nodes()")

    if ((nb_substring == 0) || (uc->nb_cplx_nodes == 0)){
        uc->size_annot_cplx_nodes = new_size;
        return;
    }

    if (new_size >= uc->size_annot_cplx_nodes) return;

    int i = 0;
    int size_line_cplx_n = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;
    int new_size_line_cplx_n = SIZE_BYTE_CPLX_N + new_size;
    int tot_size_line = nb_substring * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;

    uint8_t* cplx_n;

    if (to_realloc == 1){
        uc->suffixes = realloc(uc->suffixes, (tot_size_line + uc->nb_cplx_nodes * size_line_cplx_n) * sizeof(uint8_t));
        ASSERT_NULL_PTR(uc->suffixes, "increase_size_annot_cplx_nodes()")
    }

    cplx_n = &(uc->suffixes[tot_size_line]);

    memset(&(cplx_n[uc->nb_cplx_nodes * size_line_cplx_n]), 0,
           uc->nb_cplx_nodes * (new_size - uc->size_annot_cplx_nodes) * sizeof(uint8_t));

    for (i = uc->nb_cplx_nodes-1; i >= 0; i--){
        memcpy(&(cplx_n[i*new_size_line_cplx_n]), &(cplx_n[i*size_line_cplx_n]), size_line_cplx_n * sizeof(uint8_t));
        memset(&(cplx_n[i*size_line_cplx_n]), 0, size_line_cplx_n * sizeof(uint8_t));
    }

    uc->size_annot_cplx_nodes = new_size;

    return;
}

void decrease_size_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int new_size){

    ASSERT_NULL_PTR(uc,"decrease_size_annot_cplx_nodes()")

    if ((nb_substring == 0) || (uc->nb_cplx_nodes == 0)){
        uc->size_annot_cplx_nodes = new_size;
        return;
    }

    if (new_size >= uc->size_annot_cplx_nodes) return;

    int i = 0;
    int size_line_cplx_n = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;
    int new_size_line_cplx_n = SIZE_BYTE_CPLX_N + new_size;
    int tot_size_line = nb_substring * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;

    uint8_t* cplx_n = &(uc->suffixes[tot_size_line]);

    for (i = 0; i < uc->nb_cplx_nodes; i++){
        memmove(&(cplx_n[i*new_size_line_cplx_n]), &(cplx_n[i*size_line_cplx_n]), size_line_cplx_n * sizeof(uint8_t));
        memset(&(cplx_n[i*size_line_cplx_n]), 0, size_line_cplx_n * sizeof(uint8_t));
    }

    uc->suffixes = realloc(uc->suffixes, (tot_size_line + uc->nb_cplx_nodes * size_line_cplx_n) * sizeof(uint8_t));
    ASSERT_NULL_PTR(uc->suffixes, "decrease_size_annot_cplx_nodes()")

    uc->size_annot_cplx_nodes = new_size;

    return;
}

void create_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring){
    if (nb_substring == 0) return;
    ASSERT_NULL_PTR(uc,"create_annot_cplx_nodes()")
    ASSERT_NULL_PTR(uc->suffixes,"create_annot_cplx_nodes()")

    int delta;
    int size_line_cplx_n;
    int tot_new_size_line;

    int i = 0;
    int old_pos = 0;
    int max_size_annot = 0;
    int it_annot_extend = 0;
    int nb_possible_cplx_n = 0;
    int size_line = size_substring + uc->size_annot;

    uint8_t* new_tab_suffixes;
    uint8_t* extend_annot;
    uint8_t possible_cplx_n[nb_substring];

    uint8_t** ext_annot = get_extend_annots(uc, size_substring, nb_substring, 0, nb_substring-1);
    UC_SIZE_ANNOT_T *size_annot = min_size_per_sub(uc->suffixes, nb_substring, size_substring, uc->size_annot);
    ASSERT_NULL_PTR(size_annot,"create_annot_cplx_nodes()")

    for (i=0; i<nb_substring; i++){
        if ((ext_annot != NULL) && (ext_annot[i] != NULL) && (ext_annot[i][0] != 0)){
            nb_possible_cplx_n++;
            possible_cplx_n[i] = 1;
            max_size_annot = uc->size_annot + 1;
        }
        else if (size_annot[i] != 0){
            nb_possible_cplx_n++;
            possible_cplx_n[i] = 1;
            max_size_annot = MAX(max_size_annot, size_annot[i]);
        }
        else possible_cplx_n[i] = 0;
    }

    if ((nb_possible_cplx_n * (max_size_annot + SIZE_BYTE_CPLX_N)) >
        (nb_substring * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT)) return;

    tot_new_size_line = nb_substring * size_substring;
    size_line_cplx_n = max_size_annot + SIZE_BYTE_CPLX_N;

    new_tab_suffixes = calloc(tot_new_size_line + nb_possible_cplx_n * size_line_cplx_n, sizeof(uint8_t));
    ASSERT_NULL_PTR(new_tab_suffixes,"create_annot_cplx_nodes()")

    extend_annot = &(new_tab_suffixes[tot_new_size_line]);

    for (i=0; i<nb_substring; i++){
        memcpy(&(new_tab_suffixes[i*size_substring]), &(uc->suffixes[i*size_line]), size_substring*sizeof(uint8_t));

        if (possible_cplx_n[i] == 1){

            delta = i - old_pos;

            extend_annot[it_annot_extend] = delta >> SIZE_BITS_UINT_8T;
            extend_annot[it_annot_extend+1] = delta & 0xff;
            memcpy(&(extend_annot[it_annot_extend+2]), &(uc->suffixes[i * size_line + size_substring]), size_annot[i] * sizeof(uint8_t));

            if (size_annot[i] == uc->size_annot){
                if ((ext_annot != NULL) && (ext_annot[i] != NULL) && (ext_annot[i][0] != 0))
                    extend_annot[it_annot_extend + uc->size_annot] = ext_annot[i][0];
            }

            old_pos = i;
            it_annot_extend += size_line_cplx_n;
        }
    }

    uc->size_annot = 0;
    uc->size_annot_cplx_nodes = max_size_annot;
    uc->nb_cplx_nodes = nb_possible_cplx_n;

    free(uc->suffixes);
    uc->suffixes = new_tab_suffixes;

    free(size_annot);
    if (ext_annot != NULL) free(ext_annot);

    return;
}

/*void create_annot_cplx_nodes_marked(UC* uc, int size_substring, int nb_substring){

    if (nb_substring == 0) return;
    ASSERT_NULL_PTR(uc,"create_annot_cplx_nodes()")
    ASSERT_NULL_PTR(uc->suffixes,"create_annot_cplx_nodes()")

    int delta;
    int size_line_cplx_n;
    int tot_new_size_line;

    int i = 0;
    int old_pos = 0;
    int max_size_annot = 0;
    int it_annot_extend = 0;
    int nb_possible_cplx_n = 0;
    int size_line = size_substring + uc->size_annot;

    int tot = nb_substring * size_line
                + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT
                + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes);

    uint8_t possible_cplx_n[nb_substring];

    uint8_t* new_tab_suffixes;
    uint8_t* extend_annot;

    uint8_t** ext_annot = get_extend_annots(uc, size_substring, nb_substring, 0, nb_substring-1);

    UC_SIZE_ANNOT_T *size_annot = min_size_per_sub(uc->suffixes, nb_substring, size_substring, uc->size_annot);
    ASSERT_NULL_PTR(size_annot,"create_annot_cplx_nodes()")

    memset(possible_cplx_n, 0, nb_substring * sizeof(uint8_t));

    for (i=0; i<nb_substring; i++){

        if (get_mark_UC_4states_bis(uc, i, tot) != 0){

            if ((ext_annot != NULL) && (ext_annot[i] != NULL) && (ext_annot[i][0] != 0)){
                nb_possible_cplx_n++;
                possible_cplx_n[i] = 1;
                max_size_annot = uc->size_annot + 1;
            }
            else if (size_annot[i] != 0){
                nb_possible_cplx_n++;
                possible_cplx_n[i] = 1;
                max_size_annot = MAX(max_size_annot, size_annot[i]);
            }
        }
    }

    if ((nb_possible_cplx_n * (max_size_annot + SIZE_BYTE_CPLX_N)) >
        (nb_substring * uc->size_annot + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT)) goto CREATE_ANNOT_CPLX_NODES_MARKED;

    tot_new_size_line = nb_substring * size_substring;
    size_line_cplx_n = max_size_annot + SIZE_BYTE_CPLX_N;

    new_tab_suffixes = calloc(tot_new_size_line + nb_possible_cplx_n * size_line_cplx_n, sizeof(uint8_t));
    ASSERT_NULL_PTR(new_tab_suffixes,"create_annot_cplx_nodes()")

    extend_annot = &(new_tab_suffixes[tot_new_size_line]);

    for (i=0; i<nb_substring; i++){

        memcpy(&(new_tab_suffixes[i*size_substring]), &(uc->suffixes[i*size_line]), size_substring*sizeof(uint8_t));

        if (possible_cplx_n[i] == 1){

            delta = i - old_pos;

            extend_annot[it_annot_extend] = delta >> SIZE_BITS_UINT_8T;
            extend_annot[it_annot_extend+1] = delta & 0xff;
            memcpy(&(extend_annot[it_annot_extend+2]), &(uc->suffixes[i * size_line + size_substring]), size_annot[i] * sizeof(uint8_t));

            if (size_annot[i] == uc->size_annot){
                if ((ext_annot != NULL) && (ext_annot[i] != NULL) && (ext_annot[i][0] != 0))
                    extend_annot[it_annot_extend + uc->size_annot] = ext_annot[i][0];
            }

            old_pos = i;
            it_annot_extend += size_line_cplx_n;
        }
    }

    uc->size_annot = 0;
    uc->nb_extended_annot = 0;
    uc->nb_cplx_nodes = nb_possible_cplx_n;
    uc->size_annot_cplx_nodes = max_size_annot;

    free(uc->suffixes);
    uc->suffixes = new_tab_suffixes;

    CREATE_ANNOT_CPLX_NODES_MARKED: free(size_annot);
    if (ext_annot != NULL) free(ext_annot);

    return;
}*/
