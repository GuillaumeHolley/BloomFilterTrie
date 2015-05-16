#include "annotation_special_nodes.h"

extern uint8_t* min_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot);

void shift_annot_cplx_nodes(UC* uc, int size_substring, int nb_substring, int pos_insert){

    ASSERT_NULL_PTR(uc,"shift_annot_cplx_nodes()")

    if (uc->nb_cplx_nodes != 0){

        uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot) + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]);
        int pos = 0;
        int i = 0;
        int size_line = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;

        while (1){
            pos += ((((uint16_t)extend_annot[i * size_line]) << SIZE_CELL) | ((uint16_t)extend_annot[i * size_line + 1]));
            if (pos < pos_insert){
                i++;
                if (i == uc->nb_cplx_nodes) break;
            }
            else break;
        }

        if (i < uc->nb_cplx_nodes){
            i *= size_line;
            uint16_t delta = ((((uint16_t)extend_annot[i]) << SIZE_CELL) | ((uint16_t)extend_annot[i+1]));
            delta++;
            extend_annot[i] = delta >> SIZE_CELL;
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
            pos += ((((uint16_t)extend_annot[i * size_line]) << SIZE_CELL) | ((uint16_t)extend_annot[i * size_line + 1]));
            if (pos < pos_insert){
                i++;
                if (i == uc->nb_cplx_nodes) break;
            }
            else break;
        }

        delta = pos_insert - pos;

        if (i == 0){
            memmove(&(extend_annot[size_line]), extend_annot, uc->nb_cplx_nodes * size_line * sizeof(uint8_t));

            extend_annot[0] = pos_insert >> SIZE_CELL;
            extend_annot[1] = pos_insert & 0xff;
            memcpy(&(extend_annot[SIZE_BYTE_CPLX_N]), annot, size_annot * sizeof(uint8_t));

            delta = pos - pos_insert;
        }
        else if (i == uc->nb_cplx_nodes) memcpy(&(extend_annot[i * size_line + SIZE_BYTE_CPLX_N]), annot, size_annot * sizeof(uint8_t));
        else {
            sum = i * size_line;
            memmove(&(extend_annot[sum + size_line]), &(extend_annot[sum]), (uc->nb_cplx_nodes - i) * size_line * sizeof(uint8_t));

            delta = pos_insert - old_pos;

            extend_annot[sum] = delta >> SIZE_CELL;
            extend_annot[sum+1] = delta & 0xff;
            memcpy(&(extend_annot[sum + SIZE_BYTE_CPLX_N]), annot, size_annot * sizeof(uint8_t));

            delta = pos - pos_insert;
        }

        if (i != uc->nb_cplx_nodes){
            if (shift_or_not == 1) delta++;
        }
        else i--;

        extend_annot[(i+1) * size_line] = delta >> SIZE_CELL;
        extend_annot[(i+1) * size_line + 1] = delta & 0xff;

        uc->nb_cplx_nodes++;
    }
    else{
        extend_annot[0] = pos_insert >> SIZE_CELL;
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
            pos += ((((uint16_t)extend_annot[i]) << SIZE_CELL) | ((uint16_t)extend_annot[i+1]));

            if (pos > tmp_pos_sub_end){
                if (delete_ext_sub_array == 1){
                    pos = ((((uint16_t)extend_annot[i]) << SIZE_CELL) | ((uint16_t)extend_annot[i+1])) - (pos_sub_end - pos_sub_start + 1 - deleted);
                    extend_annot[i] = pos >> SIZE_CELL;
                    extend_annot[i+1] = pos & 0xff;
                }

                break;
            }

            if (it_pos <= pos){

                shift = pos - it_pos + 1;

                if (i < size_delta_list - size_line_cplx_n){
                    memmove(&(extend_annot[i]), &(extend_annot[i + size_line_cplx_n]), (size_delta_list - i - size_line_cplx_n) * sizeof(uint8_t));
                    pos += (((((uint16_t)extend_annot[i]) << SIZE_CELL) | ((uint16_t)extend_annot[i+1])));
                    if (delete_ext_sub_array == 1) pos -= shift;
                    extend_annot[i] = (pos - old_pos) >> SIZE_CELL;
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
            pos += ((((uint16_t)extend_annot[i]) << SIZE_CELL) | ((uint16_t)extend_annot[i+1]));
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
            pos += ((((uint16_t)extend_annot[i]) << SIZE_CELL) | ((uint16_t)extend_annot[i+1]));

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

    int old_size_line = size_substring + uc->size_annot;
    int new_size_line = old_size_line + uc->size_annot_cplx_nodes;
    int tot_size_line = nb_substring * old_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT;
    uint8_t* new_suffixes;
    int i = 0;

    new_suffixes = calloc(nb_substring * new_size_line,sizeof(uint8_t));
    ASSERT_NULL_PTR(new_suffixes,"recopy_back_annot_cplx_nodes()")

    for (i=0; i<nb_substring; i++)
        memcpy(&(new_suffixes[i*new_size_line]), &(uc->suffixes[i*old_size_line]), old_size_line*sizeof(uint8_t));

    if (uc->nb_cplx_nodes != 0){

        int pos = 0;
        int size_line_cplx_nodes = SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes;

        for (i=0; i < uc->nb_cplx_nodes * size_line_cplx_nodes; i += size_line_cplx_nodes){
            pos += ((((uint16_t)uc->suffixes[tot_size_line + i]) << SIZE_CELL) | ((uint16_t)uc->suffixes[tot_size_line + i + 1]));

            memcpy(&(new_suffixes[pos * new_size_line + old_size_line]),
                   &(uc->suffixes[tot_size_line + i + SIZE_BYTE_CPLX_N]),
                   uc->size_annot_cplx_nodes * sizeof(uint8_t*));
        }

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
    uint8_t* size_annot = min_size_per_sub(uc->suffixes, nb_substring, size_substring, uc->size_annot);
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

            extend_annot[it_annot_extend] = delta >> SIZE_CELL;
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
