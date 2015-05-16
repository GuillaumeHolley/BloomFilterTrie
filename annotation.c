#include "annotation.h"

extern int modify_annot_bis(uint8_t** current_annot, uint8_t* annot_sup, int* it_annot, int* size_current_annot, uint16_t id_genome, uint8_t flag, uint8_t flag_ext);
extern uint8_t* extract_from_annotation_array_elem(annotation_array_elem* annot_sorted, uint32_t position, int* size_annot);

void shift_extended_annot(UC* uc, int size_substring, int nb_substring, int pos_insert){

    ASSERT_NULL_PTR(uc,"shift_extended_annot()")

    if (uc->nb_extended_annot != 0){

        uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot)]);
        int pos = 0;
        int i = 0;

        while (1){

            pos += ((((uint16_t)extend_annot[i*SIZE_BYTE_EXT_ANNOT]) << SIZE_CELL) | ((uint16_t)extend_annot[i*SIZE_BYTE_EXT_ANNOT+1]));
            if (pos < pos_insert){
                i++;
                if (i == uc->nb_extended_annot) break;
            }
            else break;
        }

        if (i < uc->nb_extended_annot){
            i *= SIZE_BYTE_EXT_ANNOT;
            uint16_t delta = ((((uint16_t)extend_annot[i]) << SIZE_CELL) | ((uint16_t)extend_annot[i+1]));
            delta++;
            extend_annot[i] = delta >> SIZE_CELL;
            extend_annot[i+1] = delta & 0xff;
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
        int old_pos;
        int delta;
        int sum;

        while (1){
            old_pos = pos;
            pos += ((((uint16_t)extend_annot[i*SIZE_BYTE_EXT_ANNOT]) << SIZE_CELL) | ((uint16_t)extend_annot[i*SIZE_BYTE_EXT_ANNOT+1]));
            if (pos < pos_insert){
                i++;
                if (i == uc->nb_extended_annot) break;
            }
            else break;
        }

        delta = pos_insert - pos;

        if (i == 0){
            memmove(&(extend_annot[SIZE_BYTE_EXT_ANNOT]), &(extend_annot[0]), (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + tot) * sizeof(uint8_t));

            extend_annot[0] = pos_insert >> SIZE_CELL;
            extend_annot[1] = pos_insert & 0xff;
            extend_annot[2] = annot;

            delta = pos - pos_insert;
        }
        else if (i == uc->nb_extended_annot) extend_annot[(i+1)*SIZE_BYTE_EXT_ANNOT-1] = annot;
        else {
            sum = i*SIZE_BYTE_EXT_ANNOT;
            memmove(&(extend_annot[sum+SIZE_BYTE_EXT_ANNOT]), &(extend_annot[sum]), ((uc->nb_extended_annot - i) * SIZE_BYTE_EXT_ANNOT + tot) * sizeof(uint8_t));

            delta = pos_insert - old_pos;

            extend_annot[sum] = delta >> SIZE_CELL;
            extend_annot[sum+1] = delta & 0xff;
            extend_annot[sum+2] = annot;

            delta = pos - pos_insert;
        }

        if (i != uc->nb_extended_annot){
            if (shift_or_not == 1) delta++;
        }
        else i--;

        extend_annot[(i+1)*SIZE_BYTE_EXT_ANNOT] = delta >> SIZE_CELL;
        extend_annot[(i+1)*SIZE_BYTE_EXT_ANNOT+1] = delta & 0xff;

        uc->nb_extended_annot++;
    }
    else{
        extend_annot[0] = pos_insert >> SIZE_CELL;
        extend_annot[1] = pos_insert & 0xff;
        extend_annot[2] = annot;
        uc->nb_extended_annot = 1;
    }

    return;
}

void delete_extend_annots(UC* uc, int size_substring, int nb_substring, int pos_sub_start, int pos_sub_end, int delete_sub, int delete_ext_sub_array, int realloc_table){

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

                if (i < size_delta_list-SIZE_BYTE_EXT_ANNOT){

                    memmove(&(extend_annot[i]), &(extend_annot[i+SIZE_BYTE_EXT_ANNOT]), (size_delta_list-i-SIZE_BYTE_EXT_ANNOT + tot)*sizeof(uint8_t));
                    pos += (((((uint16_t)extend_annot[i]) << SIZE_CELL) | ((uint16_t)extend_annot[i+1])));
                    if (delete_ext_sub_array == 1) pos -= shift;
                    extend_annot[i] = (pos - old_pos) >> SIZE_CELL;
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
        int pos = 0;
        int i = 0;

        for (i=0; i<uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT; i+=SIZE_BYTE_EXT_ANNOT){
            pos += ((((uint16_t)extend_annot[i]) << SIZE_CELL) | ((uint16_t)extend_annot[i+1]));
            if (pos == pos_substring) return &(extend_annot[i+2]);
            else if (pos > pos_substring) return NULL;
        }
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

        for (i=0; i<count_sub; i++) ptr_extend_annot[i] = NULL;

        uint8_t* extend_annot = &(uc->suffixes[nb_substring * (size_substring + uc->size_annot)]);
        int pos = 0;
        int it_pos = pos_substring_begin;
        int bool_pres_extend_annot = 0;

        for (i=0; i < uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT; i+=SIZE_BYTE_EXT_ANNOT){
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

void recopy_back_annot_extend(UC* uc, int size_substring, int nb_substring){

    ASSERT_NULL_PTR(uc,"recopy_back_annot_extend()")

    if (nb_substring == 0){
        uc->size_annot++;
        return;
    }

    int old_size_line = size_substring + uc->size_annot;
    int new_size_line = old_size_line + 1;
    int tot_size_line = nb_substring * old_size_line;
    uint8_t* new_suffixes;
    int i = 0;
    int pos = 0;

    new_suffixes = calloc(nb_substring * new_size_line + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes), sizeof(uint8_t));
    ASSERT_NULL_PTR(new_suffixes,"recopy_back_annot_extend()")

    for (i=0; i<nb_substring; i++) memcpy(&(new_suffixes[i*new_size_line]), &(uc->suffixes[i*old_size_line]), old_size_line*sizeof(uint8_t));

    if (uc->nb_extended_annot != 0){

        for (i=0; i < uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT; i+=SIZE_BYTE_EXT_ANNOT){

            pos += ((((uint16_t)uc->suffixes[tot_size_line + i]) << SIZE_CELL) | ((uint16_t)uc->suffixes[tot_size_line + i + 1]));
            new_suffixes[pos * new_size_line + old_size_line] = uc->suffixes[tot_size_line + i + SIZE_BYTE_EXT_ANNOT - 1];
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

    for (i=0; i<nb_substring; i++)
        if (uc->suffixes[(i+1) * size_line - 1] != 0) nb_possible_annot_extend++;

    if ((nb_possible_annot_extend * SIZE_BYTE_EXT_ANNOT) > nb_substring) return;

    new_size_line = size_substring + uc->size_annot - 1;
    tot_new_size_line = nb_substring * new_size_line;

    new_tab_suffixes = calloc(tot_new_size_line + nb_possible_annot_extend * SIZE_BYTE_EXT_ANNOT
                              + uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes), sizeof(uint8_t));
    ASSERT_NULL_PTR(new_tab_suffixes,"create_annot_extended()")

    extend_annot = &(new_tab_suffixes[tot_new_size_line]);

    for (i=0; i<nb_substring; i++){
        memcpy(&(new_tab_suffixes[i*new_size_line]), &(uc->suffixes[i*size_line]), new_size_line*sizeof(uint8_t));

        if (uc->suffixes[(i+1) * size_line - 1] != 0){
            delta = i - old_pos;

            extend_annot[it_annot_extend] = delta >> SIZE_CELL;
            extend_annot[it_annot_extend+1] = delta & 0xff;
            extend_annot[it_annot_extend+2] = uc->suffixes[(i+1) * size_line - 1];
            old_pos = i;
            it_annot_extend += SIZE_BYTE_EXT_ANNOT;
        }
    }

    memcpy(&(new_tab_suffixes[tot_new_size_line + it_annot_extend]), &(uc->suffixes[nb_substring * size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
           uc->nb_cplx_nodes * (SIZE_BYTE_CPLX_N + uc->size_annot_cplx_nodes) * sizeof(uint8_t));

    uc->size_annot--;
    uc->nb_extended_annot = nb_possible_annot_extend;

    free(uc->suffixes);
    uc->suffixes = new_tab_suffixes;

    return;
}

uint8_t* realloc_annotation(UC* uc, int size_substring, int nb_substring, uint8_t new_size_annotation, int new_insertion, int pos_insert_extend){

    ASSERT_NULL_PTR(uc,"realloc_annotation()")

    if (nb_substring == 0){
        uc->size_annot = new_size_annotation;
        uc->nb_extended_annot = 0;
        return NULL;
    }

    ASSERT_NULL_PTR(uc->suffixes,"realloc_annotation()")

    int i = 0;
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

        for (i=0; i<nb_substring; i++) memcpy(&(new_tab_suffixes[i*new_size_line]), &(uc->suffixes[i*old_size_line]), old_size_line*sizeof(uint8_t));

        if (uc->nb_extended_annot != 0){
            int pos = 0;
            for (i=0; i < uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT; i += SIZE_BYTE_EXT_ANNOT){
                pos += ((((uint16_t)uc->suffixes[tot_size_line + i]) << SIZE_CELL) | ((uint16_t)uc->suffixes[tot_size_line + i + 1]));
                memcpy(&(new_tab_suffixes[pos * new_size_line + old_size_line]), &(uc->suffixes[tot_size_line + i + SIZE_BYTE_EXT_ANNOT - 1]), sizeof(uint8_t));
            }
        }

        if (new_insertion == 1){
            memmove(&(new_tab_suffixes[(pos_insert_extend+1) * new_size_line]),
                    &(new_tab_suffixes[pos_insert_extend * new_size_line]),
                    (nb_substring - pos_insert_extend) * new_size_line * sizeof(uint8_t));

            new_tab_suffixes[(nb_substring+1) * new_size_line] = pos_insert_extend >> SIZE_CELL;
            new_tab_suffixes[(nb_substring+1) * new_size_line + 1] = pos_insert_extend & 0xff;
            to_return = &(new_tab_suffixes[(nb_substring+1) * new_size_line + SIZE_BYTE_EXT_ANNOT - 1]);

            memcpy(&(to_return[1]), &(uc->suffixes[tot_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]), tot_size_line_cplx * sizeof(uint8_t));
        }
        else memcpy(&(new_tab_suffixes[nb_substring * new_size_line]), &(uc->suffixes[tot_size_line + uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT]),
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
            int old_pos = pos;
            int delta;
            int sum;

            i = 0;

            while (1){
                old_pos = pos;
                pos += ((((uint16_t)extended_annot[i*SIZE_BYTE_EXT_ANNOT]) << SIZE_CELL) | ((uint16_t)extended_annot[i*SIZE_BYTE_EXT_ANNOT+1]));
                if (pos < pos_insert_extend){
                    i++;
                    if (i >= uc->nb_extended_annot) break;
                }
                else break;
            }

            sum = i*SIZE_BYTE_EXT_ANNOT;

            if (i==0){
                memmove(&(extended_annot[SIZE_BYTE_EXT_ANNOT]), extended_annot, (uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + tot_size_line_cplx) * sizeof(uint8_t));

                extended_annot[0] = pos_insert_extend >> SIZE_CELL;
                extended_annot[1] = pos_insert_extend & 0xff;

                delta = pos - pos_insert_extend + new_insertion;
            }
            else if (i != uc->nb_extended_annot){
                memmove(&(extended_annot[sum + SIZE_BYTE_EXT_ANNOT]), &(extended_annot[sum]),
                        ((uc->nb_extended_annot - i) * SIZE_BYTE_EXT_ANNOT + tot_size_line_cplx) * sizeof(uint8_t));

                delta = pos_insert_extend - old_pos;

                extended_annot[sum] = delta >> SIZE_CELL;
                extended_annot[sum+1] = delta & 0xff;

                delta = pos - pos_insert_extend + new_insertion;
            }
            else delta = pos_insert_extend - pos;

            if (i == uc->nb_extended_annot){
                extended_annot[sum] = delta >> SIZE_CELL;
                extended_annot[sum+1] = delta & 0xff;
            }
            else{
                extended_annot[sum+SIZE_BYTE_EXT_ANNOT] = delta >> SIZE_CELL;
                extended_annot[sum+SIZE_BYTE_EXT_ANNOT+1] = delta & 0xff;
            }

            to_return = &(extended_annot[sum+SIZE_BYTE_EXT_ANNOT-1]);
        }
        else {
            uc->suffixes[tot_size_line] = pos_insert_extend >> SIZE_CELL;
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

        for (i=0; i<nb_substring; i++) memcpy(&(new_tab_suffixes[i*new_size_line]), &(uc->suffixes[i*old_size_line]), new_size_line*sizeof(uint8_t));

        if (new_insertion == 1){
            memmove(&(new_tab_suffixes[(pos_insert_extend+1) * new_size_line]),
                    &(new_tab_suffixes[pos_insert_extend * new_size_line]),
                    ((nb_substring - pos_insert_extend) * new_size_line)* sizeof(uint8_t));

            memset(&(new_tab_suffixes[pos_insert_extend * new_size_line]), 0, new_size_line*sizeof(uint8_t));
        }

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

int is_genome_present(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup, uint16_t id_genome){

    ASSERT_NULL_PTR(ann_inf, "is_genome_present()")
    ASSERT_NULL_PTR(annot, "is_genome_present()")

    if (id_genome >= NB_MAX_ID_GENOMES) return 0;
    if (size_annot + size_annot_sup == 0) return 0;

    int i = 0;
    int size = 0;
    int to_return = 0;

    if ((annot[0] & 0x3) == 3){
        uint32_t position = 0;
        while ((i<size_annot) && (annot[i] != 0)){
            if (i == 0) position = annot[0] >> 2;
            else if ((annot[i] & 0x1) == 1) position |= ((uint32_t)(annot[i] >> 1)) << (6+(i-1)*7);
            else break;
            i++;
        }

        if ((i >= size_annot) && (annot_sup != NULL)){
            i = 0;
            while ((i<size_annot_sup) && (annot_sup[i] != 0)){
                if ((annot_sup[i] & 0x1) == 1) position |= ((uint32_t)(annot_sup[i] >> 1)) << (6+(i+size_annot-1)*7);
                else break;
                i++;
            }
        }

        uint8_t* annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &size);
        memcpy(ann_inf->annotation, annot_tmp, size*sizeof(uint8_t));
    }
    else{
        size = size_annot;
        memcpy(ann_inf->annotation, annot, size_annot*sizeof(uint8_t));
        if (annot_sup != NULL){
            memcpy(&(ann_inf->annotation[size_annot]), annot_sup, size_annot_sup*sizeof(uint8_t));
            size += size_annot_sup;
        }
    }

    ann_inf->current_mode = ann_inf->annotation[0] & 0x3;

    if (ann_inf->current_mode == 0){
        if ((ann_inf->annotation[(id_genome+2)/SIZE_CELL] & MASK_POWER_8[(id_genome+2)%SIZE_CELL]) != 0) to_return = 1;
    }
    else if (ann_inf->current_mode == 1){ //<Present everywhere from x to y> mode

        while ((i<size) && ((ann_inf->annotation[i] & 0x3) == 1)){

            ann_inf->id_stored[ann_inf->nb_id_stored] = ann_inf->annotation[i] >> 2;
            ann_inf->nb_id_stored++;
            i++;

            if ((i<size) && ((ann_inf->annotation[i] & 0x3) == 2)){
                ann_inf->id_stored[ann_inf->nb_id_stored-1] = (ann_inf->id_stored[ann_inf->nb_id_stored-1] << 6) | (ann_inf->annotation[i] >> 2);
                i++;
            }

            if (ann_inf->nb_id_stored == 2){
                if ((id_genome >= ann_inf->id_stored[0]) && (id_genome <= ann_inf->id_stored[1])){
                    to_return = 1;
                    break;
                }

                ann_inf->nb_id_stored = 0;
            }
        }
    }
    else if (ann_inf->current_mode == 2){

        while ((i<size) && ((ann_inf->annotation[i] & 0x3) == 2)){

            ann_inf->id_stored[0] = ann_inf->annotation[i] >> 2;
            i++;

            if ((i<size) && ((ann_inf->annotation[i] & 0x3) == 1)){
                ann_inf->id_stored[0] = (ann_inf->id_stored[0] << 6) | (ann_inf->annotation[i] >> 2);
                i++;
            }

            if (id_genome == ann_inf->id_stored[0]){
                to_return = 1;
                break;
            }
        }
    }
    else ERROR( "is_genome_present(): mode 3, should not happen" )

    reinit_annotation_inform(ann_inf);

    return to_return;
}

void compute_best_mode(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup, uint16_t id_genome2insert){

    ASSERT_NULL_PTR(ann_inf,"compute_best_mode()")

    if (id_genome2insert >= NB_MAX_ID_GENOMES)
        ERROR( "compute_best_mode(): impossible to add new genome, the data structure contains already the maximum number of genomes possible to store" )

    ann_inf->last_added = INT_MIN;

    int nb_id_stored_extent = 0;
    int size;
    int i = 0;

    if (size_annot != 0){ //If the annotation is not new

        if ((annot[0] & 0x3) == 3){

            uint32_t position = 0;
            while ((i<size_annot) && (annot[i] != 0)){
                if (i == 0) position = annot[0] >> 2;
                else if ((annot[i] & 0x1) == 1) position |= ((uint32_t)(annot[i] >> 1)) << (6+(i-1)*7);
                else break;
                i++;
            }

            if ((i >= size_annot) && (annot_sup != NULL)){
                i = 0;
                while ((i<size_annot_sup) && (annot_sup[i] != 0)){
                    if ((annot_sup[i] & 0x1) == 1) position |= ((uint32_t)(annot_sup[i] >> 1)) << (6+(i+size_annot-1)*7);
                    else break;
                    i++;
                }
            }

            uint8_t* annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &size);
            memcpy(ann_inf->annotation, annot_tmp, size*sizeof(uint8_t));
        }
        else{
            size = size_annot;
            memcpy(ann_inf->annotation, annot, size_annot*sizeof(uint8_t));
            if (annot_sup != NULL){
                memcpy(&(ann_inf->annotation[size_annot]), annot_sup, size_annot_sup*sizeof(uint8_t));
                size += size_annot_sup;
            }
        }

        i = 0;

        ann_inf->current_mode = ann_inf->annotation[0] & 0x3;

        if (ann_inf->current_mode == 0){ //Bitwize mode

            int start_run = -1;

            if (size*SIZE_CELL-2 < MASK_POWER_8[6]){

                for (i=0; i<size*SIZE_CELL-2; i++){

                    if ((ann_inf->annotation[(i+2)/SIZE_CELL] & MASK_POWER_8[(i+2)%SIZE_CELL]) != 0){

                        ann_inf->last_added = i;
                        ann_inf->nb_genome_pres++;

                        if (start_run == -1){

                            start_run = 1;
                            nb_id_stored_extent++;
                        }
                    }
                    else if (start_run != -1){

                        nb_id_stored_extent++;
                        start_run = -1;
                    }
                }

                if (start_run != -1) nb_id_stored_extent++;
            }
            else{

                for (i=2; i<MASK_POWER_8[6]+2; i++){

                    if ((ann_inf->annotation[i/SIZE_CELL] & MASK_POWER_8[i%SIZE_CELL]) != 0){

                        ann_inf->last_added = i-2;
                        ann_inf->nb_genome_pres++;

                        if (start_run == -1){
                            start_run = 1;
                            nb_id_stored_extent++;
                        }
                    }
                    else if (start_run != -1){
                        start_run = -1;
                        nb_id_stored_extent++;
                    }
                }

                for (i=MASK_POWER_8[6]+2; i<size*SIZE_CELL; i++){

                    if ((ann_inf->annotation[i/SIZE_CELL] & MASK_POWER_8[i%SIZE_CELL]) != 0){

                        ann_inf->last_added = i-2;
                        ann_inf->nb_genome_pres++;
                        ann_inf->nb_ext_pres++;

                        if (start_run == -1){
                            nb_id_stored_extent += 2;
                            start_run = 1;
                        }
                    }
                    else if (start_run != -1){
                        nb_id_stored_extent += 2;
                        start_run = -1;
                    }
                }

                if (start_run != -1){
                    if (ann_inf->last_added < MASK_POWER_8[6]) nb_id_stored_extent++;
                    else nb_id_stored_extent += 2;
                }
            }
        }
        else if (ann_inf->current_mode == 1){ //<Present everywhere from x to y> mode

            while ((i<size) && ((ann_inf->annotation[i] & 0x3) == 1)){

                ann_inf->id_stored[ann_inf->nb_id_stored] = ann_inf->annotation[i] >> 2;
                ann_inf->nb_id_stored++;
                i++;

                if ((i<size) && ((ann_inf->annotation[i] & 0x3) == 2)){

                    ann_inf->id_stored[ann_inf->nb_id_stored-1] = (ann_inf->id_stored[ann_inf->nb_id_stored-1] << 6) | (ann_inf->annotation[i] >> 2);
                    i++;
                }
            }

            for (i=0; i<ann_inf->nb_id_stored; i+=2){

                ann_inf->nb_genome_pres += ann_inf->id_stored[i+1] - ann_inf->id_stored[i] + 1;

                if ((ann_inf->id_stored[i+1] >= MASK_POWER_8[6]) && (ann_inf->id_stored[i] >= MASK_POWER_8[6])){

                    ann_inf->nb_ext_pres += ann_inf->id_stored[i+1] - ann_inf->id_stored[i] + 1;
                    nb_id_stored_extent += 4;
                }
                else if (ann_inf->id_stored[i+1] >= MASK_POWER_8[6]){

                    ann_inf->nb_ext_pres += ann_inf->id_stored[i+1] - MASK_POWER_8[6] + 1;
                    nb_id_stored_extent += 3;
                }
                else nb_id_stored_extent += 2;
            }

            ann_inf->last_added = ann_inf->id_stored[ann_inf->nb_id_stored-1];
        }
        else if (ann_inf->current_mode == 2){ //<Present nowhere except in x> mode

            while ((i<size) && ((ann_inf->annotation[i] & 0x3) == 2)){

                ann_inf->id_stored[ann_inf->nb_id_stored] = ann_inf->annotation[i] >> 2;
                ann_inf->nb_genome_pres++;
                ann_inf->nb_id_stored++;
                i++;

                if ((i<size) && ((ann_inf->annotation[i] & 0x3) == 1)){
                    ann_inf->id_stored[ann_inf->nb_id_stored-1] = (ann_inf->id_stored[ann_inf->nb_id_stored-1] << 6) | (ann_inf->annotation[i] >> 2);
                    ann_inf->nb_ext_pres++;
                    i++;
                }
            }

            for (i=0; i<ann_inf->nb_id_stored; i++){

                if (i == 0){

                    if (ann_inf->id_stored[0] < MASK_POWER_8[6]) nb_id_stored_extent++;
                    else nb_id_stored_extent += 2;
                }
                else if (ann_inf->id_stored[i] != ann_inf->id_stored[i-1]+1){

                    nb_id_stored_extent += 2;

                    if (ann_inf->id_stored[i-1] >= MASK_POWER_8[6]) nb_id_stored_extent++;
                    if (ann_inf->id_stored[i] >=  MASK_POWER_8[6]) nb_id_stored_extent++;
                }
            }

            if (ann_inf->id_stored[ann_inf->nb_id_stored-1] < MASK_POWER_8[6]) nb_id_stored_extent++;
            else nb_id_stored_extent += 2;

            ann_inf->last_added = ann_inf->id_stored[ann_inf->nb_id_stored-1];
        }
        else ERROR( "compute_best_mode(): mode 3, should not happen" )
    }

    //New computed sizes for each mode after adding the current genome ID
    int new_sizes_mode[3] = {0, 0, 0};

    int size_to_add = 1;
    if (id_genome2insert >= MASK_POWER_8[6]) size_to_add = 2;

    //Compute the new possible sizes of the annot after insertion of the genonme
    new_sizes_mode[0] = CEIL(3+id_genome2insert, SIZE_CELL);

    if (id_genome2insert == ann_inf->last_added) new_sizes_mode[1] = nb_id_stored_extent;
    else if (id_genome2insert != ann_inf->last_added+1) new_sizes_mode[1] = nb_id_stored_extent+size_to_add*2;
    else if (id_genome2insert >= MASK_POWER_8[6]){
        if (ann_inf->last_added >= MASK_POWER_8[6]) new_sizes_mode[1] = nb_id_stored_extent;
        else new_sizes_mode[1] = nb_id_stored_extent+1;
    }
    else new_sizes_mode[1] = nb_id_stored_extent;

    if (id_genome2insert == ann_inf->last_added) new_sizes_mode[2] = ann_inf->nb_genome_pres + ann_inf->nb_ext_pres;
    else new_sizes_mode[2] = ann_inf->nb_genome_pres + ann_inf->nb_ext_pres + size_to_add;

    //Choose the mode which has the smaller size
    if (new_sizes_mode[2] <= new_sizes_mode[1]){
        ann_inf->min_mode = 2;
        ann_inf->min_size = new_sizes_mode[2];
    }
    else{
        ann_inf->min_mode = 1;
        ann_inf->min_size = new_sizes_mode[1];
    }

    if (ann_inf->min_size >= new_sizes_mode[0]){
        ann_inf->min_mode = 0;
        ann_inf->min_size = new_sizes_mode[0];
    }

    if ((size_annot != 0) && (new_sizes_mode[ann_inf->current_mode] == ann_inf->min_size) && (ann_inf->min_mode != ann_inf->current_mode))
        ann_inf->min_mode = ann_inf->current_mode;

    return;
}

void modify_mode_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup, uint16_t id_genome2insert){

    ASSERT_NULL_PTR(ann_inf, "modify_mode_annotation()")
    ASSERT_NULL_PTR(annot, "modify_mode_annotation()")

    uint16_t z = 0;
    uint16_t i = 0;
    int it_annot = 0;
    int size = size_annot;
    uint8_t* current_annot = annot;
    int size_current_annot = size_annot;

    if (annot_sup != NULL) size += size_annot_sup;

    memset(annot, 0, size_annot*sizeof(uint8_t));
    if (annot_sup != NULL) memset(annot_sup, 0, size_annot_sup*sizeof(uint8_t));

    if (ann_inf->current_mode == 0){ // Current mode is 0

        if (ann_inf->min_mode == 0){ // Mode with min size is 0

            memcpy(annot, ann_inf->annotation, size_annot*sizeof(uint8_t));
            if (annot_sup != NULL) memcpy(annot_sup, &(ann_inf->annotation[size_annot]), size_annot_sup*sizeof(uint8_t));

            int pos_cell = (id_genome2insert+2)/SIZE_CELL;

            if (pos_cell < size_annot) annot[pos_cell] |= MASK_POWER_8[(id_genome2insert+2)%SIZE_CELL];
            else annot_sup[pos_cell-size_annot] |= MASK_POWER_8[(id_genome2insert+2-(size_annot*SIZE_CELL))%SIZE_CELL];
        }
        else{
            if (ann_inf->min_mode == 1){ // Mode with min size is 1

                int start_run = -1;
                int stop_run = -1;

                for (z=0; z<=ann_inf->last_added; z++){

                    if ((ann_inf->annotation[(z+2)/SIZE_CELL] & MASK_POWER_8[(z+2)%SIZE_CELL]) != 0){

                        if (start_run == -1){
                            start_run = z;
                            stop_run = z;

                            modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, start_run, 1, 2);
                        }
                        else if (z == stop_run+1) stop_run++;
                    }
                    else if (stop_run != -1){

                        modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, stop_run, 1, 2);

                        start_run = -1;
                        stop_run = -1;
                    }
                }

                if (start_run == -1){

                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                }
                else if (id_genome2insert == stop_run+1){
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                }
                else {

                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, stop_run, 1, 2);

                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                }
            }
            else { // Mode with min size is 2
                while (z < size*SIZE_CELL-2){

                    if ((ann_inf->annotation[(z+2)/SIZE_CELL] & MASK_POWER_8[(z+2)%SIZE_CELL]) != 0)
                        modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, z, 2, 1);

                    z++;
                }

                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 2, 1);
            }
        }
    }
    else if (ann_inf->current_mode == 1){ // Current mode is 1

        if (ann_inf->min_mode == 0){ // Mode with min size is 0

            int pos_cell;
            int size_annot_bit = size_annot*SIZE_CELL;

            while (i < ann_inf->nb_id_stored){

                for (z = ann_inf->id_stored[i]; z <= ann_inf->id_stored[i+1]; z++){

                    pos_cell = (z+2)/SIZE_CELL;

                    if (pos_cell < size_annot) annot[pos_cell] |= MASK_POWER_8[(z+2)%SIZE_CELL];
                    else annot_sup[pos_cell-size_annot] |= MASK_POWER_8[(z+2-size_annot_bit)%SIZE_CELL];
                }

                i += 2;
            }

            pos_cell = (id_genome2insert+2)/SIZE_CELL;

            if (pos_cell < size_annot) annot[pos_cell] |= MASK_POWER_8[(id_genome2insert+2)%SIZE_CELL];
            else annot_sup[pos_cell-size_annot] |= MASK_POWER_8[(id_genome2insert+2-size_annot_bit)%SIZE_CELL];
        }
        else if (ann_inf->min_mode == 1){ // Mode with min size is 1

            memcpy(annot, ann_inf->annotation, size_annot * sizeof(uint8_t));
            if (annot_sup != NULL) memcpy(annot_sup, &(ann_inf->annotation[size_annot]), size_annot_sup*sizeof(uint8_t));

            if (annot_sup != NULL){

                it_annot = size_annot_sup-1;
                current_annot = annot_sup;
                size_current_annot = size_annot_sup;

                while (it_annot>=0){

                    if ((current_annot[it_annot] & 0x3) == 1) goto OUT2;
                    it_annot--;
                }
            }

            it_annot = size_annot-1;
            current_annot = annot;
            size_current_annot = size_annot;

            while (it_annot>=0){

                if ((current_annot[it_annot] & 0x3) == 1) break;
                it_annot--;
            }

            if (it_annot == -1) ERROR("modify_annotation() annotation.c")

            OUT2:if (id_genome2insert != ann_inf->id_stored[ann_inf->nb_id_stored-1]){

                if (id_genome2insert != ann_inf->id_stored[ann_inf->nb_id_stored-1]+1){

                        it_annot++;

                        if ((annot[it_annot] & 0x3) == 2) it_annot++;

                        /*if (it_annot >= size_annot){
                            it_annot = 0;
                            current_annot = annot_sup;
                            size_current_annot = size_annot_sup;
                        }*/

                        modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                        modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                }
                else modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
            }
        }
        else{ // Mode with min size is 2

            while (i<ann_inf->nb_id_stored){

                for (z = ann_inf->id_stored[i]; z <= ann_inf->id_stored[i+1]; z++)
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, z, 2, 1);

                i += 2;
            }

            modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 2, 1);
        }
    }
    else { // Current mode is 2

        if (ann_inf->min_mode == 0){ // Mode with min size is 0

            int pos_cell;
            int size_annot_bit = size_annot*SIZE_CELL;

            for (z=0; z<ann_inf->nb_id_stored; z++){

                pos_cell = (ann_inf->id_stored[z]+2)/SIZE_CELL;

                if (pos_cell < size_annot) annot[pos_cell] |= MASK_POWER_8[(ann_inf->id_stored[z]+2)%SIZE_CELL];
                else annot_sup[pos_cell-size_annot] |= MASK_POWER_8[(ann_inf->id_stored[z]+2-size_annot_bit)%SIZE_CELL];
            }

            pos_cell = (id_genome2insert+2)/SIZE_CELL;

            if (pos_cell < size_annot) annot[pos_cell] |= MASK_POWER_8[(id_genome2insert+2)%SIZE_CELL];
            else annot_sup[pos_cell-size_annot] |= MASK_POWER_8[(id_genome2insert+2-size_annot_bit)%SIZE_CELL];
        }
        else if (ann_inf->min_mode == 1){ // Mode with min size is 1

            for (z=0; z<ann_inf->nb_id_stored; z++){

                if (z == 0){
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, ann_inf->id_stored[z], 1, 2);
                }
                else if (ann_inf->id_stored[z] != ann_inf->id_stored[z-1]+1){
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, ann_inf->id_stored[z-1], 1, 2);
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, ann_inf->id_stored[z], 1, 2);
                }
            }

            if (ann_inf->nb_id_stored > 0){

                if (id_genome2insert == ann_inf->id_stored[ann_inf->nb_id_stored-1]+1){

                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                }
                else{
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, ann_inf->id_stored[ann_inf->nb_id_stored-1], 1, 2);

                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                }
            }
            else{
                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 1, 2);
            }
        }
        else { // Mode with min size is 2

            memcpy(annot, ann_inf->annotation, size_annot*sizeof(uint8_t));
            if (annot_sup != NULL) memcpy(annot_sup, &(ann_inf->annotation[size_annot]), size_annot_sup*sizeof(uint8_t));

            while (it_annot < size_current_annot){
                if ((current_annot[it_annot] & 0x3) == 0) goto OUT_HERE;
                it_annot++;
            }

            if (annot_sup != NULL){

                it_annot = 0;
                current_annot = annot_sup;
                size_current_annot = size_annot_sup;

                while (it_annot < size_current_annot){

                    if ((current_annot[it_annot] & 0x3) == 0) break;
                    it_annot++;
                }
            }

            OUT_HERE:modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, 2, 1);
        }
    }

    memset(ann_inf->annotation, 0, size*sizeof(uint8_t));

    return;
}

annotation_array_elem* sort_annotations(Pvoid_t* JArray_annot, int* size_array){

    PWord_t PValue_annot;
    Word_t * PValue_sizes;
    Word_t Rc_word;
    Word_t it_sizes;

    Pvoid_t PArray_sizes = (Pvoid_t) NULL;

    annotation_array_elem* annot_list = NULL;

    uint32_t next_id = 0;
    uint32_t next_id_tmp = 0;
    uint32_t tmp = 0;

    uint32_t next_id_highest_power2 = 1;
    uint32_t next_id_nb_bytes_power2 = 1;

    uint32_t next_id_tmp_nb_bytes_power2 = 1;

    uint32_t old_id = 0;

    uint32_t i = 0;
    uint32_t size_annot_list = 0;
    uint32_t pos_annot_list = 0;
    uint32_t tot_cell_index = NB_CELL_GENOMES + CEIL(NB_CELL_GENOMES, SIZE_CELL);

    uint32_t nb_annot_with_first_size;
    uint32_t first_size;

    uint8_t first_index[tot_cell_index];
    uint8_t it_index[tot_cell_index];

    memset(first_index, 255, tot_cell_index * sizeof(uint8_t));
    memset(it_index, 255, tot_cell_index * sizeof(uint8_t));
    first_index[tot_cell_index] = '\0';
    it_index[tot_cell_index] = '\0';

    PValue_annot = (void*)1;

    while ((PValue_annot != NULL) && (it_index[0] != 0)){

        JSLL(PValue_annot, *JArray_annot, it_index);

        memcpy(first_index, it_index, tot_cell_index * sizeof(uint8_t));
        first_size = first_index[0];

        if (annot_list == NULL){

            size_annot_list = first_size;
            *size_array = first_size;

            annot_list = calloc(size_annot_list, sizeof(annotation_array_elem));
            ASSERT_NULL_PTR(annot_list, "sort_annotations()")

            for (tmp = 0; tmp < size_annot_list; tmp++){
                annot_list[tmp].last_index = -1;
                annot_list[tmp].annot_array = NULL;
            }
        }

        nb_annot_with_first_size = 0;

        while ((PValue_annot != NULL) && (it_index[0] == first_size)){

            JLI(PValue_sizes, PArray_sizes, *PValue_annot);

            (*PValue_sizes)++;
            nb_annot_with_first_size++;

            JSLP(PValue_annot, *JArray_annot, it_index);
        }


        pos_annot_list = size_annot_list-first_size;
        annot_list[pos_annot_list].annot_array = calloc(nb_annot_with_first_size * first_size, sizeof(uint8_t));
        annot_list[pos_annot_list].size_annot = first_size;

        it_sizes = -1;
        tmp = 0;

        JLL(PValue_sizes, PArray_sizes, it_sizes);

        while (PValue_sizes != NULL){

            tmp = *PValue_sizes;
            next_id_tmp = next_id + tmp + 1;

            if (next_id_tmp > next_id_highest_power2) next_id_tmp_nb_bytes_power2 = get_nb_bytes_power2(round_up_next_highest_power2(next_id_tmp));

            if (next_id_nb_bytes_power2 >= first_size) *PValue_sizes = UINT32_MAX;
            else if (next_id_tmp_nb_bytes_power2 >= first_size) *PValue_sizes = UINT32_MAX;
            else {
                *PValue_sizes = next_id;
                next_id += tmp;

                if (next_id+1 > next_id_highest_power2){
                    next_id_highest_power2 = round_up_next_highest_power2(next_id+1);
                    next_id_nb_bytes_power2 = get_nb_bytes_power2(next_id_highest_power2);
                }
            }

            next_id_tmp_nb_bytes_power2 = next_id_nb_bytes_power2;

            JLP(PValue_sizes, PArray_sizes, it_sizes);
        }

        nb_annot_with_first_size = 0;

        memcpy(it_index, first_index, tot_cell_index * sizeof(uint8_t));

        JSLL(PValue_annot, *JArray_annot, it_index);

        while ((PValue_annot != NULL) && (it_index[0] == first_size)){

            JLG(PValue_sizes, PArray_sizes, *PValue_annot);

            if (*PValue_sizes == UINT32_MAX) *PValue_annot = UINT32_MAX;
            else{
                nb_annot_with_first_size++;

                memcpy(&(annot_list[pos_annot_list].annot_array[((*PValue_sizes)-old_id)*first_size]), &(it_index[1]), first_size*sizeof(uint8_t));

                for (i=1; i< first_size + 1; i++){
                    if ((it_index[i] == 254) && ((it_index[first_size + 1 + (i-1)/(SIZE_CELL-1)] & MASK_POWER_8[(i-1) % (SIZE_CELL-1) + 1]) != 0))
                        annot_list[pos_annot_list].annot_array[((*PValue_sizes)-old_id) * first_size + i - 1] = 0;
                }

                *PValue_annot = *PValue_sizes;
                *PValue_sizes += 1;
            }

            JSLP(PValue_annot, *JArray_annot, it_index);
        }

        annot_list[pos_annot_list].last_index = next_id - 1;

        JLFA(Rc_word, PArray_sizes);
        annot_list[pos_annot_list].annot_array = realloc(annot_list[pos_annot_list].annot_array, nb_annot_with_first_size * first_size * sizeof(uint8_t));

        old_id = next_id;
    }

    printf("next_id = %d\n", next_id);

    int test = 0;
    for (i = 0; i < size_annot_list; i++){
        if (annot_list[i].last_index == -1) annot_list[i].last_index = test;
        else test = annot_list[i].last_index;
    }

    return annot_list;
}

annotation_array_elem* intersection_annotations(uint8_t* annot1, int size_annot1, uint8_t* annot_sup1, int size_annot_sup1, uint8_t* annot2, int size_annot2,
                                  uint8_t* annot_sup2, int size_annot_sup2, int id_genome_max, annotation_array_elem* annot_sorted){

    int flag1, flag2;

    int i = 0;
    int size_annot_flag0 = MAX(CEIL(id_genome_max, SIZE_CELL), 1);
    int size1 = size_annot1+size_annot_sup1;
    int size2 = size_annot2+size_annot_sup2;

    uint8_t* annot1_flag0 = calloc(size_annot_flag0, sizeof(uint8_t));
    uint8_t* annot2_flag0 = calloc(size_annot_flag0, sizeof(uint8_t));

    uint8_t* annot1_real = annot1;
    uint8_t* annot2_real = annot2;

    if (annot1 != NULL){
        if ((annot1[0] & 0x3) == 3){
            int position = 0;
            while ((i<size_annot1) && (annot1[i] != 0)){
                if (i == 0) position = annot1[0] >> 2;
                else if ((annot1[i] & 0x1) == 1) position |= ((int)(annot1[i] >> 1)) << (6+(i-1)*7);
                else break;
                i++;
            }

            if ((i >= size_annot1) && (annot_sup1 != NULL)){
                i = 0;
                while ((i<size_annot_sup1) && (annot_sup1[i] != 0)){
                    if ((annot_sup1[i] & 0x1) == 1) position |= ((int)(annot_sup1[i] >> 1)) << (6+(i+size_annot1-1)*7);
                    else break;
                    i++;
                }
            }

            annot1_real = extract_from_annotation_array_elem(annot_sorted, position, &size1);
        }
        else if (annot_sup1 != NULL){
            annot1_real = malloc(size1*sizeof(uint8_t));
            memcpy(annot1_real, annot1, size_annot1*sizeof(uint8_t));
            annot1_real[size_annot1] = annot_sup1[0];
        }
    }

    i = 0;

    if (annot2 != NULL){
        if ((annot2[0] & 0x3) == 3){
            int position = 0;
            while ((i<size_annot2) && (annot2[i] != 0)){
                if (i == 0) position = annot2[0] >> 2;
                else if ((annot2[i] & 0x1) == 1) position |= ((int)(annot2[i] >> 1)) << (6+(i-1)*7);
                else break;
                i++;
            }

            if ((i >= size_annot2) && (annot_sup2 != NULL)){
                i = 0;
                while ((i<size_annot_sup2) && (annot_sup2[i] != 0)){
                    if ((annot_sup2[i] & 0x1) == 1) position |= ((int)(annot_sup2[i] >> 1)) << (6+(i+size_annot2-1)*7);
                    else break;
                    i++;
                }
            }

            annot2_real = extract_from_annotation_array_elem(annot_sorted, position, &size2);
        }
        else if (annot_sup2 != NULL){
            annot2_real = malloc(size2*sizeof(uint8_t));
            memcpy(annot2_real, annot2, size_annot2*sizeof(uint8_t));
            annot2_real[size_annot2] = annot_sup2[0];
        }
    }

    i = 0;

    if (annot1 != NULL){
        flag1 = annot1_real[0] & 0x3;

        if (flag1 == 0){
            memcpy(annot1_flag0, annot1, size_annot1*sizeof(uint8_t));
            memcpy(&(annot1_flag0[size_annot1]), annot_sup1, size_annot_sup1*sizeof(uint8_t));
        }
        else if (flag1 == 1){
            uint16_t start, stop, z;

                while (i<size1){
                    if ((annot1[i] & 0x3) == 1){
                        start = annot1[i] >> 2;
                        i++;

                        if ((i<size1) && ((annot1[i] & 0x3) == 2)){
                            start = (start << 6) | (annot1[i] >> 2);
                            i++;
                        }
                    }
                    if ((annot1[i] & 0x3) == 1){
                        stop = annot1[i] >> 2;
                        i++;

                        if ((i<size1) && ((annot1[i] & 0x3) == 2)){
                            stop = (stop << 6) | (annot1[i] >> 2);
                            i++;
                        }
                    }

                    for (z=start; z<=stop; z++) annot1_flag0[(z+2)/SIZE_CELL] |= MASK_POWER_8[(z+2)%SIZE_CELL];
                }
        }
        else if (flag1 == 2){
            uint16_t pos;

            while (i<size1){
                if ((annot1[i] & 0x3) == 2){
                    pos = annot1[i] >> 2;
                    i++;

                    if ((i<size1) && ((annot1[i] & 0x3) == 1)){
                        pos = (pos << 6) | (annot1[i] >> 2);
                        i++;
                    }
                }

                annot1_flag0[(pos+2)/SIZE_CELL] |= MASK_POWER_8[(pos+2)%SIZE_CELL];
            }
        }
    }
    else memset(annot1_flag0, 255, size_annot_flag0 * sizeof(uint8_t));

    i = 0;

    if (annot2 != NULL){

        flag2 = annot2_real[0] & 0x3;

        if (flag2 == 0){
            memcpy(annot2_flag0, annot2, size_annot2*sizeof(uint8_t));
            memcpy(&(annot2_flag0[size_annot2]), annot_sup2, size_annot_sup2*sizeof(uint8_t));
        }
        else if (flag2 == 1){
            uint16_t start, stop, z;

                while (i<size2){
                    if ((annot2[i] & 0x3) == 1){
                        start = annot2[i] >> 2;
                        i++;

                        if ((i<size2) && ((annot2[i] & 0x3) == 2)){
                            start = (start << 6) | (annot2[i] >> 2);
                            i++;
                        }
                    }
                    if ((annot2[i] & 0x3) == 1){
                        stop = annot2[i] >> 2;
                        i++;

                        if ((i<size2) && ((annot2[i] & 0x3) == 2)){
                            stop = (stop << 6) | (annot2[i] >> 2);
                            i++;
                        }
                    }

                    for (z=start; z<=stop; z++) annot2_flag0[(z+2)/SIZE_CELL] |= MASK_POWER_8[(z+2)%SIZE_CELL];
                }
        }
        else if (flag2 == 2){
            uint16_t pos;

            while (i<size2){
                if ((annot2[i] & 0x3) == 2){
                    pos = annot2[i] >> 2;
                    i++;

                    if ((i<size2) && ((annot2[i] & 0x3) == 1)){
                        pos = (pos << 6) | (annot2[i] >> 2);
                        i++;
                    }
                }

                annot2_flag0[(pos+2)/SIZE_CELL] |= MASK_POWER_8[(pos+2)%SIZE_CELL];
            }
        }
    }
    else memset(annot2_flag0, 255, size_annot_flag0 * sizeof(uint8_t));

    for (i=0; i<size_annot_flag0; i++) annot1_flag0[i] &= annot2_flag0[i];

    free(annot2_flag0);
    if ((annot1 != NULL) && ((annot1[0] & 0x3) != 3) && (annot_sup1 != NULL)) free(annot1_real);
    if ((annot2 != NULL) && ((annot2[0] & 0x3) != 3) && (annot_sup2 != NULL)) free(annot2_real);

    annotation_array_elem* ann_arr_elem = malloc(sizeof(annotation_array_elem));
    ASSERT_NULL_PTR(ann_arr_elem, "intersection_annotations()");

    ann_arr_elem->annot_array = annot1_flag0;
    ann_arr_elem->size_annot = size_annot_flag0;

    return ann_arr_elem;
}
