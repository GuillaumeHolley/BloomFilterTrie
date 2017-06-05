#include "annotation.h"

const uint8_t MASK_POWER_8[8] = {1, 2, 4, 8, 16, 32, 64, 128};

extern inline annotation_inform* create_annotation_inform(int nb_id_genomes, bool disable_flag_0);
extern inline void reinit_annotation_inform(annotation_inform* ann_inf);
extern inline void free_annotation_inform(annotation_inform* ann_inf);
extern inline int size_annot_sub(uint8_t* annot, int size_substring, int size_annot);
extern inline int modify_annot_bis(uint8_t** current_annot, uint8_t* annot_sup, int* it_annot, int* size_current_annot,
                            uint32_t id_genome, int size_id_genome, uint8_t flag, uint8_t flag_ext);
extern inline UC_SIZE_ANNOT_T *min_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot);
extern inline int max_size_per_sub(uint8_t* annot, int nb_substrings, int size_substring, int size_annot);
extern inline uint8_t* extract_from_annotation_array_elem(annotation_array_elem* annot_sorted, uint32_t position,
                                                   int* size_annot);
extern inline int getSize_from_annotation_array_elem(annotation_array_elem* annot_sorted, uint32_t position);
extern inline double getTotalSize_annotation_array_elem(annotation_array_elem* annot_sorted, int size_array);
extern inline int getMaxSize_annotation_array_elem(annotation_array_elem* annot_sorted);
extern inline void free_annotation_array_elem(annotation_array_elem** annot_sorted, int* size_array);

int is_genome_present(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot,
                      int size_annot, uint8_t* annot_sup, int size_annot_sup, uint32_t id_genome){

    ASSERT_NULL_PTR(ann_inf, "is_genome_present()")
    ASSERT_NULL_PTR(annot, "is_genome_present()")

    if (id_genome >= NB_MAX_ID_GENOMES) return 0;
    if (size_annot + size_annot_sup == 0) return 0;

    int i = 1;
    int size = 0;
    int to_return = 0;

    if ((annot[0] & 0x3) == 3){

        uint32_t position = annot[0] >> 2;

        while ((i<size_annot) && (IS_ODD(annot[i]))){
            position |= ((uint32_t)(annot[i] >> 1)) << (6+(i-1)*7);
            i++;
        }

        if ((i >= size_annot) && (annot_sup != NULL)){
            i = 0;
            while ((i<size_annot_sup) && (IS_ODD(annot_sup[i]))){
                position |= ((uint32_t)(annot_sup[i] >> 1)) << (6+(i+size_annot-1)*7);
                i++;
            }
        }

        uint8_t* annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &size);
        memcpy(ann_inf->annotation, annot_tmp, size*sizeof(uint8_t));

        i = decomp_annotation(ann_inf, ann_inf->annotation, size, NULL, 0, 0);
        if (i) size = i;
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
    ann_inf->size_annot = size;
    ann_inf->current_mode = ann_inf->annotation[0] & 0x3;

    if (ann_inf->current_mode == 0){
        if (ann_inf->annotation[(id_genome+2)/SIZE_BITS_UINT_8T] & MASK_POWER_8[(id_genome+2)%SIZE_BITS_UINT_8T]) to_return = 1;
    }
    else if (ann_inf->current_mode == 1){ //<Present everywhere from x to y> mode

        if (ann_inf->comp_annot > 0){

            for (i = 0; i < ann_inf->nb_id_stored; i += 2){

                if ((id_genome >= ann_inf->id_stored[i]) && (id_genome <= ann_inf->id_stored[i+1])){

                    to_return = 1;
                    break;
                }
                else if (id_genome < ann_inf->id_stored[i+1]) break;
            }
        }
        else{

            while ((i<size) && (ann_inf->annotation[i] & 0x1)){

                ann_inf->id_stored[ann_inf->nb_id_stored] = ann_inf->annotation[i] >> 2;
                ann_inf->nb_id_stored++;
                i++;

                while ((i<size) && (ann_inf->annotation[i] & 0x2)){
                    ann_inf->id_stored[ann_inf->nb_id_stored-1] = (ann_inf->id_stored[ann_inf->nb_id_stored-1] << 6) | (ann_inf->annotation[i] >> 2);
                    i++;
                }

                if (ann_inf->nb_id_stored == 2){

                    if ((id_genome >= ann_inf->id_stored[0]) && (id_genome <= ann_inf->id_stored[1])){
                        to_return = 1;
                        break;
                    }
                    else if (id_genome < ann_inf->id_stored[0]) break;

                    ann_inf->nb_id_stored = 0;
                }
            }
        }
    }
    else if (ann_inf->current_mode == 2){

        if (ann_inf->comp_annot > 0){

            for (i = 0; i < ann_inf->nb_id_stored; i++){

                if (id_genome == ann_inf->id_stored[i]){

                    to_return = 1;
                    break;
                }
                else if (id_genome < ann_inf->id_stored[i]) break;
            }
        }
        else{

            while ((i<size) && (ann_inf->annotation[i] & 0x2)){

                ann_inf->id_stored[0] = ann_inf->annotation[i] >> 2;
                i++;

                while ((i<size) && (ann_inf->annotation[i] & 0x1)){

                    ann_inf->id_stored[0] = (ann_inf->id_stored[0] << 6) | (ann_inf->annotation[i] >> 2);
                    i++;
                }

                if (id_genome == ann_inf->id_stored[0]){

                    to_return = 1;
                    break;
                }
                else if (id_genome < ann_inf->id_stored[0]) break;
            }
        }
    }
    else ERROR( "is_genome_present(): mode 3, should not happen" )

    reinit_annotation_inform(ann_inf);

    return to_return;
}

int is_genome_present_from_end_annot(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot,
                      int size_annot, uint8_t* annot_sup, int size_annot_sup, uint32_t id_genome){

    ASSERT_NULL_PTR(ann_inf, "is_genome_present_from_end_annot() 1")
    ASSERT_NULL_PTR(annot, "is_genome_present_from_end_annot() 2")

    if (id_genome >= NB_MAX_ID_GENOMES) return 0;
    if (size_annot + size_annot_sup == 0) return 0;

    int j = 0;
    int size = 0;
    int to_return = 0;

    int i = 1;

    if ((annot[0] & 0x3) == 3){

        uint32_t position = annot[0] >> 2;

        while ((i<size_annot) && (IS_ODD(annot[i]))){
            position |= ((uint32_t)(annot[i] >> 1)) << (6+(i-1)*7);
            i++;
        }

        if ((i >= size_annot) && (annot_sup != NULL)){
            i = 0;
            while ((i<size_annot_sup) && (IS_ODD(annot_sup[i]))){
                position |= ((uint32_t)(annot_sup[i] >> 1)) << (6+(i+size_annot-1)*7);
                i++;
            }
        }

        uint8_t* annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &size);
        memcpy(ann_inf->annotation, annot_tmp, size*sizeof(uint8_t));

        i = decomp_annotation(ann_inf, annot_tmp, size, NULL, 0, 0);
        if (i) size = i;
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
    ann_inf->size_annot = size;
    ann_inf->current_mode = ann_inf->annotation[0] & 0x3;

    if (ann_inf->current_mode == 0){
        if (ann_inf->annotation[(id_genome+2)/SIZE_BITS_UINT_8T] & MASK_POWER_8[(id_genome+2)%SIZE_BITS_UINT_8T]) to_return = 1;
    }
    else if (ann_inf->current_mode == 1){ //<Present everywhere from x to y> mode

        if (ann_inf->comp_annot > 0){

            for (i = ann_inf->nb_id_stored - 1; i >= 0; i -= 2){

                if ((id_genome >= ann_inf->id_stored[i-1]) && (id_genome <= ann_inf->id_stored[i])){

                    to_return = 1;
                    goto END_IS_GENOME_PRES_END;
                }
                else if (id_genome > ann_inf->id_stored[i]) goto END_IS_GENOME_PRES_END;
            }
        }
        else{

            j = size-1;
            while ((j >= 0) && ((ann_inf->annotation[j] & 0x3) != 1)) j--;
            i = j;

            while ((i<size) && (ann_inf->annotation[i] & 0x1)){

                ann_inf->id_stored[ann_inf->nb_id_stored] = ann_inf->annotation[i] >> 2;
                i++;

                while ((i<size) && (ann_inf->annotation[i] & 0x2)){

                    ann_inf->id_stored[ann_inf->nb_id_stored] = (ann_inf->id_stored[ann_inf->nb_id_stored] << 6) | (ann_inf->annotation[i] >> 2);
                    i++;
                }

                ann_inf->nb_id_stored++;

                if (ann_inf->nb_id_stored == 2){

                    if ((id_genome >= ann_inf->id_stored[1]) && (id_genome <= ann_inf->id_stored[0])){

                        to_return = 1;
                        goto END_IS_GENOME_PRES_END;
                    }
                    else if (id_genome > ann_inf->id_stored[1]) goto END_IS_GENOME_PRES_END;

                    ann_inf->nb_id_stored = 0;
                }

                j--;
                while ((j >= 0) && ((ann_inf->annotation[j] & 0x3) != 1)) j--;
                i = j;
            }
        }
    }
    else if (ann_inf->current_mode == 2){

        if (ann_inf->comp_annot > 0){
            for (i = ann_inf->nb_id_stored - 1; i >= 0; i--){
                if (id_genome == ann_inf->id_stored[i]){
                    to_return = 1;
                    goto END_IS_GENOME_PRES_END;
                }
                else if (id_genome > ann_inf->id_stored[i]) goto END_IS_GENOME_PRES_END;
            }
        }
        else{

            j = size-1;
            while ((j >= 0) && ((ann_inf->annotation[j] & 0x3) != 2)) j--;
            i = j;

            while ((i<size) && (ann_inf->annotation[i] & 0x2)){

                ann_inf->id_stored[0] = ann_inf->annotation[i] >> 2;
                i++;

                while ((i<size) && (ann_inf->annotation[i] & 0x1)){

                    ann_inf->id_stored[0] = (ann_inf->id_stored[0] << 6) | (ann_inf->annotation[i] >> 2);
                    i++;
                }

                if (id_genome == ann_inf->id_stored[0]){
                    to_return = 1;
                    goto END_IS_GENOME_PRES_END;
                }
                else if (id_genome > ann_inf->id_stored[0]) goto END_IS_GENOME_PRES_END;

                j--;
                while ((j >= 0) && ((ann_inf->annotation[j] & 0x3) != 2)) j--;
                i = j;
            }
        }
    }
    else ERROR( "is_genome_present_from_end_annot(): mode 3, should not happen" )

    END_IS_GENOME_PRES_END: reinit_annotation_inform(ann_inf);

    return to_return;
}

int get_last_genome_inserted(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot,
                      int size_annot, uint8_t* annot_sup, int size_annot_sup, uint32_t* last_id_genome){

    ASSERT_NULL_PTR(ann_inf, "is_genome_present_from_end_annot() 1")
    ASSERT_NULL_PTR(annot, "is_genome_present_from_end_annot() 2")

    if (size_annot + size_annot_sup == 0) ERROR("is_genome_present_from_end_annot() 3")

    int j = 0;
    int size = 0;

    int i = 1;

    if ((annot[0] & 0x3) == 3){

        uint32_t position = annot[0] >> 2;

        while ((i<size_annot) && (IS_ODD(annot[i]))){
            position |= ((uint32_t)(annot[i] >> 1)) << (6+(i-1)*7);
            i++;
        }

        if ((i >= size_annot) && (annot_sup != NULL)){
            i = 0;
            while ((i<size_annot_sup) && (IS_ODD(annot_sup[i]))){
                position |= ((uint32_t)(annot_sup[i] >> 1)) << (6+(i+size_annot-1)*7);
                i++;
            }
        }

        uint8_t* annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &size);
        memcpy(ann_inf->annotation, annot_tmp, size*sizeof(uint8_t));

        decomp_annotation(ann_inf, annot_tmp, size, NULL, 0, 0);
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
    j = size-1;

    ann_inf->size_annot = size;
    ann_inf->current_mode = ann_inf->annotation[0] & 0x3;

    if (ann_inf->current_mode == 0){

        while ((j >= 0) && (ann_inf->annotation[j] != 0)) j--;
        *last_id_genome = MAX(0, j * SIZE_BITS_UINT_8T);
        i = j;

        for (j = SIZE_BITS_UINT_8T-1; j >= 0; j--){
            if ((ann_inf->annotation[j/SIZE_BITS_UINT_8T] & MASK_POWER_8[j%SIZE_BITS_UINT_8T]) != 0)
                *last_id_genome += j - 2;
        }
    }
    else if (ann_inf->current_mode == 1){ //<Present everywhere from x to y> mode

        if (ann_inf->comp_annot > 0) *last_id_genome = ann_inf->id_stored[ann_inf->nb_id_stored - 1];
        else{

            while ((j >= 0) && ((ann_inf->annotation[j] & 0x3) != 1)) j--;
            i = j;

            ann_inf->id_stored[0] = ann_inf->annotation[j] >> 2;
            j++;

            while ((j<size) && (ann_inf->annotation[j] & 0x2)){
                ann_inf->id_stored[0] = (ann_inf->id_stored[0] << 6) | (ann_inf->annotation[j] >> 2);
                j++;
            }

            *last_id_genome = ann_inf->id_stored[0];
        }
    }
    else if (ann_inf->current_mode == 2){

        if (ann_inf->comp_annot > 0) *last_id_genome = ann_inf->id_stored[ann_inf->nb_id_stored - 1];
        else{

            while ((j >= 0) && ((ann_inf->annotation[j] & 0x3) != 2)) j--;
            i = j;

            ann_inf->id_stored[0] = ann_inf->annotation[j] >> 2;
            j++;

            while ((j<size) && (ann_inf->annotation[j] & 0x1)){
                ann_inf->id_stored[0] = (ann_inf->id_stored[0] << 6) | (ann_inf->annotation[j] >> 2);
                j++;
            }

            *last_id_genome = ann_inf->id_stored[0];
        }
    }
    else ERROR( "is_genome_present_from_end_annot(): mode 3, should not happen" )

    reinit_annotation_inform(ann_inf);

    return i;
}

void compute_best_mode(annotation_inform* ann_inf, annotation_array_elem* annot_sorted,
                       uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup,
                       uint32_t id_genome2insert, int size_id_genome){

    ASSERT_NULL_PTR(ann_inf,"compute_best_mode()")

    int size;

    int i = 1;
    int tmp = 1, tmp2 = 1;
    int tot_size_ids_flag1 = 0;
    int tot_size_ids_flag2 = 0;
    int lim_flag0;

    int new_sizes_mode[3] = {0};

    uint32_t pow2_imin = 0;
    uint32_t imin, imax;

    bool start_run = false;

    ann_inf->last_added = INT_MIN;

    if (size_annot){ //If the annotation is not new

        if ((annot[0] & 0x3) == 3){

            uint32_t position = annot[0] >> 2;

            while ((i<size_annot) && (IS_ODD(annot[i]))){
                position |= ((uint32_t)(annot[i] >> 1)) << (6+(i-1)*7);
                i++;
            }

            if ((i >= size_annot) && (annot_sup != NULL)){
                i = 0;
                while ((i<size_annot_sup) && (IS_ODD(annot_sup[i]))){
                    position |= ((uint32_t)(annot_sup[i] >> 1)) << (6+(i+size_annot-1)*7);
                    i++;
                }
            }

            uint8_t* annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &size);
            memcpy(ann_inf->annotation, annot_tmp, size*sizeof(uint8_t));

            if ((i = decomp_annotation(ann_inf, annot_tmp, size, NULL, 0, 1))) size = i;
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

        ann_inf->size_annot = size;
        ann_inf->current_mode = ann_inf->annotation[0] & 0x3;

        if (ann_inf->current_mode == 0){ //Bitwize mode

            if (size * SIZE_BITS_UINT_8T - 2 < MASK_POWER_8[6]) lim_flag0 = size * SIZE_BITS_UINT_8T;
            else lim_flag0 = 66;

            for (i=2; i < lim_flag0; i++){

                if (ann_inf->annotation[i/SIZE_BITS_UINT_8T] & MASK_POWER_8[i%SIZE_BITS_UINT_8T]){
                    ann_inf->last_added = i-2;
                    tot_size_ids_flag2++;
                    tot_size_ids_flag1 += start_run == false;
                    start_run = true;
                }
                else {
                    tot_size_ids_flag1 += start_run == true;
                    start_run = false;
                }
            }

            if (size * SIZE_BITS_UINT_8T-2 < MASK_POWER_8[6]) tot_size_ids_flag1 += start_run == true;
            else{

                for (i=66; i<size*SIZE_BITS_UINT_8T; i++){

                    if ((ann_inf->annotation[i/SIZE_BITS_UINT_8T] & MASK_POWER_8[i%SIZE_BITS_UINT_8T]) != 0){

                        ann_inf->last_added = i-2;

                        if (ann_inf->last_added >= pow2_imin){
                            pow2_imin = round_up_next_highest_power2(ann_inf->last_added);
                            tmp2 = (tmp = get_nb_bytes_power2_annot_bis(ann_inf->last_added, pow2_imin));
                        }

                        tot_size_ids_flag2 += tmp;
                        tot_size_ids_flag1 += (start_run == false) * tmp;
                        start_run = true;
                    }
                    else {
                        if (start_run){

                            if (i-2 >= pow2_imin){
                                pow2_imin = round_up_next_highest_power2(i-2);
                                tmp2 = get_nb_bytes_power2_annot_bis(i-2, pow2_imin);
                            }

                            tot_size_ids_flag1 += tmp2;
                        }

                        start_run = false;
                    }
                }

                if (start_run){
                    if (ann_inf->last_added < 0x40) tot_size_ids_flag1++;
                    else tot_size_ids_flag1 += tmp;
                }
            }
        }
        else if (ann_inf->current_mode == 1){ //<Present everywhere from x to y> mode

            if (ann_inf->comp_annot <= 0){

                while ((i<size) && (ann_inf->annotation[i] & 0x1)){

                    ann_inf->id_stored[ann_inf->nb_id_stored] = ann_inf->annotation[i] >> 2;
                    ann_inf->size_id_stored[ann_inf->nb_id_stored] = 1;
                    i++;

                    while ((i<size) && (ann_inf->annotation[i] & 0x2)){

                        ann_inf->id_stored[ann_inf->nb_id_stored] = (ann_inf->id_stored[ann_inf->nb_id_stored] << 6) | (ann_inf->annotation[i] >> 2);
                        ann_inf->size_id_stored[ann_inf->nb_id_stored]++;
                        i++;
                    }

                    ann_inf->nb_id_stored++;
                }

                ann_inf->last_added = ann_inf->id_stored[ann_inf->nb_id_stored-1];
            }

            for (i=0; i<ann_inf->nb_id_stored; i+=2){

                tot_size_ids_flag1 += ann_inf->size_id_stored[i+1] + ann_inf->size_id_stored[i];

                if (ann_inf->size_id_stored[i] == ann_inf->size_id_stored[i+1]){
                    tot_size_ids_flag2 += (ann_inf->id_stored[i+1] - ann_inf->id_stored[i] + 1) * ann_inf->size_id_stored[i];
                }
                else {

                    imin = ann_inf->id_stored[i];
                    imax = ann_inf->id_stored[i+1];

                    while (imin <= imax){

                        pow2_imin = round_up_next_highest_power2(imin);

                        if (pow2_imin < imax) tot_size_ids_flag2 += (pow2_imin - imin) * get_nb_bytes_power2_annot_bis(imin, pow2_imin);
                        else tot_size_ids_flag2 += (imax - imin) * get_nb_bytes_power2_annot_bis(imin, pow2_imin);

                        tot_size_ids_flag2 += get_nb_bytes_power2_annot_bis(pow2_imin, pow2_imin);
                        imin = pow2_imin + 1;
                    }
                }
            }
        }
        else if (ann_inf->current_mode == 2){ //<Present nowhere except in x> mode

            if (ann_inf->comp_annot <= 0){

                while ((i<size) && (ann_inf->annotation[i] & 0x2)){

                    ann_inf->id_stored[ann_inf->nb_id_stored] = ann_inf->annotation[i] >> 2;
                    ann_inf->size_id_stored[ann_inf->nb_id_stored] = 1;
                    i++;

                    while ((i<size) && (ann_inf->annotation[i] & 0x1)){
                        ann_inf->id_stored[ann_inf->nb_id_stored] = (ann_inf->id_stored[ann_inf->nb_id_stored] << 6) | (ann_inf->annotation[i] >> 2);
                        ann_inf->size_id_stored[ann_inf->nb_id_stored]++;
                        i++;
                    }

                    ann_inf->nb_id_stored++;
                }

                ann_inf->last_added = ann_inf->id_stored[ann_inf->nb_id_stored-1];
            }

            tot_size_ids_flag2 += ann_inf->size_id_stored[0];
            tot_size_ids_flag1 += ann_inf->size_id_stored[0];

            for (i=1; i<ann_inf->nb_id_stored; i++){

                tot_size_ids_flag2 += ann_inf->size_id_stored[i];

                if (ann_inf->id_stored[i] != ann_inf->id_stored[i-1]+1)
                    tot_size_ids_flag1 += ann_inf->size_id_stored[i-1] + ann_inf->size_id_stored[i];
            }

            tot_size_ids_flag1 += ann_inf->size_id_stored[ann_inf->nb_id_stored-1];
        }
        else ERROR( "compute_best_mode(): mode 3, should not happen" )
    }

    //Compute the new possible sizes of the annot after insertion of the genome
    if (ann_inf->disabled_flags & MASK_POWER_8[0]) new_sizes_mode[0] = INT_MAX;
    else new_sizes_mode[0] = CEIL(3+id_genome2insert, SIZE_BITS_UINT_8T);

    new_sizes_mode[1] = tot_size_ids_flag1;
    new_sizes_mode[2] = tot_size_ids_flag2;

    if (id_genome2insert != ann_inf->last_added){

        if (id_genome2insert != ann_inf->last_added+1) new_sizes_mode[1] += size_id_genome * 2;
        else if (ann_inf->current_mode == 0) new_sizes_mode[1] += size_id_genome - tmp;
        else new_sizes_mode[1] += size_id_genome - ann_inf->size_id_stored[ann_inf->nb_id_stored-1];
    }

    if (id_genome2insert != ann_inf->last_added) new_sizes_mode[2] += size_id_genome;

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

void modify_mode_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot, uint8_t* annot_sup,
                            int size_annot_sup, uint32_t id_genome2insert, int size_id_genome){

    ASSERT_NULL_PTR(ann_inf, "modify_mode_annotation()")
    ASSERT_NULL_PTR(annot, "modify_mode_annotation()")

    int z = 0, i = 0, it_annot = 0;

    int pos_cell;
    int start_run, stop_run;

    int size = size_annot;
    int size_current_annot = size_annot;

    uint8_t* current_annot = annot;

    if (annot_sup != NULL) size += size_annot_sup;

    memset(annot, 0, size_annot*sizeof(uint8_t));
    if (annot_sup != NULL) memset(annot_sup, 0, size_annot_sup*sizeof(uint8_t));

    if (ann_inf->current_mode == 0){ // Current mode is 0

        if (ann_inf->min_mode == 0){ // Mode with min size is 0

            memcpy(annot, ann_inf->annotation, size_annot*sizeof(uint8_t));
            if (annot_sup != NULL) memcpy(annot_sup, &(ann_inf->annotation[size_annot]), size_annot_sup*sizeof(uint8_t));

            pos_cell = (id_genome2insert+2)/SIZE_BITS_UINT_8T;

            if (pos_cell < size_annot) annot[pos_cell] |= MASK_POWER_8[(id_genome2insert+2)%SIZE_BITS_UINT_8T];
            else annot_sup[pos_cell-size_annot] |= MASK_POWER_8[(id_genome2insert+2-(size_annot*SIZE_BITS_UINT_8T))%SIZE_BITS_UINT_8T];
        }
        else if (ann_inf->min_mode == 1){ // Mode with min size is 1

            start_run = -1;
            stop_run = -1;

            for (z=0; z<=ann_inf->last_added; z++){

                if ((ann_inf->annotation[(z+2)/SIZE_BITS_UINT_8T] & MASK_POWER_8[(z+2)%SIZE_BITS_UINT_8T]) != 0){

                    if (start_run == -1){
                        start_run = z;
                        stop_run = z;

                        modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, start_run, -1, 1, 2);
                    }
                    else if (z == stop_run+1) stop_run++;
                }
                else if (stop_run != -1){

                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, stop_run, -1, 1, 2);

                    start_run = -1;
                    stop_run = -1;
                }
            }

            if (start_run == -1){

                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
            }
            else if (id_genome2insert == stop_run+1){
                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
            }
            else {

                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, stop_run, -1, 1, 2);

                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
            }
        }
        else { // Mode with min size is 2

            while (z < size*SIZE_BITS_UINT_8T-2){
                if ((ann_inf->annotation[(z+2)/SIZE_BITS_UINT_8T] & MASK_POWER_8[(z+2)%SIZE_BITS_UINT_8T]) != 0)
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, z, -1, 2, 1);

                z++;
            }

            modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 2, 1);
        }
    }
    else if (ann_inf->current_mode == 1){ // Current mode is 1

        if (ann_inf->min_mode == 0){ // Mode with min size is 0

            int pos_cell;
            int size_annot_bit = size_annot*SIZE_BITS_UINT_8T;

            while (i < ann_inf->nb_id_stored){

                for (z = ann_inf->id_stored[i]; z <= ann_inf->id_stored[i+1]; z++){

                    pos_cell = (z+2)/SIZE_BITS_UINT_8T;

                    if (pos_cell < size_annot) annot[pos_cell] |= MASK_POWER_8[(z+2)%SIZE_BITS_UINT_8T];
                    else annot_sup[pos_cell-size_annot] |= MASK_POWER_8[(z+2-size_annot_bit)%SIZE_BITS_UINT_8T];
                }

                i += 2;
            }

            pos_cell = (id_genome2insert+2)/SIZE_BITS_UINT_8T;

            if (pos_cell < size_annot) annot[pos_cell] |= MASK_POWER_8[(id_genome2insert+2)%SIZE_BITS_UINT_8T];
            else annot_sup[pos_cell-size_annot] |= MASK_POWER_8[(id_genome2insert+2-size_annot_bit)%SIZE_BITS_UINT_8T];
        }
        else if (ann_inf->min_mode == 1){ // Mode with min size is 1

            if (ann_inf->comp_annot <= 0){
                memcpy(annot, ann_inf->annotation, size_annot * sizeof(uint8_t));
                if (annot_sup != NULL) memcpy(annot_sup, &(ann_inf->annotation[size_annot]), size_annot_sup*sizeof(uint8_t));
            }
            else{
                for (i = 0; i < ann_inf->nb_id_stored; i++){
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, ann_inf->id_stored[i],
                                     ann_inf->size_id_stored[i], 1, 2);
                }
            }

            if (annot_sup != NULL){

                it_annot = size_annot_sup-1;
                current_annot = annot_sup;
                size_current_annot = size_annot_sup;

                while ((it_annot>=0) && ((current_annot[it_annot] & 0x3) != 1)) it_annot--;
                if ((it_annot>=0) && (current_annot[it_annot] & 0x1)) goto OUT2;
            }

            it_annot = size_annot-1;
            current_annot = annot;
            size_current_annot = size_annot;

            while ((it_annot>=0) && ((current_annot[it_annot] & 0x3) != 1)) it_annot--;

            if (it_annot == -1) ERROR("modify_annotation() annotation.c")

            OUT2:if (id_genome2insert != ann_inf->id_stored[ann_inf->nb_id_stored-1]){

                if (id_genome2insert != ann_inf->id_stored[ann_inf->nb_id_stored-1]+1){

                        it_annot++;
                        while (annot[it_annot] & 0x2) it_annot++;

                        modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
                        modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
                }
                else{
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
                }
            }
        }
        else{ // Mode with min size is 2

            while (i < ann_inf->nb_id_stored){

                for (z = ann_inf->id_stored[i]; z <= ann_inf->id_stored[i+1]; z++)
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, z, -1, 2, 1);

                i += 2;
            }

            modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 2, 1);
        }
    }
    else { // Current mode is 2

        if (ann_inf->min_mode == 0){ // Mode with min size is 0

            int pos_cell;
            int size_annot_bit = size_annot*SIZE_BITS_UINT_8T;

            for (z=0; z<ann_inf->nb_id_stored; z++){

                pos_cell = (ann_inf->id_stored[z]+2)/SIZE_BITS_UINT_8T;

                if (pos_cell < size_annot) annot[pos_cell] |= MASK_POWER_8[(ann_inf->id_stored[z]+2)%SIZE_BITS_UINT_8T];
                else annot_sup[pos_cell-size_annot] |= MASK_POWER_8[(ann_inf->id_stored[z]+2-size_annot_bit)%SIZE_BITS_UINT_8T];
            }

            pos_cell = (id_genome2insert+2)/SIZE_BITS_UINT_8T;

            if (pos_cell < size_annot) annot[pos_cell] |= MASK_POWER_8[(id_genome2insert+2)%SIZE_BITS_UINT_8T];
            else annot_sup[pos_cell-size_annot] |= MASK_POWER_8[(id_genome2insert+2-size_annot_bit)%SIZE_BITS_UINT_8T];
        }
        else if (ann_inf->min_mode == 1){ // Mode with min size is 1

            for (z=0; z<ann_inf->nb_id_stored; z++){

                if (z == 0){
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, ann_inf->id_stored[z], ann_inf->size_id_stored[z], 1, 2);
                }
                else if (ann_inf->id_stored[z] != ann_inf->id_stored[z-1]+1){
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, ann_inf->id_stored[z-1], ann_inf->size_id_stored[z-1], 1, 2);
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, ann_inf->id_stored[z], ann_inf->size_id_stored[z], 1, 2);
                }
            }

            if (ann_inf->nb_id_stored > 0){

                if (id_genome2insert == ann_inf->id_stored[ann_inf->nb_id_stored-1]+1){
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
                }
                else{
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot,
                                     ann_inf->id_stored[ann_inf->nb_id_stored-1], ann_inf->size_id_stored[ann_inf->nb_id_stored-1], 1, 2);

                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
                }
            }
            else{
                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
                modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 1, 2);
            }
        }
        else { // Mode with min size is 2

            if (ann_inf->comp_annot <= 0){
                memcpy(annot, ann_inf->annotation, size_annot * sizeof(uint8_t));
                if (annot_sup != NULL) memcpy(annot_sup, &(ann_inf->annotation[size_annot]), size_annot_sup*sizeof(uint8_t));
            }
            else{
                for (i = 0; i < ann_inf->nb_id_stored; i++){
                    modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, ann_inf->id_stored[i],
                                     ann_inf->size_id_stored[i], 2, 1);
                }

                it_annot = 0;
                size_current_annot = size_annot;
                current_annot = annot;
            }

            while ((it_annot < size_current_annot) && ((current_annot[it_annot] & 0x3) != 0)) it_annot++;
            if ((it_annot < size_current_annot) && ((current_annot[it_annot] & 0x3) == 0)) goto OUT_HERE;

            if (annot_sup != NULL){

                it_annot = 0;
                current_annot = annot_sup;
                size_current_annot = size_annot_sup;

                while ((it_annot < size_current_annot) && ((current_annot[it_annot] & 0x3) != 0)) it_annot++;
            }

            OUT_HERE: modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, id_genome2insert, size_id_genome, 2, 1);
        }
    }

    memset(ann_inf->annotation, 0, ann_inf->size_annot * sizeof(uint8_t));

    return;
}

annotation_array_elem* sort_annotations(Pvoid_t* JArray_annot, int* size_array, uint32_t longest_annot){

    PWord_t PValue_annot;
    PWord_t PValue_annot_tmp;
    Word_t * PValue_sizes;
    Word_t Rc_word;
    Word_t it_sizes;

    Pvoid_t PArray_sizes = (Pvoid_t) NULL;

    annotation_array_elem* annot_list = NULL;

    int bits_per_byte_checksum = SIZE_BITS_UINT_8T-1;

    //uint32_t tot_size;

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
    uint32_t tot_cell_index = longest_annot + CEIL(longest_annot, bits_per_byte_checksum) + 4;
    uint32_t length_index;

    //uint32_t nb_annot_with_first_size;
    uint32_t real_size = 1;
    uint32_t decomp_size = 1;

    uint64_t size_annot_array;

    uint8_t it0_first_size;
    uint8_t it1_first_size;
    uint8_t it2_first_size;

    uint8_t* it_index = malloc(tot_cell_index * sizeof(uint8_t));
    ASSERT_NULL_PTR(it_index, "sort_annotations()");

    uint8_t* first_index = NULL;
    uint8_t* it_index_tmp;
    uint8_t* it_index_start;

    memset(it_index, 0xff, tot_cell_index * sizeof(uint8_t));
    it_index[tot_cell_index-1] = '\0';

    PValue_annot = (PWord_t)JArray_annot;

    while (PValue_annot != NULL){

        JSLL(PValue_annot, *JArray_annot, it_index);

        real_size = (((it_index[0]-3) << 8) | ((it_index[1]-3) << 2) | ((it_index[2]) >> 4));
        decomp_size = 0;

        it0_first_size = it_index[0];
        it1_first_size = it_index[1];
        it2_first_size = it_index[2];

        //length_index = strlen((char*)it_index);
        length_index = ((uint8_t*)memchr(it_index, '\0', tot_cell_index)) - it_index;

        if (annot_list == NULL){

            first_index = malloc((length_index + 1) * sizeof(uint8_t));
            ASSERT_NULL_PTR(first_index, "sort_annotations()");

            size_annot_list = real_size;
            *size_array = real_size;

            annot_list = calloc(size_annot_list, sizeof(annotation_array_elem));
            ASSERT_NULL_PTR(annot_list, "sort_annotations()")

            for (tmp = 0; tmp < size_annot_list; tmp++){
                annot_list[tmp].last_index = -1;
                annot_list[tmp].annot_array = NULL;
            }
        }

        memcpy(first_index, it_index, (length_index + 1) * sizeof(uint8_t));

        //nb_annot_with_first_size = 0;

        while ((PValue_annot != NULL) && (it_index[0] == it0_first_size) && (it_index[1] == it1_first_size) && (it_index[2] == it2_first_size)){

            #if defined (_WORDx64)
                JLI(PValue_sizes, PArray_sizes, *PValue_annot);
            #else
                JLI(PValue_sizes, PArray_sizes, **((uint64_t**)PValue_annot) >> 32);
                decomp_size = MAX(decomp_size, **((uint64_t**)PValue_annot) & 0xffffffff);
            #endif

            (*PValue_sizes)++;
            //nb_annot_with_first_size++;

            JSLP(PValue_annot, *JArray_annot, it_index);
        }

        it_sizes = -1;
        tmp = 0;
        size_annot_array = 0;

        JLL(PValue_sizes, PArray_sizes, it_sizes);

        while (PValue_sizes != NULL){

            tmp = *PValue_sizes;
            next_id_tmp = next_id + tmp + 1;

            if (next_id_tmp > next_id_highest_power2)
                next_id_tmp_nb_bytes_power2 = get_nb_bytes_power2_comp(round_up_next_highest_power2(next_id_tmp));

            #if defined (_WORDx64)
                decomp_size = it_sizes & 0xffffffff;
            //    tot_size = it_sizes >> 32;
            #endif

            if (next_id_tmp_nb_bytes_power2 >= decomp_size) *PValue_sizes = UINT32_MAX;
            //else if (next_id_tmp_nb_bytes_power2 * (tot_size/decomp_size) + decomp_size > tot_size) *PValue_sizes = UINT32_MAX;
            else {
                *PValue_sizes = next_id;
                next_id += tmp;
                size_annot_array += tmp * real_size;

                if (next_id+1 > next_id_highest_power2){
                    next_id_highest_power2 = round_up_next_highest_power2(next_id+1);
                    next_id_nb_bytes_power2 = get_nb_bytes_power2_comp(next_id_highest_power2);
                }
            }

            next_id_tmp_nb_bytes_power2 = next_id_nb_bytes_power2;

            JLP(PValue_sizes, PArray_sizes, it_sizes);
        }

        //nb_annot_with_first_size = 0;

        memcpy(it_index, first_index, (length_index + 1) * sizeof(uint8_t));

        pos_annot_list = size_annot_list-real_size;
        annot_list[pos_annot_list].size_annot = real_size;
        annot_list[pos_annot_list].annot_array = calloc(size_annot_array, sizeof(uint8_t));

        JSLL(PValue_annot, *JArray_annot, it_index);

        while ((PValue_annot != NULL) && (it_index[0] == it0_first_size) && (it_index[1] == it1_first_size) && (it_index[2] == it2_first_size)){

            #if defined (_WORDx64)
                PValue_annot_tmp = PValue_annot;
                JLG(PValue_sizes, PArray_sizes, *PValue_annot_tmp);
            #else
                PValue_annot_tmp = *PValue_annot;
                JLG(PValue_sizes, PArray_sizes, *PValue_annot_tmp >> 32);

                if (*PValue_sizes == UINT32_MAX) free(PValue_annot_tmp);
            #endif

            if (*PValue_sizes == UINT32_MAX){
                JSLD(PValue_annot, *JArray_annot, it_index);
                //*PValue_annot_tmp = UINT32_MAX;
            }
            else{

                /*i = -1;
                it_index_start = &(annot_list[pos_annot_list].annot_array[((*PValue_sizes)-old_id) * real_size]);
                it_index_tmp = it_index_start-1;

                memcpy(it_index_start, &(it_index[3]), real_size * sizeof(uint8_t));

                while ((it_index_tmp = memchr(it_index_tmp + 1, 254, real_size - i - 1)) != NULL){
                    i = it_index_tmp - it_index_start;
                    if ((it_index[real_size + 3 + i/(SIZE_BITS_UINT_8T-1)] & MASK_POWER_8[i % (SIZE_BITS_UINT_8T-1) + 1]) != 0)
                        *it_index_tmp = 0;
                }*/

                i = 0;
                it_index_start = &(annot_list[pos_annot_list].annot_array[((*PValue_sizes)-old_id) * real_size]);
                it_index_tmp = it_index_start + real_size;

                memcpy(it_index_start, &it_index[3], real_size * sizeof(uint8_t));

                while (it_index_start < it_index_tmp){

                    if (*it_index_start == 254){
                        if (it_index[real_size + 3 + i/bits_per_byte_checksum] & MASK_POWER_8[i % bits_per_byte_checksum + 1]) *it_index_start = 0;
                        i++;
                    }

                    it_index_start++;
                }

                *PValue_annot_tmp = *PValue_sizes;
                *PValue_sizes += 1;
            }

            JSLP(PValue_annot, *JArray_annot, it_index);
        }

        annot_list[pos_annot_list].last_index = next_id - 1;

        JLFA(Rc_word, PArray_sizes);
        //annot_list[pos_annot_list].annot_array = realloc(annot_list[pos_annot_list].annot_array, nb_annot_with_first_size * real_size * sizeof(uint8_t));

        old_id = next_id;
    }

    printf("Number of compressed annotations = %d\n", next_id);

    for (i = 0, tmp = 0; i < size_annot_list; i++){
        if (annot_list[i].last_index == -1) annot_list[i].last_index = tmp;
        else tmp = annot_list[i].last_index;
    }

    if (next_id == 0){
        free_annotation_array_elem(&annot_list, size_array);
        *size_array = 0;
    }

    free(first_index);
    free(it_index);

    return annot_list;
}

void sort_annotations2(char* filename_annot_array_elem, Pvoid_t* JArray_annot, annotation_array_elem** root_comp_set_colors,
                       int* length_root_comp_set_colors, uint32_t longest_annot){

    ASSERT_NULL_PTR(filename_annot_array_elem, "write_annotation_array_elem()\n")

    FILE* file = fopen(filename_annot_array_elem, "wb");
    ASSERT_NULL_PTR(file, "sort_annotations2()\n")

    PWord_t PValue_annot;
    PWord_t PValue_annot_tmp;
    PWord_t PValue_sizes;
    Word_t Rc_word;
    Word_t it_sizes;

    Pvoid_t PArray_sizes = (Pvoid_t) NULL;

    annotation_array_elem* annot_list = calloc(1, sizeof(annotation_array_elem));
    ASSERT_NULL_PTR(annot_list, "sort_annotations()")

    int j;
    int shift_mode_3;

    int size_annot_list = 0;
    int bits_per_byte_checksum = SIZE_BITS_UINT_8T - 1;

    uint32_t position;

    uint32_t next_id = 0;
    uint32_t next_id_tmp = 0;
    uint32_t tmp = 0;

    uint32_t next_id_highest_power2 = 1;
    uint32_t next_id_nb_bytes_power2 = 1;
    uint32_t next_id_tmp_nb_bytes_power2 = 1;

    uint32_t prev_id = 0;

    uint32_t i = 0;
    uint32_t tot_cell_index = longest_annot + CEIL(longest_annot, bits_per_byte_checksum) + 4;
    uint32_t length_index;

    uint32_t real_size = 1;
    uint32_t decomp_size = 1;

    uint64_t size_annot_array;

    uint8_t it0_first_size;
    uint8_t it1_first_size;
    uint8_t it2_first_size;

    uint8_t* it_index = malloc(tot_cell_index * sizeof(uint8_t));
    ASSERT_NULL_PTR(it_index, "sort_annotations()");

    uint8_t* first_index = calloc(tot_cell_index, sizeof(uint8_t));
    ASSERT_NULL_PTR(first_index, "sort_annotations()");

    uint8_t* it_index_tmp;
    uint8_t* it_index_start;

    memset(it_index, 0xff, tot_cell_index * sizeof(uint8_t));
    it_index[tot_cell_index-1] = '\0';

    PValue_annot = (PWord_t)JArray_annot;

    fseek(file, sizeof(int), SEEK_SET);

    while (PValue_annot != NULL){

        JSLL(PValue_annot, *JArray_annot, it_index);

        real_size = (((it_index[0]-3) << 8) | ((it_index[1]-3) << 2) | ((it_index[2]) >> 4));

        decomp_size = 0;

        it0_first_size = it_index[0];
        it1_first_size = it_index[1];
        it2_first_size = it_index[2];

        //length_index = ((uint8_t*)memchr(it_index, '\0', tot_cell_index)) - it_index;
        length_index = strlen((char*) it_index);

        for (j = 0; j < *length_root_comp_set_colors; j++){

            if ((*root_comp_set_colors)[j].size_annot > real_size){

                if ((*root_comp_set_colors)[j].annot_array != NULL){
                    free((*root_comp_set_colors)[j].annot_array);
                    (*root_comp_set_colors)[j].annot_array = NULL;
                }
            }
            else break;
        }

        memcpy(first_index, it_index, (length_index + 1) * sizeof(uint8_t));

        while ((PValue_annot != NULL) && (it_index[0] == it0_first_size) && (it_index[1] == it1_first_size) && (it_index[2] == it2_first_size)){

            #if defined (_WORDx64)
                JLI(PValue_sizes, PArray_sizes, *PValue_annot);
            #else
                JLI(PValue_sizes, PArray_sizes, **((uint64_t**)PValue_annot) >> 32);
                decomp_size = MAX(decomp_size, **((uint64_t**)PValue_annot) & 0xffffffff);
            #endif

            (*PValue_sizes)++;

            JSLP(PValue_annot, *JArray_annot, it_index);
        }

        it_sizes = -1;
        tmp = 0;
        size_annot_array = 0;

        JLL(PValue_sizes, PArray_sizes, it_sizes);

        while (PValue_sizes != NULL){

            tmp = *PValue_sizes;
            next_id_tmp = next_id + tmp + 1;

            if (next_id_tmp > next_id_highest_power2)
                next_id_tmp_nb_bytes_power2 = get_nb_bytes_power2_comp(round_up_next_highest_power2(next_id_tmp));

            #if defined (_WORDx64)
                decomp_size = it_sizes & 0xffffffff;
            #endif

            if (next_id_tmp_nb_bytes_power2 >= decomp_size) *PValue_sizes = UINT32_MAX;
            else {
                *PValue_sizes = next_id;
                next_id += tmp;
                size_annot_array += tmp * real_size;

                if (next_id+1 > next_id_highest_power2){
                    next_id_highest_power2 = round_up_next_highest_power2(next_id+1);
                    next_id_nb_bytes_power2 = get_nb_bytes_power2_comp(next_id_highest_power2);
                }
            }

            next_id_tmp_nb_bytes_power2 = next_id_nb_bytes_power2;

            JLP(PValue_sizes, PArray_sizes, it_sizes);
        }

        memcpy(it_index, first_index, (length_index + 1) * sizeof(uint8_t));

        annot_list->annot_array = calloc(size_annot_array, sizeof(uint8_t));
        ASSERT_NULL_PTR(annot_list->annot_array, "sort_annotations2()\n")

        JSLL(PValue_annot, *JArray_annot, it_index);

        while ((PValue_annot != NULL) && (it_index[0] == it0_first_size) && (it_index[1] == it1_first_size) && (it_index[2] == it2_first_size)){

            #if defined (_WORDx64)
                PValue_annot_tmp = PValue_annot;
                JLG(PValue_sizes, PArray_sizes, *PValue_annot_tmp);
            #else
                PValue_annot_tmp = *PValue_annot;
                JLG(PValue_sizes, PArray_sizes, *PValue_annot_tmp >> 32);

                if (*PValue_sizes == UINT32_MAX) free(PValue_annot_tmp);
            #endif

            if (*PValue_sizes == UINT32_MAX){
                JSLD(PValue_annot, *JArray_annot, it_index);
            }
            else{

                it_index_start = &annot_list->annot_array[(*PValue_sizes - prev_id) * real_size];

                if ((it_index[3] & 0x3) != 3){

                    i = 0;
                    it_index_tmp = it_index_start + real_size;

                    memcpy(it_index_start, &it_index[3], real_size * sizeof(uint8_t));

                    while (it_index_start < it_index_tmp){

                        if (*it_index_start == 0xfe){
                            if (it_index[real_size + 3 + i/bits_per_byte_checksum] & MASK_POWER_8[i % bits_per_byte_checksum + 1]) *it_index_start = 0;
                            i++;
                        }

                        it_index_start++;
                    }
                }
                else{

                    j = 4;
                    shift_mode_3 = 5;
                    position = it_index[3] >> 2;

                    while (j < strlen((char*) it_index)){
                        position |= ((uint32_t)(it_index[j] & 0xfe)) << shift_mode_3;
                        shift_mode_3 += 7;
                        j++;
                    }

                    it_index_tmp = extract_from_annotation_array_elem(*root_comp_set_colors, position, &shift_mode_3);
                    memcpy(it_index_start, it_index_tmp, shift_mode_3 * sizeof(uint8_t));
                }

                *PValue_annot_tmp = *PValue_sizes;
                *PValue_sizes += 1;
            }

            JSLP(PValue_annot, *JArray_annot, it_index);
        }

        prev_id = next_id;

        JLFA(Rc_word, PArray_sizes);

        if (annot_list->annot_array != NULL){

            annot_list->last_index = next_id - 1;
            annot_list->size_annot = real_size;

            if (fwrite(&annot_list->last_index, sizeof(int64_t), 1, file) != 1) ERROR("sort_annotations2() 1")
            if (fwrite(&annot_list->size_annot, sizeof(int), 1, file) != 1) ERROR("sort_annotations2() 2")

            if (fwrite(annot_list->annot_array, sizeof(uint8_t), size_annot_array, file) != size_annot_array) ERROR("sort_annotations2() 3")

            free(annot_list->annot_array);

            size_annot_list++;
        }
    }

    rewind(file);

    if (fwrite(&size_annot_list, sizeof(int), 1, file) != 1) ERROR("sort_annotations2() 4")

    fclose(file);

    printf("Number of compressed annotations = %d\n", next_id);

    free(annot_list);

    free_annotation_array_elem(root_comp_set_colors, length_root_comp_set_colors);

    free(first_index);
    free(it_index);

    return;
}

void sort_annotations3(Pvoid_t* JArray_annot, uint32_t longest_annot){

    ASSERT_NULL_PTR(JArray_annot, "write_annotation_array_elem()\n")

    PWord_t PValue_annot;
    PWord_t PValue_annot_tmp;
    PWord_t PValue_sizes;
    Word_t Rc_word;
    Word_t it_sizes;

    Pvoid_t PArray_sizes = (Pvoid_t) NULL;

    int bits_per_byte_checksum = SIZE_BITS_UINT_8T - 1;

    uint32_t next_id = 0;
    uint32_t next_id_tmp = 0;
    uint32_t tmp = 0;

    uint32_t next_id_highest_power2 = 1;
    uint32_t next_id_nb_bytes_power2 = 1;
    uint32_t next_id_tmp_nb_bytes_power2 = 1;

    uint32_t decomp_size = 1;

    uint32_t tot_cell_index = longest_annot + CEIL(longest_annot, bits_per_byte_checksum) + 4;

    uint32_t length_index;

    uint8_t it0_first_size;
    uint8_t it1_first_size;
    uint8_t it2_first_size;

    uint8_t* it_index = malloc(tot_cell_index * sizeof(uint8_t));
    ASSERT_NULL_PTR(it_index, "sort_annotations()");

    uint8_t* first_index = calloc(tot_cell_index, sizeof(uint8_t));
    ASSERT_NULL_PTR(first_index, "sort_annotations()");

    memset(it_index, 0xff, tot_cell_index * sizeof(uint8_t));
    it_index[tot_cell_index-1] = '\0';

    PValue_annot = (PWord_t)JArray_annot;

    while (PValue_annot != NULL){

        JSLL(PValue_annot, *JArray_annot, it_index);

        decomp_size = 0;

        it0_first_size = it_index[0];
        it1_first_size = it_index[1];
        it2_first_size = it_index[2];

        length_index = strlen((char*) it_index);

        memcpy(first_index, it_index, (length_index + 1) * sizeof(uint8_t));

        while ((PValue_annot != NULL) && (it_index[0] == it0_first_size) && (it_index[1] == it1_first_size) && (it_index[2] == it2_first_size)){

            #if defined (_WORDx64)
                JLI(PValue_sizes, PArray_sizes, *PValue_annot);
            #else
                JLI(PValue_sizes, PArray_sizes, **((uint64_t**)PValue_annot) >> 32);
                decomp_size = MAX(decomp_size, **((uint64_t**)PValue_annot) & 0xffffffff);
            #endif

            (*PValue_sizes)++;

            JSLP(PValue_annot, *JArray_annot, it_index);
        }

        it_sizes = -1;
        tmp = 0;

        JLL(PValue_sizes, PArray_sizes, it_sizes);

        while (PValue_sizes != NULL){

            tmp = *PValue_sizes;
            next_id_tmp = next_id + tmp + 1;

            if (next_id_tmp > next_id_highest_power2)
                next_id_tmp_nb_bytes_power2 = get_nb_bytes_power2_comp(round_up_next_highest_power2(next_id_tmp));

            #if defined (_WORDx64)
                decomp_size = it_sizes & 0xffffffff;
            #endif

            if (next_id_tmp_nb_bytes_power2 >= decomp_size) *PValue_sizes = UINT32_MAX;
            else {

                *PValue_sizes = next_id;
                next_id += tmp;

                if (next_id+1 > next_id_highest_power2){
                    next_id_highest_power2 = round_up_next_highest_power2(next_id+1);
                    next_id_nb_bytes_power2 = get_nb_bytes_power2_comp(next_id_highest_power2);
                }
            }

            next_id_tmp_nb_bytes_power2 = next_id_nb_bytes_power2;

            JLP(PValue_sizes, PArray_sizes, it_sizes);
        }

        memcpy(it_index, first_index, (length_index + 1) * sizeof(uint8_t));

        JSLL(PValue_annot, *JArray_annot, it_index);

        while ((PValue_annot != NULL) && (it_index[0] == it0_first_size) && (it_index[1] == it1_first_size) && (it_index[2] == it2_first_size)){

            #if defined (_WORDx64)
                PValue_annot_tmp = PValue_annot;
                JLG(PValue_sizes, PArray_sizes, *PValue_annot_tmp);
            #else
                PValue_annot_tmp = *PValue_annot;
                JLG(PValue_sizes, PArray_sizes, *PValue_annot_tmp >> 32);

                if (*PValue_sizes == UINT32_MAX) free(PValue_annot_tmp);
            #endif

            if (*PValue_sizes == UINT32_MAX) JSLD(PValue_annot, *JArray_annot, it_index)
            else{
                *PValue_annot_tmp = *PValue_sizes;
                *PValue_sizes += 1;
            }

            JSLP(PValue_annot, *JArray_annot, it_index);
        }

        JLFA(Rc_word, PArray_sizes);
    }

    printf("Number of compressed annotations = %d\n", next_id);

    free(first_index);
    free(it_index);

    return;
}

void replace_annots_comp(annotation_array_elem* comp_colors, Pvoid_t* JArray_annot, char* filename_new_comp_colors, uint32_t longest_annot){

    ASSERT_NULL_PTR(filename_new_comp_colors, "replace_annots_comp()\n")

    PWord_t PValue_annot;

    int j;
    int shift_mode_3;
    int bits_per_byte_checksum = SIZE_BITS_UINT_8T - 1;

    int size_new_comp_colors;
    int pos_new_comp_colors = -1;
    int size_annot_read = INT_MAX;
    int size_annot_comp = 0;

    uint8_t it0_first_size = 0;
    uint8_t it1_first_size = 0;
    uint8_t it2_first_size = 0;

    uint32_t pos = 0;
    int64_t last_id = 0;
    int64_t prev_id = -1;

    uint32_t position;

    uint32_t tot_cell_index = longest_annot + CEIL(longest_annot, bits_per_byte_checksum) + 4;

    off_t curr_pos_file = 0;

    int file = open(filename_new_comp_colors, O_RDWR);
    if (file == -1) ERROR("replace_annots_comp()\n")

    uint8_t* index = malloc(tot_cell_index * sizeof(uint8_t));
    ASSERT_NULL_PTR(index, "replace_annots_comp()\n");

    uint8_t* index_tmp;

    memset(index, 0xff, tot_cell_index * sizeof(uint8_t));
    index[tot_cell_index-1] = '\0';

    if (read(file, &size_new_comp_colors, sizeof(int)) == -1)
        ERROR("replace_annots_comp(): Could not read the file.\n");

    JSLL(PValue_annot, *JArray_annot, index);

    while (PValue_annot != NULL){

        if ((index[0] != it0_first_size) || (index[1] != it1_first_size) || (index[2] != it2_first_size)){

            size_annot_comp = (((index[0]-3) << 8) | ((index[1]-3) << 2) | ((index[2]) >> 4));

            while (size_annot_read > size_annot_comp){

                if (pos_new_comp_colors != -1){

                    if (lseek(file, (last_id - prev_id) * size_annot_read * sizeof(uint8_t), SEEK_CUR) == -1)
                        ERROR("replace_annots_comp(): Could not seek the file.\n");

                    prev_id = last_id;
                }

                pos_new_comp_colors++;

                if (pos_new_comp_colors < size_new_comp_colors){

                    if (read(file, &last_id, sizeof(int64_t)) == -1)
                        ERROR("replace_annots_comp(): Could not read the file.\n");

                    if (read(file, &size_annot_read, sizeof(int)) == -1)
                        ERROR("replace_annots_comp(): Could not read the file.\n");
                }
                else break;
            }

            if (pos_new_comp_colors >= size_new_comp_colors) break;

            if ((curr_pos_file = lseek(file, 0, SEEK_CUR)) == -1)
                ERROR("replace_annots_comp(): Could not seek the file.\n");

            it0_first_size = index[0];
            it1_first_size = index[1];
            it2_first_size = index[2];
        }

        if ((index[3] & 0x3) == 3){

            pos = (uint32_t) *PValue_annot;

            j = 4;
            shift_mode_3 = 5;
            position = index[3] >> 2;

            while (j < strlen((char*) index)){
                position |= ((uint32_t)(index[j] & 0xfe)) << shift_mode_3;
                shift_mode_3 += 7;
                j++;
            }

            index_tmp = extract_from_annotation_array_elem(comp_colors, position, &size_annot_comp);

            if (pwrite(file, index_tmp, size_annot_comp, (pos - prev_id - 1) * size_annot_comp + curr_pos_file) == -1)
                ERROR("replace_annots_comp(): Could not write at the specific offset.\n");
        }

        JSLP(PValue_annot, *JArray_annot, index);
    }

    close(file);

    free(index);

    return;
}

void write_partial_comp_set_colors(char* filename_annot_array_elem, Pvoid_t* JArray_annot, uint32_t longest_annot){

    ASSERT_NULL_PTR(filename_annot_array_elem, "write_partial_comp_set_colors()\n")

    PWord_t PValue_annot;

    uint8_t it0_first_size = 0;
    uint8_t it1_first_size = 0;
    uint8_t it2_first_size = 0;

    int64_t last_id = 0;
    int64_t prev_id = -1;

    uint32_t pos = 0;
    uint32_t size_annot = 0;

    int i = 0;
    int size_annot_array_elem = 0;
    int bits_per_byte_checksum = SIZE_BITS_UINT_8T - 1;

    uint32_t tot_cell_index = longest_annot + CEIL(longest_annot, bits_per_byte_checksum) + 4;

    off_t curr_pos_file = 0;

    int file = open(filename_annot_array_elem, O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR);
    if (file == -1) ERROR("write_partial_comp_set_colors()\n")

    uint8_t* index = malloc(tot_cell_index * sizeof(uint8_t));
    ASSERT_NULL_PTR(index, "write_partial_comp_set_colors()");

    uint8_t* index_cpy = malloc(tot_cell_index * sizeof(uint8_t));
    ASSERT_NULL_PTR(index_cpy, "write_partial_comp_set_colors()");

    uint8_t* it_index_start;
    uint8_t* it_index_end;

    memset(index, 0xff, tot_cell_index * sizeof(uint8_t));
    index[tot_cell_index-1] = '\0';

    if (lseek(file, sizeof(int), SEEK_SET) == -1)
        ERROR("write_partial_comp_set_colors(): Could not seek the file.\n");

    JSLL(PValue_annot, *JArray_annot, index);

    while (PValue_annot != NULL){

        if ((index[0] != it0_first_size) || (index[1] != it1_first_size) || (index[2] != it2_first_size)){

            if (it0_first_size || it1_first_size || it2_first_size){

                if (write(file, &last_id, sizeof(int64_t)) == -1)
                    ERROR("write_partial_comp_set_colors(): Could not write at the specific offset.\n");

                if (write(file, &size_annot, sizeof(int)) == -1)
                    ERROR("write_partial_comp_set_colors(): Could not write at the specific offset.\n");

                if (lseek(file, (last_id - prev_id) * size_annot * sizeof(uint8_t), SEEK_CUR) == -1)
                    ERROR("write_partial_comp_set_colors(): Could not seek the file.\n");

                size_annot_array_elem++;
                prev_id = last_id;
            }

            size_annot = (((index[0]-3) << 8) | ((index[1]-3) << 2) | ((index[2]) >> 4));

            it0_first_size = index[0];
            it1_first_size = index[1];
            it2_first_size = index[2];

            if ((curr_pos_file = lseek(file, 0, SEEK_CUR)) == -1)
                ERROR("write_partial_comp_set_colors(): Could not seek the file.\n");

            curr_pos_file += sizeof(int64_t) + sizeof(int);
        }

        pos = (uint32_t) *PValue_annot;
        last_id = MAX(last_id, pos);

        if ((index[3] & 0x3) != 3){

            i = 0;
            it_index_start = index_cpy;
            it_index_end = index_cpy + size_annot;

            memcpy(index_cpy, &index[3], size_annot * sizeof(uint8_t));

            while (it_index_start < it_index_end){

                if (*it_index_start == 0xfe){
                    if (index[size_annot + 3 + i/bits_per_byte_checksum] & MASK_POWER_8[i % bits_per_byte_checksum + 1]) *it_index_start = 0x0;
                    i++;
                }

                it_index_start++;
            }
        }
        else memset(index_cpy, 0, size_annot * sizeof(uint8_t));

        if (pwrite(file, index_cpy, size_annot * sizeof(uint8_t), (pos - prev_id - 1) * size_annot + curr_pos_file) == -1)
            ERROR("write_partial_comp_set_colors(): Could not write at the specific offset.\n");

        JSLP(PValue_annot, *JArray_annot, index);
    }

    if (it0_first_size || it1_first_size || it2_first_size){

        if (write(file, &last_id, sizeof(int64_t)) == -1)
            ERROR("write_partial_comp_set_colors(): Could not write at the specific offset.\n");

        if (write(file, &size_annot, sizeof(int)) == -1)
            ERROR("write_partial_comp_set_colors(): Could not write at the specific offset.\n");

        size_annot_array_elem++;
    }

    if (pwrite(file, &size_annot_array_elem, sizeof(int), 0) == -1)
        ERROR("write_partial_comp_set_colors(): Could not write at the specific offset.\n");

    close(file);

    free(index);
    free(index_cpy);

    return;
}

int comp_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup){

    ASSERT_NULL_PTR(ann_inf, "comp_annotation()")
    ASSERT_NULL_PTR(annot, "comp_annotation()")

    if (size_annot + size_annot_sup == 0) return 0;
    if (((annot[0] & 0x3) == 0) || ((annot[0] & 0x3) == 3)) return 0; //There is nothing we can do about it;

    int i = 0;
    int it_annot = 0;

    int size = size_annot;
    int size_current_annot = size_annot;
    int final_size = 0;

    uint8_t flag1 = 1;
    uint8_t flag2 = 2;

    uint32_t nb_id_stored = 0;

    uint8_t* current_annot = annot;

    memcpy(ann_inf->annotation, annot, size_annot * sizeof(uint8_t));
    memset(annot, 0, size_annot * sizeof(uint8_t));

    if (annot_sup != NULL){
        memcpy(&(ann_inf->annotation[size_annot]), annot_sup, size_annot_sup * sizeof(uint8_t));
        memset(annot_sup, 0, size_annot_sup * sizeof(uint8_t));
        size += size_annot_sup;
    }

    ann_inf->size_annot = size;

    if (ann_inf->annotation[0] & 0x2){ //<Present everywhere from x to y> mode
        flag1 = 2;
        flag2 = 1;
    }

    while ((i<size) && (ann_inf->annotation[i] & flag1)){

        ann_inf->id_stored[nb_id_stored] = ann_inf->annotation[i] >> 2;
        i++;

        while ((i<size) && (ann_inf->annotation[i] & flag2)){
            ann_inf->id_stored[nb_id_stored] = (ann_inf->id_stored[nb_id_stored] << 6) | (ann_inf->annotation[i] >> 2);
            i++;
        }

        nb_id_stored++;
    }

    for (i = nb_id_stored - 1; i > 0; i--) ann_inf->id_stored[i] -= ann_inf->id_stored[i-1];

    for (i = 0; i < nb_id_stored; i++){
        final_size += modify_annot_bis(&current_annot, annot_sup, &it_annot, &size_current_annot, ann_inf->id_stored[i],
                                       get_nb_bytes_power2_annot(ann_inf->id_stored[i]), flag1, flag2);
    }

    reinit_annotation_inform(ann_inf);

    return final_size;
}

int decomp_annotation(annotation_inform* ann_inf, uint8_t* annot, int size_annot,
                       uint8_t* annot_sup, int size_annot_sup, int get_sizes){

    ASSERT_NULL_PTR(ann_inf, "decomp_annotation() 1")
    ASSERT_NULL_PTR(annot, "decomp_annotation() 2")

    if (size_annot + size_annot_sup == 0) return 0;
    if (((annot[0] & 0x3) == 0) || ((annot[0] & 0x3) == 3)) return 0;

    int i = 0;
    int size = size_annot;
    int size_return = 0;

    uint32_t pow2 = 0;
    uint32_t tmp = 1;
    uint32_t nb_id_stored = 0;

    uint8_t flag1 = 1;
    uint8_t flag2 = 2;

    ann_inf->comp_annot = 1;
    ann_inf->current_mode = annot[0] & 0x3;

    memmove(ann_inf->annotation, annot, size_annot * sizeof(uint8_t));

    if (annot_sup != NULL){
        memmove(&(ann_inf->annotation[size_annot]), annot_sup, size_annot_sup * sizeof(uint8_t));
        size += size_annot_sup;
    }

    ann_inf->size_annot = size;

    if (ann_inf->current_mode == 2){ //<Present everywhere from x to y> mode
        flag1 = 2;
        flag2 = 1;
    }

    while ((i<size) && (ann_inf->annotation[i] & flag1)){

        ann_inf->id_stored[nb_id_stored] = ann_inf->annotation[i] >> 2;
        i++;

        while ((i<size) && (ann_inf->annotation[i] & flag2)){
            ann_inf->id_stored[nb_id_stored] = (ann_inf->id_stored[nb_id_stored] << 6)
                                                        | (ann_inf->annotation[i] >> 2);
            i++;
        }

        nb_id_stored++;
    }

    if (get_sizes > 0){

        if (ann_inf->id_stored[0] >= pow2){
            pow2 = round_up_next_highest_power2(ann_inf->id_stored[0]);
            tmp = get_nb_bytes_power2_annot_bis(ann_inf->id_stored[0], pow2);
        }

        ann_inf->size_id_stored[0] = tmp;
        size_return = tmp;

        for (i = 1; i < nb_id_stored; i++){

            ann_inf->id_stored[i] += ann_inf->id_stored[i-1];

            if (ann_inf->id_stored[i] >= pow2){
                pow2 = round_up_next_highest_power2(ann_inf->id_stored[i]);
                tmp = get_nb_bytes_power2_annot_bis(ann_inf->id_stored[i], pow2);
            }

            ann_inf->size_id_stored[i] = tmp;
            size_return += tmp;
        }
    }
    else{
        for (i = 1; i < nb_id_stored; i++) ann_inf->id_stored[i] += ann_inf->id_stored[i-1];
    }

    ann_inf->last_added = ann_inf->id_stored[nb_id_stored - 1];
    ann_inf->nb_id_stored = nb_id_stored;

    return size_return;
}

void printAnnotation_CSV(FILE* file_output, uint8_t* annot, int size_annot, uint8_t* annot_sup,
                         int size_annot_sup, uint32_t id_genome_max, annotation_array_elem* annot_sorted){

    const char zero_comma[2] = "0,";
    const char one_comma[2] = "1,";
    const char zero_end[2] = "0\n";
    const char one_end[2] = "1\n";

    uint8_t flag;

    int i = 1;
    int size_annot_flag0 = MAX(CEIL(id_genome_max+3, SIZE_BITS_UINT_8T), 1);
    int size = size_annot+size_annot_sup;

    UC_SIZE_ANNOT_T start = 0;
    UC_SIZE_ANNOT_T stop = 0;
    UC_SIZE_ANNOT_T z;

    uint8_t* annot_real = annot;

    uint8_t* annot_flag0 = calloc(size_annot_flag0, sizeof(uint8_t));
    ASSERT_NULL_PTR(annot_flag0, "printAnnotation_CSV()")

    static const char* presence_genomes[16] = {"0,0,0,0,", "1,0,0,0,", "0,1,0,0,", "1,1,0,0,",
                                        "0,0,1,0,", "1,0,1,0,", "0,1,1,0,", "1,1,1,0,",
                                        "0,0,0,1,", "1,0,0,1,", "0,1,0,1,", "1,1,0,1,",
                                        "0,0,1,1,", "1,0,1,1,", "0,1,1,1,", "1,1,1,1,"};

    if (annot != NULL){

        if ((annot[0] & 0x3) == 3){

            uint32_t position = annot[0] >> 2;

            while ((i<size_annot) && (IS_ODD(annot[i]))){
                position |= (annot[i] >> 1) << (6+(i-1)*7);
                i++;
            }

            if ((i >= size_annot) && (annot_sup != NULL)){
                i = 0;
                while ((i<size_annot_sup) && (IS_ODD(annot_sup[i]))){
                    position |= (annot_sup[i] >> 1) << (6+(i+size_annot-1)*7);
                    i++;
                }
            }

            annot_real = extract_from_annotation_array_elem(annot_sorted, position, &size);
        }
        else if (annot_sup != NULL){
            annot_real = malloc(size*sizeof(uint8_t));
            memcpy(annot_real, annot, size_annot*sizeof(uint8_t));
            annot_real[size_annot] = annot_sup[0];
        }
    }

    i = 0;

    if (annot_real != NULL){

        flag = annot_real[0] & 0x3;

        if (flag == 0) memcpy(annot_flag0, annot_real, size*sizeof(uint8_t));
        else if (flag == 1){

            flag = annot[0] & 0x3;

            while (i<size){

                if (annot_real[i] & 0x1){
                    start = annot_real[i] >> 2;
                    i++;

                    while ((i<size) && (annot_real[i] & 0x2)){
                        start = (start << 6) | (annot_real[i] >> 2);
                        i++;
                    }
                }

                if (flag == 3) start += stop;

                if (annot_real[i] & 0x1){
                    stop = annot_real[i] >> 2;
                    i++;

                    while ((i<size) && (annot_real[i] & 0x2)){
                        stop = (stop << 6) | (annot_real[i] >> 2);
                        i++;
                    }
                }

                if (flag == 3) stop += start;

                for (z=start; z<=stop; z++) annot_flag0[(z+2)/SIZE_BITS_UINT_8T] |= MASK_POWER_8[(z+2)%SIZE_BITS_UINT_8T];
            }
        }
        else if (flag == 2){

            flag = annot[0] & 0x3;

            while ((i<size) && (annot_real[i] & 0x2)){
                z = annot_real[i] >> 2;
                i++;

                while ((i<size) && (annot_real[i] & 0x1)){
                    z = (z << 6) | (annot_real[i] >> 2);
                    i++;
                }

                if (flag == 3){
                    z += start;
                    start = z;
                }

                annot_flag0[(z+2)/SIZE_BITS_UINT_8T] |= MASK_POWER_8[(z+2)%SIZE_BITS_UINT_8T];
            }
        }
    }

    if (id_genome_max <= 5){

        for (i = 2; i <= id_genome_max + 1; i++){
            if ((annot_flag0[i/SIZE_BITS_UINT_8T] & MASK_POWER_8[i%SIZE_BITS_UINT_8T]) != 0){
                fwrite(&one_comma, sizeof(char), 2, file_output);
            }
            else fwrite(&zero_comma, sizeof(char), 2, file_output);
        }
    }
    else{

        for (i = 2; i < SIZE_BITS_UINT_8T; i++){
            if ((annot_flag0[i/SIZE_BITS_UINT_8T] & MASK_POWER_8[i%SIZE_BITS_UINT_8T]) != 0){
                fwrite(&one_comma, sizeof(char), 2, file_output);
            }
            else fwrite(&zero_comma, sizeof(char), 2, file_output);
        }

        for (i=1; i<size_annot_flag0-1; i++){
            fwrite(presence_genomes[annot_flag0[i] & 0xf], sizeof(char), 8, file_output);
            fwrite(presence_genomes[annot_flag0[i] >> 4], sizeof(char), 8, file_output);
        }

        for (i = (size_annot_flag0 - 1) * SIZE_BITS_UINT_8T; i <= id_genome_max + 1; i++){
            if ((annot_flag0[i/SIZE_BITS_UINT_8T] & MASK_POWER_8[i%SIZE_BITS_UINT_8T]) != 0){
                fwrite(&one_comma, sizeof(char), 2, file_output);
            }
            else fwrite(&zero_comma, sizeof(char), 2, file_output);
        }
    }

    if ((annot_flag0[i/SIZE_BITS_UINT_8T] & MASK_POWER_8[i%SIZE_BITS_UINT_8T]) != 0){
        fwrite(&one_end, sizeof(char), 2, file_output);
    }
    else fwrite(&zero_end, sizeof(char), 2, file_output);

    if (((annot[0] & 0x3) != 3) && (annot_sup != NULL)) free(annot_real);

    free(annot_flag0);

    return;
}

void get_id_genomes_from_annot(annotation_inform* ann_inf, annotation_array_elem* annot_sorted, uint8_t* annot, int size_annot,
                               uint8_t* annot_sup, int size_annot_sup){

    ASSERT_NULL_PTR(ann_inf, "get_id_genomes_from_annot()\n")

    int size;

    int i = 1;

    if (size_annot != 0){ //If the annotation is not new

        if ((annot[0] & 0x3) == 3){

            uint32_t position = annot[0] >> 2;

            while ((i<size_annot) && (IS_ODD(annot[i]))){
                position |= ((uint32_t)(annot[i] >> 1)) << (6+(i-1)*7);
                i++;
            }

            if ((i >= size_annot) && (annot_sup != NULL)){
                i = 0;
                while ((i<size_annot_sup) && (IS_ODD(annot_sup[i]))){
                    position |= ((uint32_t)(annot_sup[i] >> 1)) << (6+(i+size_annot-1)*7);
                    i++;
                }
            }

            uint8_t* annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &size);
            memcpy(ann_inf->annotation, annot_tmp, size*sizeof(uint8_t));

            i = decomp_annotation(ann_inf, annot_tmp, size, NULL, 0, 1);
            if (i) size = i;
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

        ann_inf->size_annot = size;
        ann_inf->current_mode = ann_inf->annotation[0] & 0x3;

        if (ann_inf->current_mode == 0){ //Bitwize mode

            for (i=2; i<size*SIZE_BITS_UINT_8T; i++){

                if (ann_inf->annotation[i/SIZE_BITS_UINT_8T] & MASK_POWER_8[i%SIZE_BITS_UINT_8T]){

                    ann_inf->id_stored[ann_inf->nb_id_stored] = i-2;
                    ann_inf->nb_id_stored++;
                }
            }
        }
        else if (ann_inf->current_mode == 1){ //<Present everywhere from x to y> mode

            if (ann_inf->comp_annot <= 0){

                bool it = false;

                while ((i<size) && (ann_inf->annotation[i] & 0x1)){

                    ann_inf->id_stored[ann_inf->nb_id_stored] = ann_inf->annotation[i] >> 2;
                    i++;

                    while ((i<size) && (ann_inf->annotation[i] & 0x2)){

                        ann_inf->id_stored[ann_inf->nb_id_stored] = (ann_inf->id_stored[ann_inf->nb_id_stored] << 6) | (ann_inf->annotation[i] >> 2);
                        i++;
                    }

                    if (it){

                        uint32_t j = ann_inf->id_stored[ann_inf->nb_id_stored-1];
                        uint32_t stop = ann_inf->id_stored[ann_inf->nb_id_stored];

                        if (j != stop){
                            for (j++; j <= stop; j++){
                                ann_inf->id_stored[ann_inf->nb_id_stored] = j;
                                ann_inf->nb_id_stored++;
                            }
                        }
                    }
                    else ann_inf->nb_id_stored++;

                    it = !it;
                }
            }
            else {

                uint32_t pow2 = 0;
                uint32_t tmp = 1;

                int nb_id_stored_cpy = 0;

                uint32_t* id_stored_cpy = malloc(ann_inf->nb_id_stored * sizeof(uint32_t));
                ASSERT_NULL_PTR(id_stored_cpy, "get_id_genomes_from_annot()\n")

                memcpy(id_stored_cpy, ann_inf->id_stored, ann_inf->nb_id_stored * sizeof(uint32_t));

                for (i = 0; i < ann_inf->nb_id_stored; i++){

                    if (i%2){

                        for (uint32_t z = id_stored_cpy[i-1]+1; z <= id_stored_cpy[i]; z++){

                            ann_inf->id_stored[nb_id_stored_cpy] = z;

                            if (ann_inf->id_stored[nb_id_stored_cpy] >= pow2){
                                pow2 = round_up_next_highest_power2(ann_inf->id_stored[nb_id_stored_cpy]);
                                tmp = get_nb_bytes_power2_annot_bis(ann_inf->id_stored[nb_id_stored_cpy], pow2);
                            }

                            ann_inf->size_id_stored[nb_id_stored_cpy] = tmp;

                            nb_id_stored_cpy++;
                        }
                    }
                    else{
                        ann_inf->id_stored[nb_id_stored_cpy] = id_stored_cpy[i];

                        if (ann_inf->id_stored[nb_id_stored_cpy] >= pow2){
                            pow2 = round_up_next_highest_power2(ann_inf->id_stored[nb_id_stored_cpy]);
                            tmp = get_nb_bytes_power2_annot_bis(ann_inf->id_stored[nb_id_stored_cpy], pow2);
                        }

                        ann_inf->size_id_stored[nb_id_stored_cpy] = tmp;

                        nb_id_stored_cpy++;
                    }
                }

                ann_inf->nb_id_stored = nb_id_stored_cpy;

                free(id_stored_cpy);
            }
        }
        else if (ann_inf->current_mode == 2){ //<Present nowhere except in x> mode

            if (ann_inf->comp_annot <= 0){

                while ((i<size) && (ann_inf->annotation[i] & 0x2)){

                    ann_inf->id_stored[ann_inf->nb_id_stored] = ann_inf->annotation[i] >> 2;
                    i++;

                    while ((i<size) && (ann_inf->annotation[i] & 0x1)){
                        ann_inf->id_stored[ann_inf->nb_id_stored] = (ann_inf->id_stored[ann_inf->nb_id_stored] << 6) | (ann_inf->annotation[i] >> 2);
                        i++;
                    }

                    ann_inf->nb_id_stored++;
                }
            }
        }
        else ERROR( "get_id_genomes_from_annot(): mode 3, should not happen.\n" )
    }

    return;
}

int get_count_id_genomes_from_annot(annotation_inform* ann_inf, annotation_array_elem* annot_sorted,
                                    uint8_t* annot, int size_annot, uint8_t* annot_sup, int size_annot_sup){

    ASSERT_NULL_PTR(ann_inf,"get_count_id_genomes_from_annot()\n")

    int size;

    int i = 1;
    int count_ids = 0;

    if (size_annot != 0){ //If the annotation is not new

        if ((annot[0] & 0x3) == 3){

            uint32_t position = annot[0] >> 2;

            while ((i<size_annot) && (IS_ODD(annot[i]))){
                position |= ((uint32_t)(annot[i] >> 1)) << (6+(i-1)*7);
                i++;
            }

            if ((i >= size_annot) && (annot_sup != NULL)){
                i = 0;
                while ((i<size_annot_sup) && (IS_ODD(annot_sup[i]))){
                    position |= ((uint32_t)(annot_sup[i] >> 1)) << (6+(i+size_annot-1)*7);
                    i++;
                }
            }

            uint8_t* annot_tmp = extract_from_annotation_array_elem(annot_sorted, position, &size);
            memcpy(ann_inf->annotation, annot_tmp, size*sizeof(uint8_t));

            i = decomp_annotation(ann_inf, annot_tmp, size, NULL, 0, 1);
            if (i != 0) size = i;
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

        ann_inf->size_annot = size;
        ann_inf->current_mode = ann_inf->annotation[0] & 0x3;

        if (ann_inf->current_mode == 0){ //Bitwize mode

            for (i=2; i<size*SIZE_BITS_UINT_8T; i++)
                if (ann_inf->annotation[i/SIZE_BITS_UINT_8T] & MASK_POWER_8[i%SIZE_BITS_UINT_8T]) count_ids++;
        }
        else if (ann_inf->current_mode == 1){ //<Present everywhere from x to y> mode

            if (ann_inf->comp_annot <= 0){

                bool it = false;

                while ((i<size) && (ann_inf->annotation[i] & 0x1)){

                    ann_inf->id_stored[ann_inf->nb_id_stored] = ann_inf->annotation[i] >> 2;
                    i++;

                    while ((i<size) && (ann_inf->annotation[i] & 0x2)){

                        ann_inf->id_stored[ann_inf->nb_id_stored] = (ann_inf->id_stored[ann_inf->nb_id_stored] << 6) | (ann_inf->annotation[i] >> 2);
                        i++;
                    }

                    if (it){
                        count_ids += ann_inf->id_stored[ann_inf->nb_id_stored] - ann_inf->id_stored[ann_inf->nb_id_stored-1] + 1;
                        ann_inf->nb_id_stored = 0;
                    }
                    else ann_inf->nb_id_stored++;

                    it = !it;
                }
            }
            else {
                for (i = 0; i < ann_inf->nb_id_stored; i+=2) count_ids += ann_inf->id_stored[i+1] - ann_inf->id_stored[i] + 1;
            }
        }
        else if (ann_inf->current_mode == 2){ //<Present nowhere except in x> mode

            if (ann_inf->comp_annot <= 0){

                while ((i<size) && (ann_inf->annotation[i] & 0x2)){

                    count_ids++;
                    i++;

                    while ((i<size) && (ann_inf->annotation[i] & 0x1)) i++;
                }
            }
            else count_ids = ann_inf->nb_id_stored;
        }
        else ERROR( "get_count_id_genomes_from_annot(): mode 3, should not happen.\n" )
    }

    reinit_annotation_inform(ann_inf);

    return count_ids;
}

annotation_array_elem* cmp_annots(uint8_t* annot1, int size_annot1, uint8_t* annot_sup1, int size_annot_sup1, uint8_t* annot2, int size_annot2,
                                  uint8_t* annot_sup2, int size_annot_sup2, uint32_t id_genome_max, uint8_t (*f)(const uint8_t, const uint8_t),
                                  annotation_array_elem* annot_sorted){

    int flag1, flag2;

    int i = 1;
    int size_annot_flag0 = MAX(CEIL(id_genome_max+3, SIZE_BITS_UINT_8T), 1);
    int size1 = size_annot1 + size_annot_sup1;
    int size2 = size_annot2 + size_annot_sup2;

    uint8_t* annot1_flag0 = calloc(size_annot_flag0, sizeof(uint8_t));
    uint8_t* annot2_flag0 = calloc(size_annot_flag0, sizeof(uint8_t));

    uint8_t* annot1_real = annot1;
    uint8_t* annot2_real = annot2;

    if (annot1 != NULL){

        if ((annot1[0] & 0x3) == 3){

            uint32_t position = annot1[0] >> 2;

            while ((i<size_annot1) && (IS_ODD(annot1[i]))){
                position |= (annot1[i] >> 1) << (6+(i-1)*7);
                i++;
            }

            if ((i >= size_annot1) && (annot_sup1 != NULL)){
                i = 0;
                while ((i<size_annot_sup1) && (IS_ODD(annot_sup1[i]))){
                    position |= (annot_sup1[i] >> 1) << (6+(i+size_annot1-1)*7);
                    i++;
                }
            }

            annot1_real = extract_from_annotation_array_elem(annot_sorted, position, &size1);
        }
        else if (annot_sup1 != NULL){
            annot1_real = malloc(size1 * sizeof(uint8_t));
            memcpy(annot1_real, annot1, size_annot1 * sizeof(uint8_t));
            annot1_real[size_annot1] = annot_sup1[0];
        }
    }

    i = 1;

    if (annot2 != NULL){

        if ((annot2[0] & 0x3) == 3){

            uint32_t position = annot2[0] >> 2;

            while ((i<size_annot2) && (IS_ODD(annot2[i]))){
                position |= (annot2[i] >> 1) << (6+(i-1)*7);
                i++;
            }

            if ((i >= size_annot2) && (annot_sup2 != NULL)){
                i = 0;
                while ((i<size_annot_sup2) && (IS_ODD(annot_sup2[i]))){
                    position |= (annot_sup2[i] >> 1) << (6+(i+size_annot2-1)*7);
                    i++;
                }
            }

            annot2_real = extract_from_annotation_array_elem(annot_sorted, position, &size2);
        }
        else if (annot_sup2 != NULL){
            annot2_real = malloc(size2 * sizeof(uint8_t));
            memcpy(annot2_real, annot2, size_annot2 * sizeof(uint8_t));
            annot2_real[size_annot2] = annot_sup2[0];
        }
    }

    i = 0;

    if (annot1 != NULL){

        flag1 = annot1_real[0] & 0x3;

        if (flag1 == 0){
            memcpy(annot1_flag0, annot1_real, size_annot1 * sizeof(uint8_t));
            memcpy(&(annot1_flag0[size_annot1]), annot_sup1, size_annot_sup1 * sizeof(uint8_t));
        }
        else if (flag1 == 1){

            UC_SIZE_ANNOT_T start, stop, z;

            while ((i < size1) && (annot1_real[i] & 0x1)){
                start = annot1_real[i] >> 2;
                i++;

                while ((i<size1) && (annot1_real[i] & 0x2)){
                    start = (start << 6) | (annot1_real[i] >> 2);
                    i++;
                }

                stop = annot1_real[i] >> 2;
                i++;

                while ((i < size1) && (annot1_real[i] & 0x2)){
                    stop = (stop << 6) | (annot1_real[i] >> 2);
                    i++;
                }

                for (z=start; z<=stop; z++)
                    annot1_flag0[(z+2)/SIZE_BITS_UINT_8T] |= MASK_POWER_8[(z+2)%SIZE_BITS_UINT_8T];
            }
        }
        else if (flag1 == 2){

            UC_SIZE_ANNOT_T pos;

            while ((i<size1) && (annot1[i] & 0x2)){
                pos = annot1[i] >> 2;
                i++;

                while ((i<size1) && (annot1[i] & 0x1)){
                    pos = (pos << 6) | (annot1[i] >> 2);
                    i++;
                }

                annot1_flag0[(pos+2)/SIZE_BITS_UINT_8T] |= MASK_POWER_8[(pos+2)%SIZE_BITS_UINT_8T];
            }
        }
    }

    i = 0;

    if (annot2 != NULL){

        flag2 = annot2_real[0] & 0x3;

        if (flag2 == 0){
            memcpy(annot2_flag0, annot2_real, size_annot2 * sizeof(uint8_t));
            memcpy(&(annot2_flag0[size_annot2]), annot_sup2, size_annot_sup2 * sizeof(uint8_t));
        }
        else if (flag2 == 1){

            UC_SIZE_ANNOT_T start, stop, z;

            while ((i < size2) && (annot2_real[i] & 0x1)){
                start = annot2_real[i] >> 2;
                i++;

                while ((i < size2) && (annot2_real[i] & 0x2)){
                    start = (start << 6) | (annot2_real[i] >> 2);
                    i++;
                }

                stop = annot2_real[i] >> 2;
                i++;

                while ((i < size2) && (annot2_real[i] & 0x2)){
                    stop = (stop << 6) | (annot2_real[i] >> 2);
                    i++;
                }

                for (z=start; z<=stop; z++)
                    annot2_flag0[(z+2)/SIZE_BITS_UINT_8T] |= MASK_POWER_8[(z+2)%SIZE_BITS_UINT_8T];
            }
        }
        else if (flag2 == 2){

            UC_SIZE_ANNOT_T pos;

            while ((i < size2) && (annot2_real[i] & 0x2)){
                pos = annot2_real[i] >> 2;
                i++;

                while ((i < size2) && (annot2_real[i] & 0x1)){
                    pos = (pos << 6) | (annot2_real[i] >> 2);
                    i++;
                }

                annot2_flag0[(pos+2)/SIZE_BITS_UINT_8T] |= MASK_POWER_8[(pos+2)%SIZE_BITS_UINT_8T];
            }
        }
    }

    for (i=0; i < size_annot_flag0; i++) annot1_flag0[i] = (*f)(annot1_flag0[i], annot2_flag0[i]);

    free(annot2_flag0);
    if ((annot1 != NULL) && ((annot1[0] & 0x3) != 3) && (annot_sup1 != NULL)) free(annot1_real);
    if ((annot2 != NULL) && ((annot2[0] & 0x3) != 3) && (annot_sup2 != NULL)) free(annot2_real);

    annotation_array_elem* ann_arr_elem = malloc(sizeof(annotation_array_elem));
    ASSERT_NULL_PTR(ann_arr_elem, "intersection_annotations()");

    ann_arr_elem->annot_array = annot1_flag0;
    ann_arr_elem->size_annot = size_annot_flag0;

    return ann_arr_elem;
}
