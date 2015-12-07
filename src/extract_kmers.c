#include "./../lib/extract_kmers.h"

int extract_kmers_from_node(Node* n, int lvl_node, uint8_t* kmer, int size_kmer, int bucket, int pos_in_bucket,
                            int size_kmer_root, int compression, info_per_level* restrict info_per_lvl, FILE* file_output){

    ASSERT_NULL_PTR(n,"extract_kmers_from_node()")
    ASSERT_NULL_PTR(kmer,"extract_kmers_from_node()")
    ASSERT_NULL_PTR(info_per_lvl,"extract_kmers_from_node()")
    ASSERT_NULL_PTR(file_output,"extract_kmers_from_node()")

    CC* cc;
    UC* uc;

    int i = -1, j = 0, k = 0, count = 0, size_kmer_array = CEIL(size_kmer_root*2,SIZE_BITS_UINT_8T);

    uint16_t size_bf, nb_elt, it_filter2;

    uint64_t new_substring;

    int it_filter3, first_bit, it_bucket, last_shift, last_it_children_bucket, nb_cell_children, shifting_UC;
    int it_children_pos_bucket, it_children_bucket, it_node, it_substring, size_line, size_new_substring, size_new_substring_bytes;

    int shifting1 = (NB_CHAR_SUF_PREF*2)-SIZE_BITS_UINT_8T+pos_in_bucket;
    int shifting2 = shifting1-SIZE_BITS_UINT_8T;
    int shifting3 = SIZE_BITS_UINT_8T-shifting2;

    uint8_t mask = ~(MASK_POWER_8[shifting3]-1);
    uint8_t mask_end_kmer = MASK_POWER_16[SIZE_BITS_UINT_8T - (size_kmer_array*SIZE_BITS_UINT_8T - size_kmer_root*2)]-1;

    uint8_t s, p;
    uint8_t type;

    uint8_t kmer_tmp[size_kmer_array];
    uint8_t kmer_start[size_kmer_array];

    if ((CC*)n->CC_array != NULL){
        do {
            i++;
            cc = &(((CC*)n->CC_array)[i]);

            s = (cc->type >> 1) & 0x1f;
            p = NB_CHAR_SUF_PREF*2-s;

            it_filter2 = 0;
            it_filter3 = 0;
            it_bucket = 0;
            it_node = 0;
            it_children_pos_bucket = 0;
            it_children_bucket = 0;
            it_substring = 0;

            size_bf = cc->type >> 7;
            type = (cc->type >> 6) & 0x1;
            last_it_children_bucket = -1;

            if (s==8){
                if (info_per_lvl[lvl_node].level_min == 1){
                    for (it_filter2=0; it_filter2<MASK_POWER_16[p]; it_filter2++){
                        if ((cc->BF_filter2[size_bf+it_filter2/SIZE_BITS_UINT_8T] & (MASK_POWER_8[it_filter2%SIZE_BITS_UINT_8T])) != 0){

                            first_bit = 1;

                            while((it_filter3 < cc->nb_elem) && (((cc->extra_filter3[it_filter3/SIZE_BITS_UINT_8T] & MASK_POWER_8[it_filter3%SIZE_BITS_UINT_8T]) == 0) || (first_bit == 1))){

                                memcpy(kmer_tmp, kmer, size_kmer_array*sizeof(uint8_t));

                                new_substring = (it_filter2 << SIZE_BITS_UINT_8T) | cc->filter3[it_filter3];
                                if (compression <= 0) new_substring = (new_substring >> 2) | ((new_substring & 0x3) << 16);

                                kmer_tmp[bucket] |= new_substring >> shifting1;
                                kmer_tmp[bucket+1] = new_substring >> shifting2;
                                kmer_tmp[bucket+2] = new_substring << shifting3;
                                it_bucket = bucket+2;
                                if (shifting3 == 0) it_bucket++;

                                if (size_kmer != NB_CHAR_SUF_PREF){

                                    if ((nb_elt = getNbElts(cc, it_filter3, type)) != 0){
                                        it_children_bucket = it_filter3/info_per_lvl[lvl_node].nb_ucs_skp;

                                        if (it_children_bucket != last_it_children_bucket){

                                            it_children_pos_bucket = 0;
                                            uc = &(((UC*)cc->children)[it_children_bucket]);
                                            size_line = info_per_lvl[lvl_node].size_kmer_in_bytes_minus_1 + uc->size_annot;
                                            last_it_children_bucket = it_children_bucket;
                                        }

                                        for (j=it_children_pos_bucket*size_line; j<(it_children_pos_bucket+nb_elt)*size_line; j+=size_line){

                                            extractSuffix(kmer_tmp, size_kmer, size_kmer_array, shifting3, it_bucket,
                                                          &(uc->suffixes[j]), info_per_lvl[lvl_node].size_kmer_in_bytes_minus_1, 0);
                                            memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));

                                            //kmer_tmp is the local copy of the current k-mer, it must not be modified
                                            //You can do stuff here with kmer_start (even modify it)
                                            fwrite(kmer_start, sizeof(uint8_t), size_kmer_array, file_output);

                                            for (k=0; k<bucket+3; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);
                                            if (shifting3 != 0) kmer_tmp[bucket+2] &= mask;
                                            memset(&(kmer_tmp[bucket+3]), 0, (size_kmer_array-bucket-3)*sizeof(uint8_t));
                                        }

                                        it_children_pos_bucket += nb_elt;
                                        count += nb_elt;
                                    }
                                    else{
                                        if (shifting3 == 0){
                                            count += extract_kmers_from_node(&(cc->children_Node_container[it_node]), lvl_node-1, kmer_tmp,
                                                                             size_kmer-NB_CHAR_SUF_PREF, it_bucket, 0, size_kmer_root, compression,
                                                                             info_per_lvl, file_output);
                                        }
                                        else{
                                            count += extract_kmers_from_node(&(cc->children_Node_container[it_node]), lvl_node-1, kmer_tmp,
                                                                             size_kmer-NB_CHAR_SUF_PREF, it_bucket, shifting2, size_kmer_root,
                                                                             compression, info_per_lvl, file_output);
                                        }

                                        it_node++;
                                    }
                                }
                                else{
                                    count++;
                                    memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));
                                    for (k=0; k<size_kmer_array; k++) kmer_start[k] = reverse_word_8(kmer_start[k]);

                                    //kmer_tmp is the local copy of the current k-mer, it must not be modified
                                    //You can do stuff here with kmer_start (even modify it)
                                    fwrite(kmer_start, sizeof(uint8_t), size_kmer_array, file_output);
                                }

                                it_filter3++;
                                first_bit=0;
                            }
                        }
                    }
                }
                else{
                    int it = 0;
                    int end = 0;
                    int cpt_pv = 0;
                    int nb_skp = CEIL(cc->nb_elem, info_per_lvl[lvl_node].nb_ucs_skp);

                    nb_cell_children = info_per_lvl[lvl_node].size_kmer_in_bytes_minus_1-1;

                    for (it_filter2=0; it_filter2 < MASK_POWER_16[p]; it_filter2++){

                        if ((cc->BF_filter2[size_bf+it_filter2/SIZE_BITS_UINT_8T] & (MASK_POWER_8[it_filter2%SIZE_BITS_UINT_8T])) != 0){

                            first_bit = 1;

                            while (it_children_bucket < nb_skp){

                                if (it_children_bucket == nb_skp - 1) end = cc->nb_elem - it_children_bucket * info_per_lvl[lvl_node].nb_ucs_skp;
                                else end = info_per_lvl[lvl_node].nb_ucs_skp;

                                uc = &(((UC*)cc->children)[it_children_bucket]);
                                size_line = info_per_lvl[lvl_node].size_kmer_in_bytes_minus_1 + uc->size_annot;

                                while (it < end){

                                    memcpy(kmer_tmp, kmer, size_kmer_array*sizeof(uint8_t));

                                    new_substring = (it_filter2 << SIZE_BITS_UINT_8T) | cc->filter3[it_filter3];
                                    if (compression <= 0) new_substring = (new_substring >> 2) | ((new_substring & 0x3) << 16);
                                    kmer_tmp[bucket] |= new_substring >> shifting1;
                                    kmer_tmp[bucket+1] = new_substring >> shifting2;
                                    kmer_tmp[bucket+2] = new_substring << shifting3;
                                    it_bucket = bucket+2;
                                    if (shifting3 == 0) it_bucket++;

                                    if ((nb_elt = getNbElts(cc, it_filter3, type)) == 0){

                                        if (((cc->children_Node_container[it_node].UC_array.nb_children & 0x1) == 0) || (first_bit == 1)){

                                            first_bit=0;

                                            if (shifting3 == 0){
                                                count += extract_kmers_from_node(&(cc->children_Node_container[it_node]), lvl_node-1, kmer_tmp,
                                                                                 size_kmer-NB_CHAR_SUF_PREF, it_bucket, 0, size_kmer_root,
                                                                                 compression, info_per_lvl, file_output);
                                            }
                                            else{
                                                count += extract_kmers_from_node(&(cc->children_Node_container[it_node]), lvl_node-1, kmer_tmp,
                                                                                 size_kmer-NB_CHAR_SUF_PREF, it_bucket, shifting2, size_kmer_root,
                                                                                 compression, info_per_lvl, file_output);
                                            }

                                            it_node++;
                                        }
                                        else goto OUT_LOOP_S8;
                                    }
                                    else{
                                        if ((uc->suffixes[cpt_pv*size_line+nb_cell_children] < 0x80)  || (first_bit == 1)){

                                            first_bit=0;

                                            for (j=cpt_pv*size_line; j<(cpt_pv+nb_elt)*size_line; j+=size_line){

                                                extractSuffix(kmer_tmp, size_kmer, size_kmer_array, shifting3, it_bucket,
                                                              &(uc->suffixes[j]), info_per_lvl[lvl_node].size_kmer_in_bytes_minus_1, 0);

                                                kmer_tmp[size_kmer_array-1] &= mask_end_kmer;
                                                memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));

                                                //kmer_tmp is the local copy of the current k-mer, it must not be modified
                                                //You can do stuff here with kmer_start (even modify it)
                                                fwrite(kmer_start, sizeof(uint8_t), size_kmer_array, file_output);

                                                for (k=0; k<bucket+3; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);
                                                if (shifting3 != 0) kmer_tmp[bucket+2] &= mask;
                                                memset(&(kmer_tmp[bucket+3]), 0, (size_kmer_array-bucket-3)*sizeof(uint8_t));
                                            }

                                            cpt_pv += nb_elt;
                                            count += nb_elt;
                                        }
                                        else goto OUT_LOOP_S8;
                                    }

                                    it++;
                                    it_filter3++;
                                }

                                it = 0;
                                cpt_pv = 0;
                                it_children_bucket++;
                            }
                        }

                        OUT_LOOP_S8: continue;

                    }
                }
            }
            else {

                if (info_per_lvl[lvl_node].level_min == 1){

                    for (it_filter2=0; it_filter2<MASK_POWER_16[p]; it_filter2++){

                        if ((cc->BF_filter2[size_bf+it_filter2/SIZE_BITS_UINT_8T] & (MASK_POWER_8[it_filter2%SIZE_BITS_UINT_8T])) != 0){

                            first_bit = 1;

                            while((it_filter3 < cc->nb_elem) && (((cc->extra_filter3[it_filter3/SIZE_BITS_UINT_8T] & MASK_POWER_8[it_filter3%SIZE_BITS_UINT_8T]) == 0) || (first_bit == 1))){

                                memcpy(kmer_tmp, kmer, size_kmer_array*sizeof(uint8_t));

                                if (IS_ODD(it_filter3)) new_substring = (it_filter2 << 4) | (cc->filter3[it_filter3/2] >> 4);
                                else new_substring = (it_filter2 << 4) | (cc->filter3[it_filter3/2] & 0xf);

                                if (compression <= 0) new_substring = (new_substring >> 2) | ((new_substring & 0x3) << 16);

                                kmer_tmp[bucket] |= new_substring >> shifting1;
                                kmer_tmp[bucket+1] = new_substring >> shifting2;
                                kmer_tmp[bucket+2] = new_substring << shifting3;
                                it_bucket = bucket+2;
                                if (shifting3 == 0) it_bucket++;

                                if (size_kmer != NB_CHAR_SUF_PREF){

                                    if ((nb_elt = getNbElts(cc, it_filter3, type)) != 0){
                                        it_children_bucket = it_filter3/info_per_lvl[lvl_node].nb_ucs_skp;

                                        if (it_children_bucket != last_it_children_bucket){

                                            it_children_pos_bucket = 0;
                                            uc = &(((UC*)cc->children)[it_children_bucket]);

                                            size_line = info_per_lvl[lvl_node].size_kmer_in_bytes_minus_1 + uc->size_annot;
                                            last_it_children_bucket = it_children_bucket;
                                        }

                                        for (j=it_children_pos_bucket*size_line; j<(it_children_pos_bucket+nb_elt)*size_line; j+=size_line){

                                            extractSuffix(kmer_tmp, size_kmer, size_kmer_array, shifting3, it_bucket,
                                                          &(uc->suffixes[j]), info_per_lvl[lvl_node].size_kmer_in_bytes_minus_1, 0);
                                            memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));

                                            //kmer_tmp is the local copy of the current k-mer, it must not be modified
                                            //You can do stuff here with kmer_start (even modify it)
                                            fwrite(kmer_start, sizeof(uint8_t), size_kmer_array, file_output);

                                            for (k=0; k<bucket+3; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);
                                            if (shifting3 != 0) kmer_tmp[bucket+2] &= mask;
                                            memset(&(kmer_tmp[bucket+3]), 0, (size_kmer_array-bucket-3)*sizeof(uint8_t));
                                        }

                                        it_children_pos_bucket += nb_elt;
                                        count += nb_elt;
                                    }
                                    else{

                                        if (shifting3 == 0){
                                            count += extract_kmers_from_node(&(cc->children_Node_container[it_node]), lvl_node-1, kmer_tmp,
                                                                             size_kmer-NB_CHAR_SUF_PREF, it_bucket, 0, size_kmer_root, compression,
                                                                             info_per_lvl, file_output);
                                        }
                                        else{
                                            count += extract_kmers_from_node(&(cc->children_Node_container[it_node]), lvl_node-1, kmer_tmp,
                                                                             size_kmer-NB_CHAR_SUF_PREF, it_bucket, shifting2, size_kmer_root,
                                                                             compression, info_per_lvl, file_output);
                                        }

                                        it_node++;
                                    }
                                }
                                else{
                                    count++;
                                    memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));
                                    for (k=0; k<size_kmer_array; k++) kmer_start[k] = reverse_word_8(kmer_start[k]);
                                    //kmer_tmp is the local copy of the current k-mer, it must not be modified
                                    //You can do stuff here with kmer_start (even modify it)
                                    fwrite(kmer_start, sizeof(uint8_t), size_kmer_array, file_output);
                                }

                                it_filter3++;
                                first_bit=0;
                            }
                        }
                    }
                }
                else{
                    int it = 0;
                    int end = 0;
                    int cpt_pv = 0;
                    int nb_skp = CEIL(cc->nb_elem, info_per_lvl[lvl_node].nb_ucs_skp);

                    nb_cell_children = info_per_lvl[lvl_node].size_kmer_in_bytes_minus_1-1;

                    for (it_filter2=0; it_filter2<MASK_POWER_16[p]; it_filter2++){

                        if ((cc->BF_filter2[size_bf+it_filter2/SIZE_BITS_UINT_8T] & (MASK_POWER_8[it_filter2%SIZE_BITS_UINT_8T])) != 0){

                            first_bit = 1;

                            while (it_children_bucket < nb_skp){

                                if (it_children_bucket == nb_skp - 1) end = cc->nb_elem - it_children_bucket * info_per_lvl[lvl_node].nb_ucs_skp;
                                else end = info_per_lvl[lvl_node].nb_ucs_skp;

                                uc = &(((UC*)cc->children)[it_children_bucket]);
                                size_line = info_per_lvl[lvl_node].size_kmer_in_bytes_minus_1 + uc->size_annot;

                                while (it < end){

                                    memcpy(kmer_tmp, kmer, size_kmer_array*sizeof(uint8_t));

                                    if (IS_ODD(it_filter3)) new_substring = (it_filter2 << 4) | (cc->filter3[it_filter3/2] >> 4);
                                    else new_substring = (it_filter2 << 4) | (cc->filter3[it_filter3/2] & 0xf);

                                    if (compression <= 0) new_substring = (new_substring >> 2) | ((new_substring & 0x3) << 16);

                                    kmer_tmp[bucket] |= new_substring >> shifting1;
                                    kmer_tmp[bucket+1] = new_substring >> shifting2;
                                    kmer_tmp[bucket+2] = new_substring << shifting3;
                                    it_bucket = bucket+2;
                                    if (shifting3 == 0) it_bucket++;

                                    if ((nb_elt = getNbElts(cc, it_filter3, type)) == 0){

                                        if (((cc->children_Node_container[it_node].UC_array.nb_children & 0x1) == 0) || (first_bit == 1)){

                                            first_bit=0;

                                            if (shifting3 == 0){
                                                count += extract_kmers_from_node(&(cc->children_Node_container[it_node]), lvl_node-1, kmer_tmp,
                                                                                 size_kmer-NB_CHAR_SUF_PREF, it_bucket, 0, size_kmer_root,
                                                                                 compression, info_per_lvl, file_output);
                                            }
                                            else{
                                                count += extract_kmers_from_node(&(cc->children_Node_container[it_node]), lvl_node-1, kmer_tmp,
                                                                                 size_kmer-NB_CHAR_SUF_PREF, it_bucket, shifting2, size_kmer_root,
                                                                                 compression, info_per_lvl, file_output);
                                            }

                                            it_node++;
                                        }
                                        else goto OUT_LOOP_S4;
                                    }
                                    else{
                                        if ((uc->suffixes[cpt_pv*size_line+nb_cell_children] < 0x80)  || (first_bit == 1)){

                                            first_bit=0;

                                            for (j=cpt_pv*size_line; j<(cpt_pv+nb_elt)*size_line; j+=size_line){

                                                extractSuffix(kmer_tmp, size_kmer, size_kmer_array, shifting3, it_bucket,
                                                              &(uc->suffixes[j]), info_per_lvl[lvl_node].size_kmer_in_bytes_minus_1, 0);
                                                kmer_tmp[size_kmer_array-1] &= mask_end_kmer;
                                                memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));

                                                //kmer_tmp is the local copy of the current k-mer, it must not be modified
                                                //You can do stuff here with kmer_start (even modify it)
                                                fwrite(kmer_start, sizeof(uint8_t), size_kmer_array, file_output);

                                                for (k=0; k<bucket+3; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);
                                                if (shifting3 != 0) kmer_tmp[bucket+2] &= mask;
                                                memset(&(kmer_tmp[bucket+3]), 0, (size_kmer_array-bucket-3)*sizeof(uint8_t));
                                            }

                                            cpt_pv += nb_elt;
                                            count += nb_elt;
                                        }
                                        else goto OUT_LOOP_S4;
                                    }

                                    it++;
                                    it_filter3++;
                                }

                                it = 0;
                                cpt_pv = 0;
                                it_children_bucket++;
                            }
                        }

                        OUT_LOOP_S4: continue;
                    }
                }
            }
        }
        while ((((CC*)n->CC_array)[i].type & 0x1) == 0);
    }

    if (n->UC_array.suffixes != NULL){

        size_line = info_per_lvl[lvl_node].size_kmer_in_bytes + n->UC_array.size_annot;
        nb_elt = n->UC_array.nb_children >> 1;

        count += nb_elt;

        for (j = 0; j < nb_elt * size_line; j += size_line){

            it_substring = 0;
            it_bucket = bucket;

            shifting_UC = SIZE_BITS_UINT_8T-pos_in_bucket;

            memcpy(kmer_tmp, kmer, size_kmer_array*sizeof(uint8_t));

            while (it_substring < info_per_lvl[lvl_node].size_kmer_in_bytes){

                it_substring += sizeof(uint64_t);
                new_substring = 0;

                if (it_substring > info_per_lvl[lvl_node].size_kmer_in_bytes){

                    size_new_substring = size_kmer*2-((it_substring-sizeof(uint64_t))*SIZE_BITS_UINT_8T);
                    size_new_substring_bytes = CEIL(size_new_substring, SIZE_BITS_UINT_8T);

                    for (k=0; k<size_new_substring_bytes; k++)
                        new_substring = (new_substring << 8) | reverse_word_8(n->UC_array.suffixes[j+(it_substring-sizeof(uint64_t))+k]);

                    new_substring >>= info_per_lvl[lvl_node].size_kmer_in_bytes*SIZE_BITS_UINT_8T - size_new_substring;
                }
                else{

                    size_new_substring = sizeof(uint64_t)*SIZE_BITS_UINT_8T;
                    size_new_substring_bytes = sizeof(uint64_t);

                    for (k=0; k<size_new_substring_bytes; k++)
                        new_substring = (new_substring << 8) | reverse_word_8(n->UC_array.suffixes[j+(it_substring-sizeof(uint64_t))+k]);
                }

                if (shifting_UC != 8) size_new_substring_bytes++;

                for (k=it_bucket; k<it_bucket+size_new_substring_bytes; k++){

                    last_shift = size_new_substring-shifting_UC;

                    if (last_shift >= 0) kmer_tmp[k] |= new_substring >> last_shift;
                    else kmer_tmp[k] |= new_substring << abs(last_shift);

                    shifting_UC += SIZE_BITS_UINT_8T;
                }

                shifting_UC = SIZE_BITS_UINT_8T-pos_in_bucket;

                if ((shifting_UC != 8)/* && (it_substring >= info_per_lvl[lvl_node].size_kmer_in_bytes)*/) size_new_substring_bytes--;
                it_bucket+=size_new_substring_bytes;
            }

            for (k=0; k<size_kmer_array; k++) kmer_tmp[k] = reverse_word_8(kmer_tmp[k]);

            memcpy(kmer_start, kmer_tmp, size_kmer_array*sizeof(uint8_t));
            //kmer_tmp is the local copy of the current k-mer, it must not be modified
            //You can do stuff here with kmer_start (even modify it)
            fwrite(kmer_start, sizeof(uint8_t), size_kmer_array, file_output);
        }
    }

    return count;
}
