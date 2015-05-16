#include "UC.h"

extern void reinit_annotation_inform(annotation_inform* ann_inf);
extern UC* createUC();
extern void initiateUC(UC* uc);
extern void freeUC(UC* uc);

void insertKmer_UC(UC* restrict uc, uint8_t* restrict kmer, int id_genome, int size_kmer, annotation_inform* ann_inf, annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(uc, "insertKmer_UC()")
    if (id_genome > NB_MAX_ID_GENOMES || id_genome < 0) ERROR("insertKmer_UC(): insertion of a genome with unauthorized ID (<1 or >NB_MAX_ID_GENOMES)")

    int nb_cell = CEIL(size_kmer*2, SIZE_CELL);
    int size_line = nb_cell;
    int nb_elem = uc->nb_children >> 1;
    int beginning_cluster = uc->nb_children & 0x1;

    compute_best_mode(ann_inf, annot_sorted, NULL, 0, NULL, 0, id_genome);

    if (nb_elem == 0){

        uc->suffixes = calloc(nb_cell + ann_inf->min_size, sizeof(uint8_t));
        ASSERT_NULL_PTR(uc->suffixes, "insertKmer_UC()")

        uc->nb_children = 0x2 | beginning_cluster;

        uc->size_annot = ann_inf->min_size;
        memcpy(uc->suffixes, kmer, nb_cell);
        modify_mode_annotation(ann_inf, &(uc->suffixes[nb_cell]), uc->size_annot, NULL, 0, id_genome);
    }
    else{

        uint8_t* extend_annot = NULL;

        if (ann_inf->min_size > uc->size_annot){
            extend_annot = realloc_annotation(uc, nb_cell, nb_elem, ann_inf->min_size, 1, nb_elem);
            size_line += uc->size_annot;
        }
        else{
            size_line += uc->size_annot;
            beginning_cluster = uc->nb_extended_annot * SIZE_BYTE_EXT_ANNOT + uc->nb_cplx_nodes * (uc->size_annot_cplx_nodes + SIZE_BYTE_CPLX_N);

            uc->suffixes = realloc(uc->suffixes, ((nb_elem+1) * size_line + beginning_cluster) * sizeof(uint8_t));

            memmove(&(uc->suffixes[(nb_elem+1) * size_line]), &(uc->suffixes[nb_elem * size_line]), beginning_cluster * sizeof(uint8_t));
        }

        if (uc->suffixes != NULL){
            size_line = nb_elem * size_line;

            memcpy(&(uc->suffixes[size_line]), kmer, nb_cell);
            uc->nb_children += 2;

            memset(&(uc->suffixes[size_line + nb_cell]), 0, uc->size_annot * sizeof(uint8_t));
            modify_mode_annotation(ann_inf, &(uc->suffixes[size_line + nb_cell]), uc->size_annot, extend_annot, 1, id_genome);
        }
        else ERROR("insertKmer_UC(): out of memory\n")
    }

    reinit_annotation_inform(ann_inf);

    return;
}
