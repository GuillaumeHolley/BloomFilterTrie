#ifndef DEF_UC
#define DEF_UC

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#include "./../lib/useful_macros.h"
#include "./../lib/UC_annotation.h"

inline UC* createUC(){

    UC* uc = calloc(1,sizeof(UC));
    ASSERT_NULL_PTR(uc,"createUC()")

    uc->suffixes = NULL;

    return uc;
}

inline void initiateUC(UC* uc){

    ASSERT_NULL_PTR(uc,"initiateUC()")

    uc->nb_children &= 1;
    uc->suffixes = NULL;

    return;
}

inline void initializeUC(UC* uc){

    ASSERT_NULL_PTR(uc,"initializeUC()")

    uc->size_annot = 0;
    uc->size_annot_cplx_nodes = 0;
    uc->nb_extended_annot = 0;
    uc->nb_cplx_nodes = 0;
    uc->nb_children = 0;
    uc->suffixes = NULL;

    return;
}

inline void resetUC(UC* uc){

    ASSERT_NULL_PTR(uc,"resetUC()")

    uc->size_annot = 0;
    uc->size_annot_cplx_nodes = 0;
    uc->nb_extended_annot = 0;
    uc->nb_cplx_nodes = 0;

    uc->nb_children &= 1;

    uc->suffixes = NULL;

    return;
}

inline void freeUC(UC* uc){
    free(uc->suffixes);
    free(uc);
    return;
}

#endif
