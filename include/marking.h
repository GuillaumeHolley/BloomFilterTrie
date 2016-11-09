#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#include "Node.h"
#include "CC.h"

#define MASK_COUNT_FLAG_1 85

void create_marking_Node_4states(Node* n, int lvl_node, int size_kmer, info_per_level*  info_per_lvl);
void create_marking_CC_4states(CC* cc, int lvl_cc, int size_kmer, info_per_level*  info_per_lvl);
void create_marking_UC_4states(UC* uc, int lvl_uc, info_per_level*  info_per_lvl);

void delete_marking_Node_4states(Node* n, int lvl_node, int size_kmer, info_per_level*  info_per_lvl);
void delete_marking_CC_4states(CC* cc, int lvl_cc, int size_kmer, info_per_level*  info_per_lvl);
void delete_marking_UC_4states(UC* uc, int size_substring, int nb_children);

void mark_UC_4states(UC* uc, int size_substring, int nb_children, int position, uint8_t flag);
uint8_t get_mark_UC_4states(UC* uc, int size_substring, int nb_children, int position);
uint8_t get_mark_UC_4states_bis(UC* uc, int position, int nb_bytes_before_marking);

int count_flag_mark_UC_4states(UC* uc, int size_substring, int nb_children);
