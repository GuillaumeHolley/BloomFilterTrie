#include "Node.h"

extern inline bool are_genomes_ids_overlapping(BFT_Root* root_1, BFT_Root* root_2);
extern inline uint64_t* create_hash_v_array(int rand_seed1, int rand_seed2);
extern inline Node* createNode(void);
extern inline void initiateNode(Node* node);
extern inline resultPresence* create_resultPresence();
extern inline void initialize_resultPresence(resultPresence* res);
