#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>

#include <jemalloc/jemalloc.h>
const char* malloc_conf = "narenas:1,tcache:false,lg_dirty_mult:8,lg_chunk:22";

#include <Judy.h>

#include "UC_annotation.h"
#include "insertNode.h"
#include "branchingNode.h"
#include "deleteColorsNode.h"
#include "fasta.h"
#include "printMemory.h"
#include "replaceAnnotation.h"

#define PRINT_EVERY_X_KMERS 1000000
#define NB_FILE_2_READ 20

#define SIZE_BUFFER 4096 //Size of the buffer (in bytes) used to store kmers read from input files
#define TRESH_DEL_ANNOT 43

#include "kseq.h"
//Initialize parser for FASTA/FASTQ files
KSEQ_INIT(int, read, SIZE_BUFFER)

void insertKmers(Root* restrict root,
                 uint8_t* restrict tab_kmers,
                 int size_kmers,
                 int nb_kmers,
                 int id_genome,
                 ptrs_on_func* restrict func_on_types,
                 annotation_inform* ann_inf,
                 resultPresence* res,
                 annotation_array_elem* annot_sorted);

void insert_Genomes_from_KmerCounting(Root* tree, int binary_files, int size_kmer, annotation_array_elem** ptr_annot_sorted);
void insert_Genomes_from_FASTx(Root* tree, int size_kmer, annotation_array_elem** ptr_annot_sorted);

int get_nb_cplx_nodes_from_KmerCounting(Root* tree, char* name_file, int size_kmer, uint16_t** skip_node_root, ptrs_on_func* func_on_types,
                                        annotation_inform* ann_inf);
Root* get_nb_cplx_nodes_from_FASTx(Root* tree, int size_kmer, uint16_t** skip_node_root, annotation_array_elem** ptr_annot_sorted,
                                   ptrs_on_func* func_on_types, annotation_inform* ann_inf, resultPresence* res);

extern void freeNode(Node* restrict node);

int main()
{

    Root* tree = NULL; //Create the tree by creating an empty root node
    annotation_array_elem* annot_sorted = NULL;

    int i = 0;
    int size_kmer = 63; //Size of kmers to insert
    int binary_files = 1;

    char** name_file = malloc(NB_FILE_2_READ*sizeof(char*));
    ASSERT_NULL_PTR(name_file,"main()")

    for (i=0; i<NB_FILE_2_READ; i++){
        name_file[i] = malloc(100*sizeof(char));
        ASSERT_NULL_PTR(name_file[i],"main()")
    }

    /*name_file[0] = "test_data/ERR430991.k63_cutoff3_shuf_bin";
    name_file[1] = "test_data/ERR430992.k63_cutoff3_shuf_bin";
    name_file[2] = "test_data/ERR430993.k63_cutoff3_shuf_bin";
    name_file[3] = "test_data/ERR430994.k63_cutoff3_shuf_bin";
    name_file[4] = "test_data/ERR430995.k63_cutoff3_shuf_bin";
    name_file[5] = "test_data/ERR430996.k63_cutoff3_shuf_bin";
    name_file[6] = "test_data/ERR430997.k63_cutoff3_shuf_bin";
    name_file[7] = "test_data/ERR430998.k63_cutoff3_shuf_bin";
    name_file[8] = "test_data/ERR430999.k63_cutoff3_shuf_bin";
    name_file[9] = "test_data/ERR431000.k63_cutoff3_shuf_bin";
    name_file[10] = "test_data/ERR431001.k63_cutoff3_shuf_bin";
    name_file[11] = "test_data/ERR431002.k63_cutoff3_shuf_bin";
    name_file[12] = "test_data/ERR431003.k63_cutoff3_shuf_bin";
    name_file[13] = "test_data/ERR431004.k63_cutoff3_shuf_bin";
    name_file[14] = "test_data/ERR431005.k63_cutoff3_shuf_bin";
    name_file[15] = "test_data/ERR431006.k63_cutoff3_shuf_bin";
    name_file[16] = "test_data/ERR431007.k63_cutoff3_shuf_bin";
    name_file[17] = "test_data/ERR431008.k63_cutoff3_shuf_bin";
    name_file[18] = "test_data/ERR431009.k63_cutoff3_shuf_bin";
    name_file[19] = "test_data/ERR431010.k63_cutoff3_shuf_bin";
    name_file[20] = "test_data/ERR431011.k63_cutoff3_shuf_bin";
    name_file[21] = "test_data/ERR431012.k63_cutoff3_shuf_bin";
    name_file[22] = "test_data/ERR431013.k63_cutoff3_shuf_bin";
    name_file[23] = "test_data/ERR431014.k63_cutoff3_shuf_bin";
    name_file[24] = "test_data/ERR431015.k63_cutoff3_shuf_bin";
    name_file[25] = "test_data/ERR431016.k63_cutoff3_shuf_bin";
    name_file[26] = "test_data/ERR431017.k63_cutoff3_shuf_bin";
    name_file[27] = "test_data/ERR431018.k63_cutoff3_shuf_bin";
    name_file[28] = "test_data/ERR431019.k63_cutoff3_shuf_bin";
    name_file[29] = "test_data/ERR431020.k63_cutoff3_shuf_bin";
    name_file[30] = "test_data/ERR431021.k63_cutoff3_shuf_bin";
    name_file[31] = "test_data/ERR431022.k63_cutoff3_shuf_bin";
    name_file[32] = "test_data/ERR431023.k63_cutoff3_shuf_bin";
    name_file[33] = "test_data/ERR431024.k63_cutoff3_shuf_bin";
    name_file[34] = "test_data/ERR431025.k63_cutoff3_shuf_bin";
    name_file[35] = "test_data/ERR431026.k63_cutoff3_shuf_bin";
    name_file[36] = "test_data/ERR431027.k63_cutoff3_shuf_bin";
    name_file[37] = "test_data/ERR431028.k63_cutoff3_shuf_bin";
    name_file[38] = "test_data/ERR431029.k63_cutoff3_shuf_bin";
    name_file[39] = "test_data/ERR431030.k63_cutoff3_shuf_bin";
    name_file[40] = "test_data/ERR431031.k63_cutoff3_shuf_bin";
    name_file[41] = "test_data/ERR431032.k63_cutoff3_shuf_bin";
    name_file[42] = "test_data/ERR431033.k63_cutoff3_shuf_bin";
    name_file[43] = "test_data/ERR431034.k63_cutoff3_shuf_bin";
    name_file[44] = "test_data/ERR431035.k63_cutoff3_shuf_bin";
    name_file[45] = "test_data/ERR431036.k63_cutoff3_shuf_bin";
    name_file[46] = "test_data/ERR431037.k63_cutoff3_shuf_bin";
    name_file[47] = "test_data/ERR431038.k63_cutoff3_shuf_bin";
    name_file[48] = "test_data/ERR431039.k63_cutoff3_shuf_bin";
    name_file[49] = "test_data/ERR431040.k63_cutoff3_shuf_bin";
    name_file[50] = "test_data/ERR431041.k63_cutoff3_shuf_bin";
    name_file[51] = "test_data/ERR431042.k63_cutoff3_shuf_bin";
    name_file[52] = "test_data/ERR431043.k63_cutoff3_shuf_bin";
    name_file[53] = "test_data/ERR431044.k63_cutoff3_shuf_bin";
    name_file[54] = "test_data/ERR431045.k63_cutoff3_shuf_bin";
    name_file[55] = "test_data/ERR431046.k63_cutoff3_shuf_bin";
    name_file[56] = "test_data/ERR431047.k63_cutoff3_shuf_bin";
    name_file[57] = "test_data/ERR431048.k63_cutoff3_shuf_bin";
    name_file[58] = "test_data/ERR431049.k63_cutoff3_shuf_bin";
    name_file[59] = "test_data/ERR431050.k63_cutoff3_shuf_bin";
    name_file[60] = "test_data/ERR431051.k63_cutoff3_shuf_bin";
    name_file[61] = "test_data/ERR431052.k63_cutoff3_shuf_bin";
    name_file[62] = "test_data/ERR431053.k63_cutoff3_shuf_bin";
    name_file[63] = "test_data/ERR431054.k63_cutoff3_shuf_bin";
    name_file[64] = "test_data/ERR431055.k63_cutoff3_shuf_bin";
    name_file[65] = "test_data/ERR431056.k63_cutoff3_shuf_bin";
    name_file[66] = "test_data/ERR431057.k63_cutoff3_shuf_bin";
    name_file[67] = "test_data/ERR431058.k63_cutoff3_shuf_bin";
    name_file[68] = "test_data/ERR431059.k63_cutoff3_shuf_bin";
    name_file[69] = "test_data/ERR431060.k63_cutoff3_shuf_bin";
    name_file[70] = "test_data/ERR431061.k63_cutoff3_shuf_bin";
    name_file[71] = "test_data/ERR431062.k63_cutoff3_shuf_bin";
    name_file[72] = "test_data/ERR431063.k63_cutoff3_shuf_bin";
    name_file[73] = "test_data/ERR431064.k63_cutoff3_shuf_bin";
    name_file[74] = "test_data/ERR431065.k63_cutoff3_shuf_bin";
    name_file[75] = "test_data/ERR431066.k63_cutoff3_shuf_bin";
    name_file[76] = "test_data/ERR431067.k63_cutoff3_shuf_bin";
    name_file[77] = "test_data/ERR431068.k63_cutoff3_shuf_bin";
    name_file[78] = "test_data/ERR431069.k63_cutoff3_shuf_bin";
    name_file[79] = "test_data/ERR431070.k63_cutoff3_shuf_bin";
    name_file[80] = "test_data/ERR431071.k63_cutoff3_shuf_bin";
    name_file[81] = "test_data/ERR431072.k63_cutoff3_shuf_bin";
    name_file[82] = "test_data/ERR431073.k63_cutoff3_shuf_bin";
    name_file[83] = "test_data/ERR431074.k63_cutoff3_shuf_bin";
    name_file[84] = "test_data/ERR431075.k63_cutoff3_shuf_bin";
    name_file[85] = "test_data/ERR431076.k63_cutoff3_shuf_bin";
    name_file[86] = "test_data/ERR431077.k63_cutoff3_shuf_bin";
    name_file[87] = "test_data/ERR431078.k63_cutoff3_shuf_bin";
    name_file[88] = "test_data/ERR431079.k63_cutoff3_shuf_bin";
    name_file[89] = "test_data/ERR431080.k63_cutoff3_shuf_bin";
    name_file[90] = "test_data/ERR431081.k63_cutoff3_shuf_bin";
    name_file[91] = "test_data/ERR431082.k63_cutoff3_shuf_bin";
    name_file[92] = "test_data/ERR431083.k63_cutoff3_shuf_bin";
    name_file[93] = "test_data/ERR431084.k63_cutoff3_shuf_bin";
    name_file[94] = "test_data/ERR431085.k63_cutoff3_shuf_bin";
    name_file[95] = "test_data/ERR431086.k63_cutoff3_shuf_bin";
    name_file[96] = "test_data/ERR431087.k63_cutoff3_shuf_bin";
    name_file[97] = "test_data/ERR431088.k63_cutoff3_shuf_bin";
    name_file[98] = "test_data/ERR431089.k63_cutoff3_shuf_bin";
    name_file[99] = "test_data/ERR431090.k63_cutoff3_shuf_bin";
    name_file[100] = "test_data/ERR431091.k63_cutoff3_shuf_bin";
    name_file[101] = "test_data/ERR431092.k63_cutoff3_shuf_bin";
    name_file[102] = "test_data/ERR431093.k63_cutoff3_shuf_bin";
    name_file[103] = "test_data/ERR431094.k63_cutoff3_shuf_bin";
    name_file[104] = "test_data/ERR431095.k63_cutoff3_shuf_bin";
    name_file[105] = "test_data/ERR431096.k63_cutoff3_shuf_bin";
    name_file[106] = "test_data/ERR431097.k63_cutoff3_shuf_bin";
    name_file[107] = "test_data/ERR431098.k63_cutoff3_shuf_bin";
    name_file[108] = "test_data/ERR431099.k63_cutoff3_shuf_bin";
    name_file[109] = "test_data/ERR431100.k63_cutoff3_shuf_bin";
    name_file[110] = "test_data/ERR431101.k63_cutoff3_shuf_bin";
    name_file[111] = "test_data/ERR431102.k63_cutoff3_shuf_bin";
    name_file[112] = "test_data/ERR431103.k63_cutoff3_shuf_bin";
    name_file[113] = "test_data/ERR431104.k63_cutoff3_shuf_bin";
    name_file[114] = "test_data/ERR431105.k63_cutoff3_shuf_bin";
    name_file[115] = "test_data/ERR431106.k63_cutoff3_shuf_bin";
    name_file[116] = "test_data/ERR431107.k63_cutoff3_shuf_bin";
    name_file[117] = "test_data/ERR431108.k63_cutoff3_shuf_bin";
    name_file[118] = "test_data/ERR431109.k63_cutoff3_shuf_bin";
    name_file[119] = "test_data/ERR431110.k63_cutoff3_shuf_bin";
    name_file[120] = "test_data/ERR431111.k63_cutoff3_shuf_bin";
    name_file[121] = "test_data/ERR431112.k63_cutoff3_shuf_bin";
    name_file[122] = "test_data/ERR431113.k63_cutoff3_shuf_bin";
    name_file[123] = "test_data/ERR431114.k63_cutoff3_shuf_bin";
    name_file[124] = "test_data/ERR431115.k63_cutoff3_shuf_bin";
    name_file[125] = "test_data/ERR431116.k63_cutoff3_shuf_bin";
    name_file[126] = "test_data/ERR431117.k63_cutoff3_shuf_bin";
    name_file[127] = "test_data/ERR431118.k63_cutoff3_shuf_bin";
    name_file[128] = "test_data/ERR431119.k63_cutoff3_shuf_bin";
    name_file[129] = "test_data/ERR431120.k63_cutoff3_shuf_bin";
    name_file[130] = "test_data/ERR431121.k63_cutoff3_shuf_bin";
    name_file[131] = "test_data/ERR431122.k63_cutoff3_shuf_bin";
    name_file[132] = "test_data/ERR431123.k63_cutoff3_shuf_bin";
    name_file[133] = "test_data/ERR431124.k63_cutoff3_shuf_bin";
    name_file[134] = "test_data/ERR431125.k63_cutoff3_shuf_bin";
    name_file[135] = "test_data/ERR431126.k63_cutoff3_shuf_bin";
    name_file[136] = "test_data/ERR431127.k63_cutoff3_shuf_bin";
    name_file[137] = "test_data/ERR431128.k63_cutoff3_shuf_bin";
    name_file[138] = "test_data/ERR431129.k63_cutoff3_shuf_bin";
    name_file[139] = "test_data/ERR431130.k63_cutoff3_shuf_bin";
    name_file[140] = "test_data/ERR431131.k63_cutoff3_shuf_bin";
    name_file[141] = "test_data/ERR431132.k63_cutoff3_shuf_bin";
    name_file[142] = "test_data/ERR431133.k63_cutoff3_shuf_bin";
    name_file[143] = "test_data/ERR431134.k63_cutoff3_shuf_bin";
    name_file[144] = "test_data/ERR431135.k63_cutoff3_shuf_bin";
    name_file[145] = "test_data/ERR431136.k63_cutoff3_shuf_bin";
    name_file[146] = "test_data/ERR431137.k63_cutoff3_shuf_bin";
    name_file[147] = "test_data/ERR431138.k63_cutoff3_shuf_bin";
    name_file[148] = "test_data/ERR431139.k63_cutoff3_shuf_bin";
    name_file[149] = "test_data/ERR431140.k63_cutoff3_shuf_bin";
    name_file[150] = "test_data/ERR431141.k63_cutoff3_shuf_bin";
    name_file[151] = "test_data/ERR431142.k63_cutoff3_shuf_bin";
    name_file[152] = "test_data/ERR431143.k63_cutoff3_shuf_bin";
    name_file[153] = "test_data/ERR431144.k63_cutoff3_shuf_bin";
    name_file[154] = "test_data/ERR431145.k63_cutoff3_shuf_bin";
    name_file[155] = "test_data/ERR431146.k63_cutoff3_shuf_bin";
    name_file[156] = "test_data/ERR431147.k63_cutoff3_shuf_bin";
    name_file[157] = "test_data/ERR431148.k63_cutoff3_shuf_bin";
    name_file[158] = "test_data/ERR431149.k63_cutoff3_shuf_bin";
    name_file[159] = "test_data/ERR431150.k63_cutoff3_shuf_bin";
    name_file[160] = "test_data/ERR431151.k63_cutoff3_shuf_bin";
    name_file[161] = "test_data/ERR431152.k63_cutoff3_shuf_bin";
    name_file[162] = "test_data/ERR431153.k63_cutoff3_shuf_bin";
    name_file[163] = "test_data/ERR431154.k63_cutoff3_shuf_bin";
    name_file[164] = "test_data/ERR431155.k63_cutoff3_shuf_bin";
    name_file[165] = "test_data/ERR431156.k63_cutoff3_shuf_bin";
    name_file[166] = "test_data/ERR431157.k63_cutoff3_shuf_bin";
    name_file[167] = "test_data/ERR431158.k63_cutoff3_shuf_bin";
    name_file[168] = "test_data/ERR431159.k63_cutoff3_shuf_bin";
    name_file[169] = "test_data/ERR431160.k63_cutoff3_shuf_bin";
    name_file[170] = "test_data/ERR431161.k63_cutoff3_shuf_bin";
    name_file[171] = "test_data/ERR431162.k63_cutoff3_shuf_bin";
    name_file[172] = "test_data/ERR431163.k63_cutoff3_shuf_bin";
    name_file[173] = "test_data/ERR431164.k63_cutoff3_shuf_bin";
    name_file[174] = "test_data/ERR431165.k63_cutoff3_shuf_bin";
    name_file[175] = "test_data/ERR431166.k63_cutoff3_shuf_bin";
    name_file[176] = "test_data/ERR431167.k63_cutoff3_shuf_bin";
    name_file[177] = "test_data/ERR431168.k63_cutoff3_shuf_bin";
    name_file[178] = "test_data/ERR431169.k63_cutoff3_shuf_bin";
    name_file[179] = "test_data/ERR431170.k63_cutoff3_shuf_bin";
    name_file[180] = "test_data/ERR431171.k63_cutoff3_shuf_bin";
    name_file[181] = "test_data/ERR431172.k63_cutoff3_shuf_bin";
    name_file[182] = "test_data/ERR431173.k63_cutoff3_shuf_bin";
    name_file[183] = "test_data/ERR431174.k63_cutoff3_shuf_bin";
    name_file[184] = "test_data/ERR431175.k63_cutoff3_shuf_bin";
    name_file[185] = "test_data/ERR431176.k63_cutoff3_shuf_bin";
    name_file[186] = "test_data/ERR431177.k63_cutoff3_shuf_bin";
    name_file[187] = "test_data/ERR431178.k63_cutoff3_shuf_bin";
    name_file[188] = "test_data/ERR431179.k63_cutoff3_shuf_bin";
    name_file[189] = "test_data/ERR431180.k63_cutoff3_shuf_bin";
    name_file[190] = "test_data/ERR431181.k63_cutoff3_shuf_bin";
    name_file[191] = "test_data/ERR431182.k63_cutoff3_shuf_bin";
    name_file[192] = "test_data/ERR431183.k63_cutoff3_shuf_bin";
    name_file[193] = "test_data/ERR431184.k63_cutoff3_shuf_bin";
    name_file[194] = "test_data/ERR431185.k63_cutoff3_shuf_bin";
    name_file[195] = "test_data/ERR431185.k63_cutoff3_shuf_bin";
    name_file[196] = "test_data/ERR431186.k63_cutoff3_shuf_bin";
    name_file[197] = "test_data/ERR431187.k63_cutoff3_shuf_bin";
    name_file[198] = "test_data/ERR431188.k63_cutoff3_shuf_bin";
    name_file[199] = "test_data/ERR431189.k63_cutoff3_shuf_bin";
    name_file[200] = "test_data/ERR431190.k63_cutoff3_shuf_bin";
    name_file[201] = "test_data/ERR431191.k63_cutoff3_shuf_bin";
    name_file[202] = "test_data/ERR431192.k63_cutoff3_shuf_bin";
    name_file[203] = "test_data/ERR431193.k63_cutoff3_shuf_bin";
    name_file[204] = "test_data/ERR431194.k63_cutoff3_shuf_bin";
    name_file[205] = "test_data/ERR431195.k63_cutoff3_shuf_bin";
    name_file[206] = "test_data/ERR431196.k63_cutoff3_shuf_bin";
    name_file[207] = "test_data/ERR431197.k63_cutoff3_shuf_bin";
    name_file[208] = "test_data/ERR431198.k63_cutoff3_shuf_bin";
    name_file[209] = "test_data/ERR431199.k63_cutoff3_shuf_bin";
    name_file[210] = "test_data/ERR431200.k63_cutoff3_shuf_bin";
    name_file[211] = "test_data/ERR431201.k63_cutoff3_shuf_bin";
    name_file[212] = "test_data/ERR431202.k63_cutoff3_shuf_bin";
    name_file[213] = "test_data/ERR431203.k63_cutoff3_shuf_bin";
    name_file[214] = "test_data/ERR431204.k63_cutoff3_shuf_bin";
    name_file[215] = "test_data/ERR431205.k63_cutoff3_shuf_bin";
    name_file[216] = "test_data/ERR431206.k63_cutoff3_shuf_bin";
    name_file[217] = "test_data/ERR431207.k63_cutoff3_shuf_bin";
    name_file[218] = "test_data/ERR431208.k63_cutoff3_shuf_bin";
    name_file[219] = "test_data/ERR431209.k63_cutoff3_shuf_bin";
    name_file[220] = "test_data/ERR431210.k63_cutoff3_shuf_bin";
    name_file[221] = "test_data/ERR431211.k63_cutoff3_shuf_bin";
    name_file[222] = "test_data/ERR431212.k63_cutoff3_shuf_bin";
    name_file[223] = "test_data/ERR431213.k63_cutoff3_shuf_bin";
    name_file[224] = "test_data/ERR431214.k63_cutoff3_shuf_bin";
    name_file[225] = "test_data/ERR431215.k63_cutoff3_shuf_bin";
    name_file[226] = "test_data/ERR431216.k63_cutoff3_shuf_bin";
    name_file[227] = "test_data/ERR431217.k63_cutoff3_shuf_bin";
    name_file[228] = "test_data/ERR431218.k63_cutoff3_shuf_bin";
    name_file[229] = "test_data/ERR431219.k63_cutoff3_shuf_bin";
    name_file[230] = "test_data/ERR431220.k63_cutoff3_shuf_bin";
    name_file[231] = "test_data/ERR431221.k63_cutoff3_shuf_bin";
    name_file[232] = "test_data/ERR431222.k63_cutoff3_shuf_bin";
    name_file[233] = "test_data/ERR431223.k63_cutoff3_shuf_bin";
    name_file[234] = "test_data/ERR431224.k63_cutoff3_shuf_bin";
    name_file[235] = "test_data/ERR431225.k63_cutoff3_shuf_bin";
    name_file[236] = "test_data/ERR431226.k63_cutoff3_shuf_bin";
    name_file[237] = "test_data/ERR431227.k63_cutoff3_shuf_bin";
    name_file[238] = "test_data/ERR431228.k63_cutoff3_shuf_bin";
    name_file[239] = "test_data/ERR431229.k63_cutoff3_shuf_bin";
    name_file[240] = "test_data/ERR431230.k63_cutoff3_shuf_bin";
    name_file[241] = "test_data/ERR431231.k63_cutoff3_shuf_bin";
    name_file[242] = "test_data/ERR431232.k63_cutoff3_shuf_bin";
    name_file[243] = "test_data/ERR431233.k63_cutoff3_shuf_bin";
    name_file[244] = "test_data/ERR431234.k63_cutoff3_shuf_bin";
    name_file[245] = "test_data/ERR431235.k63_cutoff3_shuf_bin";
    name_file[246] = "test_data/ERR431236.k63_cutoff3_shuf_bin";
    name_file[247] = "test_data/ERR431237.k63_cutoff3_shuf_bin";
    name_file[248] = "test_data/ERR431238.k63_cutoff3_shuf_bin";
    name_file[249] = "test_data/ERR431239.k63_cutoff3_shuf_bin";
    name_file[250] = "test_data/ERR431240.k63_cutoff3_shuf_bin";
    name_file[251] = "test_data/ERR431241.k63_cutoff3_shuf_bin";
    name_file[252] = "test_data/ERR431242.k63_cutoff3_shuf_bin";
    name_file[253] = "test_data/ERR431243.k63_cutoff3_shuf_bin";
    name_file[254] = "test_data/ERR431244.k63_cutoff3_shuf_bin";
    name_file[255] = "test_data/ERR431245.k63_cutoff3_shuf_bin";
    name_file[256] = "test_data/ERR431246.k63_cutoff3_shuf_bin";
    name_file[257] = "test_data/ERR431247.k63_cutoff3_shuf_bin";
    name_file[258] = "test_data/ERR431248.k63_cutoff3_shuf_bin";
    name_file[259] = "test_data/ERR431249.k63_cutoff3_shuf_bin";
    name_file[260] = "test_data/ERR431250.k63_cutoff3_shuf_bin";
    name_file[261] = "test_data/ERR431251.k63_cutoff3_shuf_bin";
    name_file[262] = "test_data/ERR431252.k63_cutoff3_shuf_bin";
    name_file[263] = "test_data/ERR431253.k63_cutoff3_shuf_bin";
    name_file[264] = "test_data/ERR431254.k63_cutoff3_shuf_bin";
    name_file[265] = "test_data/ERR431255.k63_cutoff3_shuf_bin";
    name_file[266] = "test_data/ERR431256.k63_cutoff3_shuf_bin";
    name_file[267] = "test_data/ERR431257.k63_cutoff3_shuf_bin";
    name_file[268] = "test_data/ERR431258.k63_cutoff3_shuf_bin";
    name_file[269] = "test_data/ERR431259.k63_cutoff3_shuf_bin";
    name_file[270] = "test_data/ERR431260.k63_cutoff3_shuf_bin";
    name_file[271] = "test_data/ERR431261.k63_cutoff3_shuf_bin";
    name_file[272] = "test_data/ERR431262.k63_cutoff3_shuf_bin";
    name_file[273] = "test_data/ERR431263.k63_cutoff3_shuf_bin";
    name_file[274] = "test_data/ERR431264.k63_cutoff3_shuf_bin";
    name_file[275] = "test_data/ERR431265.k63_cutoff3_shuf_bin";
    name_file[276] = "test_data/ERR431266.k63_cutoff3_shuf_bin";
    name_file[277] = "test_data/ERR431267.k63_cutoff3_shuf_bin";
    name_file[278] = "test_data/ERR431268.k63_cutoff3_shuf_bin";
    name_file[279] = "test_data/ERR431269.k63_cutoff3_shuf_bin";
    name_file[280] = "test_data/ERR431270.k63_cutoff3_shuf_bin";
    name_file[281] = "test_data/ERR431271.k63_cutoff3_shuf_bin";
    name_file[282] = "test_data/ERR431272.k63_cutoff3_shuf_bin";
    name_file[283] = "test_data/ERR431273.k63_cutoff3_shuf_bin";
    name_file[284] = "test_data/ERR431274.k63_cutoff3_shuf_bin";
    name_file[285] = "test_data/ERR431275.k63_cutoff3_shuf_bin";
    name_file[286] = "test_data/ERR431276.k63_cutoff3_shuf_bin";
    name_file[287] = "test_data/ERR431277.k63_cutoff3_shuf_bin";
    name_file[288] = "test_data/ERR431278.k63_cutoff3_shuf_bin";
    name_file[289] = "test_data/ERR431279.k63_cutoff3_shuf_bin";
    name_file[290] = "test_data/ERR431280.k63_cutoff3_shuf_bin";
    name_file[291] = "test_data/ERR431281.k63_cutoff3_shuf_bin";
    name_file[292] = "test_data/ERR431282.k63_cutoff3_shuf_bin";
    name_file[293] = "test_data/ERR431283.k63_cutoff3_shuf_bin";
    name_file[294] = "test_data/ERR431284.k63_cutoff3_shuf_bin";
    name_file[295] = "test_data/ERR431285.k63_cutoff3_shuf_bin";
    name_file[296] = "test_data/ERR431286.k63_cutoff3_shuf_bin";
    name_file[297] = "test_data/ERR431287.k63_cutoff3_shuf_bin";
    name_file[298] = "test_data/ERR431288.k63_cutoff3_shuf_bin";
    name_file[299] = "test_data/ERR431289.k63_cutoff3_shuf_bin";
    name_file[300] = "test_data/ERR431290.k63_cutoff3_shuf_bin";
    name_file[301] = "test_data/ERR431291.k63_cutoff3_shuf_bin";
    name_file[302] = "test_data/ERR431292.k63_cutoff3_shuf_bin";
    name_file[303] = "test_data/ERR431293.k63_cutoff3_shuf_bin";
    name_file[304] = "test_data/ERR431294.k63_cutoff3_shuf_bin";
    name_file[305] = "test_data/ERR431295.k63_cutoff3_shuf_bin";
    name_file[306] = "test_data/ERR431296.k63_cutoff3_shuf_bin";
    name_file[307] = "test_data/ERR431297.k63_cutoff3_shuf_bin";
    name_file[308] = "test_data/ERR431298.k63_cutoff3_shuf_bin";
    name_file[309] = "test_data/ERR431299.k63_cutoff3_shuf_bin";
    name_file[310] = "test_data/ERR431300.k63_cutoff3_shuf_bin";
    name_file[311] = "test_data/ERR431301.k63_cutoff3_shuf_bin";
    name_file[312] = "test_data/ERR431302.k63_cutoff3_shuf_bin";
    name_file[313] = "test_data/ERR431303.k63_cutoff3_shuf_bin";
    name_file[314] = "test_data/ERR431304.k63_cutoff3_shuf_bin";
    name_file[315] = "test_data/ERR431305.k63_cutoff3_shuf_bin";
    name_file[316] = "test_data/ERR431306.k63_cutoff3_shuf_bin";
    name_file[317] = "test_data/ERR431307.k63_cutoff3_shuf_bin";
    name_file[318] = "test_data/ERR431308.k63_cutoff3_shuf_bin";
    name_file[319] = "test_data/ERR431309.k63_cutoff3_shuf_bin";
    name_file[320] = "test_data/ERR431310.k63_cutoff3_shuf_bin";
    name_file[321] = "test_data/ERR431311.k63_cutoff3_shuf_bin";
    name_file[322] = "test_data/ERR431312.k63_cutoff3_shuf_bin";
    name_file[323] = "test_data/ERR431313.k63_cutoff3_shuf_bin";
    name_file[324] = "test_data/ERR431314.k63_cutoff3_shuf_bin";
    name_file[325] = "test_data/ERR431315.k63_cutoff3_shuf_bin";
    name_file[326] = "test_data/ERR431316.k63_cutoff3_shuf_bin";
    name_file[327] = "test_data/ERR431317.k63_cutoff3_shuf_bin";
    name_file[328] = "test_data/ERR431318.k63_cutoff3_shuf_bin";
    name_file[329] = "test_data/ERR431319.k63_cutoff3_shuf_bin";
    name_file[330] = "test_data/ERR431320.k63_cutoff3_shuf_bin";
    name_file[331] = "test_data/ERR431321.k63_cutoff3_shuf_bin";
    name_file[332] = "test_data/ERR431322.k63_cutoff3_shuf_bin";
    name_file[333] = "test_data/ERR431323.k63_cutoff3_shuf_bin";
    name_file[334] = "test_data/ERR431324.k63_cutoff3_shuf_bin";
    name_file[335] = "test_data/ERR431325.k63_cutoff3_shuf_bin";
    name_file[336] = "test_data/ERR431326.k63_cutoff3_shuf_bin";
    name_file[337] = "test_data/ERR431327.k63_cutoff3_shuf_bin";
    name_file[338] = "test_data/ERR431328.k63_cutoff3_shuf_bin";
    name_file[339] = "test_data/ERR431329.k63_cutoff3_shuf_bin";
    name_file[340] = "test_data/ERR431330.k63_cutoff3_shuf_bin";
    name_file[341] = "test_data/ERR431331.k63_cutoff3_shuf_bin";
    name_file[342] = "test_data/ERR431332.k63_cutoff3_shuf_bin";
    name_file[343] = "test_data/ERR431333.k63_cutoff3_shuf_bin";
    name_file[344] = "test_data/ERR431334.k63_cutoff3_shuf_bin";
    name_file[345] = "test_data/ERR431335.k63_cutoff3_shuf_bin";
    name_file[346] = "test_data/ERR431336.k63_cutoff3_shuf_bin";
    name_file[347] = "test_data/ERR431337.k63_cutoff3_shuf_bin";
    name_file[348] = "test_data/ERR431338.k63_cutoff3_shuf_bin";
    name_file[349] = "test_data/ERR431339.k63_cutoff3_shuf_bin";
    name_file[350] = "test_data/ERR431340.k63_cutoff3_shuf_bin";
    name_file[351] = "test_data/ERR431341.k63_cutoff3_shuf_bin";
    name_file[352] = "test_data/ERR431342.k63_cutoff3_shuf_bin";
    name_file[353] = "test_data/ERR431343.k63_cutoff3_shuf_bin";
    name_file[354] = "test_data/ERR431344.k63_cutoff3_shuf_bin";
    name_file[355] = "test_data/ERR431345.k63_cutoff3_shuf_bin";
    name_file[356] = "test_data/ERR431346.k63_cutoff3_shuf_bin";
    name_file[357] = "test_data/ERR431347.k63_cutoff3_shuf_bin";
    name_file[358] = "test_data/ERR431348.k63_cutoff3_shuf_bin";
    name_file[359] = "test_data/ERR431349.k63_cutoff3_shuf_bin";
    name_file[360] = "test_data/ERR431350.k63_cutoff3_shuf_bin";
    name_file[361] = "test_data/ERR431351.k63_cutoff3_shuf_bin";
    name_file[362] = "test_data/ERR431352.k63_cutoff3_shuf_bin";
    name_file[363] = "test_data/ERR431353.k63_cutoff3_shuf_bin";
    name_file[364] = "test_data/ERR431354.k63_cutoff3_shuf_bin";
    name_file[365] = "test_data/ERR431355.k63_cutoff3_shuf_bin";
    name_file[366] = "test_data/ERR431356.k63_cutoff3_shuf_bin";
    name_file[367] = "test_data/ERR431357.k63_cutoff3_shuf_bin";
    name_file[368] = "test_data/ERR431358.k63_cutoff3_shuf_bin";
    name_file[369] = "test_data/ERR431359.k63_cutoff3_shuf_bin";
    name_file[370] = "test_data/ERR431360.k63_cutoff3_shuf_bin";
    name_file[371] = "test_data/ERR431361.k63_cutoff3_shuf_bin";
    name_file[372] = "test_data/ERR431362.k63_cutoff3_shuf_bin";
    name_file[373] = "test_data/ERR431363.k63_cutoff3_shuf_bin";
    name_file[374] = "test_data/ERR431364.k63_cutoff3_shuf_bin";
    name_file[375] = "test_data/ERR431365.k63_cutoff3_shuf_bin";
    name_file[376] = "test_data/ERR431366.k63_cutoff3_shuf_bin";
    name_file[377] = "test_data/ERR431367.k63_cutoff3_shuf_bin";
    name_file[378] = "test_data/ERR431368.k63_cutoff3_shuf_bin";
    name_file[379] = "test_data/ERR431369.k63_cutoff3_shuf_bin";
    name_file[380] = "test_data/ERR431370.k63_cutoff3_shuf_bin";
    name_file[381] = "test_data/ERR431371.k63_cutoff3_shuf_bin";
    name_file[382] = "test_data/ERR431372.k63_cutoff3_shuf_bin";
    name_file[383] = "test_data/ERR431373.k63_cutoff3_shuf_bin";
    name_file[384] = "test_data/ERR431374.k63_cutoff3_shuf_bin";
    name_file[385] = "test_data/ERR431375.k63_cutoff3_shuf_bin";
    name_file[386] = "test_data/ERR431376.k63_cutoff3_shuf_bin";
    name_file[387] = "test_data/ERR431377.k63_cutoff3_shuf_bin";
    name_file[388] = "test_data/ERR431378.k63_cutoff3_shuf_bin";
    name_file[389] = "test_data/ERR431379.k63_cutoff3_shuf_bin";
    name_file[390] = "test_data/ERR431380.k63_cutoff3_shuf_bin";
    name_file[391] = "test_data/ERR431381.k63_cutoff3_shuf_bin";
    name_file[392] = "test_data/ERR431382.k63_cutoff3_shuf_bin";
    name_file[393] = "test_data/ERR431383.k63_cutoff3_shuf_bin";
    name_file[394] = "test_data/ERR431384.k63_cutoff3_shuf_bin";
    name_file[395] = "test_data/ERR431386.k63_cutoff3_shuf_bin";
    name_file[396] = "test_data/ERR431387.k63_cutoff3_shuf_bin";
    name_file[397] = "test_data/ERR431388.k63_cutoff3_shuf_bin";
    name_file[398] = "test_data/ERR431389.k63_cutoff3_shuf_bin";
    name_file[399] = "test_data/ERR431390.k63_cutoff3_shuf_bin";
    name_file[400] = "test_data/ERR431391.k63_cutoff3_shuf_bin";
    name_file[401] = "test_data/ERR431392.k63_cutoff3_shuf_bin";
    name_file[402] = "test_data/ERR431393.k63_cutoff3_shuf_bin";
    name_file[403] = "test_data/ERR431394.k63_cutoff3_shuf_bin";
    name_file[404] = "test_data/ERR431395.k63_cutoff3_shuf_bin";
    name_file[405] = "test_data/ERR431396.k63_cutoff3_shuf_bin";
    name_file[406] = "test_data/ERR431397.k63_cutoff3_shuf_bin";
    name_file[407] = "test_data/ERR431398.k63_cutoff3_shuf_bin";
    name_file[408] = "test_data/ERR431399.k63_cutoff3_shuf_bin";
    name_file[409] = "test_data/ERR431400.k63_cutoff3_shuf_bin";
    name_file[410] = "test_data/ERR431401.k63_cutoff3_shuf_bin";
    name_file[411] = "test_data/ERR431402.k63_cutoff3_shuf_bin";
    name_file[412] = "test_data/ERR431403.k63_cutoff3_shuf_bin";
    name_file[413] = "test_data/ERR431404.k63_cutoff3_shuf_bin";
    name_file[414] = "test_data/ERR431405.k63_cutoff3_shuf_bin";
    name_file[415] = "test_data/ERR431406.k63_cutoff3_shuf_bin";
    name_file[416] = "test_data/ERR431407.k63_cutoff3_shuf_bin";
    name_file[417] = "test_data/ERR431408.k63_cutoff3_shuf_bin";
    name_file[418] = "test_data/ERR431409.k63_cutoff3_shuf_bin";
    name_file[419] = "test_data/ERR431410.k63_cutoff3_shuf_bin";
    name_file[420] = "test_data/ERR431411.k63_cutoff3_shuf_bin";
    name_file[421] = "test_data/ERR431412.k63_cutoff3_shuf_bin";
    name_file[422] = "test_data/ERR431413.k63_cutoff3_shuf_bin";
    name_file[423] = "test_data/ERR431414.k63_cutoff3_shuf_bin";
    name_file[424] = "test_data/ERR431415.k63_cutoff3_shuf_bin";
    name_file[425] = "test_data/ERR431416.k63_cutoff3_shuf_bin";
    name_file[426] = "test_data/ERR431417.k63_cutoff3_shuf_bin";
    name_file[427] = "test_data/ERR431418.k63_cutoff3_shuf_bin";
    name_file[428] = "test_data/ERR431419.k63_cutoff3_shuf_bin";
    name_file[429] = "test_data/ERR431420.k63_cutoff3_shuf_bin";
    name_file[430] = "test_data/ERR431421.k63_cutoff3_shuf_bin";
    name_file[431] = "test_data/ERR431422.k63_cutoff3_shuf_bin";
    name_file[432] = "test_data/ERR431423.k63_cutoff3_shuf_bin";
    name_file[433] = "test_data/ERR431424.k63_cutoff3_shuf_bin";
    name_file[434] = "test_data/ERR431425.k63_cutoff3_shuf_bin";
    name_file[435] = "test_data/ERR431426.k63_cutoff3_shuf_bin";
    name_file[436] = "test_data/ERR431427.k63_cutoff3_shuf_bin";
    name_file[437] = "test_data/ERR431428.k63_cutoff3_shuf_bin";
    name_file[438] = "test_data/ERR431429.k63_cutoff3_shuf_bin";
    name_file[439] = "test_data/ERR431430.k63_cutoff3_shuf_bin";
    name_file[440] = "test_data/ERR431431.k63_cutoff3_shuf_bin";
    name_file[441] = "test_data/ERR431432.k63_cutoff3_shuf_bin";
    name_file[442] = "test_data/ERR431433.k63_cutoff3_shuf_bin";
    name_file[443] = "test_data/ERR431434.k63_cutoff3_shuf_bin";
    name_file[444] = "test_data/ERR431435.k63_cutoff3_shuf_bin";
    name_file[445] = "test_data/ERR431436.k63_cutoff3_shuf_bin";
    name_file[446] = "test_data/ERR431437.k63_cutoff3_shuf_bin";
    name_file[447] = "test_data/ERR431438.k63_cutoff3_shuf_bin";
    name_file[448] = "test_data/ERR431439.k63_cutoff3_shuf_bin";
    name_file[449] = "test_data/ERR431440.k63_cutoff3_shuf_bin";
    name_file[450] = "test_data/ERR431441.k63_cutoff3_shuf_bin";
    name_file[451] = "test_data/ERR431442.k63_cutoff3_shuf_bin";
    name_file[452] = "test_data/ERR431443.k63_cutoff3_shuf_bin";
    name_file[453] = "test_data/ERR431444.k63_cutoff3_shuf_bin";
    name_file[454] = "test_data/ERR431445.k63_cutoff3_shuf_bin";
    name_file[455] = "test_data/ERR431446.k63_cutoff3_shuf_bin";
    name_file[456] = "test_data/ERR431447.k63_cutoff3_shuf_bin";
    name_file[457] = "test_data/ERR431448.k63_cutoff3_shuf_bin";
    name_file[458] = "test_data/ERR431449.k63_cutoff3_shuf_bin";
    name_file[459] = "test_data/ERR431450.k63_cutoff3_shuf_bin";
    name_file[460] = "test_data/ERR431451.k63_cutoff3_shuf_bin";
    name_file[461] = "test_data/ERR431452.k63_cutoff3_shuf_bin";
    name_file[462] = "test_data/ERR431453.k63_cutoff3_shuf_bin";
    name_file[463] = "test_data/ERR431454.k63_cutoff3_shuf_bin";
    name_file[464] = "test_data/ERR431455.k63_cutoff3_shuf_bin";
    name_file[465] = "test_data/ERR431456.k63_cutoff3_shuf_bin";
    name_file[466] = "test_data/ERR431457.k63_cutoff3_shuf_bin";
    name_file[467] = "test_data/ERR431458.k63_cutoff3_shuf_bin";
    name_file[468] = "test_data/ERR431459.k63_cutoff3_shuf_bin";
    name_file[469] = "test_data/ERR431460.k63_cutoff3_shuf_bin";
    name_file[470] = "test_data/ERR431461.k63_cutoff3_shuf_bin";
    name_file[471] = "test_data/ERR431462.k63_cutoff3_shuf_bin";
    name_file[472] = "test_data/ERR431463.k63_cutoff3_shuf_bin";
    name_file[473] = "test_data/ERR431464.k63_cutoff3_shuf_bin";
    name_file[474] = "test_data/unique_kmers.k63_cutoff3_shuf_bin";*/

    tree = createRoot(name_file, NB_FILE_2_READ);

    //insert_Genomes_from_FASTx(tree, size_kmer, &annot_sorted);
    insert_Genomes_from_KmerCounting(tree, binary_files, size_kmer, &annot_sorted);

    //freeNode(&(tree->node));
    //free(tree->filenames);
    //free(tree);

    return EXIT_SUCCESS;
}

/* ---------------------------------------------------------------------------------------------------------------
*  insertKmers(root, tab_kmers, size_kmers, nb_kmers, id_genome, func_on_types, ann_inf)
*  ---------------------------------------------------------------------------------------------------------------
*  Insert kmers into the tree
*  ---------------------------------------------------------------------------------------------------------------
*  root: pointer on the Root structure of the tree
*  tab_kmers: array of kmers
*  size_kmers: length k of kmers in tab_kmers
*  nb_kmers: number of kmers in tab_kmers
*  id_genome: genome identity to which belongs the kmers in tab_kmers
*  func_on_types: ptr on ptrs_on_func structure, contains information to manipulate CCs field CC->children_type
*  ann_inf: ptr on annotation_inform structure, used to make the transition between reading and modifying an annot
*  ---------------------------------------------------------------------------------------------------------------
*/
void insertKmers(Root* restrict root,
                 uint8_t* restrict tab_kmers,
                 int size_kmers,
                 int nb_kmers,
                 int id_genome,
                 ptrs_on_func* restrict func_on_types,
                 annotation_inform* ann_inf,
                 resultPresence* res,
                 annotation_array_elem* annot_sorted){

    ASSERT_NULL_PTR(root,"insertKmers()")
    ASSERT_NULL_PTR(func_on_types,"insertKmers()")
    ASSERT_NULL_PTR(ann_inf,"insertKmers()")
    ASSERT_NULL_PTR(res,"insertKmers()")

    int i = 0;
    int nb_cell = CEIL(size_kmers*2, SIZE_CELL);

    uint8_t kmer[nb_cell];

    for (i=0; i<nb_kmers; i++){
        memcpy(kmer, &(tab_kmers[i*nb_cell]), nb_cell * sizeof(uint8_t));
        insertKmer_Node(&(root->node), &(root->node), &(tab_kmers[i*nb_cell]), size_kmers, kmer, size_kmers, id_genome, func_on_types, ann_inf, res, annot_sorted);
    }
}

void insert_Genomes_from_KmerCounting(Root* tree, int binary_files, int size_kmer, annotation_array_elem** ptr_annot_sorted){

    ASSERT_NULL_PTR(tree,"insert_Genomes_from_KmerCounting()")
    ASSERT_NULL_PTR(tree->filenames,"insert_Genomes_from_KmerCounting()")
    ASSERT_NULL_PTR(ptr_annot_sorted,"insert_Genomes_from_KmerCounting()")

    FILE* file;

    Pvoid_t PJArray = (PWord_t)NULL;
    Word_t Rc_word;

    time_t start = time(NULL);
    time_t last_start = time(NULL);

    annotation_array_elem* annot_sorted = *ptr_annot_sorted;
    annotation_array_elem* old_annot_sorted = NULL;

    annotation_inform* ann_inf = calloc(1,sizeof(annotation_inform)); //Initialize structure to pass information between reading and modifying an annotation
    ASSERT_NULL_PTR(ann_inf,"insert_Genomes_from_FASTx()")

    resultPresence* res = create_resultPresence();

    ptrs_on_func* func_on_types = create_ptrs_on_func(SIZE_SEED, size_kmer);

    uint8_t* tab_kmers = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(tab_kmers,"insert_Genomes_from_KmerCounting()")

    uint8_t* kmer;

    char* line = calloc(100, sizeof(char));
    ASSERT_NULL_PTR(line,"insert_Genomes_from_KmerCounting()")

    uint16_t** skip_node_root;

    int i = 0;
    int j = 0;
    int k = 0;
    int nb_cell = CEIL(size_kmer*2, SIZE_CELL);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_cell;

    size_t return_fread;
    int size_annot_sorted;
    int size_old_annot_sorted;

    double count = 0;

    uint64_t kmers_read;

    for (i=0; i<tree->nb_genomes; i++){ //For each file in input

        kmers_read = 0;
        k = 0;
        j = 0;

        file = fopen(tree->filenames[i], "r");
        ASSERT_NULL_PTR(file,"insert_Genomes_from_KmerCounting()")

        printf("\nFile %d: %s\n\n", i, tree->filenames[i]);

        if (binary_files){

            if (fgets(line, 100, file) != NULL) k = atoi(line);
            else ERROR("Cannot read header of the file")

            if (fgets(line, 100, file) != NULL) printf("%d %d-mers in the file\n\n", atoi(line), k);
            else ERROR("Cannot read header of the file")

            //if (i != NB_FILE_2_READ-1){

                while ((!ferror(file)) && (!feof(file))){

                    return_fread = fread(tab_kmers, (size_t)nb_cell, (size_t)nb_kmer_in_buf, file);

                    insertKmers(tree, tab_kmers, size_kmer, return_fread, i, func_on_types, ann_inf, res, annot_sorted);

                    memset(tab_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));

                    if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+return_fread)%PRINT_EVERY_X_KMERS)){
                        printf("%" PRIu64 " kmers read\n", kmers_read+return_fread);
                        //break;
                    }

                    kmers_read += return_fread;
                }
            /*}
            else{

                struct timeval tval_before, tval_after, tval_result;
                gettimeofday(&tval_before, NULL);

                int nb_branching_node = 0;
                int count_branching_node = 0;

                skip_node_root = build_skip_nodes(&(tree->node), size_kmer, func_on_types);

                while ((!ferror(file)) && (!feof(file))){

                    return_fread = fread(tab_kmers, (size_t)nb_cell, (size_t)nb_kmer_in_buf, file);

                    for (k=0; k<return_fread; k++){

                        count_branching_node = 0;

                        if (isBranchingRight(&(tree->node), &(tab_kmers[k*nb_cell]), size_kmer, func_on_types, skip_node_root) > 1) count_branching_node++;
                        if (isBranchingLeft(&(tree->node), &(tab_kmers[k*nb_cell]), size_kmer, func_on_types, skip_node_root) > 1) count_branching_node++;
                        if (count_branching_node > 0) nb_branching_node++;
                    }

                    memset(tab_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));
                }

                gettimeofday(&tval_after, NULL);

                tval_result.tv_sec = tval_after.tv_sec - tval_before.tv_sec;
                tval_result.tv_usec = tval_after.tv_usec - tval_before.tv_usec;
                if (tval_result.tv_usec < 0) {
                    --tval_result.tv_sec;
                    tval_result.tv_usec += 1000000;
                }

                printf("Nb complex nodes: %d\n", nb_branching_node);
                printf("Elapsed time for querying complex nodes: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
            }*/
        }
        else {
            while (fgets(line, 100, file) != NULL){

                if (parseKmerCount(line, size_kmer, tab_kmers, k) == 1){
                    k += nb_cell;
                    j++;

                    if (j == nb_kmer_in_buf){
                        insertKmers(tree, tab_kmers, size_kmer, nb_kmer_in_buf, i, func_on_types, ann_inf, res, annot_sorted);

                        j = 0;
                        k = 0;
                        memset(tab_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));

                        if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+nb_kmer_in_buf)%PRINT_EVERY_X_KMERS)){
                            printf("%" PRIu64 " kmers read\n", kmers_read+nb_kmer_in_buf);
                        }

                        kmers_read += nb_kmer_in_buf;
                    }
                }
            }

            insertKmers(tree, tab_kmers, size_kmer, j, i, func_on_types, ann_inf, res, annot_sorted);
            kmers_read += j;

            memset(tab_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));
        }

        fclose(file);

        if ((i > 5) && (i%TRESH_DEL_ANNOT == 0)){

            load_annotation_from_Node(&(tree->node), size_kmer, func_on_types, &PJArray, annot_sorted);

            old_annot_sorted = annot_sorted;
            size_old_annot_sorted = size_annot_sorted;

            annot_sorted = sort_annotations(&PJArray, &size_annot_sorted);

            compress_annotation_from_Node(&(tree->node), size_kmer, func_on_types, &PJArray, old_annot_sorted);

            free_annotation_array_elem(old_annot_sorted, size_old_annot_sorted);

            JSLFA(Rc_word, PJArray);
        }

        //Delete unecessary annotations
        /*skip_node_root = build_skip_nodes(&(tree->node), size_kmer, func_on_types);

        printf("\nCreating marking structure\n");

        create_marking_Node_4states(&(tree->node), size_kmer, func_on_types);

        printf("\nGetting complex nodes\n");

        count = get_nb_cplx_nodes_from_KmerCounting(tree, tree->filenames[i], size_kmer, skip_node_root, &annot_sorted, func_on_types, ann_inf);

        if (count/((double)kmers_read) <= TRESH_DEL_ANNOT){

            kmer = calloc(CEIL(size_kmer*2, SIZE_CELL), sizeof(uint8_t));
            ASSERT_NULL_PTR(kmer,"insert_Genomes_from_KmerCounting()")

            count = deleteColors_from_branchingNodes(&(tree->node), &(tree->node), kmer, size_kmer, 0, 0, size_kmer, i, func_on_types, skip_node_root, ann_inf, annot_sorted);
            printf("\nNumber of annotations that could be erased = %f\n", count);

            if (((double)kmers_read-count)/((double)kmers_read) <= TRESH_DEL_ANNOT){

                count = resize_annotation_Node(&(tree->node), size_kmer, func_on_types);
                printf("\n%f unnecessary annoations has been deleted, BFT has been reallocated\n\n", count);
            }
            else{
                printf("\nFor this genome, memory cannot be optimized by deleting unnecessary annotations in graph's paths\n\n");
                delete_marking_Node_4states(&(tree->node), size_kmer, func_on_types);
            }

            free(kmer);
        }
        else {
            printf("\nFor this genome, memory cannot be optimized by deleting unnecessary annotations in graph's paths\n\n");
            delete_marking_Node_4states(&(tree->node), size_kmer, func_on_types);
        }

        free_skip_nodes(&(tree->node), skip_node_root);

        if (annot_sorted != NULL) free(annot_sorted);*/

        printf("\nElapsed time: %.2f s\n", (double)(time(NULL) - last_start));
        printf("Total elapsed time: %.2f s\n", (double)(time(NULL) - start));
        printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);
        printf("Current memory: %llu mb\n", ((unsigned long long int)getCurrentRSS())/1024);

        last_start = time(NULL);
    }

    memory_Used* mem = printMemoryUsedFromNode(&(tree->node), size_kmer, func_on_types);

    if (skip_node_root != NULL) free_skip_nodes(&(tree->node), skip_node_root);

    printMemory(mem);

    free(mem);

    free(line);
    free(tab_kmers);

    free(func_on_types);
    free(ann_inf);
    free(res);

    return;
}

void insert_Genomes_from_FASTx(Root* tree, int size_kmer, annotation_array_elem** ptr_annot_sorted){

    ASSERT_NULL_PTR(tree,"insert_Genomes_from_FASTx()")
    ASSERT_NULL_PTR(tree->filenames,"insert_Genomes_from_FASTx()")

    int i = 0;
    int size_buf_tmp = 0; //How many characters are stored in buf_tmp
    int nb_kmers_buf = 0;
    int nb_cell_kmer = CEIL(size_kmer*2, SIZE_CELL); //Size of kmers in bytes
    int size_annot_sorted;
    int size_old_annot_sorted;

    annotation_inform* ann_inf = calloc(1,sizeof(annotation_inform)); //Initialize structure to pass information between reading and modifying an annotation
    ASSERT_NULL_PTR(ann_inf,"insert_Genomes_from_FASTx()")

    resultPresence* res = create_resultPresence();

    ptrs_on_func* func_on_types = create_ptrs_on_func(SIZE_SEED, size_kmer);

    Pvoid_t PJArray = (PWord_t)NULL;
    Word_t Rc_word;

    annotation_array_elem* annot_sorted = *ptr_annot_sorted;
    annotation_array_elem* old_annot_sorted = NULL;

    char* buf_tmp = calloc((size_kmer-1)*2, sizeof(char)); //Allocate temporary buffer
    ASSERT_NULL_PTR(buf_tmp,"insert_Genomes_from_FASTx()")

    uint8_t* tab_kmers = calloc(SIZE_BUFFER*nb_cell_kmer, sizeof(uint8_t)); //Allocate buffer for kmers
    ASSERT_NULL_PTR(tab_kmers,"insert_Genomes_from_FASTx()")

    uint64_t kmers_read = 0;
    uint64_t tmp_kmers_read = 0;

    for (i=0; i<tree->nb_genomes; i++){ //For each file in input
        size_buf_tmp = 0;
        kmers_read = 0;
        tmp_kmers_read = 0;
        nb_kmers_buf = 0;

        int fp = open(tree->filenames[i], O_RDONLY); //Open it
        kseq_t *seq = kseq_init(fp); //initialize the parser for this file
        int size_seq = kseq_read(seq, -1); //Start reading file, seq contains a buffer with a part of a sequence from the file

        printf("\nFile : %s\n\n", tree->filenames[i]);

        while (size_seq > -1) { //While the end of the file is not reached

            if (size_seq > 0) size_buf_tmp = 0; //New sequence

            int current_buf_length = seq->seq.l - seq->seq.z; //Number of characters put into the seq buffer

            if (current_buf_length > 0){ //If the seq buffer is not empty

                nb_kmers_buf = MAX(current_buf_length-size_kmer+1, 0); //Number of kmers it is possible to read in seq buffer

                if (size_buf_tmp == 0){ //If the number of characters in the temporary buffer is 0
                    if (nb_kmers_buf != 0){ //If there is at least one kmer in the seq buffer
                        memcpy(buf_tmp, &(seq->seq.s[nb_kmers_buf]), size_kmer-1); //Copy the last size_kmer-1 characters of seq-buffer into buf_tmp
                        size_buf_tmp = (size_kmer-1)/* *2 */;
                    }
                    else{
                        memcpy(buf_tmp, &(seq->seq.s[0]), current_buf_length); //Copy the content of seq buffer into buf_tmp
                        size_buf_tmp = current_buf_length;
                    }
                }
                else { //If the number of characters in the temporary buffer is not 0
                    //Insertion of kmers overlapping the last buffer and the current one (they are in buf_tmp)
                    int size_to_copy = MIN(size_kmer-1, current_buf_length);
                    memcpy(&(buf_tmp[size_buf_tmp]), seq->seq.s, size_to_copy);
                    size_buf_tmp += size_to_copy;

                    int nb_kmers = size_buf_tmp - size_kmer + 1;
                    if (nb_kmers > 0){
                        parseSequenceBuffer(buf_tmp, tab_kmers, &nb_kmers, size_kmer, nb_cell_kmer); //Read buf_tmp, extract the kmers in tab_kmers
                        insertKmers(tree, tab_kmers, size_kmer, nb_kmers, i, func_on_types, ann_inf, res, annot_sorted); //Insert the kmers into the tree
                        memset(tab_kmers, 0, nb_kmers*nb_cell_kmer*sizeof(uint8_t)); //Reinit tab_kmers
                        tmp_kmers_read = nb_kmers;
                    }
                    else tmp_kmers_read = 0;

                    if (nb_kmers_buf != 0){
                        memcpy(buf_tmp, &(seq->seq.s[nb_kmers_buf]), size_kmer-1);
                        size_buf_tmp = size_kmer-1;
                    }
                    else{
                        memcpy(buf_tmp, &(seq->seq.s[0]), current_buf_length);
                        size_buf_tmp = current_buf_length;
                    }
                }

                //Extraction of buffer's kmers. Insertion in the tree.
                if (nb_kmers_buf > 0){
                    parseSequenceBuffer(seq->seq.s, tab_kmers, &nb_kmers_buf, size_kmer, nb_cell_kmer);
                    insertKmers(tree, tab_kmers, size_kmer, nb_kmers_buf, i, func_on_types, ann_inf, res, annot_sorted);
                    memset(tab_kmers, 0, nb_kmers_buf*nb_cell_kmer*sizeof(uint8_t));
                    tmp_kmers_read += nb_kmers_buf;
                }

                //Display how many kmers were read
                if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+tmp_kmers_read)%PRINT_EVERY_X_KMERS))
                    printf("%" PRIu64 " kmers read\n", kmers_read+tmp_kmers_read);

                kmers_read += tmp_kmers_read;
            }

            size_seq = kseq_read(seq, size_seq);
        }

        if ((i > 5) && (i%TRESH_DEL_ANNOT == 0)){

            load_annotation_from_Node(&(tree->node), size_kmer, func_on_types, &PJArray, annot_sorted);

            old_annot_sorted = annot_sorted;
            size_old_annot_sorted = size_annot_sorted;

            annot_sorted = sort_annotations(&PJArray, &size_annot_sorted);

            compress_annotation_from_Node(&(tree->node), size_kmer, func_on_types, &PJArray, old_annot_sorted);

            free_annotation_array_elem(old_annot_sorted, size_old_annot_sorted);

            JSLFA(Rc_word, PJArray);
        }

        kseq_destroy(seq);
        close(fp);
    }

    memory_Used* mem = printMemoryUsedFromNode(&(tree->node), size_kmer, func_on_types);

    printMemory(mem);

    free(mem);

    free(func_on_types);
    free(ann_inf);
    free(res);

    free(buf_tmp);
    free(tab_kmers);

    return;
}

/* ---------------------------------------------------------------------------------------------------------------
*  Functions get_nb_cplx_nodes_from_* are currently in beta test and should not be used.
*  ---------------------------------------------------------------------------------------------------------------
*/

int get_nb_cplx_nodes_from_KmerCounting(Root* tree, char* name_file, int size_kmer, uint16_t** skip_node_root, ptrs_on_func* func_on_types, annotation_inform* ann_inf){

    ASSERT_NULL_PTR(tree,"get_nb_cplx_nodes_from_KmerCounting()")
    ASSERT_NULL_PTR(name_file,"get_nb_cplx_nodes_from_KmerCounting()")
    ASSERT_NULL_PTR(skip_node_root,"get_nb_cplx_nodes_from_KmerCounting()")
    ASSERT_NULL_PTR(func_on_types,"get_nb_cplx_nodes_from_KmerCounting()")
    ASSERT_NULL_PTR(ann_inf,"get_nb_cplx_nodes_from_KmerCounting()")

    resultPresence* res;

    UC* uc;

    uint8_t* tab_kmers = calloc(SIZE_BUFFER, sizeof(uint8_t));
    ASSERT_NULL_PTR(tab_kmers,"get_nb_cplx_nodes_from_KmerCounting()")

    char* line = calloc(100, sizeof(char));
    ASSERT_NULL_PTR(line,"get_nb_cplx_nodes_from_KmerCounting()")

    int j = 0;
    int k = 0;
    int kmers_read = 0;
    int count_branching = 0;
    int count_without_left_n = 0;
    int degree_left = 0;

    int nb_cell = CEIL(size_kmer*2, SIZE_CELL);
    int nb_kmer_in_buf = SIZE_BUFFER/nb_cell;

    FILE* file = fopen(name_file, "r");
    ASSERT_NULL_PTR(file,"get_nb_cplx_nodes_from_KmerCounting()")

    while (fgets(line, 100, file) != NULL){

        if (parseKmerCount(line, size_kmer, tab_kmers, k) == 1){
            k += nb_cell;
            j++;

            if (j == nb_kmer_in_buf){

                for (k = 0; k < nb_kmer_in_buf * nb_cell; k += nb_cell){

                    if (isBranchingRight(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types, skip_node_root) > 1){ //Flag 0

                        count_branching++;
                        res = isKmerPresent(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types);

                        uc = (UC*)res->container;

                        if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 2);
                        else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 2);

                        free(res);
                    }
                    else if ((degree_left = isBranchingLeft(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types, skip_node_root)) != 1){

                        res = isKmerPresent(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types);

                        uc = (UC*)res->container;

                        if (degree_left == 0){

                            count_without_left_n++;

                            if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 3);
                            else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 3);
                        }
                        else{

                            count_branching++;

                            if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 2);
                            else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 2);
                        }

                        free(res);
                    }
                }

                j = 0;
                k = 0;
                memset(tab_kmers, 0, SIZE_BUFFER*sizeof(uint8_t));

                if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+nb_kmer_in_buf)%PRINT_EVERY_X_KMERS)){
                    //break;
                }

                kmers_read += nb_kmer_in_buf;
            }
        }
    }

    for (k = 0; k < j * nb_cell; k += nb_cell){

        if (isBranchingRight(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types, skip_node_root) > 1){ //Flag 0

            count_branching++;
            res = isKmerPresent(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types);

            uc = (UC*)res->container;

            if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 2);
            else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 2);

            free(res);
        }
        else if ((degree_left = isBranchingLeft(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types, skip_node_root)) != 1){

            res = isKmerPresent(&(tree->node), &(tab_kmers[k]), size_kmer, func_on_types);
            uc = (UC*)res->container;

            if (degree_left == 0){

                count_without_left_n++;

                if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 3);
                else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 3);
            }
            else{

                count_branching++;

                if (res->container_is_UC == 1) mark_UC_4states(uc, res->posFilter2, uc->nb_children >> 1, res->pos_sub_bucket, 2);
                else mark_UC_4states(uc, res->posFilter2, uc->nb_children, res->pos_sub_bucket, 2);
            }

            free(res);
        }
    }

    fclose(file);

    printf("\nNumber of branching nodes = %d\n", count_branching);
    printf("Number of nodes with 0 left neighbor and 1 right neighbor = %d\n", count_without_left_n);
    printf("Total = %d\n", count_branching + count_without_left_n);

    free(line);
    free(tab_kmers);

    return count_branching /*+ count_without_left_n*/;
}

Root* get_nb_cplx_nodes_from_FASTx(Root* tree, int size_kmer, uint16_t** skip_node_root, annotation_array_elem** ptr_annot_sorted,
                                   ptrs_on_func* func_on_types, annotation_inform* ann_inf, resultPresence* res){

    ASSERT_NULL_PTR(tree,"get_nb_cplx_nodes_from_FASTx()")
    ASSERT_NULL_PTR(tree->filenames,"get_nb_cplx_nodes_from_FASTx()")
    ASSERT_NULL_PTR(skip_node_root,"get_nb_cplx_nodes_from_FASTx()")
    ASSERT_NULL_PTR(func_on_types,"get_nb_cplx_nodes_from_FASTx()")
    ASSERT_NULL_PTR(ann_inf,"get_nb_cplx_nodes_from_FASTx()")
    ASSERT_NULL_PTR(res,"get_nb_cplx_nodes_from_FASTx()")

    Root* bft_cplx_nodes = createRoot(NULL,0);

    int i = 0, j = 0;
    int size_buf_tmp = 0; //How many characters are stored in buf_tmp
    int nb_kmers_buf = 0;
    int nb_cell_kmer = CEIL(size_kmer*2, SIZE_CELL); //Size of kmers in bytes

    char* buf_tmp = calloc((size_kmer-1)*2, sizeof(char)); //Allocate temporary buffer
    ASSERT_NULL_PTR(buf_tmp,"get_nb_cplx_nodes_from_FASTx()")

    uint8_t* tab_kmers = calloc(SIZE_BUFFER*nb_cell_kmer, sizeof(uint8_t)); //Allocate buffer for kmers
    ASSERT_NULL_PTR(tab_kmers,"get_nb_cplx_nodes_from_FASTx()")

    uint64_t kmers_read = 0;
    uint64_t tmp_kmers_read = 0;

    for (i=0; i<tree->nb_genomes; i++){ //For each file in input
        size_buf_tmp = 0;
        kmers_read = 0;
        tmp_kmers_read = 0;
        nb_kmers_buf = 0;

        int fp = open(tree->filenames[i], O_RDONLY); //Open it
        kseq_t *seq = kseq_init(fp); //initialize the parser for this file
        int size_seq = kseq_read(seq, -1); //Start reading file, seq contains a buffer with a part of a sequence from the file

        while (size_seq > -1) { //While the end of the file is not reached

            if (size_seq > 0) size_buf_tmp = 0; //New sequence

            int current_buf_length = seq->seq.l - seq->seq.z; //Number of characters put into the seq buffer

            if (current_buf_length > 0){ //If the seq buffer is not empty

                nb_kmers_buf = MAX(current_buf_length-size_kmer+1, 0); //Number of kmers it is possible to read in seq buffer

                if (size_buf_tmp == 0){ //If the number of characters in the temporary buffer is 0
                    if (nb_kmers_buf != 0){ //If there is at least one kmer in the seq buffer
                        memcpy(buf_tmp, &(seq->seq.s[nb_kmers_buf]), size_kmer-1); //Copy the last size_kmer-1 characters of seq-buffer into buf_tmp
                        size_buf_tmp = size_kmer-1;
                    }
                    else{
                        memcpy(buf_tmp, &(seq->seq.s[0]), current_buf_length); //Copy the content of seq buffer into buf_tmp
                        size_buf_tmp = current_buf_length;
                    }
                }
                else { //If the number of characters in the temporary buffer is not 0
                    //Insertion of kmers overlapping the last buffer and the current one (they are in buf_tmp)
                    int size_to_copy = MIN(size_kmer-1, current_buf_length);
                    memcpy(&(buf_tmp[size_buf_tmp]), seq->seq.s, size_to_copy);
                    size_buf_tmp += size_to_copy;

                    int nb_kmers = size_buf_tmp - size_kmer + 1;
                    if (nb_kmers > 0){
                        parseSequenceBuffer(buf_tmp, tab_kmers, &nb_kmers, size_kmer, nb_cell_kmer); //Read buf_tmp, extract the kmers in tab_kmers

                        for (j=0; j < nb_kmers * nb_cell_kmer; j += nb_cell_kmer){
                            if ((isBranchingRight(&(tree->node), &(tab_kmers[j]), size_kmer, func_on_types, skip_node_root) > 1) ||
                                (isBranchingLeft(&(tree->node), &(tab_kmers[j]), size_kmer, func_on_types, skip_node_root) != 1))
                                   insertKmers(bft_cplx_nodes, &(tab_kmers[j]), size_kmer, 1, 0, func_on_types, ann_inf, res, *ptr_annot_sorted);
                        }

                        memset(tab_kmers, 0, nb_kmers*nb_cell_kmer*sizeof(uint8_t)); //Reinit tab_kmers
                        tmp_kmers_read = nb_kmers;
                    }
                    else tmp_kmers_read = 0;

                    if (nb_kmers_buf != 0){
                        memcpy(buf_tmp, &(seq->seq.s[nb_kmers_buf]), size_kmer-1);
                        size_buf_tmp = size_kmer-1;
                    }
                    else{
                        memcpy(buf_tmp, &(seq->seq.s[0]), current_buf_length);
                        size_buf_tmp = current_buf_length;
                    }
                }

                //Extraction of buffer's kmers. Insertion in the tree.
                if (nb_kmers_buf > 0){
                    parseSequenceBuffer(seq->seq.s, tab_kmers, &nb_kmers_buf, size_kmer, nb_cell_kmer);

                    for (j=0; j< nb_kmers_buf * nb_cell_kmer; j += nb_cell_kmer){
                        if ((isBranchingRight(&(tree->node), &(tab_kmers[j]), size_kmer, func_on_types, skip_node_root) > 1) ||
                            (isBranchingLeft(&(tree->node), &(tab_kmers[j]), size_kmer, func_on_types, skip_node_root) != 1))
                                insertKmers(bft_cplx_nodes, &(tab_kmers[j]), size_kmer, 1, 0, func_on_types, ann_inf, res, *ptr_annot_sorted);
                    }

                    memset(tab_kmers, 0, nb_kmers_buf*nb_cell_kmer*sizeof(uint8_t));
                    tmp_kmers_read += nb_kmers_buf;
                }

                //Display how many kmers were read
                if ((kmers_read%PRINT_EVERY_X_KMERS) > ((kmers_read+tmp_kmers_read)%PRINT_EVERY_X_KMERS)){
                    printf("%" PRIu64 " kmers read\n", kmers_read+tmp_kmers_read);
                    //break;
                }

                kmers_read += tmp_kmers_read;
            }

            size_seq = kseq_read(seq, size_seq);
        }

        kseq_destroy(seq);
        close(fp);
    }

    free(func_on_types);
    free(buf_tmp);
    free(tab_kmers);

    return bft_cplx_nodes;
}
