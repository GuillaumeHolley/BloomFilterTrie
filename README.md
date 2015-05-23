# BloomFilterTrie

This repository contains the source code of the Bloom Filter Trie (BFT) version 0.1. A readme file and a manual will be available in few days.

To compile the code, Judy's Array and Jemalloc must be installed on your machine. For Ubuntu, they corresponds to packages libjemalloc1 and libjudydebian1.

To compile the code:
```
cd <BFT directory>
make
```

Usage:
```
./bft build k {fastx|kmers|kmers_comp} list_genome_files output_file [Options]
./bft load file_bft [-add_genomes {fastx|kmers|kmers_comp} list_genome_files output_file] [Options]

Options:
[-query_kmers {kmers|kmers_comp} list_kmer_files]
[-query_branching {kmers|kmers_comp} list_kmer_files]
```

*fastx*: files of *list_genome_files* are FASTA/FASTQ files. All *k*-mers are extracted from the files and inserted in the BFT.
*kmers*: files of *list_genome_files* are *k*-mers files. Each file contains one *k*-mer (plein text) per line, eventually followed by a count.
*kmers_comp*: files of *list_genome_files* are compressed *k*-mers files. Each file is built as the following: First line is *k* (plain text), second line is the number of *k*-mers in the file and third line is the concatenation of all compressed *k*-mers. A compressed *k*-mer is encoded with two bits per nucleotid (one byte for 4 nucleotids). A byte is always encoded from the Less Significant Bit to the Most Significant Bit. If a bytes cannot be entirely filled in with nucleotids, it is padded with 0.
Example: ACTTGTCTG -> 11110100 11011110 00000010

Command **build** creates the BFT for the files listed in *list_genome_files* and write the BFT in file *output_file*.

*k*: length of *k*-mers
*list_genome_files*: file that contains a list of files (one path and name per line) to be inserted in the BFT.
*output_file*: file where to write the BFT.

Command **load** load a BFT from file *file_bft*.

*file_bft*: file that contains a BFT

Options:

**-add_genomes** adds the genomes listed in *list_genome_files* to the BFT stored in *file_bft*, the new BFT is written in *output_file*
**-query_kmers** queries the BFT for the number of *k*-mers written in file *list_kmer_files* that are present in the BFT.
**-query_branching** queries the BFT for the number of branching *k*-mers written in file *list_kmer_files* that are present in the BFT.

New options will be available soon.

For any question, feedback or problem, please contact me at gholley[AT]cebitec.uni-bielefeld.de
