# BloomFilterTrie

This repository contains the source code of the Bloom Filter Trie (BFT) version 0.1.2.

To compile the code, Judy's Array and Jemalloc must be installed on your machine. For Ubuntu, they corresponds to packages *libjemalloc1* and *libjudydebian1*.

To compile the code:
```
cd <BFT directory>
make
```

### Usage:
```
./bft build k treshold_compression {fastx|kmers|kmers_comp} list_genome_files output_file [Options]
./bft load file_bft [-add_genomes {fastx|kmers|kmers_comp} list_genome_files output_file] [Options]

Options:
[-query_kmers {kmers|kmers_comp} list_kmer_files]
[-query_branching {kmers|kmers_comp} list_kmer_files]
```
### Commands

Command **build** creates the BFT for the files listed in *list_genome_files* and writes the BFT in file *output_file*.

* *k*: length of *k*-mers
* *treshold_compression*: number of genomes inserted that triggers a compression of the BFT's colors sets. This compression step will be triggered every *treshold_compression* genomes inserted and on the last genome inserted but does not start before insertion of 7 genomes. A *treshold_compression* equals to 0 means no compression of the BFT's colors sets. Example: if *treshold_compression* = 2 and 10 genomes have to be inserted, the compression step will be triggered after insertion of the 8th and 10th genomes.
* *list_genome_files*: file that contains a list of files (one path and name per line) to be inserted in the BFT.
* *output_file*: file where to write the BFT.

Command **load** loads a BFT from file *file_bft*.

* *file_bft*: file that contains a BFT

### Options

* **-add_genomes** adds the genomes listed in *list_genome_files* to the BFT stored in *file_bft*, the new BFT is written in *output_file*
* **-query_kmers** queries the BFT for *k*-mers written in the files of *list_kmer_files*. For each file of *list_kmer_files* is output a CSV file: columns are the genomes represented in the BFT, rows are the queried *k*-mers, the intersection of a column and a row is a binary value indicating if the *k*-mer represented by the row is present in the genome represented by the column.
* **-query_branching** queries the BFT for the number of *k*-mers written in the files of *list_kmer_files* that are branching in the colored de-Bruijn graph represented by the BFT.

New options will be available soon.

### Input file types

* *fastx*: files of *list_genome_files* are FASTA/FASTQ files. All *k*-mers are extracted from the files and inserted in the BFT.
* *kmers*: files of *list_genome_files* are *k*-mers files. Each file contains one *k*-mer (plain text) per line, eventually followed by a count.
* *kmers_comp*: files of *list_genome_files* are compressed *k*-mers files. Each file is built as the following: First line is *k* (plain text), second line is the number of *k*-mers in the file (plain text) and third line is the concatenation of all compressed *k*-mers. A compressed *k*-mer is encoded with two bits per nucleotid (one byte for 4 nucleotids) with A=00, C=01, G=10 and T=11. A byte is always encoded from the Less Significant Bit to the Most Significant Bit. If a byte cannot be entirely filled in with nucleotids, it is padded with 0s.
Example: ACTTGTCTG -> 11110100 11011110 00000010


### Contact

For any question, feedback or problem, please contact me at gholley@cebitec.uni-bielefeld.de
