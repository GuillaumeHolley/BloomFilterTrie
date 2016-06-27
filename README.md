# BloomFilterTrie

This repository contains the source code of the Bloom Filter Trie (BFT) version 0.6 for Linux and Mac OS X systems.

The BFT library depends on two external libraries that are required to compile and use the library: Judy Arrays (http://judy.sourceforge.net) and Jemalloc (http://www.canonware.com/jemalloc). Both can be downloaded and installed by following instructions on their respective websites.

Both can be easily downloaded and installed on Ubuntu/Debian via their packages:
```
sudo apt-get install libjemalloc1 libjemalloc-dev
sudo apt-get install libjudydebian1 libjudy-dev
```

For Mac OS X, Jemalloc can be downloaded and installed via Homebrew:
```
brew install jemalloc
```

The BFT library compiles with GNU GCC and G++. Compatibility with Clang LLVM is in preparation.

On Ubuntu/Debian, you can verify their presences on your system with:
```
gcc -v
g++ -v
```

If not present, they can be installed with:
```
sudo apt-get install build-essential
```

For Mac OS X, GCC and G++ can be installed via Homebrew:
```
brew install gcc-5
brew install g++-5
```

The code can be then compiled with:
```
cd <BFT directory>
make
```

## API documentation:

Documentation for the BFT library is available in the /doc/doxygen folder (HTML and Latex).

The following command regenerates the documentation:
```
cd <BFT directory>
doxygen Doxyfile
```

## Binary usage:
```
./bft build k treshold_compression {kmers|kmers_comp} list_genome_files output_file [Options]
./bft load file_bft [-add_genomes {kmers|kmers_comp} list_genome_files output_file] [Options]

Options:
[-query_kmers {kmers|kmers_comp} list_kmer_files]
[-query_branching {kmers|kmers_comp} list_kmer_files]
[-extract_kmers {kmers|kmers_comp} kmers_file]

Version:
./bft --version
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
* **-extract_kmers** extracts the *k*-mers stored in the BFT and writes them to a *k*-mers file named *kmers_file* (see below for input file types).

New options will be available soon.

### I/O file types

* *kmers*: files are *k*-mers files. Each file contains one *k*-mer (plain text) per line, eventually followed by a count.
* *kmers_comp*: files of *list_genome_files* are compressed *k*-mers files. Each file is built as the following: First line is *k* (plain text), second line is the number of *k*-mers in the file (plain text) and third line is the concatenation of all compressed *k*-mers. A compressed *k*-mer is encoded with two bits per nucleotid (one byte for 4 nucleotids) with A=00, C=01, G=10 and T=11. A byte is always encoded from the Less Significant Bit to the Most Significant Bit. If a byte cannot be entirely filled in with nucleotids, it is padded with 0s.
Example: ACTTGTCTG -> 11110100 11011110 00000010

## Citation

If you want to cite the Bloom Filter Trie, please use:
```
@inproceedings{holley2015bloom,
  title="{Bloom Filter Trie--A Data Structure for Pan-Genome Storage}",
  author={Holley, Guillaume and Wittler, Roland and Stoye, Jens},
  booktitle={Proceedings of 15th International Workshop on Algorithms in Bioinformatics},
  volume={9289},
  pages={217-230},
  year={2015},
  publisher={Springer}
}
```

## Contact

For any question, feedback or problem, please contact me at gholley{at}cebitec{dot}uni-bielefeld{dot}de
