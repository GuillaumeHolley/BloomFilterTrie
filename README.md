# BFT: Bloom Filter Trie

This repository contains the source code of the Bloom Filter Trie (BFT) library. The BFT is an alignment-free, reference-free and incremental succinct data structure for colored de Bruijn graphs. It is based on the burst trie and use Bloom filters for efficient trie and graph traversals. The data structure indexes *k*-mers and their colors based on a new representation of trie vertices that
compress and index shared substrings. A typical application of the BFT is pan-genome indexing.

## Dependencies

In order to compile and use the BFT, you need a machine running a 64 bits Linux or Mac OS operating system, POSIX compliant and accepting SSE4 instructions.

The library depends on two external libraries: Judy Arrays (http://judy.sourceforge.net) and Jemalloc (http://www.canonware.com/jemalloc).

Both can be downloaded and installed by following the instructions on their respective websites. It is however most likely that at least few of them are available via a package manager for your operating system.

If you operating system is Ubuntu or Debian:
```
sudo apt-get install libjudydebian1 libjudy-dev libjemalloc1 libjemalloc-dev
```

If you operating system is Mac OS, Jemalloc can be easily downloaded and installed via Homebrew:
```
brew install jemalloc
```

## Compilation and installation

The library compiles with GNU GCC and G++ (compatibility with Clang is in preparation). It successfully compiles and runs on Ubuntu 14.04 and 15.04.

### Linux

On Linux, you can verify the presence of gcc and g++ on your system with:
```
gcc -v
g++ -v
```

If not present (unlikely), they can be installed for Ubuntu or Debian with:
```
sudo apt-get install build-essential
```

Compiling the library should then be as simple as:
```
cd <BFT_directory>
./configure
make
make install
```

You can also install it in a specific directory (for example because you are not root on the machine you are using) with:
```
cd <BFT_directory>
./configure --prefix=<a_directory>
make
make install
```

Make sure that your environment variables are all set with <a_directory>.

### Mac OS

For Mac OS, you will need the "real" GCC and G++, not the Clang interface that is called when you use GCC or G++. Both can be installed via Homebrew:
```
brew install gcc-x
brew install g++-x
```

in which *x* is the latest major version of GCC and G++.

Compiling the library should then be as simple as:
```
cd <BFT_directory>
./configure CC=gcc-x
make
make install
```

in which *x* is the version of GCC and G++ that you installed via Homebrew or other.

To install the BFT library in a specific directory, see Linux compilation and installation.

## API Usage

Using the BFT library is very simple. Once the library is installed on your system, just use
```
#include<bft/bft.h>
```
in your C or C++ code. Then, compile your code with the flag
```
-lbft
```

## API documentation:

Documentation for the BFT library is available in the /doc/doxygen folder (HTML).

The following command regenerates the documentation:
```
cd <BFT_directory>
doxygen Doxyfile
```

The documentation contains a description of all the functions and structures of the library as well as code snippets.

## Binary usage:

Installing the BFT library also produces a binary that shows what it is possible to do with the library. Therefore, the binary can perform a limited number of operations described in the following.

```
bft build k treshold_compression {kmers|kmers_comp} list_genome_files output_file [Options]
bft load file_bft [-add_genomes {kmers|kmers_comp} list_genome_files output_file] [Options]

Options:
[-query_kmers {kmers|kmers_comp} list_kmer_files]
[-query_branching {kmers|kmers_comp} list_kmer_files]
[-extract_kmers {kmers|kmers_comp} kmers_file]

Version:
bft --version
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
@article{holley2016bloom,
  title="{Bloom Filter Trie: an alignment-free and reference-free data structure for pan-genome storage}",
  author={Holley, Guillaume and Wittler, Roland and Stoye, Jens},
  journal={Algorithm. Mol. Biol.},
  volume={11},
  number={1},
  pages={1},
  year={2016}
}
```

## Contact

For any question, feedback or problem, please contact me at gholley[At]cebitec[D0t]uni-bielefeld[D0t]de
