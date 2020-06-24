Contents
====

[![Build Status](https://travis-ci.org/chungongyu/siga.svg?branch=master)](https://travis-ci.org/chungongyu/siga)
[![GitHub release](https://img.shields.io/github/release/chungongyu/siga.svg)](https://github.com/chungongyu/siga/releases)
[![GitHub license](https://img.shields.io/github/license/chungongyu/siga.svg)](https://github.com/chungongyu/siga)
[![GitHub issues](https://img.shields.io/github/issues/chungongyu/siga.svg)](https://github.com/chungongyu/siga/issues)

* [Overview](#overview)
* [Quick Start](#quick-start)
	* [Dependencies](#dependencies)
	* [Compiling SIGA](#compiling-siga)
	* [Installing SIGA](#installing-siga)
	* [Running SIGA](#running-siga)
* [Benchmarks](#benchmarks)
* [Citation](#citation)
* [Publications](#publications)
* [FAQ](#faq)
* [Support](#support)
* [Authors](#authors)

Overview
========

SIGA is an open-source _de_ novo assembly toolkit containing various assembly pipelines which are totally compatible with [SGA](https://github.com/jts/sga) file format. Similar to [SGA](https://github.com/jts/sga), it is designed as a [modular set of programs](#running-siga), together which form an assembly pipeline.

Quick Start
===========

### Dependencies

* c++ compiler that supports [OpenMP](http://www.openmp.org) and [c++ 11](https://en.wikipedia.org/wiki/C%2B%2B11) such as [gcc](http://gcc.gnu.org).
* [autoconf](http://www.gnu.org/software/autoconf)
* [automake](http://www.gnu.org/software/automake)
* [boost](https://www.boost.org/)
* [log4cxx](https://logging.apache.org/log4cxx)
* [rapidjson](https://github.com/Tencent/rapidjson)
* [gperftools](https://github.com/gperftools/gperftools) (optional but suggested)

Dependencies may be installed using the package manager [Homebrew](https://homebrew.sh) on macOS and [Linxubrew](http://linuxbrew.sh) on Linux and Windows, using [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/).

### Compiling SIGA

If you cloned the repository from github, run `autogen.sh` from the root directory 
to generate the configure file:

```
./autogen.sh
```

If the dependencies have been installed in standard locations (like `/usr/local`) you
can run `configure` without any parameters then run `make`:

```
./configure
make
```

### Installing SIGA

Running make install will install siga into `/usr/local/bin/` by default. To specify the install
location use the --prefix option to configure:

```
./configure --prefix=/home/chaunceyyu/ && make && make install
```

This command will copy siga to `/home/chaunceyyu/bin/siga`

### Running SIGA

SIGA consists of a number of subprograms, together which form the assembly pipeline. The subprograms
can also be used to perform other interesting tasks, like read error correction or removing PCR duplicates.
Each program and subprogram will print a brief description and its usage instructions if the -h or --help 
flag is used.

To get a listing of all subprograms, run `siga --help`.

Examples of an SIGA assembly are provided in the [examples](https://github.com/chungongyu/siga/tree/master/examples) directory. It is suggested to look at these examples to become familiar with the flow of data through the program.

The major subprograms are:

```
siga preprocess
```

Prepare reads for assembly. It can perform optional quality filtering/trimming. By default
it will discard reads that have uncalled bases ('N' or '.'). It is mandatory to run this command 
on real data. 

If your reads are paired, the --pe-mode option should be specified. The paired reads can be input in two 
files (siga preprocess READS1 READS2) where the first read in READS1 is paired with the first read on READS2 
and so on. Alternatively, they can be specified in a single file where the two reads are expected to appear 
in consecutive records. By default, output is written to stdout.

```
siga index READS
```

Build the FM-index for READS, which is a fasta or fastq file. 
This program is threaded (-t N).

```
siga correct READS
```

Perform error correction on READS file. Overlap and kmer-based correction algorithms
are implemented. By default, a k-mer based correction is performed. 

Many options exist for this program, see --help. 

```
siga overlap -m N READS
```

Find overlaps between reads to construct the string graph. The -m parameter specifies
the minimum length of the overlaps to find. By default only non-transitive (irreducible) edges are output and edges
between identical sequences. If all overlaps between reads are desired, the --exhaustive option can be specified.
This program is threaded. The output file is READS.asqg.gz by default.

```
siga assemble READS.asqg.gz
```

Assemble takes the output of the overlap step and constructs contigs. The output is in contigs.fa by default. Options
exist for cleaning the graph before assembly which will substantially increase assembly continuity. 
See the --cut-terminal, --bubble, --resolve-small options.

Detail usage information for each command is printed from the --help option. For example, this command will print the 
options for the index subprogram:

```
siga index --help
```

Benchmarks
========

Comming soon ...

Citation
========

## [SIGA](https://github.com/chungongyu/siga)

Chungong Yu, Yu Lin, Guozheng Wei, Bing Wang, Yanbo Li and Dongbo Bu. **SIGA** : An open-source Sensitive and Intelligent *de* novo Genome Assembler.

Publications
============

None yet!   ^_~

FAQ
====

1. **Where can I get further help or advice?**

	See [Support](#support) or the [siga wiki](https://github.com/chungongyu/siga/wiki).

2. **What parameters should I tune to improve my assembly?**
    
	See [here](https://github.com/chungongyu/siga/wiki/SIGA-parameter-tuning)

Support
=======

[Create a new issue](https://github.com/chungongyu/siga/issues) on GitHub.

Contact [siga@ict.ac.cn](mailto:siga@ict.ac.cn)

Authors
=======

Written by **[Chungong Yu](http://bioinfo.ict.ac.cn/~yuchungong)**.<br>
The algorithms were developed by **[Chungong Yu](http://bioinfo.ict.ac.cn/~yuchungong)**, **[Yu Lin](https://cecs.anu.edu.au/people/yu-lin)**, **[Guozheng Wei](http://bioinfo.ict.ac.cn/~weiguozheng)**, **[Bing Wang](https://jp.linkedin.com/in/bing-wang-155128151)**,
Yanbo Li, **[Bin Huang](https://github.com/huangbinapple)**, **[Shiwei Sun](http://bioinfo.ict.ac.cn/~dwsun)** and **[Dongbo Bu](http://bioinfo.ict.ac.cn/~dbu)**.
