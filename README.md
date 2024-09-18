# dbghaplo - long-read haplotype phasing for diverse small sequences

## Introduction

**dbghaplo** is a method that separates long-read sequencing of a mixture of sequences into individual haplotypes. This is called "phasing" or "haplotyping".

dbghaplo works on both PacBio and Nanopore. It is a "local haplotyping" method, so it works best when the sequence-of-interest is approximately the size of the reads. 

Example use cases include:

* mixed viral long-read samples (e.g. co-infections)
* deconvolving amplicon/enriched sequencing of specific genes
* haplotyping small sections of multi-strain bacterial communities

### Why dbghaplo?

Similar tools exist for "viral quasispecies" detection. dbghaplo was developed to fill the following gaps:

* **Speed**: dbghaplo can haplotype 

## Install + Quick start 

#### Option 1 - bioconda

```sh
conda install -c bioconda floria
```

#### Option 2 - compile from scratch

1. [rust](https://www.rust-lang.org/tools/install) **version > 1.63.0** and associated tools such as cargo are required and assumed to be in PATH.
2. [cmake](https://cmake.org/download/) **version > 3.12** is required. It's sufficient to download the binary from the link and do `PATH="/path/to/cmake-3.xx.x-linux-x86_64/bin/:$PATH"` before installation. 
3. make 
4. GCC 

If you're using an **x86-64 architecture with SSE instructions (most linux systems)**: 

```sh
git clone https://github.com/bluenote-1577/dbghaplo
cd dbghaplo

cargo install --path . 
dbghap -h # binary is available in PATH
```

If you're using an **ARM architecture with NEON instructions** (e.g. Mac M1): 

```sh

# If using ARM architecture with NEON instructions
cargo install --path . --root ~/.cargo --features=neon --no-default-features
dbghap -h # binary is available in PATH

```

#### Option 3 - precompiled static binary on **x86-64-linux**

The static binary is only for x86-64 linux with SSE instructions currently. 

```sh
wget https://github.com/bluenote-1577/floria/releases/download/latest/floria
chmod +x floria
./floria -h
```

### Quick Start after install 

```sh
dbghap -b mapped_reads.bam -r reference.fa -v variants.vcf
```
