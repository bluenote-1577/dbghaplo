# dbghaplo - long-read haplotype phasing for diverse small sequences

## Introduction

**dbghaplo** is a method that separates long-read sequencing (Nanopore or PacBio) of a mixture of sequences into individual haplotypes. This is called "phasing" or "haplotyping".

dbghaplo is a "local haplotyping" method, so it works best when the sequence-of-interest is approximately the size of the reads. 

Example use cases:

* mixed viral long-read samples (e.g. co-infections)
* deconvolving amplicon/enriched sequencing of specific genes
* haplotyping small sections of multi-strain bacterial communities

![dbghap-fig1-github](https://github.com/user-attachments/assets/c0a82bb5-7feb-4d13-ab59-04da2bce52b3)

![dbghap-github-fig2](https://github.com/user-attachments/assets/34cb8bcf-8f23-47e4-b2f6-8515a21d3cf4)

### Why dbghaplo?

Similar tools exist for "viral quasispecies" detection. dbghaplo was developed to fill the following gaps:

* **Speed and low-memory** - dbghaplo can haplotype > 30,000x coverage of ~1kb gene in a few minutes on a laptop.
* **High heterogeneity** - dbghaplo uses a de Bruijn Graph approach, which scales up to very diverse samples.
* **Ease-of-use** - conda installable, engineered in rust, simple command line. 

## Install + Quick start 

#### Option 1 - bioconda

```sh
conda install -c bioconda dbghaplo
```

#### Option 2 - compile from scratch

1. [rust](https://www.rust-lang.org/tools/install) **version > 1.63.0** and associated tools such as cargo are required and assumed to be in PATH.
2. [cmake](https://cmake.org/download/) **version > 3.12** is required. It's sufficient to download the binary from the link and do `PATH="/path/to/cmake-3.xx.x-linux-x86_64/bin/:$PATH"` before installation. 
3. make 

Optional:

1. minimap2
2. lofreq
3. tabix (bcftools)

If you're using an **x86-64 architecture with SSE instructions (most linux systems)**: 

```sh
git clone https://github.com/bluenote-1577/dbghaplo
cd dbghaplo

cargo install --path . 
dbghaplo -h # binary is available in PATH
```

If you're using an **ARM architecture with NEON instructions** (e.g. Mac M1): 

```sh

# If using ARM architecture with NEON instructions
cargo install --path . --root ~/.cargo --features=neon --no-default-features
dbghaplo -h # binary is available in PATH

```

#### Option 3 - precompiled static binary on **x86-64-linux**

The static binary is only for x86-64 linux with SSE instructions currently. 

```sh
wget https://github.com/bluenote-1577/floria/releases/download/latest/dbghaplo
chmod +x dbghaplo
./dbghaplo -h
```

### Quick Start after install 

```sh
```
