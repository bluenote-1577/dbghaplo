# dbghaplo - long-read haplotype phasing for diverse small sequences

## Introduction

### Inputs

## Install + Quick start 

#### Option 1 - compile from scratch

A relatively recent standard toolchain is needed.

1. [rust](https://www.rust-lang.org/tools/install) **version > 1.63.0** and associated tools such as cargo are required and assumed to be in PATH.
2. [cmake](https://cmake.org/download/) **version > 3.12** is required. It's sufficient to download the binary from the link and do `PATH="/path/to/cmake-3.xx.x-linux-x86_64/bin/:$PATH"` before installation. 
3. make 
4. GCC 

If you're using an **x86-64 architecture with SSE instructions (most linux systems)**: 

```sh
git clone https://github.com/bluenote-1577/floria
cd floria

cargo install --path . 
floria -h # binary is available in PATH
```

If you're using an **ARM architecture with NEON instructions** (e.g. Mac M1): 

```sh

# If using ARM architecture with NEON instructions
cargo install --path . --root ~/.cargo --features=neon --no-default-features
floria -h # binary is available in PATH

```

#### Option 2 - bioconda

```sh
conda install -c bioconda floria
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
