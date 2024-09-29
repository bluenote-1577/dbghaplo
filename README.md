# dbghaplo - long-read haplotypes from mixtures of "small" sequences

**dbghaplo** is a method that separates long reads (Nanopore or PacBio) of a mixture of sequences into groups with similar alleles. This is called "phasing" or "haplotyping". 

dbghaplo is a "local haplotyping" method, so it works best when the sequence-of-interest is approximately the size of the reads. 

### Example use cases:

* mixed viral long-read samples (e.g. co-infections)
* deconvolving amplicon/enriched sequencing of specific genes
* haplotyping small sections of multi-strain bacterial communities

<p align="center">
  <img width="600" height="200" src="https://github.com/user-attachments/assets/c0a82bb5-7feb-4d13-ab59-04da2bce52b3", caption="asdf">
</p>
<p align="center">
  <i>
High-depth, heterogeneous sequencing that spans a 1kb gene.
  </i>
</p>

<p align="center">
  <img width="600" height="200" src="https://github.com/user-attachments/assets/34cb8bcf-8f23-47e4-b2f6-8515a21d3cf4", caption="asdf">
</p>
<p align="center">
  <i>
Separated groups ("haplotypes") after running dbghaplo.
  </i>
</p>

### Why dbghaplo?

Similar tools exist for detection of similar haplotypes in mixtures. dbghaplo was developed to fill the following gaps:

* **Speed and low-memory** - dbghaplo scales approximately linearly with sequencing depth and # of SNPs. > 30,000x coverage genes can be haplotyped in minutes. 
* **High heterogeneity and coverage** - dbghaplo uses a de Bruijn Graph approach, which works with very diverse samples (> 10 haplotypes)
* **Ease-of-use + interpretable outputs** - conda installable, engineered in rust, simple command line. Outputs are easy to interpret (haplotagged BAM or MSA). 

## Install + Quick start 

#### Option 1 - bioconda

```sh
conda install -c bioconda dbghaplo
```

#### Option 2 - compile from scratch

1. [rust](https://www.rust-lang.org/tools/install) **version >= 1.70.0** and associated tools such as cargo are required and assumed to be in PATH.
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
