# 🦠🧬OTU Calculation Using Abundance Greedy Clustering

This repository contains a Python program to calculate OTUs (Operational Taxonomic Units) from a bacterial sequencing dataset. The focus is on processing "mock" sequencing data, with eight bacterial species expected in the analysis.

The program performs the following:
- Full-length sequence dereplication
- Clustering using the **Abundance Greedy Clustering** algorithm

## 🔄Installation

To set up the environment and install the required dependencies, use the following commands:

```bash
conda env create -f environment.yml
conda activate abundance-greedy-clustering
```

## 🧑‍💻️Usage

First, clone the repository and navigate to the project folder:

```bash
git clone git@github.com:yourusername/abundance-greedy-clustering.git
cd abundance-greedy-clustering
```

The program processes sequences in FASTA format for OTU calculation and accepts the following arguments:

- `-i`, `--amplicon_file`: Path to the input FASTA file
- `-s`, `--minseqlen`: Minimum sequence length (optional, default: 400)
- `-m`, `--mincount`: Minimum sequence occurrence count (optional, default: 10)
- `-c`, `--chunk_size`: Chunk size for sequence partitioning (optional, default: 100)
- `-k`, `--kmer_size`: K-mer size for sequence analysis (optional, default: 8)
- `-o`, `--output_file`: Path to the output file where the calculated OTUs in FASTA format will be saved

To run the program, execute the following command:

```bash
python3 agc/agc.py -i data/amplicon.fasta.gz -o output/OTU.fasta
```

## ⚙️Testing

To run unit tests and measure code coverage, use:

```bash
pytest --cov=agc -v -s --ignore=tests/test_chimera_removal.py 
```

## 🎁Example Usage

To run the program and calculate OTUs, execute the following command:

```bash
python3 agc/agc.py -i data/amplicon.fasta.gz
```

This will output the OTU sequences in a file named `OTU.fasta`.

### Verifying Results with `vsearch`

To assess the quality of the OTUs generated, you can use `vsearch` to compare them against a reference 16S rRNA database (e.g., `mock_16S.fasta`).

Align OTUs against the reference database using the `usearch_global` function in `vsearch`:

   ```bash
   vsearch --usearch_global OTU.fasta --db data/mock_16S.fasta --id 0.8 --blast6out results.tsv
   ```

### Example `vsearch` output:

```bash
vsearch v2.29.0_linux_x86_64, 7.6GB RAM, 12 cores
https://github.com/torognes/vsearch

Reading file data/mock_16S.fasta 100%  
15480 nt in 10 seqs, min 1526, max 1568, avg 1548
Masking 100% 
Counting k-mers 100% 
Creating k-mer index 100% 
Searching 100%  
Matching unique query sequences: 116 of 117 (99.15%)
```

This output demonstrates that the OTUs generated are well-matched against the reference sequences.

## ✉️Contact

For questions or support, please contact [nadiajuckova@gmail.com](mailto:nadiajuckova@gmail.com).
