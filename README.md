# OTU Calculation with Abundance Greedy Clustering

This repository contains a Python program designed to calculate OTUs (Operational Taxonomic Units) from a "mock" sequencing dataset, focusing on bacterial sequences. Eight species are expected in the analysis.

The program includes functionality for full-length sequence dereplication, chimeric sequence detection, and clustering based on an Abundance Greedy Clustering algorithm.

## Dependency Installation

To set up the environment, run the following commands:
```bash
conda env create -f environment.yml
conda activate abundance-greedy-clustering
```

## Usage

The program processes sequences in FASTA format and performs dereplication, chimeric detection, and clustering. It accepts the following arguments:

- `-i`, `--amplicon_file`: Input file containing sequences in FASTA format
- `-s`, `--minseqlen`: Minimum sequence length (optional, default: 400)
- `-m`, `--mincount`: Minimum sequence count (optional, default: 10)
- `-c`, `--chunk_size`: Sequence partition size (optional, default: 100)
- `-k`, `--kmer_size`: K-mer length (optional, default: 8)
- `-o`, `--output_file`: Output file to save OTUs in FASTA format

## Testing

To run tests, use the command:
```
pytest --cov=agc
```

## Contact

For any questions, please contact nadiajuckova@gmail.com.
