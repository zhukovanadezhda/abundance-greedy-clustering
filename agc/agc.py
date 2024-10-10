#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
from collections import Counter
import gzip
from pathlib import Path
import sys
import textwrap
from typing import Iterator, List
import nwalign3 as nw


__author__ = "Nadezhda Zhukova"
__copyright__ = "Université Paris Cité"
__credits__ = ["Nadezhda Zhukova"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Nadezhda Zhukova"
__email__ = "nadiajuckova@gmail.com"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__,
                                     usage=f"{sys.argv[0]} -h")
    parser.add_argument('-i',
                        '-amplicon_file',
                        dest='amplicon_file',
                        type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s',
                        '-minseqlen',
                        dest='minseqlen',
                        type=int,
                        default = 400,
                        help=("Minimum sequence length for dereplication "
                              "(default 400)"))
    parser.add_argument('-m',
                        '-mincount',
                        dest='mincount',
                        type=int,
                        default = 10,
                        help="Minimum count for dereplication (default 10)")
    parser.add_argument('-o',
                        '-output_file',
                        dest='output_file',
                        type=Path,
                        default=Path("OTU.fasta"),
                        help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, 'rt') as fasta_file:
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                # If there is a sequence ready and its length is valid, yield
                if sequence and len("".join(sequence)) >= minseqlen:
                    yield "".join(sequence)
                # Reset sequence buffer for the next entry
                sequence = []
            else:
                # Continue collecting the sequence
                sequence.append(line)

        # Yield the last sequence if it meets the length requirement
        if sequence and len("".join(sequence)) >= minseqlen:
            yield "".join(sequence)


def dereplication_fulllength(amplicon_file: Path,
                             minseqlen: int,
                             mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of 
             sequence with a count >= mincount and a length >= minseqlen.
    """
    # Get all sequences that satisfy the length requirement
    sequence_counter = Counter(read_fasta(amplicon_file, minseqlen))

    # Filter out sequences that do not meet the mincount requirement
    # Sort by occurrence (most common first)
    for sequence, count in sequence_counter.most_common():
        if count >= mincount:
            yield [sequence, count]


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences 
                            in the format ["SEQUENCE1", "SEQUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    seq1, seq2 = alignment_list
    matches = sum(res1 == res2 for res1, res2 in zip(seq1, seq2))
    total_length = len(seq1)

    # Calculate identity as the percentage of matches over the total length
    return (matches / total_length) * 100 if total_length > 0 else 0.0


def abundance_greedy_clustering(amplicon_file: Path,
                                minseqlen: int,
                                mincount: int,
                                chunk_size: int,
                                kmer_size: int) -> List:
    """Compute an abundance greedy clustering in sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    # List to store OTUs
    otu_list = []

    # Get sequences in descending order of abundance
    for sequence, count in dereplication_fulllength(amplicon_file,
                                                    minseqlen,
                                                    mincount):
        is_otu = True

        # Compare with existing OTUs
        for otu, _ in otu_list:
            # Align the sequence with the OTU
            alignment = [sequence, otu]

            # Calculate identity
            identity = get_identity(alignment)

            # If the identity is greater than 97%, don't add it as an OTU
            if identity >= 97:
                is_otu = False
                break

        # If the sequence is not similar to existing OTUs, add it to the OTU
        if is_otu:
            otu_list.append([sequence, count])

    return otu_list


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with output_file.open('w') as f:
        for idx, (sequence, count) in enumerate(OTU_list, 1):
            # Write the header in the format ">OTU_{index} occurrence:{count}"
            f.write(f">OTU_{idx} occurrence:{count}\n")
            # Use textwrap to wrap the sequence to 80 characters per line
            wrapped_sequence = textwrap.fill(sequence, width=80)
            f.write(f"{wrapped_sequence}\n")


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    otu_list = abundance_greedy_clustering(args.amplicon_file,
                                           args.minseqlen,
                                           args.mincount,
                                           100,
                                           8)
    write_OTU(otu_list, args.output_file)


if __name__ == '__main__':
    main()
