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
import os
from pathlib import Path
import statistics
import sys
import textwrap
from typing import Iterator, Dict, List
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
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
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

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    seq1, seq2 = alignment_list
    matches = sum(res1 == res2 
                  for res1, res2 
                  in zip(seq1, seq2) 
                  if res1 != '-' and res2 != '-')
    length = sum(res1 != '-' and res2 != '-' 
                 for res1, res2 
                 in zip(seq1, seq2))
    return (matches / length) * 100 if length > 0 else 0.0


def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    pass


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    pass


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici



if __name__ == '__main__':
    main()
