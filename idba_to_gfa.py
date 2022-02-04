#!/usr/bin/env python3
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/IDBA-to-GFA

This file is part of IDBA to GFA. IDBA to GFA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. IDBA to GFA is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with IDBA to GFA. If
not, see <http://www.gnu.org/licenses/>.
"""



import argparse
import collections
import subprocess
import os
import shutil
import sys


def main():
    args = get_arguments()
    if shutil.which(args.print_graph) is None:
        sys.exit('Error: could not find print_graph')

    # Run the IDBA print_graph tool to make a graph.
    temp_fasta = 'idba_graph_to_gfa_temp.fasta'
    connection_str = subprocess.check_output([args.print_graph, '-k', str(args.kmer), '--max_length', '1000000000', args.idba_assembly, temp_fasta]).decode()
    depths = load_depths(args.idba_assembly)
    sequences = load_fasta(temp_fasta)
    connections = load_connections(connection_str, sequences, args.kmer)
    os.remove(temp_fasta)

    # Print the GFA segment lines.
    for seg_num in sorted(sequences.keys()):
        if seg_num in depths :
            print('\t'.join(['S', str(seg_num), sequences[seg_num], "RC:"+depths[seg_num]]))
        else :
            print('\t'.join(['S', str(seg_num), sequences[seg_num]]))

    # Print the GFA link lines.
    overlap_str = str(args.kmer - 1) + 'M'
    for start_seg in sorted(connections.keys(), key=lambda x: (abs(x), -x)):
        start_strand = '+' if start_seg > 0 else '-'
        end_segs = sorted(connections[start_seg])
        for end_seg in end_segs:
            end_strand = '+' if end_seg > 0 else '-'
            print('\t'.join(['L', str(abs(start_seg)), start_strand, str(abs(end_seg)), end_strand, overlap_str]))


def get_arguments():
    parser = argparse.ArgumentParser(description='IDBA to GFA: a tool for converting IDBA assemblies to GFA graphs',
                                     epilog='output: GFA to stdout')
    parser.add_argument('idba_assembly', type=str,
                        help='Assembly from IDBA, e.g. scaffold.fa')
    parser.add_argument('kmer', type=int,
                        help='assembly k-mer size (e.g. 100)')
    parser.add_argument('--print_graph', type=str, default='print_graph',
                        help="Location of IDBA's print_graph tool (required if print_graph is not in PATH)")

    return parser.parse_args()

def load_depths(filename):
    fasta_depths = {}
    
    with open(filename, 'rt') as fasta_file:
        
        for line in fasta_file:
            
            if line[0] == '>':  # Header line = start of new contig
                if 'count' in line.split()[-1].split('_')[1] :
                    name = int(line[1:].split()[0].split('_')[1])
                    fasta_depths[name] = line.split()[-1].split('_')[2]
    return fasta_depths
    
def load_fasta(filename):
    fasta_seqs = {}
    
    with open(filename, 'rt') as fasta_file:
        name = None
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name is not None:
                    fasta_seqs[name] = sequence
                    sequence = ''
                name = int(line[1:].split('_')[1])
                
            else:
                sequence += line
        if name:
            fasta_seqs[name] = sequence
    return fasta_seqs

def load_connections(connection_str, sequences, kmer_size):
    connections = collections.defaultdict(set)
    k_minus_1 = kmer_size - 1

    connection_lines = connection_str.split('\n')

    for line in connection_lines:
        if line.count(' ') != 1 or line.count('_') != 6:
            continue

        part_1, part_2 = line.strip().split()
        part_1_parts = part_1.split('_')
        part_2_parts = part_2.split('_')

        segment_1_num = int(part_1_parts[0])
        segment_1_strand = part_1_parts[1]
        segment_1_seq = sequences[segment_1_num]
        if segment_1_strand == '1':
            segment_1_seq = reverse_complement(segment_1_seq)

        segment_2_num = int(part_2_parts[0])
        segment_2_strand = part_2_parts[1]
        segment_2_seq = sequences[segment_2_num]
        if segment_2_strand == '1':
            segment_2_seq = reverse_complement(segment_2_seq)

        segment_1_signed_num = segment_1_num if segment_1_strand == '0' else -segment_1_num
        segment_2_signed_num = segment_2_num if segment_2_strand == '0' else -segment_2_num

        # Figure out which way the connection goes.
        segment_1_seq_start = segment_1_seq[:k_minus_1]
        segment_1_seq_end = segment_1_seq[-k_minus_1:]
        segment_2_seq_start = segment_2_seq[:k_minus_1]
        segment_2_seq_end = segment_2_seq[-k_minus_1:]
        one_to_two = (segment_1_seq_end == segment_2_seq_start)
        two_to_one = (segment_2_seq_end == segment_1_seq_start)

        if one_to_two:
            connections[segment_1_signed_num].add(segment_2_signed_num)
            connections[-segment_2_signed_num].add(-segment_1_signed_num)
        if two_to_one:
            connections[segment_2_signed_num].add(segment_1_signed_num)
            connections[-segment_1_signed_num].add(-segment_2_signed_num)

    return connections


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N',
                 'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
                 'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
                 'd': 'h', 'h': 'd', 'n': 'n',
                 '.': '.', '-': '-', '?': '?'}


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


if __name__ == '__main__':
    main()
