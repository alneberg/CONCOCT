#!/usr/bin/env python
"""extract_fasta_bins.py 

Extract a fasta file for each cluster from a concoct result file.
"""

import argparse
import sys
import os
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

def main(args):
    all_seqs = {}
    for i, seq in enumerate(SeqIO.parse(args.fasta_file, "fasta")):
        all_seqs[seq.id] = seq
    # Contig ids consisting of digits only should still be read as strings
    df = pd.read_csv(args.cluster_file, dtype={'contig_id': str})
    try:
        assert df.columns[0] == 'contig_id'
        assert df.columns[1] == 'cluster_id'
    except AssertionError:
        sys.stderr.write("ERROR! Header line was not 'contig_id, cluster_id', please adjust your input file. Exiting!\n")
        sys.exit(-1)

    missing = [contig_id for contig_id in df['contig_id'] if contig_id not in all_seqs]
    if missing:
        sys.stderr.write("ERROR! {} out of {} contigs in the cluster file are not present in the fasta file, for example '{}'.\n".format(
            len(missing), len(df), missing[0]))
        if 'concoct_part_' in missing[0]:
            sys.stderr.write("The cluster file seems to contain cut up contigs. Either run merge_cutup_clustering.py on the "
                "cluster file first or use the cut up fasta file as input.\n")
        sys.exit(-1)

    cluster_to_contigs = defaultdict(list)
    for contig_id, cluster_id in zip(df['contig_id'], df['cluster_id']):
        cluster_to_contigs[cluster_id].append(contig_id)
    
    for cluster_id, contig_ids in cluster_to_contigs.items():
        output_file = os.path.join(args.output_path, "{0}.fa".format(cluster_id))
        seqs = [all_seqs[contig_id] for contig_id in contig_ids] 
        with open(output_file, 'w') as ofh:
            SeqIO.write(seqs, ofh, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_file", help="Input Fasta file.")
    parser.add_argument("cluster_file", help="Concoct output cluster file")
    parser.add_argument("--output_path", default=os.getcwd(), 
            help="Directory where files will be printed")
    args = parser.parse_args()

    main(args)
