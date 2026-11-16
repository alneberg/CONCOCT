#!/usr/bin/env python
"""
A script to generate the input coverage table for CONCOCT using a BEDFile.

Output is written to stdout. The BEDFile defines the regions used as 
subcontigs for concoct. This makes it possible to get the coverage for 
subcontigs without specifically mapping reads against the subcontigs.

@author: inodb, alneberg
"""
from __future__ import print_function
import sys
import os
import argparse
import subprocess
import glob
from io import StringIO
from signal import signal, SIGPIPE, SIG_DFL
import pandas as pd

def check_bed_file_for_errors(bedfile):
    with open(bedfile) as ifh:
        for line in ifh:
            line = line.strip()
            original_id, _, _, cutup_id = line.split('\t')
            try:
                assert 'concoct_part_' not in original_id
                assert 'concoct_part_' in cutup_id
            except AssertionError:
                sys.stderr.write(("ERROR! Something is wrong with the line:\n'{}'\n"
                "Perhaps 'concoct_part_' is misplaced or missing? Exiting!\n").format(line))
                sys.exit(-1)

def exit_with_error(msg):
    sys.stderr.write("ERROR! {}\n".format(msg))
    sys.exit(-1)

def run_samtools(args):
    """Runs samtools with args and returns the returncode and stdout as a string."""
    try:
        p = subprocess.Popen(["samtools"] + args, stdout=subprocess.PIPE)
    except OSError:
        exit_with_error("Could not run samtools. Please make sure samtools is installed and available in your PATH.")
    out, _ = p.communicate()
    return p.returncode, out.decode('utf-8')

def get_bam_contigs(bamfile):
    """Returns the ids of the contigs in the header of bamfile, in order."""
    returncode, out = run_samtools(["view", "-H", bamfile])
    if returncode != 0:
        exit_with_error("Could not read the header of '{}', see the samtools error message above.".format(bamfile))
    contigs = []
    for line in out.splitlines():
        if line.startswith('@SQ'):
            for field in line.split('\t')[1:]:
                if field.startswith('SN:'):
                    contigs.append(field[3:])
    return contigs

def check_bam_files_for_errors(bedfile, bamfiles):
    """samtools bedcov looks up the contigs using the header of the first
    bamfile only. Contigs missing from the header are skipped and bamfiles
    with a different contig order silently get the coverage of the wrong
    contigs, so make sure this is not the case before running it."""
    with open(bedfile) as ifh:
        bed_contigs = set(line.split('\t')[0] for line in ifh if line.strip())

    first_bam_contigs = None
    for bamfile in bamfiles:
        bam_contigs = get_bam_contigs(bamfile)
        missing = bed_contigs.difference(bam_contigs)
        if missing:
            exit_with_error(("{} out of {} contigs in '{}' are not present in the header of '{}', for example '{}'. "
                "Were the reads mapped against the same fasta file that was given to cut_up_fasta.py?").format(
                    len(missing), len(bed_contigs), bedfile, bamfile, sorted(missing)[0]))
        if first_bam_contigs is None:
            first_bam_contigs = bam_contigs
        elif bam_contigs != first_bam_contigs:
            exit_with_error(("The contigs in the headers of '{}' and '{}' are not identical and in the same order. "
                "All bamfiles need to be mapped against the same fasta file.").format(bamfiles[0], bamfile))

def has_index(bamfile):
    candidates = [bamfile + '.bai', bamfile + '.csi', os.path.splitext(bamfile)[0] + '.bai']
    return any(os.path.isfile(c) for c in candidates)

def generate_input_table(bedfile, bamfiles, samplenames=None):
    """Reads input files into dictionaries then prints everything in the table
    format required for running CONCOCT."""

    returncode, out = run_samtools(["bedcov", bedfile] + bamfiles)

    if returncode != 0 or not out.strip():
        for bamfile in bamfiles:
            if not has_index(bamfile):
                sys.stderr.write(("No index was found for '{0}'. The bamfiles need to be sorted and indexed, "
                    "the index is created with 'samtools index {0}' and is expected to "
                    "be named '{0}.bai'.\n").format(bamfile))
        exit_with_error("Running 'samtools bedcov' failed or gave no output, see any samtools error message above.")
    else:
        # Header
        if samplenames == None:
            # Use index if no sample names given in header
            col_names = [os.path.splitext(os.path.basename(bamfile))[0] for bamfile in bamfiles]
        else:
            # Use given sample names in header
            col_names = samplenames
        header=["cov_mean_sample_{}".format(n) for n in col_names]

        # Content
        fh = StringIO(out)
        df = pd.read_csv(fh, header=None, sep='\t')
        avg_coverage_depth = df[df.columns[4:]].divide((df[2]-df[1]), axis=0)
        avg_coverage_depth.index = df[3]
        avg_coverage_depth.columns = header 
        avg_coverage_depth.to_csv(sys.stdout, index_label='contig', sep='\t', float_format='%.3f')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bedfile", help="Contigs BEDFile with four columns representing: 'Contig ID, Start Position, "
    "End Position and SubContig ID' respectively. The Subcontig ID must contain the pattern "
    "'concoct_part_[0-9]*' while the contigs which are not cutup cannot contain this pattern. "
    "This file can be generated by the cut_up_fasta.py script.")
    parser.add_argument("bamfiles", nargs='+', help="BAM files with mappings to the original contigs.")
    parser.add_argument("--samplenames", default=None, help="File with sample names, one line each. Should be same nr "
    "of bamfiles. Default sample names used are the file names of the bamfiles, excluding the file extension.")
    args = parser.parse_args()

    # Get sample names
    if args.samplenames != None:
        samplenames = [ s.strip() for s in open(args.samplenames).readlines() if s.strip() ]
        if len(samplenames) != len(args.bamfiles):
            raise Exception("Nr of names ({0}) in samplenames should be equal to nr of given "
                            "bamfiles ({1})".format(len(samplenames), len(args.bamfiles)))
    else:
        samplenames=None

    # ignore broken pipe error when piping output
    # http://newbebweb.blogspot.pt/2012/02/python-head-ioerror-errno-32-broken.html
    signal(SIGPIPE,SIG_DFL)
    check_bed_file_for_errors(args.bedfile)
    check_bam_files_for_errors(args.bedfile, args.bamfiles)
    generate_input_table(args.bedfile, args.bamfiles, samplenames=samplenames)
