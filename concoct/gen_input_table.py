"""
@author: inodb, alneberg
"""
import sys
import os
import subprocess
import errno
from signal import signal, SIGPIPE, SIG_DFL
import pysam
import pandas as p
from Bio import SeqIO


def get_bedcov_dict(bedcoverage):
    """Uses the BEDTools genomeCoverageBed histogram output to determine mean
    coverage and percentage covered for each contig.
    
    Returns dict with fasta id as key and percentage covered and cov_mean as
    keys for the inner dictionary."""
    out_dict = {}

    # Check if given argument is a file, otherwise use the content of the
    # variable
    if os.path.isfile(bedcoverage):
        fh = open(bedcoverage)
    else:
        fh = bedcoverage.split('\n')[:-1]

    for line in fh:
        cols = line.split()

        if cols[0] == 'genome':
            continue
        
        if cols[0] not in out_dict:
            out_dict[cols[0]] = 0

        if int(cols[1]) != 0:
            out_dict[cols[0]] += int(cols[1]) * float(cols[4])
    
    return out_dict

def get_bam_dict(bamfile):
    """Uses samtools pileup engine to calculate average coverage for each contig"""
    output_dict = {}
    samfile = pysam.AlignmentFile(bamfile, "rb")

    for ref, rlen in zip(samfile.references, samfile.lengths):
        iterator = samfile.pileup(reference = ref)
        avg_covs = sum(x.n for x in iterator) / float(rlen)
        output_dict[ref] = avg_covs

    return output_dict

def print_sample_columns(t):
    sys.stdout.write(("\tcov_mean_sample_%s" * len(t)) % t)

def print_input_table(covdicts):
    """Write output""" 
    df = p.DataFrame.from_dict(covdicts)
    df.fillna(0.0)
    df.to_csv(sys.stdout, sep='\t')

def print_input_table_old(fastadict, bedcovdicts, samplenames=None):
    """Writes the input table for Probin to stdout. See hackathon google
    docs."""

    # Header
    sys.stdout.write("contig")
    if samplenames == None:
        # Use index if no sample names given in header
        print_sample_columns(tuple(range(len(bedcovdicts))))
    else:
        # Use given sample names in header
        assert(len(samplenames) == len(bedcovdicts))
        print_sample_columns(tuple(samplenames))
    sys.stdout.write("\n")

    # Content
    for acc in fastadict:
        # fasta stats
        sys.stdout.write("%s\t%s"  %
            (
                acc,
                fastadict[acc]['length']
            )
        )

        # Print mean
        for bcd in bedcovdicts:
            try:
                # Print cov mean
                sys.stdout.write("\t%f" % (bcd[acc]["cov_mean"]))
            except KeyError:
                # No reads mapped to this contig
                sys.stdout.write("\t0")

        sys.stdout.write("\n")

def generate_input_table_old(fastafile, bamfiles, samplenames=None, isbedfiles=False):
    """Reads input files into dictionaries then prints everything in the table
    format required for running CONCOCT."""
    bedcovdicts = []
    
    # Determine coverage information from bam file using BEDTools
    for i, bf in enumerate(bamfiles):
        if isbedfiles == False:
            p = subprocess.Popen(["genomeCoverageBed", "-ibam", bf], stdout=subprocess.PIPE)
            out, err = p.communicate()
            if p.returncode != 0:
                sys.stderr.write(out)
                raise Exception('Error with genomeCoverageBed')
            else:
                bedcovdicts.append(get_bedcov_dict(out))
        else:
            bedcovdicts.append(get_bedcov_dict(bf))
                
    print_input_table(get_gc_and_len_dict(fastafile), bedcovdicts, samplenames=samplenames)

def generate_input_table(bamfiles, samplenames=None, isbedfiles=False):
    covdicts = {}

    if samplenames is None:
        samplenames = [os.path.basename(bf) for bf in bamfiles]
    
        if len(samplenames) != len(set(samplenames)):
            raise Exception("No sample names were given, but file names are not unique")
    else:
        if len(samplenames) != len(set(samplenames)):
            raise Exception("Sample names given are not unique")

    for bamfile, sample_name in zip(bamfiles, samplenames):
        if isbedfiles == False:
            covdicts[sample_name] = get_bam_dict(bamfile)
        else:
            covdicts[sample_name] = get_bedcov_dict(bamfile)

    print_input_table(covdicts)

