"""
@author: inodb
"""
import sys
import os
import subprocess
import errno
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO

def get_gc_and_len_dict(fastafile):
    """Creates a dictionary with the fasta id as key and GC and length as keys
    for the inner dictionary."""
    out_dict = {}

    for rec in SeqIO.parse(fastafile, "fasta"):
        out_dict[rec.id] = {}
        out_dict[rec.id]["length"] = len(rec.seq)

    return out_dict


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

        try:
            d = out_dict[cols[0]]
        except KeyError:
            d = {}
            out_dict[cols[0]] = d

        if int(cols[1]) == 0:
            d["percentage_covered"] = 100 - float(cols[4]) * 100.0
        else:
            d["cov_mean"] = d.get("cov_mean", 0) + int(cols[1]) * float(cols[4])

    return out_dict

def print_sample_columns(t):
    sys.stdout.write(("\tcov_mean_sample_%s" * len(t)) % t)


def print_input_table(fastadict, bedcovdicts, samplenames=None):
    """Writes the input table for Probin to stdout. See hackathon google
    docs."""

    # Header
    sys.stdout.write("contig\tlength")
    if samplenames == None:
        # Use index if no sample names given in header
        print_sample_columns(tuple(range(len(bedcovdicts))))
    else:
        # Use given sample names in header
        assert(len(samplenames) == len(bedcovdicts))
        print_sample_columns(tuple(samplenames))
    sys.stdout.write("\n")

    # Content
    assert(len(fastadict) > 0)
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

def generate_input_table(fastafile, bamfiles, samplenames=None, isbedfiles=False):
    """Reads input files into dictionaries then prints everything in the table
    format required for running ProBin."""
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
