import os
import sys
from random import randint
from argparse import ArgumentParser, ArgumentTypeError
from concoct.entry_points import main, parse

def set_random_state(seed):
    ERROR="'{0}' should be converatable to integer".format(seed)
    try:
        seed = int(seed)
        if seed < 0:
            raise ArgumentTypeError("'" + seed + "' should be >= 0")
        elif seed == 0:
            seed = randint(2,10000)
        return seed
    except ValueError as e:
        raise ArgumentTypeError(ERROR)

def get_version():
  from concoct import __version__
  return '%(prog)s {version}'.format(version=__version__)

def arguments():
    parser = ArgumentParser()

    subparsers = parser.add_subparsers()

    ## Parse parser
    ###################
    parser_parse = subparsers.add_parser('parse', help='Construct the CONCOCT input file from bam or BED files.')
    parser_parse.set_defaults(func=parse)
    parser_parse.add_argument("fastafile", help="Contigs fasta file")
    parser_parse.add_argument("bamfiles", nargs='+', help="BAM files with mappings to contigs")
    parser_parse.add_argument("--samplenames", default=None, help="File with sample names, one line each. Should be same nr as bamfiles.")
    parser_parse.add_argument("--isbedfiles", action='store_true',
        help="The bamfiles argument are outputs of genomeCoverageBed, not the actual bam file. Skips running genomeCoverageBed from within this script.")

    ## Cluster parser
    ###################
    parser_cluster = subparsers.add_parser('cluster', help='Run the clustering algorithm, this is the main method of CONCOCT.')
    parser_cluster.set_defaults(func=main)
    #Input files
    parser_cluster.add_argument('--coverage_file',
        help=("specify the coverage file, containing a table where each row "
              "correspond to a contig, and each column correspond to a sample. "
              "The values are the average coverage for this contig in that sample. "
              "All values are separated with tabs."))
    parser_cluster.add_argument('--composition_file',
        help=("specify the composition file, containing sequences in fasta format. "
              "It is named the composition file since it is used to calculate the "
              "kmer composition (the genomic signature) of each contig."))

    #Handle cluster number parsing
    parser_cluster.add_argument('-c', '--clusters', default=400, type=int,
      help='specify maximal number of clusters for VGMM, default 400.')
    #Kmer length, kmer count threshold and read length
    parser_cluster.add_argument('-k','--kmer_length', type=int, default=4,
        help='specify kmer length, default 4.')
    parser_cluster.add_argument('-l','--length_threshold', type=int, default=1000,
        help=("specify the sequence length threshold, contigs shorter than this "
              "value will not be included. Defaults to 1000."))
    parser_cluster.add_argument('-r','--read_length', type=int, default=100,
        help='specify read length for coverage, default 100')
    #Joined PCA
    parser_cluster.add_argument('--total_percentage_pca', default=90, type=int,
                        help=('The percentage of variance explained'
                              ' by the principal components for the'
                              ' combined data.'))
    #Output
    parser_cluster.add_argument('-b', '--basename', default=os.curdir,
      help=("Specify the basename for files or directory where output"
            "will be placed. Path to existing directory or basename"
            "with a trailing '/' will be interpreted as a directory."
            "If not provided, current directory will be used."))
    parser_cluster.add_argument('-s','--seed',type=set_random_state, default=set_random_state(1),
      help=('Specify an integer to use as seed for clustering. '
            '0 gives a random seed, 1 is the default seed and '
            'any other positive integer can be used. Other values '
            'give ArgumentTypeError.'))
    parser_cluster.add_argument('-i','--iterations',type=int, default=500,
      help=('Specify maximum number of iterations for the VBGMM. '
            'Default value is 500'))
    parser_cluster.add_argument('-e','--epsilon',type=float, default=1.0e-6,
      help=('Specify the epsilon for VBGMM. '
            'Default value is 1.0e-6'))
    parser_cluster.add_argument('--no_cov_normalization', default=False, action="store_true",
      help=("By default the coverage is normalized with regards to samples, "
            "then normalized with regards of contigs and finally log transformed. "
            "By setting this flag you skip the normalization and only do log "
            "transorm of the coverage."))
    parser_cluster.add_argument('--no_total_coverage', default=False, action="store_true",
      help=("By default, the total coverage is added as a new column in the coverage "
            "data matrix, independently of coverage normalization but previous to "
            "log transformation. Use this tag to escape this behaviour."))
    parser_cluster.add_argument('-o','--converge_out', default=False, action="store_true",
      help=('Write convergence info to files.'))


    parser.add_argument('-v','--version', action='version',
      version=get_version())

    args  = parser.parse_args()

    if args.func == main:
        check_main_args(args, parser)

    return args

def check_main_args(args, parser):
    # Make sure at least one input file was given
    if not (args.coverage_file or args.composition_file): 
        parser.error("No input data supplied, add file(s) using --coverage_file <cov_file> and/or "
                     "--composition_file <comp_file>")
    
    # PCA would interpret 100/100.0 as 1 dimension and not all
    if args.total_percentage_pca == 100:
        args.pca_components = "All"
    else:
        args.pca_components = args.total_percentage_pca/100.0


