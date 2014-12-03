
Usage
=====

CONCOCT uses different subcommands to accomplish different things. To find out which subcommands are available, run ``concoct -h`` on the command line::

	usage: concoct [-h] [-v] {parse,cluster} ...
	
	positional arguments:
	  {parse,cluster}
	    parse          Construct the CONCOCT input file from bam or BED files.
	    cluster        Run the clustering algorithm, this is the main method of
	                   CONCOCT.
	
	optional arguments:
	  -h, --help       show this help message and exit
	  -v, --version    show program's version number and exit


Parse
-----
The `parse` subcommand is used to generate the input table. This functionality was previously included in the script 'gen_input_table.py'. To see the options run ``concoct parse -h``::

	usage: concoct parse [-h] [--samplenames SAMPLENAMES] [--isbedfiles]
	                     fastafile bamfiles [bamfiles ...]
	
	positional arguments:
	  fastafile             Contigs fasta file
	  bamfiles              BAM files with mappings to contigs
	
	optional arguments:
	  -h, --help            show this help message and exit
	  --samplenames SAMPLENAMES
	                        File with sample names, one line each. Should be same
	                        nr as bamfiles.
	  --isbedfiles          The bamfiles argument are outputs of
	                        genomeCoverageBed, not the actual bam file. Skips
	                        running genomeCoverageBed from within this script.


Cluster
------
For the subcommand `cluster` concoct has several command line options to control the clustering, here is a complete documentation of these. These can also be viewed by typing ``concoct cluster -h`` on the command line.::

	usage: concoct cluster [-h] [--coverage_file COVERAGE_FILE]
	                       [--composition_file COMPOSITION_FILE] [-c CLUSTERS]
	                       [-k KMER_LENGTH] [-l LENGTH_THRESHOLD] [-r READ_LENGTH]
	                       [--total_percentage_pca TOTAL_PERCENTAGE_PCA]
	                       [-b BASENAME] [-s SEED] [-i ITERATIONS] [-e EPSILON]
	                       [--no_cov_normalization] [--no_total_coverage] [-o]
	
	optional arguments:
	  -h, --help            show this help message and exit
	  --coverage_file COVERAGE_FILE
	                        specify the coverage file, containing a table where
	                        each row correspond to a contig, and each column
	                        correspond to a sample. The values are the average
	                        coverage for this contig in that sample. All values
	                        are separated with tabs.
	  --composition_file COMPOSITION_FILE
	                        specify the composition file, containing sequences in
	                        fasta format. It is named the composition file since
	                        it is used to calculate the kmer composition (the
	                        genomic signature) of each contig.
	  -c CLUSTERS, --clusters CLUSTERS
	                        specify maximal number of clusters for VGMM, default
	                        400.
	  -k KMER_LENGTH, --kmer_length KMER_LENGTH
	                        specify kmer length, default 4.
	  -l LENGTH_THRESHOLD, --length_threshold LENGTH_THRESHOLD
	                        specify the sequence length threshold, contigs shorter
	                        than this value will not be included. Defaults to
	                        1000.
	  -r READ_LENGTH, --read_length READ_LENGTH
	                        specify read length for coverage, default 100
	  --total_percentage_pca TOTAL_PERCENTAGE_PCA
	                        The percentage of variance explained by the principal
	                        components for the combined data.
	  -b BASENAME, --basename BASENAME
	                        Specify the basename for files or directory where
	                        outputwill be placed. Path to existing directory or
	                        basenamewith a trailing '/' will be interpreted as a
	                        directory.If not provided, current directory will be
	                        used.
	  -s SEED, --seed SEED  Specify an integer to use as seed for clustering. 0
	                        gives a random seed, 1 is the default seed and any
	                        other positive integer can be used. Other values give
	                        ArgumentTypeError.
	  -i ITERATIONS, --iterations ITERATIONS
	                        Specify maximum number of iterations for the VBGMM.
	                        Default value is 500
	  -e EPSILON, --epsilon EPSILON
	                        Specify the epsilon for VBGMM. Default value is 1.0e-6
	  --no_cov_normalization
	                        By default the coverage is normalized with regards to
	                        samples, then normalized with regards of contigs and
	                        finally log transformed. By setting this flag you skip
	                        the normalization and only do log transorm of the
	                        coverage.
	  --no_total_coverage   By default, the total coverage is added as a new
	                        column in the coverage data matrix, independently of
	                        coverage normalization but previous to log
	                        transformation. Use this tag to escape this behaviour.
	  -o, --converge_out    Write convergence info to files.

