from __future__ import division

import sys
import logging

import vbgmm

from concoct.output import Output
from concoct.cluster import cluster
from concoct.input import load_data
from concoct.transform import perform_pca

def main(args):

    # Initialize output handling
    Output(args.basename,args)

    composition, cov, cov_range = load_data(args)

    # If there are zero or one contig that exceed the filter, do not continue
    if len(composition) < 2:
        logging.error('Not enough contigs pass the threshold filter. Exiting!')
        sys.exit(-1)
    
    if cov is not None:
        joined = composition.join(cov.ix[:,cov_range[0]:cov_range[1]],how="inner")
    else:
        joined = composition

    # Fix special case in pca_components
    if args.pca_components == "All":
        args.pca_components = joined.shape[1]

    #PCA on the contigs that have kmer count greater than length_threshold
    transform_filter, pca = perform_pca(
        joined,
        args.pca_components
        )

    logging.info('Performed PCA, resulted in %s dimensions' % transform_filter.shape[1])
    
    Output.write_original_data(
        joined,
        args.length_threshold
        )

    Output.write_pca(
        transform_filter,
        args.length_threshold,
        joined.index,
        )

    Output.write_pca_components(
        pca.components_,
        args.length_threshold
        )

    logging.info('PCA transformed data.')

    logging.info('Will call vbgmm with parameters: %s, %s, %s' % (Output.CONCOCT_PATH, args.clusters, args.length_threshold))

    vbgmm.fit(Output.CONCOCT_PATH, args.clusters, args.length_threshold,args.seed,args.iterations,args.epsilon,args.converge_out)

    logging.info("CONCOCT Finished")


