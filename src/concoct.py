#!/usr/bin/env python
from __future__ import division

import sys
import os
import re

import pandas as p
import numpy as np

from itertools import tee, izip, product
from argparse import ArgumentParser,ArgumentTypeError
from datetime import datetime

from Bio import SeqIO

from sklearn.preprocessing import scale
from sklearn.mixture import GMM
from sklearn.decomposition import PCA



class Output(object):
    """
    Class to print out result information to their files
    """
    CONCOCT_PATH = None
    DT = None
    BIC_FILE = None
    ARGS_FILE = None
    PCA_FILE_BASE = None
    CLUSTERING_FILE_BASE = None
    MEANS_FILE_BASE = None
    VARIANCE_FILE_BASE = None
    RESPONSIBILITY_FILE_BASE = None
    @classmethod
    def __init__(self,outdir,args):
        """
        Set output parameters and create output folders and bic.csv and args.txt
        """
        self.DT = datetime.now().strftime("%y%m%d_%H%M")
        self.CONCOCT_PATH = os.path.join(outdir,"concoct_{0}".format(self.DT))
        os.makedirs(self.CONCOCT_PATH)
        print >> sys.stderr, "Results created in folder {0}".format(self.CONCOCT_PATH)
        self.BIC_FILE = os.path.join(self.CONCOCT_PATH,"bic.csv")
        self.ARGS_FILE = os.path.join(self.CONCOCT_PATH,"args.txt")
        self.ORIGINAL_FILE_BASE = os.path.join(self.CONCOCT_PATH,"original_data_gt{0}_{1}")
        self.PCA_FILE_BASE = os.path.join(self.CONCOCT_PATH,"PCA_gt{0}_{1}.csv")
        self.CLUSTERING_FILE_BASE = os.path.join(self.CONCOCT_PATH,"clustering_{0}.csv")
        self.MEANS_FILE_BASE = os.path.join(self.CONCOCT_PATH,"means_{0}_gt{1}.csv")
        self.VARIANCE_FILE_BASE = os.path.join(self.CONCOCT_PATH,"variance_{0}_gt{1}.csv")
        self.RESPONSIBILITY_FILE_BASE = os.path.join(self.CONCOCT_PATH,"responsibility_{0}_gt{1}.csv")
        #Write header to bic.csv
        with open(self.BIC_FILE,"a+") as fh:
            print >> fh, "cluster_count,bic_value"
        with open(self.ARGS_FILE,"w+") as fh:
            print >> fh, args
    
    @classmethod
    def write_pca(self,transform,threshold,prefix):
        np.savetxt(self.PCA_FILE_BASE.format(threshold,prefix),transform)
    
    @classmethod
    def write_original_data(self,original,threshold,prefix):
        original.to_csv(self.ORIGINAL_FILE_BASE.format(threshold,prefix))

    
    @classmethod
    def write_clustering(self,dataframe,threshold_filter,threshold,c):
        dataframe.clustering.to_csv(self.CLUSTERING_FILE_BASE.format("{0}_full".format(c)))
        dataframe[threshold_filter].clustering.to_csv(self.CLUSTERING_FILE_BASE.format("{0}_gt{1}_filtered".format(c,threshold)))
        
    @classmethod
    def write_bic(self,bic,c):
        with open(self.BIC_FILE,"a+") as fh:
            print >> fh, "{0},{1}".format(bic,c)
    
    @classmethod
    def write_cluster_means(self,means,threshold,c):
        np.savetxt(self.MEANS_FILE_BASE.format(c,threshold),means)
            
    @classmethod
    def write_cluster_variance(self,means,threshold,c):
        np.savetxt(self.VARIANCE_FILE_BASE.format(c,threshold),means)

    @classmethod
    def write_cluster_responsibilities(self,res,threshold,c):
        np.savetxt(self.RESPONSIBILITY_FILE_BASE.format(c,threshold),res)        

def cluster(comp_file, cov_file, kmer_len, threshold, read_length, clusters_range, cov_range, split_pca, inits, iters, outdir, args=None):
    Output(outdir,args)
    #Composition
    #Generate kmer dictionary
    feature_mapping, nr_features = generate_feature_mapping(kmer_len)
    #Count lines in composition file
    count_re = re.compile("^>")
    seq_count = 0
    with open(comp_file) as fh:
        for line in fh:
            if re.match(count_re,line):
                seq_count += 1

    #Initialize with ones since we do pseudo count, we have i contigs as rows and j features as columns
    composition = np.ones((seq_count,nr_features))
    
    
    contigs_id = []
    for i,seq in enumerate(SeqIO.parse(comp_file,"fasta")):
        contigs_id.append(seq.id)
        for kmer_tuple in window(seq.seq.tostring().upper(),kmer_len):
            composition[i,feature_mapping["".join(kmer_tuple)]] += 1
    composition = p.DataFrame(composition,index=contigs_id,dtype=float)
    #Select contigs to cluster on
    threshold_filter = composition.sum(axis=1) > threshold
    
    #log(p_ij) = log[(X_ij +1) / rowSum(X_ij+1)]
    composition = np.log(composition.divide(composition.sum(axis=1),axis=0))
    
    #Coverage import, file has header and contig ids as index
    cov = p.read_table(cov_file,header=0,index_col=0)
    #TODO: Here we expect the data to be coverage not read counts. The commented section calculates log coverage from read counts.
    ###log(q_ij) = log[(Y_ij + 1).R/L_i]) where L_i is the length of contig i and R is the read length.
    ##cov.ix[:,cov_range[0]:cov_range[1]] = np.log((cov.ix[:,cov_range[0]:cov_range[1]] + 1).mul((read_length/cov.length)))
    cov.ix[:,cov_range[0]:cov_range[1]] = np.log((cov.ix[:,cov_range[0]:cov_range[1]] + 0.01))
    #cov = scale(cov.ix[:,cov_range[0]:cov_range[1]])

    if split_pca:
        raise NotImplementedError("Not implemented yet to run seperate PCA")
    else:
        joined = composition.join(cov.ix[:,cov_range[0]:cov_range[1]],how="inner")
        Output.write_original_data(joined,threshold,"joined")
        #PCA on the contigs that have kmer count greater than threshold
        pca = PCA(n_components=0.9).fit(joined[threshold_filter])
        transform_filter = pca.transform(joined[threshold_filter])
        Output.write_pca(pca.transform(transform_filter),threshold,"joined_filtered")
    
    cv_type='full'
    for c in clusters_range:
        #Run GMM on the pca transform of contigs with kmer count greater than threshold
        gmm = GMM(n_components=c, covariance_type=cv_type, n_init=inits,n_iter=iters).fit(transform_filter)
        print >> sys.stderr, "Convergence for cluster number {0}: {1}".format(c,gmm.converged_)
        #Classify all datapoints based on the clustering of filtered contigs
        joined["clustering"] = gmm.predict(pca.transform(joined))
        Output.write_clustering(joined,threshold_filter,threshold,c)
        Output.write_bic(gmm.bic(pca.transform(transform_filter)),c)
        Output.write_cluster_means(pca.inverse_transform(gmm.means_),threshold,c)
        Output.write_cluster_variance(pca.inverse_transform(gmm.covars_),threshold,c)
        Output.write_cluster_responsibilities(pca.inverse_transform(gmm.predict_proba(transform_filter)),threshold,c)
        
def window(seq,n):
    els = tee(seq,n)
    for i,el in enumerate(els):
        for _ in xrange(i):
            next(el, None)
    return izip(*els)

def generate_feature_mapping(kmer_len):
    BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
    kmer_hash = {}
    counter = 0
    for kmer in product("ATGC",repeat=kmer_len):
        kmer = ''.join(kmer)
        if kmer not in kmer_hash:
            kmer_hash[kmer] = counter
            rev_compl = ''.join([BASE_COMPLEMENT[x] for x in reversed(kmer)])
            kmer_hash[rev_compl] = counter
            counter += 1
    return kmer_hash, counter+1

def parse_cluster_list(cc):
    ERROR="'" + cc + "' is not a valid range of number. Expected forms like '20,100,2'."
    try:
        first, last, step = map(int,cc.split(","))
    except ValueError as e:
        raise ArgumentTypeError(ERROR)
    except Exception as e:
        raise ArgumentTypeError(ERROR)
    return xrange(first, last+1, step)
    
def parse_coverage_columns(cov_string):
    ERROR="'" + cov_string + "' is not valid. Expected 'first_column_name,last_column_name'."
    try:
        cov = cov_string.split(",")
    except ValueError as e:
        raise ArgumentTypeError(ERROR)
    if not len(cov) == 2:
        raise ArgumentTypeError(ERROR)
    return cov

def parse_taxonomy_cluster_list(tax_file):
    raise NotImplementedError("This functionality has not been added yet. Please use -c and specify range")

def arguments():
    parser = ArgumentParser()

    #Input files
    parser.add_argument('coverage_file',
        help='specify the coverage file')
    parser.add_argument('composition_file',
        help='specify the composition file')

    #Handle cluster number parsing
    cluster_count = parser.add_mutually_exclusive_group()
    cluster_count.add_argument('-c', '--clusters', default=range(20,101,2), type=parse_cluster_list,
        help='specify range of clusters to try out on format first,last,step. \
              default 20,100,2.')
    cluster_count.add_argument('-t', '--taxonomy_file', type=parse_taxonomy_cluster_list,
        help='specify a taxonomy file to estimate species number from (X). \
              Will use range X*0.5,X*1.5,2')
    #Columns in coverage file to use
    parser.add_argument('-n','--coverage_file_column_names', type=parse_coverage_columns, default=None,
        help='specify the first and last column names for continuous coverage \
              range of read counts as first,last')   

    #Kmer length, kmer count threshold and read length
    parser.add_argument('-k','--kmer_length', type=int, default=4,
        help='specify kmer length, defaults to tetramer')
    parser.add_argument('-l','--limit_kmer_count', type=int, default=1000,
        help='specify the kmer count for threshold in running PCA on \
              composition contigs, default 1000')
    parser.add_argument('-r','--read_length', type=int, default=100,
        help='specify read length for coverage, default 100')
    #Joined PCA or seperate PCA
    parser.add_argument('-s','--split_pca', default=False, action="store_true",
        help='specify this flag to first do PCA for the composition \
              and using that component number that explaines 90 percent \
              of variance for the coverage as well. Default join composition \
              and coverage before PCA.')
    #Clustering Parameters
    parser.add_argument('-e', '--executions',type=int, default=5,
        help='How often to initialize each cluster count. default 5 times')
    parser.add_argument('-i', '--iterations',type=int, default=100,
        help='Maximum number of iterations if convergance not achieved')
    #Output
    parser.add_argument('-o', '--outdir', default=os.curdir,
        help='specify the output directory, if not provided current directory \
              used. All files will be created in:\
              folder/CONCOCT_YYMMDD_HHMM')
    
    
    return parser.parse_args()

        
if __name__=="__main__":
    args = arguments()
    results = cluster(args.composition_file, 
                      args.coverage_file,
                      args.kmer_length, 
                      args.limit_kmer_count, 
                      args.read_length, 
                      args.clusters, 
                      args.coverage_file_column_names, 
                      args.split_pca, 
                      args.executions, 
                      args.iterations, 
                      args.outdir, 
                      args)
