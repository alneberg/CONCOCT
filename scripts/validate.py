#!/usr/bin/env python
"""
Calculates validation statistics for a predicted clustering compared
to a reference clustering.
"""

from sklearn import metrics
import pandas as pd
import numpy as np
import argparse
import sys

def int_labels(label_ser):
    return_labels = np.empty_like(label_ser.values)
    for i, label in enumerate(label_ser.unique()):
        return_labels[label_ser == label] = i
    return return_labels


def clustering_statistics(ref, pred):
    stats = {}
    stats['AdjRand'] = metrics.adjusted_rand_score(ref, pred)
    stats['Completeness'] = metrics.completeness_score(ref, pred)
    stats['Homogenity'] = metrics.homogeneity_score(ref, pred)
    return pd.Series(stats)
        
def main(args):
    ref = pd.Series.from_csv(args.reference_clustering_file) 
    pred = pd.Series.from_csv(args.predicted_clustering_file)
    pred_intersection = pred.ix[ref.index.intersection(pred.index)]
    ref_intersection = ref.ix[pred_intersection.index]
    stats = clustering_statistics(int_labels(ref_intersection), 
                                  int_labels(pred_intersection))
    stats.to_csv(sys.stdout, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_clustering_file", 
        help="CSV file containing reference clustering for some contigs.")
    parser.add_argument("predicted_clustering_file",
        help="CSV file containing predicted clustering for each contig.")
    args = parser.parse_args()
    main(args)
        

