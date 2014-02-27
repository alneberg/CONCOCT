#!/usr/bin/env python
"""
Calculates validation statistics for a predicted clustering compared
to a reference clustering.
"""

from sklearn.metrics import adjusted_rand_score
import pandas as pd
import numpy as np

def int_labels_from_names(label_ser):
    int_labels = np.empty_like(label_ser.values)
    for label, i in enumerate(label_ser.uniq):
        int_labels[label_ser == label] = i
    return int_labels

        
def main(args):
    ref = pd.Series.from_csv(args.reference_clustering_file) 
    pred = pd.Series.from_csv(args.reference_clustering_file)
    print adjusted_rand_score(int_labels(ref), int_labels(pred))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_clustering_file", type=file,
        help="CSV file containing reference clustering for some contigs.")
    parser.add_argument("predicted_clustering_file", type=file,
        help="CSV file containing predicted clustering for each contig.")
    args = parser.parse_args()
    main(args)
        

