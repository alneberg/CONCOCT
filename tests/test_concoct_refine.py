#!/usr/bin/env python
import unittest
from os.path import isdir
import os
import subprocess
import pandas as p

file_path = os.path.realpath(__file__)
test_dir_path = os.path.dirname(file_path)
data_path = os.path.abspath(os.path.join(test_dir_path,"test_data"))
tmp_dir_path = os.path.join(test_dir_path, 'tmp_output')

class TestCMD(unittest.TestCase):
    def setUp(self):
        """Create temporary dir if necessary,
        otherwise clear contents of it"""
        if not isdir(tmp_dir_path):
            os.mkdir(tmp_dir_path)
        self.tearDown()

        # Use a concoct run as input
        subprocess.check_call(
            "concoct --coverage_file {0}/coverage --composition_file {0}/composition.fa "
            "--basename {1}/ -c 10 2> /dev/null".format(data_path, tmp_dir_path),
            shell=True)
        self.cluster_file = os.path.join(tmp_dir_path, 'clustering_gt1000.csv')
        self.original_data = os.path.join(tmp_dir_path, 'original_data_gt1000.csv')
        self.clustering = p.read_csv(self.cluster_file, index_col=0)

        # Pretend the largest cluster has two copies of each single copy gene
        self.largest_cluster = self.clustering.cluster_id.value_counts().index[0]
        nr_clusters = self.clustering.cluster_id.max() + 1
        scg_freq = p.DataFrame(1, index=range(nr_clusters), columns=['COG1', 'COG2', 'COG3'])
        scg_freq.loc[self.largest_cluster] = 2
        self.scg_file = os.path.join(tmp_dir_path, 'scg_freq.csv')
        scg_freq.to_csv(self.scg_file)

    def tearDown(self):
        """remove temporary output files"""
        for f in os.listdir(tmp_dir_path):
            os.remove(os.path.join(tmp_dir_path, f))
        assert os.listdir(tmp_dir_path) == []

    def run_command(self, cluster_file):
        # The result is written to the working directory
        p_refine = subprocess.Popen(["concoct_refine", cluster_file, self.original_data, self.scg_file],
                                    cwd=tmp_dir_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, err = p_refine.communicate()
        self.c = p_refine.returncode
        self.err = err.decode('utf-8')

    def check_refined_clustering(self):
        self.assertEqual(self.c, 0, msg=self.err)
        refined = p.read_csv(os.path.join(tmp_dir_path, 'clustering_refine.csv'), index_col=0)
        self.assertEqual(list(refined.index), list(self.clustering.index))
        refined = refined.iloc[:, 0]

        in_largest = self.clustering.cluster_id == self.largest_cluster
        # Only the cluster with duplicated single copy genes is touched
        self.assertTrue((refined[~in_largest] == self.clustering.cluster_id[~in_largest]).all())
        # and that one is split into several clusters
        self.assertTrue(refined[in_largest].nunique() > 1)
        self.assertTrue(refined.nunique() > self.clustering.cluster_id.nunique())

    def test_cluster_file_with_header(self):
        self.run_command(self.cluster_file)
        self.check_refined_clustering()

    def test_cluster_file_without_header(self):
        no_header = os.path.join(tmp_dir_path, 'clustering_no_header.csv')
        self.clustering.to_csv(no_header, header=False)
        self.run_command(no_header)
        self.check_refined_clustering()
