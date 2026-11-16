#!/usr/bin/env python
import unittest
from os.path import isdir
import os
import subprocess
from Bio import SeqIO

file_path = os.path.realpath(__file__)
test_dir_path = os.path.dirname(file_path)
tmp_dir_path = os.path.join(test_dir_path, 'tmp_output')
script_path = os.path.join(test_dir_path, '..', 'scripts', 'extract_fasta_bins.py')

class TestCMD(unittest.TestCase):
    def setUp(self):
        """Create temporary dir if necessary,
        otherwise clear contents of it"""
        if not isdir(tmp_dir_path):
            os.mkdir(tmp_dir_path)
        self.tearDown()
        self.fasta_file = os.path.join(tmp_dir_path, 'contigs.fa')
        self.cluster_file = os.path.join(tmp_dir_path, 'clustering_merged.csv')

    def tearDown(self):
        """remove temporary output files"""
        for f in os.listdir(tmp_dir_path):
            os.remove(os.path.join(tmp_dir_path, f))
        assert os.listdir(tmp_dir_path) == []

    def write_input(self, contig_ids, clustering):
        with open(self.fasta_file, 'w') as fh:
            for contig_id in contig_ids:
                fh.write(">{}\nACGTACGT\n".format(contig_id))
        with open(self.cluster_file, 'w') as fh:
            fh.write("contig_id,cluster_id\n")
            for contig_id, cluster_id in clustering:
                fh.write("{},{}\n".format(contig_id, cluster_id))

    def run_command(self):
        call = ["python", script_path, self.fasta_file, self.cluster_file, "--output_path", tmp_dir_path]
        p = subprocess.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, err = p.communicate()
        self.c = p.returncode
        self.err = err.decode('utf-8')

    def bin_content(self, cluster_id):
        with open(os.path.join(tmp_dir_path, "{}.fa".format(cluster_id))) as fh:
            return [rec.id for rec in SeqIO.parse(fh, "fasta")]

    def test_basic(self):
        self.write_input(['contig_a', 'contig_b', 'contig_c', 'unclustered'],
                         [('contig_a', 0), ('contig_b', 1), ('contig_c', 0)])
        self.run_command()
        self.assertEqual(self.c, 0, msg=self.err)
        self.assertEqual(self.bin_content(0), ['contig_a', 'contig_c'])
        self.assertEqual(self.bin_content(1), ['contig_b'])

    def test_numeric_contig_ids(self):
        self.write_input(['1', '2', '0003'], [('1', 0), ('2', 1), ('0003', 0)])
        self.run_command()
        self.assertEqual(self.c, 0, msg=self.err)
        self.assertEqual(self.bin_content(0), ['1', '0003'])
        self.assertEqual(self.bin_content(1), ['2'])

    def test_cutup_clustering_with_original_contigs(self):
        self.write_input(['contig_a', 'contig_b'],
                         [('contig_a.concoct_part_0', 0), ('contig_b.concoct_part_0', 1)])
        self.run_command()
        self.assertNotEqual(self.c, 0)
        self.assertTrue("2 out of 2 contigs" in self.err, msg=self.err)
        self.assertTrue("merge_cutup_clustering.py" in self.err, msg=self.err)
        self.assertFalse("Traceback" in self.err, msg=self.err)
