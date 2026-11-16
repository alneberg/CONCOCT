#!/usr/bin/env python
import unittest
from os.path import isdir
import os
import shutil
import subprocess
import pandas as pd

file_path = os.path.realpath(__file__)
test_dir_path = os.path.dirname(file_path)
data_path = os.path.abspath(os.path.join(test_dir_path,"test_data","map"))
tmp_dir_path = os.path.join(test_dir_path, 'tmp_output')
script_path = os.path.join(test_dir_path, '..', 'scripts', 'concoct_coverage_table.py')

ten_reads_bam = os.path.join(data_path, 'ten_reads', 'ten_reads_two_contigs-s.bam')
twenty_reads_bam = os.path.join(data_path, 'twenty_reads', 'twenty_reads_two_contigs-s.bam')

class TestCMD(unittest.TestCase):
    def setUp(self):
        """Create temporary dir if necessary,
        otherwise clear contents of it"""
        if not isdir(tmp_dir_path):
            os.mkdir(tmp_dir_path)
        self.tearDown()
        self.bedfile = os.path.join(tmp_dir_path, 'contigs.bed')
        self.write_bedfile([('contig-75000034', 1615), ('contig-21000001', 9998)])

    def tearDown(self):
        """remove temporary output files"""
        for f in os.listdir(tmp_dir_path):
            os.remove(os.path.join(tmp_dir_path, f))
        assert os.listdir(tmp_dir_path) == []

    def write_bedfile(self, contigs):
        with open(self.bedfile, 'w') as fh:
            for contig_id, length in contigs:
                fh.write("{0}\t0\t{1}\t{0}.concoct_part_0\n".format(contig_id, length))

    def run_command(self, bam_files, tags=[]):
        self.output_file = os.path.join(tmp_dir_path, 'coverage_table.tsv')
        call = ["python", script_path, self.bedfile] + bam_files + tags
        with open(self.output_file, 'w') as ofh:
            p = subprocess.Popen(call, stdout=ofh, stderr=subprocess.PIPE)
            _, err = p.communicate()
        self.c = p.returncode
        self.err = err.decode('utf-8')

    def test_basic(self):
        self.run_command([ten_reads_bam, twenty_reads_bam])
        self.assertEqual(self.c, 0,
                     msg = "Command exited with nonzero status")

        df = pd.read_csv(self.output_file, sep='\t', index_col=0)
        self.assertEqual(list(df.columns), ['cov_mean_sample_ten_reads_two_contigs-s',
                                            'cov_mean_sample_twenty_reads_two_contigs-s'])
        ten, twenty = df.columns
        self.assertAlmostEqual(df[ten].loc['contig-75000034.concoct_part_0'], 10*100.0/1615, 3)
        self.assertAlmostEqual(df[ten].loc['contig-21000001.concoct_part_0'], 10*100.0/9998, 3)
        self.assertAlmostEqual(df[twenty].loc['contig-75000034.concoct_part_0'], 20*100.0/1615, 3)
        self.assertAlmostEqual(df[twenty].loc['contig-21000001.concoct_part_0'], 20*100.0/9998, 3)

    def test_samplenames(self):
        # No newline after the last sample name
        samplenames = os.path.join(tmp_dir_path, 'sample_names')
        with open(samplenames, 'w') as fh:
            fh.write("first\nsecond")
        self.run_command([ten_reads_bam, twenty_reads_bam], tags=["--samplenames", samplenames])
        self.assertEqual(self.c, 0,
                     msg = "Command exited with nonzero status")
        df = pd.read_csv(self.output_file, sep='\t', index_col=0)
        self.assertEqual(list(df.columns), ['cov_mean_sample_first', 'cov_mean_sample_second'])

    def test_missing_index(self):
        bam_file = os.path.join(tmp_dir_path, 'no_index.bam')
        shutil.copy(ten_reads_bam, bam_file)
        self.run_command([bam_file])
        self.assertNotEqual(self.c, 0)
        self.assertTrue("No index was found for '{}'".format(bam_file) in self.err, msg=self.err)
        self.assertFalse("Traceback" in self.err, msg=self.err)

    def test_missing_bamfile(self):
        self.run_command([os.path.join(tmp_dir_path, 'does_not_exist.bam')])
        self.assertNotEqual(self.c, 0)
        self.assertTrue("Could not read the header" in self.err, msg=self.err)
        self.assertFalse("Traceback" in self.err, msg=self.err)

    def test_contig_missing_from_bamfile(self):
        self.write_bedfile([('contig-75000034', 1615), ('other_contig', 9998)])
        self.run_command([ten_reads_bam])
        self.assertNotEqual(self.c, 0)
        self.assertTrue("1 out of 2 contigs" in self.err, msg=self.err)
        self.assertTrue("'other_contig'" in self.err, msg=self.err)
        self.assertFalse("Traceback" in self.err, msg=self.err)

    def test_bamfiles_with_different_contig_order(self):
        """samtools bedcov would silently report the coverage for the wrong contig"""
        header = subprocess.check_output(["samtools", "view", "-H", ten_reads_bam]).decode('utf-8')
        sq_lines = [l for l in header.splitlines() if l.startswith('@SQ')]
        other_lines = [l for l in header.splitlines() if not l.startswith('@SQ')]
        self.assertEqual(len(sq_lines), 2)

        new_header = os.path.join(tmp_dir_path, 'header.sam')
        with open(new_header, 'w') as fh:
            fh.write("\n".join(other_lines[:1] + sq_lines[::-1] + other_lines[1:]) + "\n")
        bam_file = os.path.join(tmp_dir_path, 'reordered.bam')
        with open(bam_file, 'wb') as fh:
            subprocess.check_call(["samtools", "reheader", new_header, ten_reads_bam], stdout=fh)
        subprocess.check_call(["samtools", "index", bam_file])

        self.run_command([ten_reads_bam, bam_file])
        self.assertNotEqual(self.c, 0)
        self.assertTrue("not identical and in the same order" in self.err, msg=self.err)
