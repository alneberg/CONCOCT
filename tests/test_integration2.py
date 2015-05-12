#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal, nottest, assert_false
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
import pandas as p

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))
test_dir_path = os.path.dirname(file_path)
tmp_dir_path = test_dir_path + '/nose_tmp_output'
tmp_basename_dir = tmp_dir_path + '/1'
tmp_basename_file = tmp_dir_path + '/file'

CWD = os.getcwd()

class TestCMD(object):
    def setUp(self):
        """Create temporary dir if necessary,
        otherwise clear contents of it"""
        if not isdir(tmp_dir_path):
            os.mkdir(tmp_dir_path)
        self.tearDown()
        os.mkdir(tmp_basename_dir)
        os.chdir(test_dir_path)

    def tearDown(self):
        """remove temporary output files"""
        for d in os.listdir(tmp_dir_path):
            d_path = os.path.join(tmp_dir_path,d)
            try:
                os.remove(d_path)
            except:
                for f in os.listdir(d_path):
                    f_path = os.path.join(d_path,f)
                    os.remove(f_path)
                os.rmdir(d_path)
        assert os.listdir(tmp_dir_path) == []


    def run_command(self,cov_file='coverage',comp_file='composition.fa',
                basename='nose_tmp_output/1'):
        print(sys.executable)
        print(sys.path)
        print(subprocess.check_output(['which', 'concoct']))
        call_string = "concoct --coverage_file test_data/{0} --composition_file test_data/{1} --basename {2} -c 10 --no_total_coverage".format(cov_file,comp_file,basename)
        self.c = 0 # Exit code
        print(call_string)
        try:
            self.op = subprocess.check_output(
                call_string,
                shell=True)
        except subprocess.CalledProcessError as exc:
            self.c = exc.returncode
        print(self.op)

    def test_big_file_validation(self):
        outdir = os.path.join(tmp_dir_path, 'large_contigs/')
        self.run_command(cov_file='large_contigs/coverage_table.tsv', 
                         comp_file='large_contigs/contigs.fa',
                         basename= outdir)

        L = listdir(outdir)
        assert_true(len(L) == 26,
                    msg = "Wrong number of output files, observed {0}".format(L))
        
