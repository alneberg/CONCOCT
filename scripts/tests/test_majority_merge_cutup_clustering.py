import sys
import os
from os.path import join as ospj
from nose.tools import assert_equal, ok_
import pandas as pd
import collections
import subprocess 
import concoct.utils.dir_utils as dir_utils

if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

FILE_PATH = os.path.realpath(__file__)
TEST_DIR_PATH = os.path.dirname(FILE_PATH)
DATA_PATH = os.path.abspath(ospj(TEST_DIR_PATH, "test_data", "majority_merge"))
SCRIPT_PATH = ospj(TEST_DIR_PATH, '..')

# Add script dir to python path to import functions
sys.path.append(SCRIPT_PATH)
from extract_scg_bins import get_approved_bins, sum_bases_in_bins, \
    get_winning_bins, write_approved_bins, main

CWD = os.getcwd()


# DATA
"""
# Basic case where majority is in one bin
contig1.0, 1
contig1.1, 0
contig1,2, 0
contig2.0, 1
contig2.1, 1
contig2.2, 0

# Basic case where no majority exists
contig1.0, 1
contig1.1, 0
contig2.0, 1
contig2.1, 0
contig2.2, 0
contig2.3, 1
"""

class TestDnaDiff(object):
    def run_majority_merge_cutup_clustering(self, input_path):
        """Run the majority merge cutup clustering script on the given path
        and returns the captured stdout as a dataframe."""
        call_string = ["python", "../majority_merge_cutup_clustering.py", input_path]
        output = subprocess.check_output(call_string)
        return pd.DataFrame.from_csv(StringIO(output))

    def test_simplest_case(self):
        """Test simplest case of majority merge cutup clustering"""
        # Run script and capture stdout 
        df = self.run_majority_merge_cutup_clustering(ospj(DATA_PATH, "basic_case_all_in_same.tsv"))
        # assert correct clustering
        assert_equal(2, len(df.index))

