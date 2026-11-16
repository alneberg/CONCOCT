#!/usr/bin/env python
import unittest
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
import pandas as p

file_path = os.path.realpath(__file__)
test_dir_path = os.path.dirname(file_path)
tmp_dir_path = test_dir_path + '/tmp_output'
tmp_basename_dir = tmp_dir_path + '/1'
tmp_basename_dir2 = tmp_dir_path + '/2'
tmp_basename_file = tmp_dir_path + '/file'

CWD = os.getcwd()

class TestCMD(unittest.TestCase):
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
                    tags=[],basename='tmp_output/1'):
        call_string = "concoct --coverage_file test_data/{0} --composition_file test_data/{1} --basename {2} -c 10 --no_total_coverage 2> /dev/null".format(cov_file,comp_file,basename)
        for tag in tags:
            call_string += " " + tag
        self.c = 0 # Exit code
        try:
            self.op = subprocess.check_output(
                call_string,
                shell=True)
        except subprocess.CalledProcessError as exc:
            self.c = exc.returncode

    def file_len(self,fh):
        i=0
        with open(fh) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def md5sum(self,fh):
        infile = open("filename", 'rb')
        content = infile.read()
        infile.close()
        m = hashlib.md5()
        m.update(content)
        return m.hexdigest()

    def test_no_errors(self):
        self.run_command()
        self.assertEqual(self.c,0,
                     msg = "Command exited with nonzero status")

    def test_directory_creation(self):
        self.run_command()
        self.assertTrue(isdir(tmp_basename_dir),
                    msg = "Temporary directory not created")
        m_time_first = os.path.getmtime(tmp_basename_dir+'/clustering_gt1000.csv')

        # Rerun the concoct and see that the directory is overwritten
        self.run_command()
        m_time_second = os.path.getmtime(tmp_basename_dir+'/clustering_gt1000.csv')
        self.assertTrue(m_time_first != m_time_second,
                     msg = "basename dir is not overwritten")
        L = listdir(tmp_dir_path)
        self.assertTrue(len(L) == 1,
                    msg = "Multiple output directories or files was created")

        # File creation
        self.run_command(basename=tmp_basename_file)
        self.assertTrue(isfile(tmp_basename_file+'_clustering_gt1000.csv'),
                    msg = "Clustering file is not created, when file is used as basename")
        L = listdir(tmp_basename_dir)
        self.assertTrue(len(L) == 6,
                    msg = "Wrong number of output files, observed {0}".format(L))

    def test_prior_to_clustering(self):
        self.run_command()
        d_p = os.path.join(tmp_basename_dir)
        self.assertTrue(isfile(d_p+ '/args.txt'),
                           msg="Args file is not created")
        self.assertTrue(isfile(d_p+ '/log.txt'),
                           msg="Log file is not created")
        self.assertTrue(isfile(d_p+ '/original_data_gt1000.csv'),
                           msg="Original data file is not created")
        self.assertTrue(isfile(d_p+ '/PCA_transformed_data_gt1000.csv'),
                           msg="PCA transformed data file is not created")


    def test_output_files_creation(self):
        # dir as basename
        self.run_command()
        d_p = os.path.join(tmp_basename_dir)
        self.assertTrue(
            isfile(d_p+ '/clustering_gt1000.csv'),
            msg='Large contigs clustering file is not created'
            )
        self.assertTrue(
            isfile(d_p+ '/PCA_transformed_data_gt1000.csv'),
            msg='PCA file is not created'
            )
        self.assertTrue(
            isfile(d_p+ '/original_data_gt1000.csv'),
            msg='Original data file is not created'
            )
        self.assertTrue(
            isfile(d_p+ '/log.txt'),
            msg='Log file is not created'
            )
        # dir as file
        self.run_command(basename=tmp_basename_file)
        d_p = tmp_basename_file +'_'
        self.assertTrue(
            isfile(d_p+ 'clustering_gt1000.csv'),
            msg='Large contigs clustering file is not created'
            )
        self.assertTrue(
            isfile(d_p+ 'PCA_transformed_data_gt1000.csv'),
            msg='PCA file is not created'
            )
        self.assertTrue(
            isfile(d_p+ 'original_data_gt1000.csv'),
            msg='Original data file is not created'
            )
        self.assertTrue(
            isfile(d_p+ 'log.txt'),
            msg='Log file is not created'
            )

    def test_threshold_functionality(self):
        self.run_command()
        d_p = tmp_basename_dir
        od_1 = d_p+'/original_data_gt1000.csv'
        clust_gt_1 = d_p+'/clustering_gt1000.csv'
        odl_1 = self.file_len(od_1)
        clust_gtl_1= self.file_len(clust_gt_1)

        self.run_command(comp_file='composition_some_shortened.fa',
                         basename=tmp_basename_dir2+'/')
        d_p2 = tmp_basename_dir2
        od_2 = d_p2+'/original_data_gt1000.csv'
        clust_gt_2 = d_p2+'/clustering_gt1000.csv'
        odl_2 = self.file_len(od_2)
        clust_gtl_2= self.file_len(clust_gt_2)

        self.assertTrue(odl_1!=odl_2,
                    msg='Original data have the same lengths')
        self.assertTrue(clust_gtl_1!=clust_gtl_2,
                    msg='Filtered clustering files have the same lengths')

    def test_logging(self):
        self.run_command()
        with open(tmp_basename_dir+'/log.txt','r') as log:
            log_content = log.read()
            self.assertTrue(len(log_content)>10,
                        "Log content is too small")
            pca_report = [row for row in log_content.split('\n') if 'Performed PCA, resulted in ' in row][0]
            pca_dimensions_log = int(pca_report.split()[-2])
            with open(tmp_basename_dir+'/PCA_transformed_data_gt1000.csv', 'r') as pca_comps:
                header = pca_comps.readlines()[0]
                header = header.strip()
                last_dim = int(header.split(',')[-1])
                pca_dimensions = last_dim + 1
            self.assertEqual(pca_dimensions, pca_dimensions_log)

    def test_seed(self):
        #Test default behaviour, seed = 11
        self.run_command()
        first_time = os.path.getmtime(tmp_basename_dir+'/clustering_gt1000.csv')
        with open(tmp_basename_dir+'/clustering_gt1000.csv','r') as clustering:
            first_file=clustering.read()

        self.run_command()
        second_time = os.path.getmtime(tmp_basename_dir+'/clustering_gt1000.csv')
        with open(tmp_basename_dir+'/clustering_gt1000.csv','r') as clustering:
            second_file=clustering.read()
        self.assertTrue(not (first_time==second_time),
                    msg='clustering_gt1000.csv did not change')
        self.assertTrue(first_file == second_file,
                    msg='Clustering outcomes were not the same with same seeds')

        #Should be equal to both above since default seed is 1
        self.run_command(tags=["--seed","1"])
        first_time = os.path.getmtime(tmp_basename_dir+'/clustering_gt1000.csv')
        with open(tmp_basename_dir+'/clustering_gt1000.csv','r') as clustering:
            first_file=clustering.read()
        self.assertTrue(not (first_time==second_time),
                    msg='clustering_gt1000.csv did not change')
        self.assertTrue(first_file == second_file,
                    msg='Clustering outcomes were not the same with same seeds')

        #Test that 0 gives different seed
        self.run_command(tags=['--seed','0'])
        first_time = os.path.getmtime(tmp_basename_dir+'/clustering_gt1000.csv')
        with open(tmp_basename_dir+'/clustering_gt1000.csv','r') as clustering:
            first_file=clustering.read()


        #Should give different clustering
        self.run_command(tags=['--seed','0'])
        second_time = os.path.getmtime(tmp_basename_dir+'/clustering_gt1000.csv')
        with open(tmp_basename_dir+'/clustering_gt1000.csv','r') as clustering:
            second_file=clustering.read()
        self.assertTrue(not (first_time==second_time),
                    msg='clustering_gt1000.csv did not change')
        self.assertTrue(not (first_file == second_file),
                    msg='Clustering outcomes were the same with random seeds')


        #Test that two differnet seeds give different clustering
        #Should give clustering 2
        self.run_command(tags=['--seed','2'])
        first_time = os.path.getmtime(tmp_basename_dir+'/clustering_gt1000.csv')
        with open(tmp_basename_dir+'/clustering_gt1000.csv','r') as clustering:
            first_file=clustering.read()

        #Should give clustering 3
        self.run_command(tags=['--seed','3'])
        second_time = os.path.getmtime(tmp_basename_dir+'/clustering_gt1000.csv')
        with open(tmp_basename_dir+'/clustering_gt1000.csv','r') as clustering:
            second_file=clustering.read()
        self.assertTrue(not (first_time==second_time),
                    msg='clustering_gt1000.csv did not change')
        self.assertTrue(not (first_file == second_file),
                    msg='Clustering outcomes were the same with different seeds')

    def test_log_coverage(self):
        self.run_command()
        original_coverage_data_path = os.path.join(tmp_basename_dir,'original_data_gt1000.csv')
        df = p.io.parsers.read_csv(original_coverage_data_path,index_col=0,sep=',')

        true_pseudo_cov = -1.3143
        calc_pseudo_cov = df.sample_1.iloc[0]
        self.assertAlmostEqual(true_pseudo_cov,calc_pseudo_cov,places=4)

    def test_log_coverage_no_cov_normalization(self):
        self.run_command(tags=["--no_cov_normalization"])
        original_coverage_data_path = os.path.join(tmp_basename_dir,'original_data_gt1000.csv')
        df = p.io.parsers.read_csv(original_coverage_data_path,index_col=0,sep=',')

        true_pseudo_cov = -1.8107
        calc_pseudo_cov = df.sample_1.iloc[0]
        self.assertAlmostEqual(true_pseudo_cov,calc_pseudo_cov,places=4)


    def test_big_file_validation(self):
        """ Run Validate.pl on the result files after running a larger input
        file and make sure the statistics are good enough. """
        self.run_command(cov_file='large_contigs/coverage_table.tsv',
                         comp_file='large_contigs/contigs.fa',
                         basename=os.path.join(tmp_dir_path, 'large_contigs/'))

        validate_path = os.path.join(test_dir_path, '..', 'scripts', 'Validate.pl')
        clustering_reference = os.path.join(test_dir_path, 'test_data', 'large_contigs',
                                            'clustering_gt1000_taxassign.csv')
        clustering_file = os.path.join(tmp_dir_path,'large_contigs',
                                       'clustering_gt1000.csv')

        self.assertTrue(isfile(validate_path))
        self.assertTrue(isfile(clustering_reference))
        self.assertTrue(isfile(clustering_file))
        validate_so = subprocess.check_output(['perl', validate_path,
                                               '--sfile={}'.format(clustering_reference),
                                               '--cfile={}'.format(clustering_file) ])
        print("Results for large clustering file: ")
        print(validate_so)

        headers = validate_so.split(b'\n')[0].split(b'\t')
        stats = validate_so.split(b'\n')[1].split(b'\t')
        stats_dict = dict(list(zip(headers, stats)))

        self.assertTrue(float(stats_dict[b'AdjRand']) > 0.85,
                    msg=("Insufficient adjusted rand index "
                         "reached, requires > 0.85"))
        self.assertTrue(float(stats_dict[b'Prec.']) > 0.95,
                    msg=("Insufficient precision reached, "
                         "requires > /0.95"))
        self.assertTrue(float(stats_dict[b'Rec.']) > 0.90,
                    msg=("Insufficient recall reached, "
                         "requires > 0.90"))

        conf_file = os.path.join(test_dir_path, 'Conf.csv')
        if isfile(conf_file):
            os.remove(conf_file)

    def test_numeric_contig_ids(self):
        """Contig ids consisting of digits only should not be read as integers"""
        cov_file = os.path.join(tmp_dir_path, 'numeric_coverage')
        comp_file = os.path.join(tmp_dir_path, 'numeric_composition.fa')
        new_ids = {}
        with open(os.path.join(test_dir_path, 'test_data', 'composition.fa')) as ifh, open(comp_file, 'w') as ofh:
            for line in ifh:
                if line.startswith('>'):
                    contig_id = line[1:].split()[0]
                    new_ids[contig_id] = "{:04d}".format(len(new_ids))
                    line = ">{}\n".format(new_ids[contig_id])
                ofh.write(line)
        with open(os.path.join(test_dir_path, 'test_data', 'coverage')) as ifh, open(cov_file, 'w') as ofh:
            ofh.write(ifh.readline())
            for line in ifh:
                contig_id, rest = line.split('\t', 1)
                ofh.write("{}\t{}".format(new_ids[contig_id], rest))

        self.run_command()
        expected = p.read_csv(tmp_basename_dir+'/clustering_gt1000.csv')

        self.run_command(cov_file='../tmp_output/numeric_coverage',
                         comp_file='../tmp_output/numeric_composition.fa',
                         basename=tmp_basename_dir2+'/')
        self.assertEqual(self.c, 0,
                     msg = "Command exited with nonzero status")
        observed = p.read_csv(tmp_basename_dir2+'/clustering_gt1000.csv', dtype={'contig_id': str})
        self.assertEqual(list(observed.contig_id), [new_ids[c] for c in expected.contig_id])
        self.assertEqual(list(observed.cluster_id), list(expected.cluster_id))

    def test_one_contig_threshold(self):
        """Make sure we don't execute clustering of 0 or 1 contig"""
        # Make sure the error code is not set before running command
        self.assertFalse(hasattr(self,"c"))
        # Longest contig is 33356 so we put the threshold just below
        self.run_command(tags=["--length_threshold 33350"])
        # The command should have failed with code 255
        self.assertTrue(hasattr(self,"c"))
        self.assertEqual(self.c,255)
