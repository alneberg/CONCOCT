#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os
from distutils.core import Extension

version = '0.1'

module1 = Extension('vbgmm',
        libraries =['gsl',  'gslcblas'],
        sources = ['c-concoct/vbgmmmodule.c'])

setup(name='concoct',
      version=version,
      description="Clustering cONtigs with COverage and ComposiTion",
      long_description="""\
To be done""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python Scilifelab Metagenomics Binning Clustering Contig',
      author='Brynjar Smari Bjarnason, Johannes Alneberg, Christopher Quince, Anders Andersson, Ino de Bruijn',
      author_email='binni@binnisb.com',
      maintainer='Johannes Alneberg',
      maintainer_email='johannes.alneberg@scilifelab.se',
      url='www.github.com/BinPro/CONCOCT',
      license='FreeBSD',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts=["bin/concoct"],
      include_package_data=True,
      zip_safe=False,
      ext_modules = [module1],
      install_requires=['argparse==1.2.1',
                        'numpy==1.7.1',
                        'scipy==0.13.0'],
                        'pandas==0.12.0',
                        'scikit-learn==0.14.1',
                        'biopython==1.62',
                        'nose==1.3.0'],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )

