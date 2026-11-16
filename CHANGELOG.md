# Changelog

A list of changes per version.
When changing something in the develop branch, it should be added here.

## [Unreleased]

### `Changed`

 - Python 2 is no longer supported, Python 3.10 or later is required
 - The tests are now run with pytest instead of the unmaintained nose
 - Continuous integration is moved from Travis CI to GitHub Actions
 - The build requirements (cython and numpy) are declared in `pyproject.toml` so that `pip install .` works in a clean environment
 - [#118](https://github.com/BinPro/CONCOCT/issues/118) - `concoct` exits with an error message when no composition file is given, the coverage file is optional
 - `concoct` no longer prints the debug line with the number of contigs and threads to stdout

### `Added`

 - [#284](https://github.com/BinPro/CONCOCT/pull/284) - The `--iterations` argument for `concoct_refine`
 - [#318](https://github.com/BinPro/CONCOCT/issues/318) - Documentation of the output files

### `Fixed`

 - [#321](https://github.com/BinPro/CONCOCT/issues/321) - `TypeError: Feature names are only supported if all input features have string names` with scikit-learn >= 1.2
 - `ModuleNotFoundError: No module named 'pkg_resources'` with recent setuptools versions
 - `nose` is no longer an install requirement
 - `extract_scg_bins.py` works with pandas >= 1.0
 - The Dockerfile builds again
 - [#300](https://github.com/BinPro/CONCOCT/issues/300) - `concoct_coverage_table.py` no longer hides the samtools error behind `TypeError: write() argument must be str, not bytes`, and points out missing bam indices
 - [#307](https://github.com/BinPro/CONCOCT/issues/307) - `concoct_coverage_table.py` exits with an explanation when contigs in the bedfile are missing from a bamfile
 - `concoct_coverage_table.py` refuses bamfiles with differing contig order in their headers, for which `samtools bedcov` silently reports the coverage of the wrong contigs
 - `concoct_coverage_table.py` handles a samplenames file without a trailing newline
 - [#150](https://github.com/BinPro/CONCOCT/issues/150), [#281](https://github.com/BinPro/CONCOCT/issues/281) - Contig ids consisting of digits only are now handled by `concoct` and `extract_fasta_bins.py`
 - [#312](https://github.com/BinPro/CONCOCT/issues/312) - `extract_fasta_bins.py` explains what is wrong when the contigs of the cluster file are not found in the fasta file
 - `cut_up_fasta.py` wrote the wrong coordinates to the bedfile when a non zero overlap was used
 - [#292](https://github.com/BinPro/CONCOCT/issues/292), [#330](https://github.com/BinPro/CONCOCT/issues/330) - `concoct` used more cores than given with `--threads`
 - [#316](https://github.com/BinPro/CONCOCT/issues/316), [#329](https://github.com/BinPro/CONCOCT/issues/329) - Endless `OpenBLAS Warning : Detect OpenMP Loop and this application may hang`, which also made the clustering many times slower
 - [#324](https://github.com/BinPro/CONCOCT/issues/324) - `concoct_refine` failed with `fit() takes exactly 5 positional arguments (4 given)`
 - `concoct_refine` accepts the cluster file as written by `concoct`, with a header line
 - `concoct` failed with pandas >= 3 when the coverage file only contained integer values

## [v1.1.0] 2019-08-02

### `Changed`

 - [#236](https://github.com/BinPro/CONCOCT/pull/236) - Always add suffix to contigs at cutup, even when they are not cut.
 - [#254](https://github.com/BinPro/CONCOCT/pull/254) - Slight cleanup of concoct refine
 - [#258](https://github.com/BinPro/CONCOCT/pull/258) - New suffices (.concoct_part_XX) are now used for contig parts 
 - [#261](https://github.com/BinPro/CONCOCT/pull/261) - Epsilon argument removed as it was not working and is not very useful
 - [#262](https://github.com/BinPro/CONCOCT/pull/262) - Rewrote documentation, including installation instructions
 - [#264](https://github.com/BinPro/CONCOCT/pull/264) - `concoct_part_` suffix is enforced in subcontig for coverage script 
 - [#264](https://github.com/BinPro/CONCOCT/pull/264) - Header line is enforced for input for `merge_cutup_clustering.py` and `extract_fasta_bins.py`
 - [#267](https://github.com/BinPro/CONCOCT/pull/267) - Updated documentation

### `Added`

 - [#253](https://github.com/BinPro/CONCOCT/pull/253) - A dockerfile useful to test the conda installation
 - [#258](https://github.com/BinPro/CONCOCT/pull/258) - Tests for all fundamental scripts, including a new integration test data repository
 - [#259](https://github.com/BinPro/CONCOCT/pull/259) - This changelog
 - [#262](https://github.com/BinPro/CONCOCT/pull/262) - Added documentation for the core scripts used with concoct
 - [#265](https://github.com/BinPro/CONCOCT/pull/265) - A warning is now printed when concoct runs in single threaded mode

### `Fixed`

 - [#230](https://github.com/BinPro/CONCOCT/pull/230) - Enable at least single threaded installation on Mac OSX
 - [#231](https://github.com/BinPro/CONCOCT/pull/231) - Replace pandas .ix with .loc to fix deprecation warnings
 - [#246](https://github.com/BinPro/CONCOCT/pull/246) - Limit some dependency version numbers for python 2
 - [#254](https://github.com/BinPro/CONCOCT/pull/254) - Concoct refine now works with python 3
 - [#258](https://github.com/BinPro/CONCOCT/pull/258) - Seed tests now working again
 - [#260](https://github.com/BinPro/CONCOCT/pull/260) - Fix the dockerfile build by adding integration test data
