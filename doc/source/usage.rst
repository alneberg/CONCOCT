
Basic Usage
===========

This guide assumes you have your original contigs assembled into a file ``original_contigs.fa`` and that you have mapped reads from several samples to these contigs into ``.bam`` files.
Note that the assembly can be constructed using either one single sample or several (usually all) samples.
In either case, all sample reads should be mapped against the assembly to achieve the best binning performance.

The next step is then to cut contigs into smaller parts::

    cut_up_fasta.py original_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa


Generate table with coverage depth information per sample and subcontig.
This step assumes the directory 'mapping' contains sorted and indexed bam files where each sample has been mapped against the original contigs::

    concoct_coverage_table.py contigs_10K.bed mapping/Sample*.sorted.bam > coverage_table.tsv

The bam files need to be sorted and indexed, which for each file is done with::

    samtools sort mapping/Sample1.bam -o mapping/Sample1.sorted.bam
    samtools index mapping/Sample1.sorted.bam

The index needs to be located next to the bam file with ``.bai`` added to the complete file name, in this case ``mapping/Sample1.sorted.bam.bai``.


Run concoct::

    concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/


Merge subcontig clustering into original contig clustering::

    merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv


Extract bins as individual FASTA::

    mkdir concoct_output/fasta_bins
    extract_fasta_bins.py original_contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins

These bins should now be evaluated and filtered for completeness and contamination using for example `CheckM <https://ecogenomics.github.io/CheckM/>`__ or `BUSCO <https://busco.ezlab.org/>`__.


Output Files
------------

With ``-b concoct_output/`` the following files are created in the directory ``concoct_output``.
If the basename is not an existing directory and does not end with a ``/``, it is instead used as a prefix for the file names.
The ``1000`` in the file names is the length threshold given with ``--length_threshold``.

``clustering_gt1000.csv``
    The result of the clustering. One line per contig with the contig id and the id of the cluster it was assigned to.

``original_data_gt1000.csv``
    The input data for all contigs passing the length threshold, as it looks before the PCA.
    The columns are the log transformed kmer frequencies followed by the log transformed coverage for each sample.
    This file can be big and is not written if ``--no_original_data`` is used.

``PCA_transformed_data_gt1000.csv``
    The data which the clustering is performed on, the contigs in the coordinates of the principal components.

``PCA_components_data_gt1000.csv``
    The principal components, one per line, as weights for the columns of the original data.

``args.txt``
    The arguments the run was started with.

``log.txt``
    Timestamps for the steps of the run.

While running, the progress of the clustering is written to stderr with one line per iteration,
containing the iteration number, the variational lower bound and how much it changed since the previous iteration.
