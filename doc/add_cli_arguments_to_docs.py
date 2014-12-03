#!/usr/bin/env python
import os
import subprocess
from argparse import ArgumentParser

TEMPLATE = """
Usage
=====

CONCOCT uses different subcommands to accomplish different things. To find out which subcommands are available, run ``concoct -h`` on the command line::

"""

PARSE_TEMPLATE = """
Parse
-----
The `parse` subcommand is used to generate the input table. This functionality was previously included in the script 'gen_input_table.py'. To see the options run ``concoct parse -h``::

"""

CLUSTERING_TEMPLATE = """
Cluster
------
For the subcommand `cluster` concoct has several command line options to control the clustering, here is a complete documentation of these. These can also be viewed by typing ``concoct cluster -h`` on the command line.::

"""

def indent(s):
    return "\n".join(map(lambda w: '\t'+w, s.splitlines()))


def help_doc_rst(script):
    """ Fetch the --help info from the script, outputs rst formatted string """
    sp = subprocess.Popen(script + ["-h"],
            stdout=subprocess.PIPE)
    stdout,stderr = sp.communicate()

    # Add help message to template
    return "{0}\n\n".format(indent(stdout))

def main(main_template, parse_template, clustering_template):
    file_path = os.path.dirname(os.path.realpath(__file__))
    main_template += help_doc_rst(["concoct"])

    main_template += parse_template + help_doc_rst(["concoct", "parse"])
    main_template += clustering_template + help_doc_rst(["concoct", "cluster"])

    # Print the help message to a sphinx (restructured text) markup file
    docs_path = os.path.join(file_path, 'source', 'usage.rst')
    with open(docs_path,'w') as doc_f:
        doc_f.write(main_template)


if __name__ == "__main__":
    # Argumentparser only to add --help option
    args = ArgumentParser(description=("Generates basic command line documentation for "
                                       "concoct. This is used since the hassle of getting "
                                       "readthedocs.org to install concoct within a virtualenv "
                                       "felt like the worse option.")).parse_args()

    main(TEMPLATE, PARSE_TEMPLATE, CLUSTERING_TEMPLATE)

    
