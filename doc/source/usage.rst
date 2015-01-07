
Usage
=====

CONCOCT uses different subcommands to accomplish different things. To find out which subcommands are available, run ``concoct -h`` on the command line:

.. program-output:: (echo 'import conf'; cat ../../concoct/parser.py; echo 'args=arguments()') | python - --help
   :shell:


Parse
-----
The `parse` subcommand is used to generate the input table. This functionality was previously included in the script 'gen_input_table.py'. To see the options run ``concoct parse -h``:

.. program-output:: (echo 'import conf'; cat ../../concoct/parser.py; echo 'args=arguments()') | python - parse --help
   :shell:


Cluster
------
For the subcommand `cluster` concoct has several command line options to control the clustering, here is a complete documentation of these. These can also be viewed by typing ``concoct cluster -h`` on the command line.:

.. program-output:: (echo 'import conf'; cat ../../concoct/parser.py; echo 'args=arguments()') | python - cluster --help
   :shell:

