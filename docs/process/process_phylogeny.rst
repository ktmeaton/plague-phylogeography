Phylogeny
***************************

IQ-TREE
-------

Maximum likelihood tree search and model selection, iqtree phylogeny.

========================================= =========================== ===========================
Input                                     Type                        Description
========================================= =========================== ===========================
ch_snippy_core_filter_iqtree              fasta                       Multi fasta of filtered core genome sites from process :ref:`snippy_multi_filter<Snippy Multi Filter>`.
========================================= =========================== ===========================

========================================= =========================== ===========================
Output                                    Type                        Description
========================================= =========================== ===========================
ch_iqtree_treefile_augur_refine           newick                      Newick treefile phylogeny with branch supports for process :ref:`augur_refine<Augur Refine>`.
========================================= =========================== ===========================

========================================================= =========================== ===========================
Publish                                                   Type                        Description
========================================================= =========================== ===========================
iqtree.core-filter\*_bootstrap.treefile                   newick                      Newick treefile phylogeny with branch supports.
iqtree\*                                                  misc                        All default output of iqtree.
========================================================= =========================== ===========================

**Shell script**::

      # Remember to change outgroup here later
      # A thorough tree search for model selection can be done with -m MF -mtree
      iqtree \
        -s ${snippy_core_filter_aln} \
        -m MFP \
        -nt ${task.cpus} \
        -o ${params.iqtree_outgroup} \
        -seed ${params.iqtree_rng} \
        -pre iqtree.core-filter${params.snippy_multi_missing_data_text}_bootstrap \
        -bb 1000 \
        -alrt 1000 \
        2>&1 | tee iqtree.core-filter${params.snippy_multi_missing_data_text}_bootstrap.output
