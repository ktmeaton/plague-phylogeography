# v0.2.6

1. Put project specific results in different directories:
    - ```project/test```
    - ```project/main```
    - ```project/denmark```
    - ```project/megapestis```
1. Added 2020 Latvia samples to database, mark as low coverage for now.
1. Move log directory to within results.
1. Add jupyter notebooks to snakemake pipeline.
1. Remove xml tree output from notebooks.
1. Create scripts to backup, restore, and clean projects.
1. Remove plot_table rules and script.
1. Only include visual files into the report.
1. Put config files inside the associated project directory.
1. Update function ```identify_local_sample``` to use database.
1. Environment addition: ugur, cartopy, ffpmeg, snp-dists, bokeh.
1. Create Pairwise SNP matrix.
1. Remove results folder from docs.
1. Create new output folder collect for detect snp density.
1. Move project specific results to a new repository.
1. Remove US Kim strain, mark as laboratory manipulation.
1. Rerun projects
