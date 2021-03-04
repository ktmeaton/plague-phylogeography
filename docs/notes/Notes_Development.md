# Development

1. Added 2020 Latvia samples to database.
1. Put project specific results in different directories:
    - ```project/main```
    - ```project/denmark```
    - ```project/megapestis```
1. Move log directory to within results.
1. Add parse_tree notebook to snakemake pipeline.
1. Add clock_model notebook to snakemake pipeline.
1. Rename 'clock_model' output directory to 'clock'
1. Remove xml tree output from notebooks.
1. Create scripts to clean and backup results.
1. Remove plot_table rules and script.
1. Remove almost all attached files from report.
1. Put config files inside the associated project directory.
1. Update function ```identify_local_sample``` to use database.

## TO-DO

1. Maybe write a function to condense notebook import/exporting.
