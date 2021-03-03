# Development

1. Added 2020 Latvia samples to database.
1. Added several results directories to version control.
    - snippy_multi
    - detect*
    - iqtree
    - multiqc
1. Decided not to upload qualimap, too many small files!
1. Put the version controlled results in different directories:
    - ```project/main```
    - ```project/denmark```
    - ```project/megapestis```
1. Move log directory to within results.
1. Add parse_tree notebook to snakemake pipeline.
1. Add clock_model notebook to snakemake pipeline.
1. Rename 'clock_model' output directory to 'clock'
1. Remove xml tree output.
1. Create scripts to clean and backup results.
1. Remove plot_table rules and script.
1. Remove almost all attached files from report.

## TO-DO

1. Check if multiqc still runs without eager output dir...
1. Maybe write a function to condense exporting.
