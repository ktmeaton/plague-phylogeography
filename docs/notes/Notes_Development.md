# Development

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
1. Add augur and cartopy to environment.

## TO-DO

1. Maybe write a function to condense notebook import/exporting.
1. Rethink the use of large dill objects that can't be version controlled...

## GFF INDEX

```yaml
col1: id (col 9, parse)
col2: loci (col 1)
col3: type (col 3)
col4: start (col 4)
col5: end (col 5)
```

```bash
cd results/data/reference/GCA_000009065.1_ASM906v1_genomic/
echo -e "id\tloci\ttype\tstart\tend" > GCA_000009065.1_ASM906v1_genomic.gff.index
grep -v "#" GCA_000009065.1_ASM906v1_genomic.gff | while read line;
do
    id=`echo -e "$line" | cut -f 9 | cut -d ";" -f 1 | sed 's/ID=//g'`;
    echo -e "$id\t$line" | cut -f 1,2,4,5,6 >> GCA_000009065.1_ASM906v1_genomic.gff.index;
done
```
