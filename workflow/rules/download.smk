# -----------------------------------------------------------------------------#
#                                Data Download                                 #
# -----------------------------------------------------------------------------#
rule download_fna:
  """
  Download fasta files, by seaching for sample name matches.
  """
    input:
        "results/sqlite_import/{download_dir}.txt"
    output:
        "results/{download_dir}/{sample}.fna"
    run:
        for file in input:
            with open(file) as temp_file:
                file_parse = [line.rstrip() for line in temp_file if wildcards.sample in line]
            if len(file_parse):
                match = file_parse[0]
        shell("wget --quiet -O - {match} | gunzip -c > {output}")
