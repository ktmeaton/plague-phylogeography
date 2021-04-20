import random # generate a random seed for iqtree

# -----------------------------------------------------------------------------#
#                         Maximum Likelihood Phylogeny                         #
# -----------------------------------------------------------------------------#
rule iqtree:
    """
    Construct a maximum likelihood phylogeny.
    """
    input:
        constant_sites = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/full/snippy-multi.constant_sites.txt",
        aln        = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/snippy-multi.snps.aln",
    output:
        tree           = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.nex",
        iqtree         = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.iqtree",
        log            = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.log",
        outgroup       = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.filter-taxa.txt",
    params:
        #seed = random.randint(0, 99999999),
        seed           = config["iqtree_seed"],
        prefix         = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree",
        outgroup       = config["iqtree_outgroup"],
        model          = config["iqtree_model"],
        runs           = config["iqtree_runs"],
        other          = config["iqtree_other"],
    resources:
        load           = 100,
        time_min       = 600,
      	cpus           = workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb         = workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
    shell:
        """
        echo {params.outgroup} | tr "," "\n" >> {output.outgroup};
        iqtree \
            -s {input.aln} \
		    {params.model} \
            --threads-max {resources.cpus} \
            -nt {resources.cpus} \
            -o {params.outgroup} \
            -seed {params.seed} \
            --runs {params.runs} \
            -fconst `cat {input.constant_sites}` \
            {params.other} \
            -redo \
            -pre {params.prefix} > {output.log};

        {scripts_dir}/newick2nexus.py {params.prefix}.treefile {output.tree}
        """

rule filter_taxa:
    """
    Remove taxa from an alignment based on a tree.
    """
    input:
        tree           = results_dir + "/{rule}/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/{rule}.nex",
        tsv            = results_dir + "/metadata/{reads_origin}/metadata.tsv",
        taxa           = results_dir + "/{rule}/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/{rule}.filter-taxa.txt",
        aln            = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/snippy-multi.snps.aln",
    output:
        nex            = results_dir + "/{rule}/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/{rule}.filter.nex",
        nwk            = results_dir + "/{rule}/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/{rule}.filter.nwk",
        tsv            = results_dir + "/{rule}/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/metadata.tsv",
        aln            =  results_dir + "/{rule}/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/{rule}.filter.aln",
    params:
        taxa           = config["iqtree_outgroup"],
        outdir         = results_dir + "/{rule}/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/",
    shell:
        """
        workflow/scripts/filter_alignment.py \
            --tree {input.tree} \
            --aln {input.aln} \
            --outdir {params.outdir} \
            --metadata {input.tsv} \
            --prune-tips {input.taxa}
        """


rule lsd:
    """
    Estimate a time-scaled phylogeny using LSD2 in IQTREE.
    """
    input:
        tsv            = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/metadata.tsv",
        tree           = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.filter.nwk",
        aln            = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.filter.aln",
        constant_sites = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/full/snippy-multi.constant_sites.txt",
    output:
        timetree       = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.nex",
        dates          = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.dates.txt",
        log            = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.timetree.log",
        taxa           = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.filter-taxa.txt",
    params:
        seed           = config["iqtree_seed"],
        prefix         = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd",
    shell:
        """
        cut -f 1,4 {input.tsv}  | tail -n+2 | sed 's/\[\|\]//g' > {output.dates};
        echo -e "Reference\t"{config[reference_date_bp]} >> {output.dates}
        iqtree \
            -s {input.aln} \
		    {config[iqtree_model]} \
            --date {output.dates} \
            --prefix {params.prefix} \
            -fconst `cat {input.constant_sites}` \
            --date-ci 100 \
            --date-outlier 3 \
            -redo \
            -te {input.tree} > {output.log};

        grep -A 1 "outliers" {output.log} | tail -n 1 | tr " " "\n" | tail -n+2 > {output.taxa};
        cp {params.prefix}.timetree.nex {output.timetree}
        """

rule beast_geo:
    """
    Continuous phylogeography with BEAST
    """
    input:
        tsv      = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/metadata.tsv",
        dates    = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.dates.txt",
        tree     = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.filter.nex",
    output:
        latlon   = results_dir + "/beast/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/beast.latlon.txt",
        tree     = results_dir + "/beast/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/beast.nex",
        dates    = results_dir + "/beast/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/beast.dates.txt",

    shell:
        """
        echo -e "traits\tlat\tlon" > {output.latlon};
        echo -e "Reference\t"{config[reference_lat]}"\t"{config[reference_lon]} >> {output.latlon};
        echo -e "Reference\t"{config[reference_date_bp]} > {output.dates}

        tail -n+2 {input.tsv} | while read line; \
        do
            sample=`echo "$line" | cut -f 1`;
            lat=`echo "$line" | cut -f 9`;
            lon=`echo "$line" | cut -f 10`;
            date=`echo "$line" | cut -f 4 | sed "s/\[\|\]//g" | tr ":" "\n" | awk '{{sum+=$0}}END{{print 0-sum/NR}}'`;
            if [[ $lat == "NA" ]]; then
                lat=`echo "$line" | cut -f 7`;
                lon=`echo "$line" | cut -f 8`;
            fi;
            echo -e $sample"\t"$lat"\t"$lon >> {output.latlon};
            echo -e $sample"\t"$date >> {output.dates};
        done

        {scripts_dir}/multi2bi.py --tree {input.tree} --out {output.tree}
        """
