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
        aln            = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/snippy-multi.snps.aln",
    output:
        nwk            = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.treefile",
        nex            = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.nex",
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
        {scripts_dir}/newick2nexus.py {output.nwk} {output.nex}
        """

rule iqtree_filter:
    """
    Filter IQTREE output to remove the output and sync alignment and metadata.
    """
    input:
        aln    = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/snippy-multi.snps.aln",
        tsv    = results_dir + "/metadata/{reads_origin}/metadata.tsv",
        tree   = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.treefile",
        taxa   = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.filter-taxa.txt",
    output:
        taxa_aln    = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/filter-taxa/snippy-multi.snps.aln",
        taxa_tsv    = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/filter-taxa/metadata.tsv",
        taxa_tree   = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/filter-taxa/iqtree.treefile",
        sites_aln   = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/filter-sites/snippy-multi.snps.aln",
        sites_log   = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/filter-sites/snippy-multi.snps.log",
    params:
        taxa_outdir = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/filter-taxa/",
        sites_outdir = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/filter-sites/",
        keep_singleton = config["snippy_keep_singleton"],
    shell:
        """
        python3 {scripts_dir}/filter_taxa.py \
            --tree {input.tree} \
            --aln {input.aln} \
            --outdir {params.taxa_outdir} \
            --metadata {input.tsv} \
            --prune-tips {input.taxa};

        python3 {scripts_dir}/filter_sites.py \
            --fasta {output.taxa_aln} \
            --missing 100 \
            {params.keep_singleton} \
            --output {output.sites_aln} \
            --log {output.sites_log};
        """

rule lsd:
    """
    Estimate a time-scaled phylogeny using LSD2 in IQTREE.
    """
    input:
        tsv            = results_dir + "/metadata/{reads_origin}/metadata.tsv",
        tree           = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/iqtree.treefile",
        aln            = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/snippy-multi.snps.aln",
        constant_sites = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/full/snippy-multi.constant_sites.txt",
    output:
        timetree       = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.nex",
        dates          = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.dates.txt",
        outgroups      = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.outgroups.txt",
        log            = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.log",
        taxa           = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.filter-taxa.txt",
    params:
        seed           = config["iqtree_seed"],
        prefix         = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd",
        outgroup       = config["iqtree_outgroup"],
    shell:
        """
        cut -f 1,4 {input.tsv} | grep -v "NA" | wc -l | cut -d " " -f 1 > {output.dates};
        cut -f 1,4 {input.tsv} | grep -v "NA" | tail -n+2 | sed 's/\[/b(/g'  | sed 's/\]/)/g' | sed 's/:/,/g' >> {output.dates};
        echo -e "Reference\t"{config[reference_date_bp]} >> {output.dates};
        outgroups=`echo {params.outgroup} | tr ',' '\n'`;
        echo -e "$outgroups" | wc -l > {output.outgroups};
        echo -e "$outgroups" >> {output.outgroups};
				constant_sites=`awk -F "," '{{print ($1 + $2 + $3 + $4)}}' {input.constant_sites}`;

        lsd2 \
            -i {input.tree} \
            -s ${{constant_sites}} \
            -o {params.prefix} \
            -f 100 \
			-l '-1' \
			-q 0.2 \
            -r k \
            -d {output.dates} \
            -g {output.outgroups} \
            -G > {output.log};

        mv {params.prefix}.nexus {params.prefix}.divtree.nex
        mv {params.prefix}.nwk {params.prefix}.divtree.nwk
        mv {params.prefix}.date.nexus {params.prefix}.nex
        if  [[ `grep -A 1 "outliers" {output.log}` ]]; then
            grep -A 1 "outliers" {output.log} | tail -n 1 | tr " " "\n" | tail -n+2 > {output.taxa};
        else
            touch {output.taxa}
        fi
        """

rule beast:
    """
    Prepare input files for beast1 and beast2
    """
    input:
        tsv      = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/filter-taxa/metadata.tsv",
        divtree  = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/filter-taxa/iqtree.treefile",
        timetree = results_dir + "/lsd/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/lsd.nex",
        aln      = results_dir + "/iqtree/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/filter-sites/snippy-multi.snps.aln",
        constant_sites = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/full/snippy-multi.constant_sites.txt",
    output:
        latlon   = results_dir + "/beast/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/beast.latlon.txt",
        timetree = results_dir + "/beast/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/beast.timetree.nwk",
        divtree  = results_dir + "/beast/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/beast.divtree.nwk",
        dates    = results_dir + "/beast/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/beast.dates.txt",
        aln      = results_dir + "/beast/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/beast.fasta",
        constant_sites = results_dir + "/beast/{reads_origin}/{locus_name}/{prune}/filter{missing_data}/beast.constant-sites.txt",

    shell:
        """
        tail -n+2 {input.tsv} | cut -f 1,19,20 > {output.dates};
        cut -f 1,21,22 {input.tsv} > {output.latlon};
        cp {input.aln} {output.aln};
        cp {input.divtree} {output.divtree};
        cp {input.timetree} {output.timetree};
        cp {input.constant_sites} {output.constant_sites};
        """
