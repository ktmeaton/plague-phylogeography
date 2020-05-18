# Treetime Notes

Note: I edited treetime's treeanc.py file to force verbosity to the highest level (6).

## Parameters of Interest

```bash
--clock-filter int
--tip-slack int
--covariation
--gtr
--gtr-params
--branch-length-mode auto,input,joint,marginal
--confidence
--resolve-polytomies
--relax
```

## Default

```bash
treetime \
  --tree  ../iqtree/iqtree.core-filter0_bootstrap.treefile \
  --dates metadata_timetree_edit.tsv \
  --aln ../snippy_multi/snippy-core.full_CHROM.fasta \
  --reroot GCA_000016445.1_ASM1644v1_genomic \
  --outdir default1/
```

## Augur vs. TreeTime

Under some parameters, TreeTime will default to excluding samples that fall samples that fall 3 deviations/intervals away from the root-to-tip regression line. Samples with HIGH substitions are projected too far in the future (ex. in the year 3000/4000).

Augur refine disables this parameter by default, allowing samples with overly high/low substitions to be part of the optimization process.

### Augur Refine Mimic

```bash
treetime \
  --tree  ../iqtree/iqtree.core-filter0_bootstrap.treefile \
  --dates metadata_timetree_edit.tsv \
  --aln ../snippy_multi/snippy-core.full_CHROM.fasta \
  --reroot GCA_000016445.1_ASM1644v1_genomic \
  --coalescent opt \
  --covariation \
  --branch-length-mode joint \
  --gtr JC69 \
  --outdir augur_mimic_1
```

```bash
tt.run(infer_gtr=infer_gtr, root=reroot, Tc=Tc, time_marginal=marginal,
           branch_length_mode=branch_length_inference, resolve_polytomies=resolve_polytomies,
           max_iter=max_iter, fixed_pi=fixed_pi, fixed_clock_rate=clock_rate,
           vary_rate=vary_rate, use_covariation=covariance, **kwarks)

time_marginal is a boolean

treetime \
  --tree  ../iqtree/iqtree.core-filter0_bootstrap.treefile \
  --dates metadata_timetree_edit.tsv \
  --aln ../snippy_multi/snippy-core.full_CHROM.fasta \
  --reroot GCA_000016445.1_ASM1644v1_genomic \
  --gtr infer \
  --coalescent opt \
  --branch-length-mode auto \
  --max-iter 2 \
  --covariation \
  --outdir augur_mimic_2
```

self <treetime.treetime.TreeTime object at 0x7fcc22693e90>
root ['GCA_000016445.1_ASM1644v1_genomic']
infer_gtr True
relaxed_clock None
n_iqd 3
resolve_polytomies True
max_iter 2
Tc opt
fixed_clock_rate None
time_marginal False
sequence_marginal False
branch_length_mode auto
vary_rate False
use_covariation True
kwargs {'reconstruct_tip_states': False, 'fixed_pi': None, 'n_points': 20}

Disable the clock filter

```bash
treetime \
  --tree  ../iqtree/iqtree.core-filter0_bootstrap.treefile \
  --dates metadata_timetree_edit.tsv \
  --aln ../snippy_multi/snippy-core.full_CHROM.fasta \
  --reroot GCA_000016445.1_ASM1644v1_genomic \
  --gtr infer \
  --coalescent opt \
  --branch-length-mode auto \
  --max-iter 2 \
  --covariation \
  --clock-filter 0 \
  --outdir augur_mimic_3
```

## Augur Refine Actual

```bash
augur refine \
  --tree  ../iqtree/iqtree.core-filter0_bootstrap.treefile \
  --alignment ../snippy_multi/snippy-core.full_CHROM.fasta \
  --metadata metadata_timetree_edit.tsv \
  --timetree \
  --root GCA_000016445.1_ASM1644v1_genomic \
  --coalescent opt \
  --output-tree augur_1/tree.nwk \
  --output-node-data augur_1/branch_lengths.json;
```

self <treetime.treetime.TreeTime object at 0x7f29914d8d10>
root GCA_000016445.1_ASM1644v1_genomic
infer_gtr True
relaxed_clock None
n_iqd None
resolve_polytomies True
max_iter 2
Tc opt
fixed_clock_rate None
time_marginal False
sequence_marginal False
branch_length_mode auto
vary_rate False
use_covariation True
kwargs {'fixed_pi': None}

## Relaxed Clock

--relax
Strength of the gaussian priors on branch specific rate deviation and the coupling of parent and offspring rates can be specified.
e.g. as –relax 1.0 0.5.
Values around 1.0 correspond to weak priors.
Larger values constrain rate deviations more strongly.
Coupling 0 (–relax 1.0 0) corresponds to an un-correlated clock.

```bash
treetime \
  --tree  ../iqtree/iqtree.core-filter0_bootstrap.treefile \
  --dates metadata_timetree_edit.tsv \
  --aln ../snippy_multi/snippy-core.full_CHROM.fasta \
  --reroot GCA_000016445.1_ASM1644v1_genomic \
  --relax 1.0 0 \
  --outdir relax_uncorrelated/
```

## Model

--gtr custom
--gtr-params kappa=0.02 pi=0.263,0.237,0.238,0.262 mu=1.00000,1.91199,0.42141,0.42141,1.91199,1.00000
