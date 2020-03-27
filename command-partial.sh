nextflow run pipeline.nf \
			  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
			  --max_datasets 10 \
			  -with-trace \
			  -with-timeline \
			  -with-dag pipeline.pdf \
			  -with-report \
			  -name PlaguePhyloPartial2 \
				-resume
