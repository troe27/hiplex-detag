#!/bin/bash
## print rulegraph #
snakemake  -s ./genotype_pinus.snek --configfile config/config.yaml  --rerun-incomplete --use-conda --rulegraph | dot | display

#snakemake  -s ./genotype_pinus.snek --configfile config/config.yaml  --rerun-incomplete --use-conda
