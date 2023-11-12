# Monopogen SNV Calling Commands

## Preparation Steps (only needs to be run once)

1. `get_all_imputation.sh`: downloads imputation panel from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/

## Analysis (run for each dataset)

### Preprocessing Step 

1. `submit_index`: uses samtools to index scATAC.bam file
2. `submit_preprocess`: runs Monopogen preprocessing steps

### Calling Step

1. `submit_calling`
