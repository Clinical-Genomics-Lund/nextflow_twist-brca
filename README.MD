## Nextflow pipeline for TWIST-Ärftlig Cancer

Still under development!

## Deploy
bash deploy_hopper.sh (to /fs1/pipelines/twist-brca/)

To run nextflow with singularity on hopper you need first module load.

`module load Java singularity nextflow`

the run command should look like this.

`nextflow run main.nf --csv input.csv -with-singularity /fs1/resources/containers/container_twist-brca.sif -resume --reads 1,1,0 --manta --melt`


## CSV-file
needs to have the following headers
`clarity_sample_id,id,assay,sex,diagnosis,phenotype,group,clarity_pool_id,platform,read1,read2`

example-files can be located /fs1/pipelines/twist-brca/test_files/ any csv that should work.


## Known Issues

Melt will fail if it cannot find any mobile elements at all. As such, errors from this process is ignored and no vcf will be produced.
