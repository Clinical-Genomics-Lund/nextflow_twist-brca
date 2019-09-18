#!/usr/bin/env nextflow



// global variables /////////////////
/////////////////////////////////////
outdir = params.outdir
genome_file = params.genome_file
cadd = params.cadd
vep_fast = params.vep_fast
maxentscan = params.maxentscan
vep_cache = params.vep_cache
gnomad = params.gnomad
gerp = params.gerp
phylop =  params.phylop
phastcons = params.phastcons
snpsift = params.snpsift
clinvar = params.clinvar
swegen = params.swegen
spidex = params.spidex
sentieon_model = params.sentieon_model
////////////////////////////////////////

//csv = file(params.csv)
//mode = csv.countLines() > 2 ? "paired" : "single"
//println(mode)

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple( row.type ) }
    .set { type_check }

type_check
    .collect()
    .set {type_checker}
process tumor_normal_check {

    input:
    val(type) from type_checker.view()
    output:
    set mode, file("hej") into mode
    script:
    
    mode = type.size() > 1 ? "paired" : "single"
    println(mode)
    if(mode == "single") {
    mode = type == "tumor" ? "tumor" : "blood"
    }
    """
    echo $mode > hej
    """


}

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple( row.type, row.id, row.sex, row.read1, row.read2) }
    .set { fastq }


// process sentieon_bwa {
//     cpus 56
//     input:
//     set val(type), val(id), val(sex), r1, r2 from fastq
//     output:
//     set type, id, file("${type}_${id}_bwa.sort.bam"), file("${type}_${id}_bwa.sort.bam.bai") into bam
//     """
//     sentieon-bwa mem -M \\
//     -R '@RG\\tID:${id}_${type}\\tSM:${id}_${type}\\tPL:illumina' \\
//     -t ${task.cpus} \\
//     $genome_file \\
//     ${r1} ${r2} | \\
//     sentieon util sort \\
//     -r $genome_file \\
//     -o ${type}_${id}_bwa.sort.bam \\
//     -t ${task.cpus} \\
//     --sam2bam -i -
//     """
// }

