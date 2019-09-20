#!/usr/bin/env nextflow



// global variables /////////////////
/////////////////////////////////////
outdir = params.outdir
targets = params.targets
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
brca_adapt = params.brca_adapt
known1_indels = params.known1
known2_indels = params.known2
dbsnp = params.dbsnp
////////////////////////////////////////

csv = file(params.csv)
mode = csv.countLines() > 2 ? "paired" : "single"
println(mode)

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple( row.group, row.type, row.id, row.read1, row.read2 ) }
    .set { fastq }

// trimming dual-duplex umis from both reads. 
process trimmomatic {
    cpus 56
    input:
    set group, val(type), val(id), r1, r2 from fastq.view()
    output:
    set group, val(type), val(id), file("${type}_R1_umitrimmed.fq.gz"),file("${type}_R2_umitrimmed.fq.gz") into fastq_trimmed
    """
    trimmomatic -Xmx12g PE -phred33 -threads ${task.cpus} \\
    $r1 $r2 \\
    ${type}_R1_umitrimmed.fq.gz /dev/null \\
    ${type}_R2_umitrimmed.fq.gz /dev/null \\
    HEADCROP:5 MINLEN:30
    """
}

// aligning with sentieon bwa
process sentieon_bwa {
    cpus 56
    input:
    set group, val(type), id, file(r1), file(r2) from fastq_trimmed
    output:
    set group, type, id, file("${type}_${id}.bwa.sort.bam"), file("${type}_${id}.bwa.sort.bam.bai") into bam
    """
    sentieon bwa mem -M \\
    -R '@RG\\tID:${id}_${type}\\tSM:${id}_${type}\\tPL:illumina' \\
    -t ${task.cpus} \\
    $genome_file \\
    ${r1} ${r2} | \\
    sentieon util sort \\
    -r $genome_file \\
    -o ${type}_${id}.bwa.sort.bam \\
    -t ${task.cpus} \\
    --sam2bam -i -
    """
}

bam
    .into{ bam_qc; bam_post }

// locus collector + dedup of bam
process locus_collector_dedup {
    cpus 16
    input:
    set group, val(type), val(id), file(bam), file(bai) from bam_post
    output:
    set group, val(type), val(id), file("${type}_${id}.bwa.sort.dedup.bam"), file("${type}_${id}.bwa.sort.dedup.bam.bai") into dedup_bam
    """
    sentieon driver -t ${task.cpus} -i $bam --algo LocusCollector --fun score_info ${type}_score.gz
    sentieon driver -t ${task.cpus} -i $bam --algo Dedup --rmdup --score_info ${type}_score.gz --metrics ${type}_metrics.txt ${type}_${id}.bwa.sort.dedup.bam
    """
}

// indel realignment
process indel_realign {
    cpus 16
    input:
    set group, val(type), val(id), file(bam), file(bai) from dedup_bam
    output:
    set group, val(type), val(id), file("${type}_${id}.bwa.sort.dedup.realigned.bam"), file("${type}_${id}.bwa.sort.dedup.realigned.bam.bai") into realigned_bam
    """
    sentieon driver -t ${task.cpus} -r $genome_file -i $bam --algo Realigner -k $known1_indels -k $known2_indels ${type}_${id}.bwa.sort.dedup.realigned.bam   
    """
}

// calculate base quality score recalibration
process bqsr {
    cpus 16
    input:
    set group, val(type), val(id), file(bam), file(bai) from realigned_bam
    output:
    set group, val(type), val(id), file(bam), file(bai), file("${type}_${id}.bqsr") into bqsr_table
    """
    sentieon driver -t ${task.cpus} -r $genome_file -i $bam --algo QualCal -k $known1_indels -k $known2_indels ${type}_${id}.bqsr
    """
}

// apply bqsr to bam and plot before and after
process recal_bam {
    cpus 16
    input:
    set group, val(type), val(id), file(bam), file(bai), file(score) from bqsr_table
    output:
    set group, val(type), val(id), file("${type}_${id}.bwa.sort.dedup.realigned.recal.bam"), file("${type}_${id}.bwa.sort.dedup.realigned.recal.bam.bai"), file("${type}_${id}.bqsr.pdf") into bam_done
    """
    sentieon driver -t ${task.cpus} -r $genome_file -i $bam -q $score \\
    --algo QualCal -k $known1_indels -k $known2_indels ${type}_${id}.recal_data.post --algo ReadWriter ${type}_${id}.bwa.sort.dedup.realigned.recal.bam
    sentieon driver -t ${task.cpus} --algo QualCal --plot \\
    --before $score --after ${type}_${id}.recal_data.post RECAL_RESULT.CSV
    sentieon plot QualCal -o ${type}_${id}.bqsr.pdf RECAL_RESULT.CSV
    """
}

// Collect various QC data
// process sentieon_qc {
//     cpus 56
//     memory '64 GB'
//     input:
// 	set type, id, file(bam), file(bai) from bam_qc
//     output:
//     file("*.txt")
//     """
//     sentieon driver -r $genome_file -t ${task.cpus} -i ${bam} --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt \
//     --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt \
//     --algo HsMetricAlgo --targets_list $targets --baits_list $targets hs_metrics.txt --algo CoverageMetrics cov_metrics.txt
//     """

// }


if(mode == "paired") {
    bam_done
        .groupTuple()
        .into{ bam_tnscope; bam_freebayes; bam_melt; bam_manta; }
}
else {
    bam_done
        .into{ bam_tnscope; bam_freebayes; bam_melt; bam_manta; }
}

process tnscope {
    cpus 16
    input:
    set group, type, id, file(bam), file(bai), file(plot) from bam_tnscope
    output:
    set group, file("${group}_tnscope.vcf"), file("${group}_tnscope.vcf.idx") into tnscope_vcf
    script:
    if(mode == "paired") { 
    tumor_index = bam.findIndexOf{ it ==~ /tumor_.+/ }
    tumor = bam[tumor_index]
    tumor_id = id[tumor_index]
    normal_index = bam.findIndexOf{ it ==~ /normal_.+/ }
    normal = bam[normal_index]
    normal_id = id[normal_index]
    """
    sentieon driver -t ${task.cpus} -r $genome_file -i $tumor -i $normal --algo TNscope \\
    --tumor_sample ${tumor_id}_tumor --normal_sample ${normal_id}_normal --dbsnp $dbsnp ${group}_tnscope.vcf
    """
    }
    else {
    """
    sentieon driver -t ${task.cpus} -r $genome_file -i $bam --algo TNscope --tumor_sample $id --dbsnp $dbsnp ${group}_tnscope.vcf
    """
    }
}

process freebayes{
    cpus 16
    input:
    set group, type, id, file(bam), file(bai), file(plot) from bam_freebayes.view()
    output:
    set group, file("${group}_freebayes.vcf") into freebayes_vcf
    script:
    if(mode == "paired") { 
    tumor_index = bam.findIndexOf{ it ==~ /tumor_.+/ }
    tumor = bam[tumor_index]
    tumor_id = id[tumor_index]
    normal_index = bam.findIndexOf{ it ==~ /normal_.+/ }
    normal = bam[normal_index]
    normal_id = id[normal_index]
    """
    echo "$tumor_id:$tumor   $normal_id:$normal" > ${group}_freebayes.vcf
    """
    }
    else {
    """
    ${group}_freebayes.vcf
    """
    }
}

// process manta {

// }

// process melt {

// }

// process VEP {

// }