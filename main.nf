#!/usr/bin/env nextflow

// BED AND TARGETS //
targets = params.targets
targets_zip = params.targets_zip
targets_brca12 = params.targets_brca12
target_intervals = params.target_intervals
target_intervals_big = params.target_intervals_big
bed_melt = params.bed_melt

// FASTA //
genome_file = params.genome_file

// VEP FASTA AND ANNOTATIONS DBS //
cadd = params.cadd
vep_fasta = params.vep_fasta
vep_cache = params.vep_cache
gnomad = params.gnomad

// ADAPTORS FOR TRIMMOMATIC //
brca_adapt = params.brca_adapt

// ANNOTATION DBS //
known1_indels = params.known1
known2_indels = params.known2
dbsnp = params.dbsnp
cosmic = params.cosmic
swegen = params.swegen

// PIPELINE CONFIGS //
mei_list = params.mei_list
outdir = params.outdir

////////////////////////////////////


csv = file(params.csv)
mode = csv.countLines() > 2 ? "paired" : "single"
println(mode)

// Split bed file in to smaller parts to be used for parallel variant calling
Channel
	.fromPath("${params.targets}")
	.ifEmpty { exit 1, "Regions bed file not found: ${params.targets}" }
	.splitText( by: 170, file: 'bedpart.bed' )
	.set { beds_freebayes; }

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple( row.group, row.type, row.id, row.read1, row.read2 ) }
	.set { fastq }

// Trimming of reads, umis to header and remove 2 overhang + minlen 30
process trimmomatic {
	cpus 10

	input:
		set group, val(type), val(id), r1, r2 from fastq

	output:
		set group, val(type), val(id), file("${type}_R1_a_q_u_trimmed.fq.gz"),file("${type}_R2_a_q_u_trimmed.fq.gz") into fastq_trimmed

	script:

		if( params.trimmer == "trimmo" ) {
			"""
			trimmomatic -Xmx12g PE -phred33 -threads ${task.cpus} \\
				$r1 $r2 \\
				${type}_R1_a_trimmed.fq.gz /dev/null \\
				${type}_R2_a_trimmed.fq.gz /dev/null \\
				ILLUMINACLIP:$brca_adapt:3:12:7:1:true MINLEN:30
			trimmomatic -Xmx12G PE -phred33 -threads ${task.cpus} \\
				${type}_R1_a_trimmed.fq.gz ${type}_R2_a_trimmed.fq.gz \\
				${type}_R1_a_q_trimmed.fq.gz /dev/null \\
				${type}_R2_a_q_trimmed.fq.gz /dev/null \\
				MAXINFO:30:0.25 MINLEN:30
			trimmomatic -Xmx12G PE -phred33 -threads ${task.cpus} \\
				${type}_R1_a_q_trimmed.fq.gz ${type}_R2_a_q_trimmed.fq.gz \\
				${type}_R1_a_q_u_trimmed.fq.gz /dev/null \\
				${type}_R2_a_q_u_trimmed.fq.gz /dev/null \\
				HEADCROP:5
			"""
		}
		else if ( params.trimmer == "fastp" ) {
			"""
			fastp -i $r1 -I $r2 --stdout \\
				-U --umi_loc=per_read --umi_len=3 \\
				-w ${task.cpus} \\
			| fastp --stdin --interleaved_in -f 2 -F 2 \\
				-o ${type}_R1_a_q_u_trimmed.fq.gz \\
				-O ${type}_R2_a_q_u_trimmed.fq.gz \\
				-l 30 \\
				-w ${task.cpus}
			"""
		}
}

// aligning with sentieon bwa
process sentieon_bwa {
	cpus 40

	input:
		set group, val(type), id, file(r1), file(r2) from fastq_trimmed

	output:
		set group, type, id, file("${type}_${id}.bwa.sort.bam"), file("${type}_${id}.bwa.sort.bam.bai") into bams, qc_bam

	"""
	sentieon bwa mem -M \\
		-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
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

// locus collector + dedup of bam
process locus_collector_dedup {
	cpus 16

	input:
		set group, val(type), val(id), file(bam), file(bai) from bams

	output:
		set group, val(type), val(id), file("${type}_${id}.bwa.sort.dedup.bam"), file("${type}_${id}.bwa.sort.dedup.bam.bai") into dedup_bam
		set group, val(type), val(id), file(bam), file(bai), file("dedup_metrics.txt") into bam_dedup_stats

	"""
	sentieon driver -t ${task.cpus} -i $bam \\
		--algo LocusCollector --fun score_info ${type}_score.gz
	sentieon driver -t ${task.cpus} -i $bam \\
		--algo Dedup --rmdup --score_info ${type}_score.gz \\
		--metrics dedup_metrics.txt ${type}_${id}.bwa.sort.dedup.bam
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
	sentieon driver \\
		-t ${task.cpus} \\
		-r $genome_file \\
		-i $bam --algo Realigner \\
		-k $known1_indels \\
		-k $known2_indels \\
		${type}_${id}.bwa.sort.dedup.realigned.bam   
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
	sentieon driver \\
		-t ${task.cpus} \\
		-r $genome_file \\
		-i $bam \\
		--algo QualCal -k $known1_indels -k $known2_indels ${type}_${id}.bqsr
	"""
}

// apply bqsr to bam and plot before and after
process recal_bam {
	cpus 16
	publishDir "${outdir}/bam/twist-brca/", mode: 'copy', overwrite: 'true'

	input:
		set group, val(type), val(id), file(bam), file(bai), file(score) from bqsr_table

	output:
		set group, val(type), val(id), file("${type}_${id}.bwa.sort.dedup.realigned.recal.bam"), file("${type}_${id}.bwa.sort.dedup.realigned.recal.bam.bai"), file("${type}_${id}.bqsr.pdf") into bam_done

	"""
	sentieon driver -t ${task.cpus} -r $genome_file -i $bam -q $score \\
		--algo QualCal -k $known1_indels -k $known2_indels ${type}_${id}.recal_data.post \\
		--algo ReadWriter ${type}_${id}.bwa.sort.dedup.realigned.recal.bam
	sentieon driver -t ${task.cpus} \\
		--algo QualCal --plot --before $score --after ${type}_${id}.recal_data.post RECAL_RESULT.CSV
	sentieon plot QualCal -o ${type}_${id}.bqsr.pdf RECAL_RESULT.CSV
	"""
}

//Collect various QC data MOVE QC_SENTIEON.PL TO CONTAINER TODO!
process sentieon_qc {
	cpus 40
	memory '64 GB'
	publishDir "${outdir}/postmap/twist-brca/", mode: 'copy', overwrite: 'true'

	input:
		set group, type, id, file(bam), file(bai), file(dedup) from bam_dedup_stats

	output:
		set group, type, id, file("${id}.QC") into qc_val

	"""
	sentieon driver \\
		-r $genome_file --interval $targets \\
		-t ${task.cpus} -i ${bam} \\
		--algo MeanQualityByCycle mq_metrics.txt \\
		--algo QualDistribution qd_metrics.txt \\
		--algo GCBias --summary gc_summary.txt gc_metrics.txt \\
		--algo AlignmentStat aln_metrics.txt \\
		--algo InsertSizeMetricAlgo is_metrics.txt \\
		--algo CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt 
	sentieon driver \\
		-r $genome_file \\
		-t ${task.cpus} -i ${bam} \\
		--algo HsMetricAlgo --targets_list $target_intervals_big --baits_list $target_intervals_big hs_metrics.txt
	qc_sentieon.pl $id panel > ${id}.QC
	"""

}


if(mode == "paired") {
	bam_done
		.groupTuple()
		.into{ bam_tnscope; bam_freebayes; bam_melt; bam_manta; }
}
else {
	bam_done
		.into{ bam_tnscope; bam_freebayes; bam_melt; bam_manta; }
}
// VARIANT CALLING
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
process tnscope {
	cpus 16

	input:
		set group, type, id, file(bam), file(bai), file(plot) from bam_tnscope

	output:
		set group, file("${group}_tnscope.vcf") into tnscope_vcf

	script:
		if(mode == "paired") { 
			tumor_index = bam.findIndexOf{ it ==~ /tumor_.+/ }
			tumor = bam[tumor_index]
			tumor_id = id[tumor_index]
			normal_index = bam.findIndexOf{ it ==~ /normal_.+/ }
			normal = bam[normal_index]
			normal_id = id[normal_index]

			"""
			sentieon driver \\
				-t ${task.cpus} \\
				-r $genome_file \\
				-i $tumor \\
				-i $normal \\
				--interval $targets \\
				--algo TNscope \\
				--tumor_sample ${tumor_id}_tumor \\
				--normal_sample ${normal_id}_normal \\
				--cosmic $cosmic \\
				--dbsnp $dbsnp \\
				--disable_detector sv \\
				${group}_tnscope.vcf
			#bedtools intersect -a ${group}_tnscope.vcf.raw -b $targets -header > ${group}_tnscope.vcf
			"""
		}
		else {
			"""
			sentieon driver \\
				-t ${task.cpus} \\
				-r $genome_file \\
				-i $bam \\
				--interval $targets \\
				--algo TNscope \\
				--tumor_sample $id \\
				--cosmic $cosmic \\
				--dbsnp $dbsnp \\
				--disable_detector sv \\
				${group}_tnscope.vcf
			#bedtools intersect -a ${group}_tnscope.vcf.raw -b $targets -header >  ${group}_tnscope.vcf
			"""
		}
}

process freebayes{
	cpus 1
	
	input:
		set group, type, id, file(bam), file(bai), file(plot) from bam_freebayes
		each file(bed) from beds_freebayes;
	output:
		set group, file("${group}_${bed}_freebayes.vcf") into freebayes_vcf

	script:
		if(mode == "paired") { 
			tumor_index = bam.findIndexOf{ it ==~ /tumor_.+/ }
			tumor = bam[tumor_index]
			tumor_id = id[tumor_index]
			normal_index = bam.findIndexOf{ it ==~ /normal_.+/ }
			normal = bam[normal_index]
			normal_id = id[normal_index]
			
			"""
			freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete -F 0.03 $tumor $normal > ${group}__${bed}_freebayes.vcf
			vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ${group}__${bed}_freebayes.vcf > ${group}_${bed}_freebayes.vcf
			"""
		}
		else {
			"""
			freebayes -f $genome_file -C 2 -F 0.01 --pooled-continuous --genotype-qualities -t $bed $bam > ${group}__${bed}_freebayes.vcf
			vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ${group}__${bed}_freebayes.vcf > ${group}_${bed}_freebayes.vcf
			"""
		}
}

process concatenate_vcfs {
	input:
		set group, file(vcfs) from freebayes_vcf.groupTuple()

	output:
		set group, file("${group}_freebayes_concat.vcf") into freebayes_vcf_concat

	"""
	vcf-concat $vcfs | vcf-sort -c > ${group}_freebayes_concat.vcf
	"""
}

// MANTA SINGLE AND PAIRED
process manta {
	cpus 16
	
	input:
		set group, type, id, file(bam), file(bai), file(plot) from bam_manta

	output:
		set group, file("${group}_manta.vcf") into manta_vcf

	when:
		params.manta
	
	script:
		if(mode == "paired") { 
			tumor_index = bam.findIndexOf{ it ==~ /tumor_.+/ }
			tumor = bam[tumor_index]
			tumor_id = id[tumor_index]
			normal_index = bam.findIndexOf{ it ==~ /normal_.+/ }
			normal = bam[normal_index]
			normal_id = id[normal_index]

			"""
			configManta.py \\
				--tumorBam $tumor \\
				--normalBam $normal \\
				--reference $genome_file \\
				--exome \\
				--callRegions $targets_zip \\
				--generateEvidenceBam \\
				--runDir .
			python runWorkflow.py -m local -j ${task.cpus}
			#filter_manta_paired.pl results/variants/somaticSV.vcf.gz > ${group}_manta.vcf
			mv results/variants/tumorSV.vcf.gz ${group}_manta.vcf.gz
			gunzip ${group}_manta.vcf.gz
			"""
		}
		else {
			"""
			configManta.py \\
				--tumorBam $bam \\
				--reference $genome_file \\
				--exome \\
				--callRegions $targets_zip \\
				--generateEvidenceBam \\
				--runDir .
			python runWorkflow.py -m local -j ${task.cpus}
			#filter_manta.pl results/variants/tumorSV.vcf.gz > ${group}_manta.vcf
			mv results/variants/tumorSV.vcf.gz ${group}_manta.vcf.gz
			gunzip ${group}_manta.vcf.gz
			"""
		}
}

qc_val
	.groupTuple()
	.set{ qc_tables }
process melt {
	cpus 16
	errorStrategy 'ignore'

	input:
		set group, type, id, file(bam), file(bai), file(plot) from bam_melt
		set group2, type2, id2, qc from qc_tables

	when:
		params.melt

	output:
		set group, file("${group}_melt.vcf") into melt_vcf

	script:
		if(mode == "paired") { 
			normal_index = bam.findIndexOf{ it ==~ /normal_.+/ }
			normal = bam[normal_index]
			normal_id = id[normal_index]
			qc_index = type2.findIndexOf{ it ==~ /normal/ }
			qc = qc[qc_index]
		}
		else {
			qc = qc[0]
		}
		// Collect qc-data if possible from normal sample, if only tumor; tumor
		qc.readLines().each{
			if (it =~ /\"(ins_size_dev)\" : \"(\S+)\"/) {
				ins_dev = it =~ /\"(ins_size_dev)\" : \"(\S+)\"/
			}
			if (it =~ /\"(mean_coverage)\" : \"(\S+)\"/) {
				coverage = it =~ /\"(mean_coverage)\" : \"(\S+)\"/
			}
			if (it =~ /\"(ins_size)\" : \"(\S+)\"/) {
				ins_size = it =~ /\"(ins_size)\" : \"(\S+)\"/
			}
		}
		// might need to be defined for -resume to work "def INS_SIZE" and so on....
		INS_SIZE = ins_size[0][2]
		MEAN_DEPTH = coverage[0][2]
		COV_DEV = ins_dev[0][2]

	"""
	java -jar  /opt/MELT.jar Single \\
		-bamfile $bam \\
		-r 150 \\
		-h $genome_file \\
		-n $bed_melt \\
		-z 50000 \\
		-d 50 -t $mei_list \\
		-w . \\
		-b 1/2/3/4/5/6/7/8/9/10/11/12/14/15/16/18/19/20/21/22 \\
		-c $MEAN_DEPTH \\
		-cov $COV_DEV \\
		-e $INS_SIZE
	mv ALU.final_comp.vcf ${group}_melt.vcf
	"""

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

tnscope_vcf
	.mix(freebayes_vcf_concat, manta_vcf, melt_vcf)
	.set{ vcfs }


// Post-processing of variant callers
process normalize {
	cpus 16

	input:
		set group, file(vcf) from vcfs

	output:
		set group, file("${vcf}.norm") into norm_vcf


	"""
	vcfbreakmulti $vcf > ${group}.multibreak
	bcftools norm -m-both -c w -O v -f $genome_file -o ${vcf}.norm ${group}.multibreak
	rm ${group}.multibreak
	"""
}

process annotate_vep {
	container = '/fs1/resources/containers/container_VEP.sif'
	cpus 40

	input:
		set group, file(vcf) from norm_vcf

	output:
		set group, file("${vcf}.vep") into vep

	"""
	vep \\
		-i ${vcf} \\
		-o ${vcf}.vep \\
		--offline \\
		--merged \\
		--everything \\
		--vcf \\
		--format vcf  \\
		--no_stats \\
		--fork ${task.cpus} \\
		--force_overwrite \\
		--plugin CADD,$cadd \\
		--fasta $vep_fasta \\
		--dir_cache $vep_cache \\
		--dir_plugins $vep_cache/Plugins \\
		--distance 200 \\
		-cache \\
		-custom $gnomad
	"""
}

process bgzip_index {
	cpus 16
	publishDir "${outdir}/vcf/twist-brca/", mode: 'copy', overwrite: 'true'

	input:
		set group, file(vcf) from vep

	output:
		set group, file("${vcf}.gz"), file("${vcf}.gz.tbi") into vcf_done

	"""
	bgzip -@ ${task.cpus} $vcf -f
	tabix ${vcf}.gz -f
	"""
}

process aggregate_vcf {
	publishDir "${outdir}/vcf/twist-brca/", mode: 'copy', overwrite: 'true'

    input:
    	set group, file(vcf), file(tbi) from vcf_done.groupTuple()

    output:
		file("${group}_aggregated.vcf")

    script:
		manta_index = vcf.findIndexOf{ it ==~ /.+_manta\..+/ }
		freebayes_index = vcf.findIndexOf{ it ==~ /.+_freebayes_.+/ }
		tnscope_index = vcf.findIndexOf{ it ==~ /.+_tnscope\..+/ }
		freebayes = vcf[freebayes_index]
		tnscope = vcf[tnscope_index]
		manta = vcf[manta_index]

    """
    aggregate_vcf.pl --freebayes $freebayes --tnscope $tnscope --manta $manta --base freebayes > ${group}_aggregated.vcf
    """
}


