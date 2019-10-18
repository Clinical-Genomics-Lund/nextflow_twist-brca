#!/usr/bin/env nextflow



// global variables /////////////////
/////////////////////////////////////
outdir = params.outdir
targets = params.targets
targets_zip = params.targets_zip
targets_brca12 = params.targets_brca12
target_intervals = params.target_intervals
bed_melt = params.bed_melt
mei_list = params.mei_list
genome_file = params.genome_file
genome_file_fai = params.genome_file_fai
cadd = params.cadd
vep_fasta = params.vep_fasta
vep_cache = params.vep_cache
gnomad = params.gnomad
swegen = params.swegen
sentieon_model = params.sentieon_model
brca_adapt = params.brca_adapt
known1_indels = params.known1
known2_indels = params.known2
dbsnp = params.dbsnp
cosmic = params.cosmic
////////////////////////////////////////

csv = file(params.csv)
mode = csv.countLines() > 2 ? "paired" : "single"
println(mode)

// Split bed file in to smaller parts to be used for parallel variant calling
Channel
	.fromPath("${params.targets}")
	.ifEmpty { exit 1, "Regions bed file not found: ${params.targets}" }
	.splitText( by: 500, file: 'bedpart.bed' )
	.set { beds_freebayes; }

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple( row.group, row.type, row.id, row.read1, row.read2 ) }
	.set { fastq }

process fq2sam {
	cpus 40

	input:
		set group, val(type), val(id), r1, r2 from fastq

	output:
		set group, val(type), val(id), file("${id}_${type}.unaligned") into unaligned1

	"""
	picard -Xmx60g -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=. FastqToSam \\
		O=${id}_${type}.unaligned \\
		F1=$r1 \\
		F2=$r2 \\
		SM=$id \\
		LB=$id \\
		PU=CMD \\
		PL=illumina 
	"""
}

process ExtractUmisFromBam {
	cpus 40

	input:
		set group, val(type), val(id), file(sam) from unaligned1

	output:
		set group, val(type), val(id), file("${sam}.umi") into umi_unaligned

	"""
	fgbio --tmp-dir=. ExtractUmisFromBam \\
		--input=$sam \\
		--output=${sam}.umi \\
		--read-structure=3M2S146T 3M2S146T \\
		--molecular-index-tags=ZA ZB \\
		--single-tag=RX
	"""

}

process MarkIlluminaAdapters {
	cpus 40

	input:
		set group, val(type), val(id), file(sam) from umi_unaligned

	output:
		set group, val(type), val(id), file("${sam}.marked") into markedadapt

	"""
	picard -Xmx60g -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=. MarkIlluminaAdapters \\
		I=$sam \\
		O=${sam}.marked M=adapterMetrics.txt
	"""
}

process SamToFastq {
	cpus 40

	input:
		set group, val(type), val(id), file(sam) from markedadapt

	output:
		set group, val(type), val(id), file("${sam}.aligned") into fastq2

	"""
	picard -Xmx60g -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=. SamToFastq \\
		I=$sam \\
		CLIPPING_ATTRIBUTE=XT \\
		CLIPPING_ACTION=X \\
		CLIPPING_MIN_LENGTH=36 \\
		NON_PF=true \\
		F=/dev/stdout \\
		INTERLEAVE=true \\
	| sentieon bwa mem -M -t ${task.cpus} $genome_file -p /dev/stdin \\
	| picard -Xmx60g -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=. MergeBamAlignment \\
		UNMAPPED=$sam \\
		ALIGNED=/dev/stdin \\
		O=${sam}.aligned \\
		R=$genome_file \\
		CLIP_ADAPTERS=false \\
		VALIDATION_STRINGENCY=SILENT \\
		CREATE_INDEX=true \\
		EXPECTED_ORIENTATIONS=FR \\
		MAX_GAPS=-1 \\
		SO=coordinate \\
		ALIGNER_PROPER_PAIR_FLAGS=false
	"""
}


process GroupReadsByUmi {
	cpus 40

	input:
		set group, val(type), val(id), file(bam) from fastq2

	output:
		set group, val(type), val(id), file("${bam}.grouped") into grouped_umi

	"""
	fgbio --tmp-dir=. -Xmx60g GroupReadsByUmi \\
		--strategy=paired \\
		--input=$bam \\
		--output=${bam}.grouped \\
		--raw-tag=RX \\
		--assign-tag=MI \\
		--min-map-q=10 \\
		--edits=1
	"""
}
grouped_umi
	.into { grouped_umi1; grouped_umi2 }
process CallDuplexConsensusReads {
	cpus 20
	memory '120 GB'

	input:
		set group, val(type), val(id), file(bam) from grouped_umi1

	output:
		set group, val(type), val(id), file("${bam}.consensus") into consensus

	"""
	fgbio --tmp-dir=/tmp/ -Xmx60g CallDuplexConsensusReads \\
		--input=$bam \\
		--output=${bam}.consensus \\
		--threads=${task.cpus} \\
		--min-reads=3 2 1
	"""
}

process CollectDuplexSeqMetrics {
	cpus 40

	input:
		set group, val(type), val(id), file(bam) from grouped_umi2

	output:
		set group, val(type), val(id), file("${id}_metrics*") into seqmetrics

	"""
	fgbio --tmp-dir=. -Xmx60g CollectDuplexSeqMetrics \\
		--input=$bam \\
		--output=${id}_metrics
	"""
}


process SamToFastq2 {
	cpus 40

	input:
		set group, val(type), val(id), file(bam) from consensus

	output:
		set group, val(type), val(id), file("pool1_consensus.aligned.bam"), file("pool1_consensus.aligned.bai") into bams, qc_bam

	"""
	picard -Xmx60g -XX:ParallelGCThreads=56 -Djava.io.tmpdir=. SamToFastq \\
		VALIDATION_STRINGENCY=SILENT \\
		I=$bam \\
		F=/dev/stdout \\
		INTERLEAVE=true \\
		INCLUDE_NON_PF_READS=true \\
		CREATE_INDEX=true \\
		| sentieon bwa mem -K 1000000 -p -t 56 $genome_file /dev/stdin > tmp.aligned.bam
	picard -Xmx60g -XX:ParallelGCThreads=56 -Djava.io.tmpdir=. MergeBamAlignment \\
		VALIDATION_STRINGENCY=SILENT \\
		UNMAPPED=$bam \\
		ALIGNED=tmp.aligned.bam \\
		OUTPUT=pool1_consensus.aligned.bam \\
		R=$genome_file \\
		CLIP_ADAPTERS=false \\
		CREATE_INDEX=true \\
		ORIENTATIONS=FR \\
		MAX_GAPS=-1 \\
		SO=coordinate \\
		ALIGNER_PROPER_PAIR_FLAGS=false \\
		ATTRIBUTES_TO_RETAIN=X0 \\
		ATTRIBUTES_TO_RETAIN=ZS \\
		ATTRIBUTES_TO_RETAIN=ZI \\
		ATTRIBUTES_TO_RETAIN=ZM \\
		ATTRIBUTES_TO_RETAIN=ZC \\
		ATTRIBUTES_TO_RETAIN=ZN \\
		ATTRIBUTES_TO_REVERSE=ad \\
		ATTRIBUTES_TO_REVERSE=bd \\
		ATTRIBUTES_TO_REVERSE=cd \\
		ATTRIBUTES_TO_REVERSE=ae \\
		ATTRIBUTES_TO_REVERSE=be \\
		ATTRIBUTES_TO_REVERSE=ce
	"""
}


// // locus collector + dedup of bam
// process locus_collector_dedup {
// 	cpus 16

// 	input:
// 		set group, val(type), val(id), file(bam), file(bai) from bams

// 	output:
// 		set group, val(type), val(id), file("${type}_${id}.bwa.sort.dedup.bam"), file("${type}_${id}.bwa.sort.dedup.bam.bai") into dedup_bam
// 		set group, val(type), val(id), file(bam), file(bai), file("dedup_metrics.txt") into bam_dedup_stats

// 	"""
// 	sentieon driver -t ${task.cpus} -i $bam \\
// 		--algo LocusCollector --fun score_info ${type}_score.gz
// 	sentieon driver -t ${task.cpus} -i $bam \\
// 		--algo Dedup --rmdup --score_info ${type}_score.gz \\
// 		--metrics dedup_metrics.txt ${type}_${id}.bwa.sort.dedup.bam
// 	"""
// }

// indel realignment
process indel_realign {
	cpus 16

	input:
		set group, val(type), val(id), file(bam), file(bai) from bams

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
		set group, type, id, file(bam), file(bai) from qc_bam

	output:
		set group, type, id, file("${id}.QC") into qc_val

	"""
	sentieon driver \\
		-r $genome_file --interval $target_intervals \\
		-t ${task.cpus} -i ${bam} \\
		--algo MeanQualityByCycle mq_metrics.txt \\
		--algo QualDistribution qd_metrics.txt \\
		--algo GCBias --summary gc_summary.txt gc_metrics.txt \\
		--algo AlignmentStat aln_metrics.txt \\
		--algo InsertSizeMetricAlgo is_metrics.txt \\
		--algo HsMetricAlgo --targets_list $target_intervals --baits_list $target_intervals hs_metrics.txt \\
		--algo CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt 
	/fs1/pipelines/wgs_germline/annotation/qc_sentieon.pl $id panel > ${id}.QC
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
				--algo TNscope \\
				--tumor_sample ${tumor_id}_tumor \\
				--normal_sample ${normal_id}_normal \\
				--cosmic $cosmic \\
				--dbsnp $dbsnp \\
				${group}_tnscope.vcf
			#/opt/bin/filter_mutect.pl ${group}_tnscope.vcf.raw > ${group}_tnscope.vcf
			"""
		}
		else {
			"""
			sentieon driver \\
				-t ${task.cpus} \\
				-r $genome_file \\
				-i $bam \\
				--algo TNscope \\
				--tumor_sample $id \\
				--cosmic $cosmic \\
				--dbsnp $dbsnp \\
				${group}_tnscope.vcf
			#/opt/bin/filter_mutect_single.pl ${group}_tnscope.vcf.raw > ${group}_tnscope.vcf
			"""
		}
}

process freebayes{
	cpus 16
	
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
			freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete -F 0.03 $tumor $normal > ${group}_${bed}_freebayes.vcf
			#vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ${group}__freebayes.vcf > ${group}_${bed}_freebayes.vcf
			"""
		}
		else {
			"""
			freebayes -f $genome_file -C 2 -F 0.01 --pooled-continuous --genotype-qualities -t $bed $bam > ${group}_${bed}_freebayes.vcf
			#vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ${group}__freebayes.vcf > ${group}_${bed}_freebayes.vcf
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
			#/opt/bin/filter_manta.pl results/variants/somaticSV.vcf.gz > ${group}_manta.vcf
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
			#/opt/bin/filter_manta.pl results/variants/tumorSV.vcf.gz > ${group}_manta.vcf
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

// process aggregate_vcf {
//     input:
//     set group, file(vcf), file(tbi) from vcf_done.groupTuple()
//     output:
//     script:
//     manta_index = vcf.findIndexOf{ it ==~ /manta/ }
//     freebayes_index = vcf.findIndexOf{ it ==~ /freebayes/ }
//     tnscope_index = vcf.findIndexOf{ it ==~ /tnscope/ }
//     freebayes = vcf[freebayes_index]
//     mutect = vcf[tnscope_index]
//     manta = vcf[manta_index]
//     """
//     /opt/bin/aggregate_vcf.pl --freebayes $freebayes --mutect $mutect --manta $manta --base freebayes > ${group}_aggregated.vcf
//     """
// }


