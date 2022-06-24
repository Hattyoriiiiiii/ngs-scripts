#!/bin/bash
# usage : sh scripts/run.sh > log/run.log 2>&1 &

#----------------------------------------------------------------------------------------------------
############## config ##############
#----------------------------------------------------------------------------------------------------

# common
exp="fastp"
SAMPLES="others/sample_tf.txt"
cores=32

# fastqc : 010
INPUT='data/*.fastq.gz'
OUTDIR_FASTQC_BEFORE='01_qc/before'
REPODIR_FASTQC_BEFORE='01_qc/before/multiqc_before'

if [ ! -d ${OUTDIR_FASTQC_BEFORE} ]; then mkdir -p ${OUTDIR_FASTQC_BEFORE}; fi
if [ ! -d ${REPODIR_FASTQC_BEFORE} ]; then mkdir -p ${REPODIR_FASTQC_BEFORE}; fi

# trimming : 020
OUTDIR_TRIM='02_trim/fastp'
REPODIR_TRIM='01_qc/fastp'

if [ ! -d ${OUTDIR_TRIM} ]; then mkdir ${OUTDIR_TRIM}; fi
if [ ! -d ${REPODIR_TRIM} ]; then mkdir ${REPODIR_TRIM}; fi


OUTDIR_FASTQC_AFTER="01_qc/after"
REPODIR_FASTQC_AFTER="01_qc/after/multiqc_after"

if [ ! -d ${OUTDIR_FASTQC_AFTER} ]; then mkdir -p ${OUTDIR_FASTQC_AFTER}; fi
if [ ! -d ${REPODIR_FASTQC_AFTER} ]; then mkdir -p ${REPODIR_FASTQC_AFTER}; fi

# mapping : 030
INDIR_ALIGN='02_trim/fastp'
OUTDIR_ALIGN='03_align/fastp'
REPODIR_ALIGN='summary/bowtie2/fastp'

if [ ! -d ${OUTDIR_ALIGN} ]; then mkdir -p ${OUTDIR_ALIGN}; fi
if [ ! -d ${REPODIR_ALIGN} ]; then mkdir -p ${REPODIR_ALIGN}; fi

# filtering : 040
INDIR='03_align/fastp'
OUTDIR_FILTER='04_filtering/fastp'
REPODIR_FILTER='summary/picard/fastp'

if [ ! -d ${OUTDIR_FILTER} ]; then mkdir -p ${OUTDIR_FILTER}; fi
if [ ! -d ${REPODIR_FILTER} ]; then mkdir -p ${REPODIR_FILTER}; fi

# # peak calling
# INDIR_PEAKCALL="04_filtering/${exp}"
# OUTDIR_PEAKCALL="05_peaks/fastp/macs2"
# REPODIR_PEAKCALL="summary/macs2/${exp}"
# pval=1e-3
# pval_name="p0001"

# if [ ! -d ${OUTDIR_PEAKCALL} ]; then mkdir -p ${OUTDIR_PEAKCALL}; fi
# if [ ! -d ${REPODIR_PEAKCALL} ]; then mkdir -p ${REPODIR_PEAKCALL}; fi

#----------------------------------------------------------------------------------------------------
############## QC before trim ##############
#----------------------------------------------------------------------------------------------------

# 010: fastqc before
while read sample
do
    fastqc \
        -t ${cores} \
        --nogroup \
        -f fastq \
        ${sample} \
        -o ${OUTDIR_FASTQC_BEFORE} \
		data/${sample}_1.fastq.gz \
		data/${sample}_2.fastq.gz
done < ${SAMPLES}

multiqc ${OUTDIR_FASTQC_BEFORE} -o ${REPODIR_FASTQC_BEFORE}
mv ${OUTDIR_FASTQC_BEFORE}/*.html ${REPODIR_FASTQC_BEFORE}/

#----------------------------------------------------------------------------------------------------
############## trimming ##############
#----------------------------------------------------------------------------------------------------

while read sample
do
    echo ${sample}_1 ${sample}_2
    fastp \
        -i data/${sample}_1.fastq.gz \
        -I data/${sample}_2.fastq.gz \
        -o ${OUTDIR_TRIM}/${sample}_1.fastp.fastq.gz \
        -O ${OUTDIR_TRIM}/${sample}_2.fastp.fastq.gz \
        -h ${REPODIR_TRIM}/${sample}.html \
        -j ${REPODIR_TRIM}/${sample}_fastp.json \
        --thread 32 \
        -q 30 \
        --detect_adapter_for_pe
done < ${SAMPLES}

multiqc ${REPODIR_TRIM} -o ${REPODIR_TRIM}

#----------------------------------------------------------------------------------------------------
############## QC after trim ##############
#----------------------------------------------------------------------------------------------------

while read sample
do
    fastqc \
	    -t ${cores} \
        --nogroup \
        -f fastq \
        ${sample} \
        -o ${OUTDIR_FASTQC_AFTER} \
		${OUTDIR_TRIM}/${sample}_1.fastp.fastq.gz \
		${OUTDIR_TRIM}/${sample}_2.fastp.fastq.gz
done < ${SAMPLES}

multiqc ${OUTDIR_FASTQC_AFTER} -o ${REPODIR_FASTQC_AFTER}
mv ${OUTDIR_FASTQC_AFTER}/*.html ${REPODIR_FASTQC_AFTER}/

#----------------------------------------------------------------------------------------------------
############## mapping ##############
#----------------------------------------------------------------------------------------------------

while read sample; do
    bowtie2 \
        -p ${cores} \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        --phred33 \
        -x ~/indexes/mm10 \
        -1 ${INDIR_ALIGN}/${sample}_1.fastp.fastq.gz \
        -2 ${INDIR_ALIGN}/${sample}_2.fastp.fastq.gz \
        2> ${REPODIR_ALIGN}/${sample}_bowtie2.txt | \

    samtools view \
        -bhS \
        -@ ${cores} \
        -F 0x4 | \

    samtools sort \
        -@ ${cores} - > ${OUTDIR_ALIGN}/${sample}.sorted.bam

    samtools index ${OUTDIR_ALIGN}/${sample}.sorted.bam
done < ${SAMPLES}

#----------------------------------------------------------------------------------------------------
############## filtering ##############
#----------------------------------------------------------------------------------------------------

while read sample
do
	samtools view \
		-bh \
		-@ ${cores} \
		-f 0x2 \
		${OUTDIR_ALIGN}/${sample}.sorted.bam \
		> ${OUTDIR_FILTER}/${sample}.sorted.concordant.bam

	samtools index ${OUTDIR_FILTER}/${sample}.sorted.concordant.bam

	# mark duplicates
	picard MarkDuplicates \
		I=${OUTDIR_FILTER}/${sample}.sorted.concordant.bam \
		O=${OUTDIR_FILTER}/${sample}.sorted.concordant.dupMark.bam \
		M=${REPODIR_FILTER}/${sample}_picard.dupMark.txt

	# remove duplicates
	picard MarkDuplicates \
		I=${OUTDIR_FILTER}/${sample}.sorted.concordant.bam \
		O=${OUTDIR_FILTER}/${sample}.sorted.concordant.rmdup.bam \
		M=${REPODIR_FILTER}/${sample}_picard.rmdup.txt \
		REMOVE_DUPLICATES=true
	
	samtools view \
		-@ ${cores} \
		-h ${OUTDIR_FILTER}/${sample}.sorted.concordant.rmdup.bam | grep -v chrM | \
	
	samtools sort \
		-@ ${cores} \
		-O bam \
		-o ${OUTDIR_FILTER}/${sample}.sorted.concordant.rmdup.rmChrM.bam

	samtools index ${OUTDIR_FILTER}/${sample}.sorted.concordant.rmdup.rmChrM.bam

	# mkdir -p bigWig
	# bamCoverage \
	# 	-b ${OUTDIR_FILTER}/${sample}.sorted.uniq.rmdup.rmChrM.bam \
	# 	-o bigWig/${sample}.sorted.uniq.rmdup.rmChrM.bw \
	# 	--binSize 20 \
	# 	--normalizeUsing BPM \
	# 	--smoothLength 60 \
	# 	--extendReads 150 \
	# 	--centerReads \
	# 	-p 16

	# bamCoverage \
	# 	--bam ${OUTDIR_FILTER}/${sample}.sorted.uniq.rmdup.rmChrM.bam \
	# 	-o ${OUTDIR_FILTER}/${sample}.sorted.uniq.rmdup.rmChrM.bw \
	# 	--numberOfProcessors $cores \
	# 	--binSize 10 \
	# 	--normalizeUsing CPM \
	# 	--minMappingQuality 10 \
	# 	--extendReads
done < ${SAMPLES}


# #----------------------------------------------------------------------------------------------------
# ########## Peak Call ##########
# #----------------------------------------------------------------------------------------------------

# printf "\n <<< ########## Peak Calling ########## >>> \n\n"

# macs2 callpeak \
# 	-t ${INDIR_PEAKCALL}/Sox30.sorted.uniq.rmdup.rmChrM.bam \
# 	-c ${INDIR_PEAKCALL}/Input.sorted.uniq.rmdup.rmChrM.bam \
# 	--outdir ${OUTDIR_PEAKCALL} \
# 	-f BAMPE \
# 	-g 'mm' \
# 	-B \
# 	-n Sox30 \
# 	-p ${pval} \
# 	--cutoff-analysis
	
# 	# calculate the number of peaks
# 	wc -l ${OUTDIR_PEAKCALL}/Sox30_peaks.narrowPeak > QualityAssessment_${pval_name}/Sox30_num_peaks.txt &


# # sort peaks by -log10(p-value)
# # sort -k8,8nr ${OUTDIR_PEAKCALL}/${sample}_peaks.narrowPeak > ${OUTDIR_PEAKCALL}/${sample}_peaks.narrowPeak
