#!/usr/bin/bash
# usage : scripts/030_mapping_pe.sh > log/030_mapping_pe.log 2>&1 &

############## config ##############

cores=32
exp="fastp"

# fastqc : 010
INPUT="data/*.fastq.gz"
OUTDIR_FASTQC_BEFORE="01_qc/before"
REPODIR_FASTQC_BEFORE="01_qc/before/multiqc_before"

if [ ! -d ${OUTDIR_FASTQC_BEFORE} ]; then mkdir -p ${OUTDIR_FASTQC_BEFORE}; fi
if [ ! -d ${REPODIR_FASTQC_BEFORE} ]; then mkdir -p ${REPODIR_FASTQC_BEFORE}; fi

# trimming : 020
OUTDIR_TRIM="02_trim/${exp}"
REPODIR_TRIM="01_qc/${exp}""

if [ ! -d ${OUTDIR_TRIM} ]; then mkdir ${OUTDIR_TRIM}; fi
if [ ! -d ${REPODIR_TRIM} ]; then mkdir ${REPODIR_TRIM}; fi

# mapping : 030
INDIR_ALIGN="02_trim/${exp}"
OUTDIR_ALIGN="03_align/${exp}"
REPODIR_ALIGN='summary/bowtie2/${exp}"

if [ ! -d ${OUTDIR_ALIGN} ]; then mkdir -p ${OUTDIR_ALIGN}; fi
if [ ! -d ${REPODIR_ALIGN} ]; then mkdir -p ${REPODIR_ALIGN}; fi

# filtering : 040
INDIR="03_align/${exp}"
OUTDIR="04_filtering/${exp}"
REPODIR="summary/picard/${exp}"

if [ ! -d ${OUTDIR_FILTER} ]; then mkdir -p ${OUTDIR_FILTER}; fi
if [ ! -d ${REPODIR_FILTER} ]; then mkdir -p ${REPODIR_FILTER}; fi

# common
SAMPLES='others/sample.txt'



############## QC before trim ##############

# 010: fastqc before
for sample in `ls ${INPUT}`; do
    fastqc \
        -t ${cores} \
        --nogroup \
        -f fastq \
        ${sample} \
        -o ${OUTDIR_FASTQC_BEFORE}
done

multiqc ${OUTDIR_FASTQC_BEFORE} -o ${REPORT_FASTQC_BEFORE}
mv ${OUTDIR_FASTQC_BEFORE}/*.html ${REPODIR_FASTQC_BEFORE}/

############## trimming ##############

while read sample
do
    echo ${sample}_1 ${sample}_2
    fastp \
        -i data/${sample}_1.fastq.gz \
        -I data/${sample}_2.fastq.gz \
        -o ${OUTDIR_TRIM}/${sample}_1.${exp}.fastq.gz \
        -O ${OUTDIR_TRIM}/${sample}_2.${exp}.fastq.gz \
        -h ${REPODIR_TRIM}/${sample}.html \
        -j ${REPODIR_TRIM}/${sample}_${exp}.json \
        -w 16 \
        -3 \
        -q 30 \
        --detect_adapter_for_pe
done < ${SAMPLES}

multiqc ${REPODIR_TRIM} -o ${REPODIR_TRIM}

############## QC after trim ##############

INPUT="02_trim/${exp}/*.fastq.gz"
OUTDIR="01_qc/after"
REPODIR="01_qc/after/multiqc_after"

for sample in `ls ${INPUT}`
do
    fastqc \
	    -t ${cores} \
        --nogroup \
        -f fastq \
        ${INPUT} \
        -o ${OUTDIR}
done

multiqc ${OUTDIR} -o ${REPODIR}
mv ${OUTDIR}/*.html ${REPODIR}/

############## mapping ##############

while read sample; do
bowtie2 \
	-p ${cores} \
	--very-sensitive \
	--no-mixed \
	--no-discordant \
	--phred33 \
	-I 10 -X 700 \
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


############## filtering ##############

while read sample
do
	samtools view \
		-bh \
		-@ ${cores} \
		-f 0x2 \
		${INDIR}/${sample}.sorted.bam \
		> ${OUTDIR}/${sample}.sorted.concordant.bam

	samtools index ${OUTDIR}/${sample}.sorted.concordant.bam

	# mark duplicates
	picard MarkDuplicates \
		I=${OUTDIR}/${sample}.sorted.concordant.bam \
		O=${OUTDIR}/${sample}.sorted.concordant.dupMark.bam \
		M=${REPODIR}/${sample}_picard.dupMark.txt

	# remove duplicates
	picard MarkDuplicates \
		I=${OUTDIR}/${sample}.sorted.concordant.bam \
		O=${OUTDIR}/${sample}.sorted.concordant.rmdup.bam \
		M=${REPODIR}/${sample}_picard.rmdup.txt \
		REMOVE_DUPLICATES=true
	
	samtools view \
		-@ ${cores} \
		-h ${OUTDIR}/${sample}.sorted.concordant.rmdup.bam | grep -v chrM | \
	
	samtools sort \
		-@ ${cores} \
		-O bam \
		-o ${OUTDIR}/${sample}.sorted.concordant.rmdup.rmChrM.bam

	samtools index ${OUTDIR}/${sample}.sorted.concordant.rmdup.rmChrM.bam

	picard CollectInsertSizeMetrics \
		INPUT=${OUTDIR}/${sample}.sorted.concordant.rmdup.bam \
		OUTPUT=${REPODIR}/${sample}_insert_size_metrics.txt \
		HISTOGRAM_FILE=${REPODIR}/picard/fastp/${sample}_hist.pdf \
		MINIMUM_PCT=0
done < ${SAMPLES}
