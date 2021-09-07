#!/usr/bin/bash
# usage: sh scripts/010_fastqc_bf.sh > log/010_fastqc_bf.log 2>&1 &

cores=32
SAMPLES="others/sample.txt"
exp="fastp"

########## FastQC - Before trimming ##########
printf "\n <<< ########## FastQC - Before trimming ########## >>> \n\n"

INPUT_FASTQC_BEFORE="data/*.fastq.gz"
OUTDIR_FASTQC_BEFORE="01_qc/before"
REPODIR_FASTQC_BEFORE="01_qc/before/multiqc_before"

if [ ! -d ${OUTDIR_FASTQC_BEFORE} ]; then mkdir -p ${OUTDIR_FASTQC_BEFORE}; fi
if [ ! -d ${REPODIR_FASTQC_BEFORE} ]; then mkdir -p ${REPODIR_FASTQC_BEFORE}; fi

# 010: fastqc before
for sample in `ls ${INPUT_FASTQC_BEFORE}`; do
    fastqc \
        -t ${cores} \
        --nogroup \
        -f fastq \
        ${sample} \
        -o ${OUTDIR_FASTQC_BEFORE}
done

multiqc ${OUTDIR_FASTQC_BEFORE} -o ${REPODIR_FASTQC_BEFORE}
mv ${OUTDIR_FASTQC_BEFORE}/*.html ${REPODIR_FASTQC_BEFORE}/


########## Read Trimming ##########
printf "\n <<< ########## Read Trimming ########## >>> \n\n"

OUTDIR_TRIM="02_trim/${exp}"
REPODIR_TRIM="01_qc/${exp}"

if [ ! -d ${OUTDIR_TRIM} ]; then mkdir ${OUTDIR_TRIM}; fi
if [ ! -d ${REPODIR_TRIM} ]; then mkdir ${REPODIR_TRIM}; fi


while read sample
do
    echo ${sample}_1 ${sample}_2
    fastp \
        -i data/${sample}.fastq.gz \
        -o ${OUTDIR_TRIM}/${sample}.${exp}.fastq.gz \
        -h ${REPODIR_TRIM}/${sample}.html \
        -j ${REPODIR_TRIM}/${sample}_${exp}.json \
        --thread 16 \
        --length_required 36 \
        -q 30
done < ${SAMPLES}

multiqc ${REPODIR_TRIM} -o ${REPODIR_TRIM}


########## FastQC - After trimming ##########
printf "\n <<< ########## FastQC - After trimming ########## >>> \n\n"

INPUT_FASTQC_AFTER="02_trim/${exp}/*.fastq.gz"
OUTDIR_FASTQC_AFTER="01_qc/after"
REPODIR_FASTQC_AFTER="01_qc/after/multiqc_after"

if [ ! -d ${OUTDIR_FASTQC_AFTER} ]; then mkdir -p ${OUTDIR_FASTQC_AFTER}; fi
if [ ! -d ${REPODIR_FASTQC_AFTER} ]; then mkdir -p ${REPODIR_FASTQC_AFTER}; fi

for sample in `ls ${INPUT_FASTQC_AFTER}`
do
    fastqc \
	    -t ${cores} \
        --nogroup \
        -f fastq \
        ${INPUT_FASTQC_AFTER} \
        -o ${OUTDIR_FASTQC_AFTER}
done

multiqc ${OUTDIR_FASTQC_AFTER} -o ${REPODIR_FASTQC_AFTER}
mv ${OUTDIR_FASTQC_AFTER}/*.html ${REPODIR_FASTQC_AFTER}/


########## Mapping ##########
printf "\n <<< ########## Mapping ########## >>> \n\n"

# INDIR_ALIGN <- OUTDIR_TRIM
OUTDIR_ALIGN="03_align/${exp}"
REPODIR_ALIGN="summary/bowtie2"

if [ ! -d ${OUTDIR_ALIGN} ]; then mkdir -p ${OUTDIR_ALIGN}; fi
if [ ! -d ${REPODIR_ALIGN} ]; then mkdir -p ${REPODIR_ALIGN}; fi

while read sample; do
	bowtie2 \
		--end-to-end \
		-p ${cores} \
		--very-sensitive \
		--phred33 \
		-X 1000 \
		-x ~/indexes/mm10 \
		-U ${OUTDIR_TRIM}/${sample}.${exp}.fastq.gz \
		2> ${REPODIR_ALIGN}/${sample}_bowtie2.txt | \

	samtools view \
		-bhS \
		-@ ${cores} \
		-F 0x4 -q 42 | \

	samtools sort \
		-@ ${cores} - > ${OUTDIR_ALIGN}/${sample}.sorted.uniq.bam

	samtools index ${OUTDIR_ALIGN}/${sample}.sorted.uniq.bam
done < ${SAMPLES}


########## Filtering ##########
printf "\n <<< ########## Filtering ########## >>> \n\n"

# INDIR_FILTER=OUTDIR_ALIGN
OUTDIR_FILTER="04_filtering/${exp}"
REPODIR_FILTER="summary/picard/${exp}"

if [ ! -d ${OUTDIR_FILTER} ]; then mkdir -p ${OUTDIR_FILTER}; fi
if [ ! -d ${REPODIR_FILTER} ]; then mkdir -p ${REPODIR_FILTER}; fi

while read sample
do

	# mark duplicates
	picard MarkDuplicates \
		I=${OUTDIR_ALIGN}/${sample}.sorted.uniq.bam \
		O=${OUTDIR_FILTER}/${sample}.sorted.uniq.dupMark.bam \
		M=${REPODIR_FILTER}/${sample}_picard.dupMark.txt

	# remove duplicates
	picard MarkDuplicates \
		I=${OUTDIR_ALIGN}/${sample}.sorted.uniq.bam \
		O=${OUTDIR_FILTER}/${sample}.sorted.uniq.rmdup.bam \
		M=${REPODIR_FILTER}/${sample}_picard.rmdup.txt \
		REMOVE_DUPLICATES=true
	
	samtools view \
		-@ ${cores} \
		-h ${OUTDIR_FILTER}/${sample}.sorted.uniq.rmdup.bam | grep -v "chrM" | \

	samtools view \
		-@ ${cores} \
		-hb \
		> ${OUTDIR_FILTER}/${sample}.sorted.uniq.rmdup.rmChrM.bam

	# To use downstream analysis
	samtools index ${OUTDIR_FILTER}/${sample}.sorted.uniq.rmdup.rmChrM.bam

	picard CollectInsertSizeMetrics \
		INPUT=${OUTDIR_FILTER}/${sample}.sorted.uniq.rmdup.bam \
		OUTPUT=${REPODIR_FILTER}/${sample}_insert_size_metrics.txt \
		HISTOGRAM_FILE=${REPODIR_FILTER}/${sample}_hist.pdf \
		MINIMUM_PCT=0
	
	bamCoverage \
		--bam ${OUTDIR_FILTER}/${sample}.sorted.uniq.rmdup.rmChrM.bam \
		-o ${OUTDIR_FILTER}/${sample}.sorted.uniq.rmdup.rmChrM.bw \
		--numberOfProcessors $cores \
		--binSize 10 \
		--normalizeUsing CPM \
		--minMappingQuality 10 \
		--extendReads

done < ${SAMPLES}
