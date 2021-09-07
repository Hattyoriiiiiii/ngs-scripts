#!/usr/bin/bash
# usage : sh scripts/050_peakcalling.sh > log/050_peakcalling.log 2>&1 &

cores=32
SAMPLES="others/sample.txt"
exp="fastp"

########## Peak Call ##########

INDIR_PEAKCALL="04_filtering/${exp}"
OUTDIR_PEAKCALL="05_peaks/macs2"
REPODIR_PEAKCALL="summary/macs2"

if [ ! -d ${OUTDIR_PEAKCALL} ]; then mkdir -p ${OUTDIR_PEAKCALL}; fi
if [ ! -d ${REPODIR_PEAKCALL} ]; then mkdir -p ${REPODIR_PEAKCALL}; fi

while read sample
do
    macs2 callpeak \
        -t ${INDIR_PEAKCALL}/${sample}.sorted.uniq.rmdup.bam \
        --outdir ${OUTDIR_PEAKCALL}\
        -f BAM \
        --nomodel \
        --shift -100 \
        --extsize 200 \
        -g 'mm' \
        -B \
        -n ${sample} \
        -q 0.01
	
	# calculate the number of peaks
	wc -l ${OUTDIR_PEAKCALL}/${sample}_peaks.narrowPeak > QualityAssessment/${sample}_num_peaks.txt &

done < ${SAMPLES}


########## IDR ##########

OUTDIR_IDR="QualityAssessment/idr"

if [ ! -d ${OUTDIR_IDR} ]; then 
mkdir -p ${OUTDIR_IDR}
fi

for sample in `cat others/sample.txt | sed "s/_rep[12]//g" | uniq`
do

        # sort peaks by -log10(p-value)
        sort -k8,8nr ${OUTDIR_PEAKCALL}/${sample}_rep1_peaks.narrowPeak > ${OUTDIR_PEAKCALL}/${sample}_rep1_peaks.sorted.narrowPeak &
        sort -k8,8nr ${OUTDIR_PEAKCALL}/${sample}_rep2_peaks.narrowPeak > ${OUTDIR_PEAKCALL}/${sample}_rep2_peaks.sorted.narrowPeak && \

        idr --samples ${OUTDIR_PEAKCALL}/${sample}_rep1_peaks.sorted.narrowPeak ${OUTDIR_PEAKCALL}/${sample}_rep2_peaks.sorted.narrowPeak \
                --input-file-type narrowPeak \
                --rank p.value \
                --output-file ${OUTDIR_IDR}/${sample}-idr \
                --plot \
                --log-output-file log/${sample}.idr.log &
done