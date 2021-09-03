# ChIP-seq (single-end paired-end) の解析のワークフロー

---

- [Read mapping](#Read mapping)

---

## ChIP-seq analysis

### 考慮しておきたいこと
- 抗体の特異性
- S/N比

- read数

CUT&Tag, CUT&RUN, ATAC
- fragment size distribution

### Read mapping
- Bowtie2, BWA
- Duplicateのフィルタリングを行ったあとの残ったnonredundant readsを解析に用いる。

### Filtering
- retrotransposons, pseudogenes, paralogous genesなどのrepetitive genomic regionsに関心があるときはmuti-mapped readsを許容する

#### Blacklisted regions
にmappingされたfragmentsは下流の解析から除く

### Peak calling
- relaxed thresholdで多くのpeakをcallする
    - この時、TPだけではなく多くのノイズを含んでいることが予想される
- IDRを用いてsensitivityを上げる
    - TN増やしてFPを減らす

### Data Quality Assesment
- **Mapping ratio**
    - read qualityと目的のgenomic DNAからのreadがシーケンスされているかを反映する
    - HiSeq2500などでは、80%を越えるべきである
    - IgGなどnon-DNA-binding proteinでは60%を下回ることがある
- **Read depth**
    - nonredundant mapped readsの数
    - 十分な量のread depthはゲノムサイズや用いた抗体のS/N比に依存する
    - ENCODEでは、ヒトでsharp-mode peaksサンプルを解析するには少なくとも10 Mのuniquely mapped readsが必要とされている
    - S/N比が弱いbroad histone marksではpeak callingにさらにリードを要する
        - ヒトで40 M以上
- **Library complexity**
- **The normalized strand coefficent**

There are four steps in ChIP quality control:

1. Sample correlation clustering: Clustering of the pair-wise correlations between genome-wide signal profiles.
2. Data visualization in a genomic browser.
3. Average fragment length determination: Determining whether the ChIP was enriched for fragments of a certain length.
4. Visualization of GC bias. Here we will plot the ChIP enrichment versus the average GC content in the corresponding genomic bin.

---

![](https://i.imgur.com/r7Mx3kG.png)

- QC metrics
    - [ ] NSC > 1.1
    - [ ] RSC > 1
    - [ ] SSD
    - [ ] FRiP > 5%
    - [ ] REGI
    - [ ] RiBL

https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/lectures/ChIP-seq_workflow_scope.pdf

---


https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html

Irreproducible Discovery Rate (IDR) - 1つの実験における2つの生物学的複製間の一貫性を測定することで、ハイスループット実験の再現性を評価します。ChIP-seqやATAC-seqの評価に用いられる。

複製されたピークセットを操作して、個々の複製/擬似複製ピークセットにおけるこれらのピークのランクの一貫性を比較する統計的手順です。高いランクの一貫性を持つピークは保持されます。IDRは、「保守的」な出力ピークセットとなる一組の真の複製に含まれるピークに対して、または「最適」な出力ピークセットとなる一組の擬似複製に含まれるピークに対して操作することができる。保守的なピークセットのピークは、真の生物学的複製の間で再現可能な事象を表し、真の生物学的および技術的ノイズを考慮した、信頼性の高いピークと解釈することができます。最適なピークセットのピークは、再現性のあるイベントを表し、読み取りのサンプリングノイズを考慮した、信頼度の高いピークと解釈することができます。最適セットは、特にレプリケートの一方が他方よりもデータ品質が低い場合に、より感度が高くなります。
自己矛盾比率は、1つのデータセット内での一貫性を測定します。
レスキュー比は、1つの実験内の複製が比較できない場合に、データセット間の一貫性を測定します。

 IDRピークセットは信頼性の高いpeaksであるってわけだ。
 sensitivityが高い、つまりは陽性と判断したピークは実際に陽性である可能性が高い。特に片方のreplicateの品質がもう片方よりも低い時に。



https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/12_annotation_functional_analysis.html



- inputは1つ以上のpeak callセット。


1. Typically, you will want to annotate and characterize your most confident set of peak calls for each sample group. Thus, rather than annotating individual replicates, you will want to use your IDR peakset so you have only one peakset per sample group.
	- それぞれのreplicateをannotationするのではなくて、IDR peaksetを使用。だからサンプルグループごとに1つのpeaksetが存在する。
2. The peakcalls we have in our working directory are generated for a subset of the full dataset. The resulting target gene set that we get back from annotating will not be comprehensive enough to give us meaningful results in the downstream functional analysis.

These were obtained post-IDR analysis, (i.e. concordant peaks between replicates) and are provided in BED format which is optimal input for the ChIPseeker package.

---

### Reference
Nakato, Sakata, 2021
compgenomer
https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/lectures/ChIP-seq_workflow_scope.pdf

#### コード付き
https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html

#### IDR
https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html

#### Rのコード
https://hbctraining.github.io/Intro-to-ChIPseq/lessons/06_combine_chipQC_and_metrics.html

#### 日本語
http://catway.jp/bioinformatics/ChIP-Seq/motifdist.html
- 2つのTF ChIP-seq の共通結合領域を探す
- 2つのTF ChIP-seq の共通結合領域をカウントしベン図を描く
- ChIPpeakAnnoが多い

#### 理論的な話 (galaxy)
https://galaxyproject.org/tutorials/chip/


<br>

---

<br>

## テスト



### test data
- single-end vs paired-end
    - sharp peaks vs broad peaks
    
- single-end
    - Stra8, Meiosin


### ディレクトリ構成

```shell
ChIPseq_pe/
├── 01_qc
│   ├── after
│   │   └── multiqc_after
│   └── before
│       └── multiqc_before
├── 02_trim
├── 03_align
├── 04_filtering
├── 05_peaks
├── REPORTS_ChIP_PE_TEST
│   ├── log -> ../log
│   ├── multiqc_after -> ../01_qc/after/multiqc_after
│   ├── multiqc_before -> ../01_qc/before/multiqc_before
│   ├── others -> ../others
│   ├── scripts -> ../scripts
│   └── summary -> ../summary
├── data
├── log
├── others
├── scripts
└── summary
```

- `others/sample.txt`に入れる内容

```sh
data/RIPA_GFP-ChIP.fastq.gz
data/RIPA_Meiosin_C-ChIP.fastq.gz
data/RIPA_Meiosin_N-ChIP.fastq.gz
data/RIPA_input.fastq.gz

```

```txt
HSF5_PS_34_rep1
HSF5_PS_34_rep2
HSF5_PS_42_rep1
HSF5_PS_42_rep2
HSF5_RS_34_rep2
HSF5_RS_42_rep1
HSF5_RS_42_rep2
HSF5_w_34_cell
HSF5_w_34_rep1
HSF5_w_34_rep2
HSF5_w_42_cell
HSF5_w_42_rep1
HSF5_w_42_rep2

```

## Dockerfile


```
mkdir -p 01_qc/{before/multiqc_before,after/multiqc_after} 02_trim 03_align 04_filtering 05_peaks

mkdir REPORTS_ChIP_PE_TEST
cd REPORTS_ChIP_PE_TEST/
ln -s ../log log
ln -s ../01_qc/before/multiqc_before multiqc_before
ln -s ../01_qc/after/multiqc_after multiqc_after
ln -s ../others others
ln -s ../scripts scripts
ln -s ../summary summary
```

```
docker build -t hattyoriiiiiii/bioinformatics-bulk:version1.0 .
```


### `calculate_FRiPscore.sh`

```
#!/usr/bin/bash

########## config ##########

INDIR="04_filtering/fastp"
SAMPLES="others/sample.txt"


########## FRiP score ##########

for sample in `ls 05_peaks/*_peaks.narrowPeak | xargs basename -s _peaks.narrowPeak | cut -f 1,2 -d '_' | uniq`

while read sample; do
    echo ${sample}
    
    total_reads=$(samtools view -c ${INDIR}/${sample}.sorted.mapped.rmdup.bam)

    reads_in_peaks=$(bedtools sort -i 05_peaks/${sample}_peaks.narrowPeak | \
        bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
        -a ${INDIR}/${sample}.sorted.mapped.rmdup.bam -b stdin -ubam | samtools view -c)

    FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
    echo ${FRiP}
done < ${SAMPLES}

```


## Downstream Analysis

### Coverage

```
bamCoverage \
    --bam HSF5_RS_42_rep1.sorted.concordant.rmdup.bam \
    -o test01_HSF5_RS_42_rep1.sorted.concordant.rmdup.bw \
    --numberOfProcessors $cores \
    --binSize 10 \
    --normalizeUsing CPM \
    --effectiveGenomeSize $effect_genome_size \
    --minMappingQuality 10 \
    --extendReads
```

### motif探索

```
mkdir 06_motif
```


### Summary

`Rscript scripts/log_summary.R && sh scripts/weblogfmt.sh summary/log_summary.txt > summary/weblog.html`