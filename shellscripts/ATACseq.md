
## Reference

https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/

## Peak Call

### For single-end data
```
--nomodel --shift -100 --extsize 200
```

### For paired-end data

```
# -f BAMPE, use paired-end information
# --keep-dup all, keep all duplicate reads.
macs2 callpeak -f BAMPE -g hs --keep-dup all --cutoff-analysis -n sample1 \
  -t sample1.shifted.bam --outdir macs2/sample1 2> macs2.log
```

### Creating browser tracks

`bamCoverage`

- binサイズを小さくするとカバレッジtrackの解像度が高くなるが、ファイルサイズも大きくなる
- 1x 正規化(RPGC)ではリファレンスゲノムのマッピング可能な領域のゲノムサイズ(the effective genome size)の値を入力とする。
- 

> The 1x normalization (RPGC) requires the input of a value for the effective genome size, which is the mappable part of the reference genome. Of course, this value is species-specific. The command line help of this tool offers suggestions for a number of model species.


```
"""
using additional options (smaller bin size for higher resolution, normalizing coverage to 1x mouse genome size)
"""

cores=16
effect_genome_size=2494787188  # https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html


bamCoverage \
    --bam HSF5_RS_42_rep1.sorted.concordant.rmdup.bam \
    -o test01_HSF5_RS_42_rep1.sorted.concordant.rmdup.bw \
    --numberOfProcessors $cores \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize $effect_genome_size \
    --minMappingQuality 10 \
    --extendReads
```

- test01
```
for i 
bamCoverage \
    --bam HSF5_RS_42_rep1.sorted.concordant.rmdup.bam \
    -o test01_HSF5_RS_42_rep1.sorted.concordant.rmdup.bw \
    --numberOfProcessors $cores \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize $effect_genome_size \
    --minMappingQuality 10 \
    --extendReads


```

```
Due to filtering, 48.11754050368845% of the aforementioned alignments will be used 387121.0109651348
normalization: 1x (effective genome size 2494787188)
bamFilesList: ['HSF5_RS_42_rep1.sorted.concordant.rmdup.bam']
binLength: 10
numberOfSamples: None
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 331
numberOfProcessors: 16
verbose: False
region: None
bedFile: None
minMappingQuality: 10
ignoreDuplicates: False
chrsToSkip: []
stepSize: 10
center_read: False
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: None
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 1324

Due to filtering, 27.95391449819985% of the aforementioned alignments will be used 123670.35405319599
normalization: 1x (effective genome size 2494787188)
bamFilesList: ['HSF5_RS_42_rep2.sorted.concordant.rmdup.bam']
binLength: 10
numberOfSamples: None
blackListFileName: None
skipZeroOverZero: False
bed_and_bin: False
genomeChunkSize: None
defaultFragmentLength: 285
numberOfProcessors: 16
verbose: False
region: None
bedFile: None
minMappingQuality: 10
ignoreDuplicates: False
chrsToSkip: []
stepSize: 10
center_read: False
samFlag_include: None
samFlag_exclude: None
minFragmentLength: 0
maxFragmentLength: 0
zerosToNans: False
smoothLength: None
save_data: False
out_file_for_raw_data: None
maxPairedFragmentLength: 1140
```


わからないからとりあえずこれ。

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