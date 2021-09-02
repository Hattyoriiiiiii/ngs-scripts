FROM jupyter/datascience-notebook:latest

LABEL maintainer "Tatsuya Hattori <www.tatsuya92@gmail.com>"


# install required R packages
RUN R -e "install.packages('BiocManager', repos='https://cran.ism.ac.jp/')" 

# RNA-seq
RUN R -e "BiocManager::install(c('edgeR', 'DESeq2', 'limma', 'topGO', 'org.Mm.eg.db', 'ReactomePA', 'pathview', 'clusterProfiler'))"

# ChIP-seq/ATAC-seq
RUN R -e "BiocManager::install(c('MotifDb', 'TFBSTools', 'rGADEM', 'JASPAR2018'))"

# others
RUN R -e "BiocManager::install(c('qvalue', 'plot3D', 'ggplot2', 'pheatmap', 'cowplot', 'cluster', 'NbClust', 'fastICA', 'NMF', 'matrixStats', 'Rtsne', 'mosaic', 'knitr', 'genomation', 'ggbio', 'Gviz', 'RUVSeq', 'gProfileR', 'ggfortify', 'corrplot', 'gage', 'EDASeq', 'citr', 'formatR', 'svglite', 'ShortRead', 'QuasR', 'methylKit', 'FactoMineR', 'iClusterPlus', 'enrichR', 'caret', 'xgboost','glmnet', 'DALEX', 'kernlab', 'pROC', 'nnet', 'RANN', 'ranger', 'GenomeInfoDb', 'GenomicRanges', 'GenomicAlignments', 'ComplexHeatmap', 'circlize', 'rtracklayer', 'tidyr', 'AnnotationHub', 'GenomicFeatures', 'normr'))"

# reference
RUN R -e "options(timeout=6000); BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38', 'BSgenome.Mmusculus.UCSC.mm10', 'BSgenome.Scerevisiae.UCSC.sacCer3', 'BSgenome.Dmelanogaster.UCSC.dm6'))"

# read processing and quality check
RUN R -e "BiocManager::install(c('Rqc', 'fastqcr', 'QuasR', 'ShortRead', 'ATACseqQC'))" && \
    R -e "BiocManager::install(c('EnsDb.Hsapiens.v75','EnsDb.Hsapiens.v79', 'EnsDb.Hsapiens.v86', 'EnsDb.Mmusculus.v79'))"

# Common R Packages for data science
RUN R -e "install.packages(c('data.table', 'dplyr', 'tidyr', 'stringr', 'tidyverse', 'readr', 'caret', 'Matrix'), repos='https://cran.ism.ac.jp/')" && \
    R -e "install.packages(c('prcomp', 'randomForest'), repos='https://cran.ism.ac.jp/')"

# R for jupyterbook
RUN R -e "install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest'), repos='https://cran.ism.ac.jp/')" && \
    R -e "remotes::install_github('IRkernel/IRkernel')" && \
    R -e "IRkernel::installspec()"

# csaw
RUN R -e "BiocManager::install(c('csaw', 'chipseqDBData', 'edgeR', 'Rsubread', 'chromVAR', 'motifmatchr', 'SummarizedExperiment', 'BiocParallel'))"

RUN R -e "BiocManager::install(c('GO.db', 'reactome.db', 'org.Mm.eg.db', 'TxDb.Mmusculus.UCSC.mm10.knownGene', 'JASPAR2020'))"

# ggpubr, ATACseqQC