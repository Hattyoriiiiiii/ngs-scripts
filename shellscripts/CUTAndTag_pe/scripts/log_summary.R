#!/usr/bin/Rscript

suppressPackageStartupMessages(library(dplyr))
suppressMessages(library(stringr))
suppressWarnings(library(ggplot2))


projPath <- getwd()
projName <- projPath %>% basename()
seqType <- strsplit(projName, "_")[[1]][1]

# summaryディレクトリへのパス
alignResPath <- paste0(projPath, "/summary/bowtie2/")
dupRes <- paste0(projPath, "/summary/picard/")


summaryBowtie2 <- function(alignResPath) {

    # sampleのprefix
    sampleList <- alignResPath %>% list.files() %>% lapply(function(x) strsplit(x, "_bowtie2.txt")[[1]][1]) %>% unlist()

    ########## factorの取得 ##########
    # factor_celltype_condition_replicate: conditionはwtとかでも書く
    factor_target <- sampleList %>% lapply(function(x) strsplit(x, "_")[[1]][1]) %>% unlist() %>% unique()
    factor_cell <- sampleList %>% lapply(function(x) strsplit(x, "_")[[1]][2]) %>% unlist() %>% unique()
    factor_condition <- sampleList %>% lapply(function(x) strsplit(x, "_")[[1]][3]) %>% unlist() %>% unique()
    factor_replicate <- sampleList %>% lapply(function(x) strsplit(x, "_")[[1]][4]) %>% unlist() %>% unique()
    factor_experiment <- sampleList %>% lapply(function(x) paste(strsplit(x, "_")[[1]][2], strsplit(x, "_")[[1]][3], sep="_")) %>% unlist() %>% unique()

    # print(factor_target)
    # print(factor_cell)
    # print(factor_condition)
    # print(factor_replicate)
    # print(factor_experiment)


    ########## 集計 ##########

    alignResult <- c()
    for (sample in sampleList) {
        alignRes <- read.table(paste0(alignResPath, sample, "_bowtie2.txt"), header=FALSE, fill=TRUE)
        alignRate <- substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
        alignUniq <- substr(alignRes$V2[4], 2, nchar(as.character(alignRes$V2[4]))-2)

        tmp <- sub(paste0(factor_target, "_"), replacement="", sample) %>% strsplit("_")
        cell_type <- tmp[[1]][1]
        condition <- tmp[[1]][2]
        replicate <- tmp[[1]][3]
        experiment <- paste0(tmp[[1]][1], "_", tmp[[1]][2])
        
        alignResult <- data.frame(Cell = cell_type,
                                Condition = condition,
                                Replicate = replicate,
                                Experiment = experiment,
                                SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric,
                                MappedFragNum_mm10 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric,
                                AlignmentRate_mm10 = alignRate %>% as.numeric,
                                AlignmentUniquely = alignUniq %>% as.numeric) %>% rbind(alignResult, .)
    }

    alignResult$Cell = factor(alignResult$Cell, levels=factor_cell)
    alignResult$Condition = factor(alignResult$Condition, levels=factor_condition)
    alignResult$Replicate = factor(alignResult$Replicate, levels=factor_replicate)
    alignResult$Experiment = factor(alignResult$Experiment, levels=factor_experiment)

    alignResult <- alignResult %>% mutate(AlignmentRate_mm10 = paste0(AlignmentRate_mm10, "%"))
    return(alignResult)
}


alignResult <- summaryBowtie2(
    alignResPath=alignResPath
)

write.table(
    alignResult, 
    "./summary/log_summary.txt", 
    quote=FALSE, 
    row.names=FALSE,
    col.names=FALSE, 
    sep="\t"
    )


########## main実行部 ##########

if (seqType == "CUTAndTag") {
    print("CUT&Tag-seq")
    # factor_celltype_condition_replicate
    # CUTAndTag_summary()
} else if (seqType == "RNA") {
    print("RNA-seq")
    # RNA_summary()
} else if (seqType == "ATAC") {
    print("ATAC-seq")
    # ATAC_summary()
} else if (seqType == "ChIP") {
    print("ChIP-seq")
    # ChIP_summary()
} else {
    print("Invalid seqType")
}


