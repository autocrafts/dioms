library(SingleCellExperiment)
# 定义要加载的包向量
BiocManager::install(c('scuttle', 'scran', 'scater', 'uwot', 'rtracklayer'))
getwd()
setwd("D:/REXC/ex")
#第一步标准化数据，为数据建立类似与Read10×的一个对象
#通过BioFileCache 在ebi下载count_Calero文件，
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
calero.counts <- bfcrpath(bfc, file.path("https://www.ebi.ac.uk/biostudies", 
                                         "files/E-MTAB-5522/counts_Calero_20160113.tsv"))
  #EBI数据库 → 58012ba2759_counts_Calero_20160113.tsv（原始数据）
  
 # BiocFileCache.sqlite（记录缓存元数据） + LOCK文件（保证操作安全）                                      "files/E-MTAB-5522/counts_Calero_20160113.tsv"))
 #读取数据，得到基因表达计数矩阵
 mat <- read.delim(calero.counts, header = TRUE, row.names = 1, check.names = FALSE)
 #header = TRUE:第一行是列名——细胞条形码，cell barcodes。  row.names = 1  第一列为为行名——基因标识符，如基因名或 Ensembl ID

 #筛选出 ERCC--外部 RNA 对照联盟  ^ 是正则表达式中的 “起始位置” 符号
 spike.mat <- mat[grepl("^ERCC-", rownames(mat)), ]
 #ERCC通常会被用来进行数据质量评估，标准化校正，批次效应去除。
 
 #保留小数基因，mat[..., ]列不变，保留细胞id
 mat <- mat[grepl("^ENSMUSG", rownames(mat)),] 
 #[col, row] 
 gene.length <- mat[,1]
 #在 R 中，[ ] 索引里的负号（-）是一种排除语法，用于从数据结构中移除指定的元素（行或列）
 mat <- as.matrix(mat[, -1]) #保留行，清除第一列
 #构建 SingleCellExperiment 对象 assays 接受一个命名列表，用于存储不同类型的表达矩阵
sce <- SingleCellExperiment(assays = list(counts = mat))
  #assays(sce)：存储表达数据的容器                       
 #colData(sce)：存储细胞水平的元数据，目前为空
  #rowData(sce)：存储基因水平的元数据，目前仅包含基因名
 #dimnames(sce)：记录基因名（行名）和细胞名（列名），与 mat 的行名和列名一致。


#添加更多 assays
#对象中的原始表达计数计数数据counts进行对数标准化，存储在sce@assays@data@listData[["logcounts"]]
sce <- logNormCounts(sce)

#dim() 函数返回该矩阵的维度：第一个值：矩阵的行数基因数， 第二个值：矩阵的列数细胞数
dim(logcounts(sce))
 

#演示在assays中对counts出来，并添加在assays中
counts_100 <- counts(sce) +100
assay(sce, "counts_100") <- counts_100  #注意并非assays

#删除刚生成的counts_100
assays(sce) <- assays(sce)[1:2]
names(assays(sce))   

#处理元数据

