library(SingleCellExperiment)
library(scuttle)
# 定义要加载的包向量
# assays （主数据）、 
# colData （单元元数据）、 
# rowData / rowRanges （特征元数据）
# metadata 槽（ SingleCellExperiment 其他）
#reducedDims 插槽专门设计用于存储通过 PCA 和  t-SNE 等方法获得的主数据的降维表示

BiocManager::install(c('scuttle', 'scran', 'scater', 'uwot', 'rtracklayer'))
getwd()
setwd("D:/REXC/ex")
#第一步标准化数据，为数据建立类似与Read10×的一个对象
#通过BioFileCache 在ebi下载count_Calero文件，
library(BiocFileCache)
#install.packages("XML")
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
 #colData(sce)：存储细胞水平的元数据，目前为有细胞id但是没有更多信息？
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
#下载col的元数据
lun.sdrf <- bfcrpath(bfc, file.path("https://www.ebi.ac.uk/arrayexpress/files",
                                    "E-MTAB-5522/E-MTAB-5522.sdrf.txt"))
#读取细胞元数据
coldata <- read.delim(lun.sdrf, check.names=FALSE)

#筛选counts_Calero的数据，先提取列用于判断条件，再基于条件筛选行并保留所有列
#前面的，表示，筛选Derived Array Data File这一列，保留所有行，后面的，表示，保留符合条件所有列，包括了样本 ID、处理组等其他所有列
coldata <- coldata[coldata[, "Derived Array Data File"] == "counts_Calero_20160113.tsv", ]

#将colota数据转换为 SingleCellExperiment 对象的 colData，添加细胞元数据
sce <- SingleCellExperiment(assays = list(counts = mat), colData = coldata)
colnames(colData(sce))
#计算每个细胞的质量控制（QC）指标，并添加到sce的colData中
sce <- scuttle::addPerCellQC(sce)
colnames(colData(sce))

#为sce添加rowData(基因元数据)
rowData(sce)$Length <- gene.length

rowRanges(sce)
#准备填充rowRanges

# 下载数据rowRanges数据，包含了基因的位置数据的gtf文件
mm10.gtf <- bfcrpath(bfc, file.path("http://ftp.ensembl.org/pub/release-82",
                                    "gtf/mus_musculus/Mus_musculus.GRCm38.82.gtf.gz"))
library(XML)
library(restfulr)
#读取gtf文件
gene.data <- rtracklayer::import(mm10.gtf)

#从 gene.data 中筛选出 type 列的值为 "gene" 的记录，仅保留 “基因” 类型的注释信息，过滤掉其他类型（如转录本、外显子等）的记录。
gene.data <- gene.data[gene.data$type == "gene"]
#将 gene.data 对象的元素名称（names）设置为该对象中 gene_id 列的值，
#本质上是为每个个基因注释记录分配一个唯一的标识符（基因 ID），方便后续按基因 ID 快速索引或匹配
names(gene.data)
#names其实不是gene.data的某个子对象，而是类似与字典的索引功能
#以把 gene.data 想象成一个贴着标签的文件柜：
# 每个文件（元素）里面装着基因的详细信息（位置、类型等，对应 gene.data 中的实际数据）。
# names(gene.data) 就是每个文件的标签（贴在文件外面），标签内容是 gene_id。
# 当你用 gene.data["ENSMUSG00000000001"] 查找时，就像直接按标签找文件，无需打开每个文件查看里面的 gene_id 内容
names(gene.data) <- gene.data$gene_id

#在gene.data获取这些 “gene_” 开头列的位置索引，而非直接提取列名和数据
is.gene.related <- grep("gene_", colnames(mcols(gene.data)))
#筛选
mcols(gene.data) <- mcols(gene.data)[, is.gene.related]
colnames(mcols(gene.data))
#将gene.data通过基因id与rowRang关联 rownames(),实际上就是为rowranges添加细胞的gene.data的数据
rowRanges(sce) <- gene.data[rownames(sce)]
rowRanges(sce)[1:10]

#添加其他元数据
#在metad中添加数据的方法1
my_genes <- c("gene_1", "gene_5")
metadata(sce) <- list(favorite_genes = my_genes)

#2
your_genes <- c("gene_4", "gene_8")
metadata(sce)$your_genes <- your_genes
#筛选列，细胞
wt.only <- sce[, sce$phenotype == "wild type phenotype"]
ncol(counts(wt.only))
colData(wt.only)

# 筛选行 基因
coding.only <- sce[rowData(sce)$gene_biotype == "protein_coding",]
nrow(counts(coding.only))     
rowData(coding.only)

#合并多个对象的列
sce2 <-  cbind(sce, sce)
colData(sce2)

#合并多个对象的行
sce2 <- rbind(sce, sce)
nrow(counts(sce2))


sce <- scater::logNormCounts(sce)
sce <- scater::runPCA(sce)
#dim() 函数返回矩阵的维度，格式为 c(细胞数, 主成分数)
dim(reducedDim(sce, "PCA"))
#单细胞数据的表达矩阵通常是高维的（比如 10000 个基因 × 5000 个细胞），每个基因都是一个 “维度”。但这些基因的表达并非完全独立（比如功能相关的基因表达模式相似），PCA 的作用就是：
# 找到少数几个 “综合维度”（主成分），用它们来 “概括” 原始数据的主要信息；
# 每个主成分都是多个基因表达量的加权组合（权重反映基因对该主成分的贡献度）。
# 例如：
# PC1（第一主成分）：是能解释原始数据中最大变异的维度，代表数据中最显著的表达模式（比如细胞间最主要的差异来源，可能是细胞类型差异）；
# PC2（第二主成分）：解释第二大变异的维度，且与 PC1 不相关（代表另一类独立的表达模式）；
# 以此类推，PC3、PC4…… 解释的变异程度依次递减。
# 在单细胞分析中的意义：
# 降维简化分析：用少数主成分（比如前 20 个）替代上万基因，大幅降低计算复杂度，同时保留核心信息。
# 去除噪音：PCA 会优先保留变异大的信号（生物学差异），过滤掉变异小的噪音（技术误差）。
# 可视化基础：后续的 UMAP、t-SNE 等降维可视化通常基于 PCA 结果（而非原始基因），让细胞分群更清晰（比如 PC1-PC2 散点图可初步展示细胞聚类）。