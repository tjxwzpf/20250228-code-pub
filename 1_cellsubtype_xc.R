# 发现关键细胞亚群

library(Seurat)
gc()

# 
# # 读取数据
# SLE1=Read10X('./SLE1/')
# SLE2=Read10X('./SLE2/')
# CT1=Read10X('./CT1/')
# CT3=Read10X('./CT3/')
# CT3=CT3[[1]]
# OLD1=Read10X('./OLD1/')
# OLD1=OLD1[[1]]
# OLD2=Read10X('./OLD2/')
# OLD2=OLD2[[1]]
# 
# rm(scRNA)
# 
# 
# ## 建立seurat对象
# SLE1=CreateSeuratObject(SLE1,project = 'SLE1')
# SLE2=CreateSeuratObject(SLE2,project = 'SLE2')
# CT1=CreateSeuratObject(CT1,project = 'CT1')
# CT3=CreateSeuratObject(CT3,project = 'CT3')
# OLD1=CreateSeuratObject(OLD1,project = 'OLD1')
# OLD2=CreateSeuratObject(OLD2,project = 'OLD2')

library(Seurat)


# # scRNA <- readRDS('./GSE157007 衰老外周血文献/GSE157007_RAW/scelist_9个样本_6衰老3年轻.RDS')
# scRNA <- readRDS('./scRNA_anno_正确origident_eqtl_2vs2.RDS')
# ### 合并scRNA1,scRNA2,scRNA3
# ## 没有正常样本则不加scRNA1
# ### 更多样本可以自行加入
# #scRNA <- merge(scRNA1, y = c(scRNA2,scRNA3))
# scRNA[1]
# scRNA[-1]
# scRNA <- merge(scRNA[[1]], y = scRNA[-1])
# # scRNA <- merge(SLE1,y=c(SLE2,CT1,CT3,OLD1,OLD2))
# 
# # metadata为样本信息，我们需要定义分组
# scRNA@meta.data$tissue_type=scRNA@meta.data$orig.ident
# table(scRNA$orig.ident)
# 
# scRNA@meta.data$tissue_type <- gsub("F012","ctrl",scRNA@meta.data$tissue_type)
# scRNA@meta.data$tissue_type <- gsub("F013","ctrl",scRNA@meta.data$tissue_type)
# scRNA@meta.data$tissue_type <- gsub("F014","ctrl",scRNA@meta.data$tissue_type)
# scRNA@meta.data$tissue_type <- gsub("F020","old",scRNA@meta.data$tissue_type)
# scRNA@meta.data$tissue_type <- gsub("F021","old",scRNA@meta.data$tissue_type)
# scRNA@meta.data$tissue_type <- gsub("F023","old",scRNA@meta.data$tissue_type)
# scRNA@meta.data$tissue_type <- gsub("OH14","old",scRNA@meta.data$tissue_type)
# scRNA@meta.data$tissue_type <- gsub("OH15","old",scRNA@meta.data$tissue_type)
# scRNA@meta.data$tissue_type <- gsub("OH17","old",scRNA@meta.data$tissue_type)
# table(scRNA$tissue_type)
# saveRDS(scRNA,'./GSE157007 衰老外周血文献/GSE157007 衰老外周血 3young6old.RDS')
# 
# 
# #衰老外周血,
# # scRNA <- readRDS('./GSE157007 衰老外周血文献/GSE157007 衰老外周血 3young6old.RDS')
# scRNA <- readRDS('../../eqtl_心衰/GSE157007 衰老外周血文献/GSE157007 衰老外周血 3young6old.RDS')
# table(scRNA$tissue_type)
# table(scRNA$orig.ident)
# 
# # scRNA <- subset(scRNA, tissue_type %in% c("ctrl"))
# 
# 
# # #东方医院17样本
# # scRNA <- readRDS('./sce_data.merge_17pbmc_东方医院GSE213516.RDS')
# # table(scRNA$tissue_type)
# # 
# # scRNA$tissue_type <- "ctrl"
# # table(scRNA$tissue_type)
# 
# 
# 
# 
# ########################################## sepsis五个样本              ############################              
# ############################ ##########################################
# 
# # fs=list.files('./','txt.gz')
# fs=list.files('../0_Cellranger/',)
# fs
# fs <- fs[c(1,3,5,7,9)]
# fs
# 
# samples <- fs
# 
# sce <- list()
# 
# for (samples in fs) {
#   #pro=samples[1]
#   # folder=file.path(dir ,pro,paste0('/outs/filtered_feature_bc_matrix') )
#   # aFB=read.table(paste0('./',samples))
#   
#   # print(pro)
#   # print(folder)
#   # print(list.files(folder))
#   sce[[samples]]=CreateSeuratObject(counts = Read10X(paste0('../0_Cellranger/',samples)),
#                                     min.cells = 5,
#                                     min.features = 500,
#                                     project = samples )
# }
# sce
# 
# 
# # YoungLiver=CreateSeuratObject(YoungLiver,project='YoungLiver')
# # OldLiver=CreateSeuratObject(OldLiver,project='OldLiver')
# # sham3=CreateSeuratObject(sham3,project='sham3')
# # MCAO1=CreateSeuratObject(MCAO1,project='MCAO1')
# # MCAO2=CreateSeuratObject(MCAO2,project='MCAO2')
# # MCAO3=CreateSeuratObject(MCAO3,project='MCAO3')
# 
# # sce <- sceList
# 
# sce[[1]]
# # ### merge一键融合，从15-20个样本
# pbmc <- merge(sce[[1]],sce[-1])
# 
# pbmc
# table(pbmc@meta.data$orig.ident)
# # pbmc=merge(YoungLiver,c(OldLiver,sham3,MCAO1,MCAO2,MCAO3))
# # pbmc=merge(YoungLiver,c(OldLiver))
# table(pbmc$orig.ident)
# 
# # pbmc <- ss2
# #MSCs样本的切割
# a <- pbmc$orig.ident
# a
# table(a)
# 
# 
# library(stringr)
# a <- gsub("A--","A",a)
# a <- gsub("AS-","AS",a)
# a <- gsub("S-","S",a)
# a <- gsub("5d","7d",a)
# a
# table(a)
# result=str_split(a,'-',simplify = T)
# dim(result)
# result=as.data.frame(result)
# # result$gene <- paste0(gs$gene)
# rownames(result) <- rownames(pbmc@meta.data)
# dim(result)
# colnames(result) <- c("Patient","Time")
# 
# # metadata <- result
# head(metadata)
# pbmc@meta.data$Patient=result$Patient
# pbmc@meta.data$Time=result$Time
# # pbmc@meta.data$Sample=metadata$Sample
# # pbmc@meta.data$Capture=metadata$Capture
# 
# 
# # pbmc@meta.data$sample <- NULL
# table(pbmc@meta.data$Patient)
# table(pbmc@meta.data$Time)
# # table(pbmc@meta.data$sample)
# table(pbmc@meta.data$Capture)
# table(pbmc@meta.data$Sample)
# 
# colnames(pbmc)
# 
# ##########################################      再次更改         
# ############################               ############################ ##########################################
# #MSCs样本的切割
# a <- pbmc$orig.ident
# a
# table(a)
# 
# 
# library(stringr)
# a <- gsub("A--","",a)
# a <- gsub("X","S-",a)
# # a <- gsub("S-","S",a)
# a <- gsub("5d","7d",a)
# a
# table(a)
# table(a)
# result=str_split(a,'-',simplify = T)
# dim(result)
# result=as.data.frame(result)
# 
# rownames(result) <- rownames(pbmc@meta.data)
# dim(result)
# result$Patient <- paste0(result$V1,"-",result$V2)
# table(result$Patient)
# table(pbmc$Patient)
# colnames(result)[3] <- "Time"
# 
# table(result$Time)
# table(pbmc$Time)
# # colnames(result) <- c("Patient","Time")
# 
# # metadata <- result
# head(metadata)
# pbmc@meta.data$Patient=result$Patient
# pbmc@meta.data$orig.ident=a
# pbmc@meta.data$Time=result$Time
# # pbmc@meta.data$Time=result$
# table(pbmc$orig.ident)
# 
# pbmc$tissue_type <- "sepsis"
# 
# Idents(pbmc) <- "orig.ident"
# 
# saveRDS(pbmc,file = './scRNA_anno_正确origident5samples.RDS')
# 
# 
# 
# ##########################################   开始整合两个数据            ############################              
# ############################ ##########################################
# 
# 
# #衰老外周血, 这里不适用
# # scRNA <- readRDS('./GSE157007 衰老外周血文献/GSE157007 衰老外周血 3young6old.RDS')
# scRNA <- readRDS('../../eqtl_心衰/GSE157007 衰老外周血文献/GSE157007 衰老外周血 3young6old.RDS')
# # scRNA <- merge(scRNA[[1]],scRNA[-1]) #
# 不需要从merge这一步开始，其实是一样的，都还没经过过滤
# 
# scRNA
# table(scRNA$tissue_type)
# table(scRNA$orig.ident)
# scRNA
# DimPlot(scRNA)
# 
# Idents(scRNA) <- "orig.ident"
# 
# 
# 
# table(scRNA$tissue_type)
# scRNA$tissue_type <- gsub("ctrl", "young", scRNA$tissue_type)
# table(scRNA$tissue_type)
# 
# scRNA2 <- readRDS('./scRNA_anno_正确origident5samples.RDS')
# # Idents(scRNA2) <- "orig.ident"
# # saveRDS(scRNA2,file = './scRNA_anno_正确origident5samples.RDS')
# 
# 
# 
# table(scRNA2$tissue_type)
# 
# table(scRNA2$orig.ident)
# # scRNA2$tissue_type <- "ICM_HF"
# table(scRNA2$tissue_type)
# table(scRNA$tissue_type)
# 
# table(Idents(scRNA))
# table(Idents(scRNA2))
# 
# scRNA <- merge(scRNA,scRNA2)
# scRNA
# 
# table(scRNA$tissue_type)
# table(scRNA$orig.ident)
# table(scRNA$Patient)







MSCs治疗心衰 PBMC 从这里开始
##########################################               ############################              
############################ ##########################################

scRNA <- readRDS('./scRNA_anno_正确origident_eqtl_2vs2.RDS')












#install.packages('stringr')
# 去除数字,获得干净的分组
# scRNA@meta.data$tissue_type=stringr::str_remove(scRNA@meta.data$tissue_type,'[0-9]')

# 质控
##计算质控指标
#计算细胞中线粒体核糖体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)
library(ggplot2)
Idents(scRNA) <- "orig.ident"
col.num <- length(levels(scRNA@active.ident))
# 过滤前
table(Idents(scRNA))

P_pre1 <- VlnPlot(scRNA,features = c("nFeature_RNA", "nCount_RNA", 
                                     "percent.mt","percent.HB"), 
        cols =rainbow(col.num), 
        pt.size = 0, #不需要显示点，可以设置pt.size = 0
        ncol = 4) 
  # theme(axis.title.x=element_blank(), 
  #       axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

P_pre1

##设置质控标准，比较随意
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=200
maxGene=4000
pctMT=10

##数据质控
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(scRNA@active.ident))
p_after1 <- VlnPlot(scRNA,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        cols =rainbow(col.num), 
        pt.size = 0, 
        ncol = 4) 
  # theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

p_after1

P_pre1+p_after1

library(patchwork) 
pdf(file="./figures/day1_质控.pdf", width=13.5, height=6)
print(P_pre1/p_after1)
dev.off()



scRNA

table(scRNA$orig.ident)

# 标准化
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

#降维聚类#########################
library(Seurat)
library(tidyverse)
library(patchwork)


scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")

### harmony去批次
#BiocManager::
library(harmony)
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")
DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident")

ElbowPlot(scRNA_harmony,reduction = 'harmony')

# 一定要指定“harmony”！
scRNA <- FindNeighbors(scRNA_harmony, dims = 1:30, reduction = "harmony")
scRNA <- FindClusters(scRNA)
scRNA <- RunUMAP(scRNA, dims = 1:30,reduction = 'harmony')


# scRNA <- pbmc
# 去批次成功
dev.off()

DimPlot(scRNA,split.by = 'tissue_type')

scRNA
table(scRNA$orig.ident)
table(scRNA$tissue_type)
saveRDS(scRNA, file = "scRNA_anno.RDS")


rm(scRNA_harmony)
gc()

# 
# scRNA <- readRDS('scRNA_anno.RDS')
# #BiocManager
# library(SingleR)
# # 人用下面
# refdata <- SingleR::HumanPrimaryCellAtlasData()
# 
# # 鼠用下面
# #refdata <- SingleR::MouseRNAseqData()
# 
# library(Seurat)
# testdata <- GetAssayData(scRNA, slot="data")
# clusters <- scRNA@meta.data$seurat_clusters
# cellpred <- SingleR(test = testdata, ref = refdata,
#                     labels =refdata$label.main,
#                     method = "cluster", clusters = clusters, 
#                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
# 
# celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
# 
# 
# scRNA@meta.data$celltype = "NA"
# for(i in 1:nrow(celltype)){
#   scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
# 
# 
# Idents(scRNA)=scRNA$celltype
# 
# library(RColorBrewer)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
# 
# 
# scRNA$predicted.celltype.l2
# 
# table(scRNA$celltype)
# 
# scRNA$celltype <- gsub("HSC_-G-CSF","Monocyte",scRNA$celltype)
# 
# scRNA$celltype 
# scRNA$celltype_singleR <- scRNA$celltype 
# scRNA$celltype 
# 
# 
# col_vector <- readRDS('col_vector用于eqtl.RDS')
# library(patchwork) 
# dir.create("./figures")
# pdf(file="./figures/day1_dimplot_singleR.pdf", width=10, height=3.5)
# 
# DimPlot(scRNA,split.by = 'tissue_type',
#         # group.by = "predicted.celltype.l2",
#         group.by = "celltype",
#         label = T, repel = T,
#         cols = col_vector)
# # print(p_top)
# dev.off()
# 
# 
# # library(Azimuth)
# 
# ## 保存数据
# saveRDS(scRNA,file = 'scRNA_anno.RDS')


##########################################      azimuth注释         ############################              
############################ ##########################################

# 现在的服务器用hdf5r
dyn.load("/home/wangxc/miniconda3/envs/scVelo/lib/libhdf5.so.200")
dyn.load("/usr/bin/R-4.1.1/lib64/R/library/hdf5r/libs/libhdf5_hl.so.200")
library(hdf5r)
library("Azimuth")



library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)


scRNA1 <- readRDS('scRNA_anno.RDS')
# scRNA1 <- scRNA
table(scRNA1$orig.ident)
scRNA1

gc()
# The RunAzimuth function can take a Seurat object as input
scRNA1Azi <- RunAzimuth(scRNA1, reference = "pbmcref")

# `RunAzimuth`
# ?RunAzimuth
scRNA1Azi

# ?RunAzimuth
colnames(scRNA1Azi@meta.data)

p1 <- DimPlot(scRNA1Azi, group.by = "predicted.celltype.l2", 
              reduction = "ref.umap",repel = T,
              label = TRUE, label.size = 3) + NoLegend()\

dev.off()
p1

scRNA1Azi$celltype <- scRNA1Azi$predicted.celltype.l2
p2 <- DimPlot(scRNA1Azi, group.by = 'celltype',reduction = "ref.umap")
p2
# p2 <- DimPlot(scRNA1Azi, group.by = 'orig.ident',reduction = "ref.umap")
p1 + p2

p1 <- DimPlot(scRNA1Azi, group.by = "predicted.celltype.l2", 
              # reduction = "ref.umap",
              repel = T,
              label = TRUE, label.size = 3) + NoLegend()
p1
p2 <- DimPlot(scRNA1Azi, group.by = 'celltype')
p1 + p2

DimPlot(scRNA1Azi, group.by = 'celltype',split.by = "tissue_type",
        repel = T, label = T)

dev.off()
gc()
p3 <- DimPlot(scRNA1Azi, group.by = 'predicted.celltype.l2',
              reduction = "ref.umap",
              label = T, repel = T,
              split.by = 'tissue_type') + NoLegend()
p3
gc()

# library(RColorBrewer)
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))
# mycolors <- getPalette(25)
# 
# col_vector <- sample(getPalette(30))
# 
# sample_color <- col_vector[1:10]
# 
# saveRDS(col_vector, file = "col_vector用于eqtl.RDS")

# scRNA1Azi <- scRNA
scRNA1Azi$tissue_type <- factor(scRNA1Azi$tissue_type, levels = c("before", "after"))

library(ggplot2)
pdf(file="./figures/day1_dimplot.pdf", width=7, height=3)

DimPlot(scRNA1Azi,split.by = 'tissue_type',
        group.by = "predicted.celltype.l2",
        # group.by = "celltype",
        reduction = "ref.umap",
        label = F, repel = T,
        cols = col_vector)+theme_bw()
# print(p_top)
dev.off()

pdf(file="./figures/day1_dimplot2.pdf", width=5, height=5)

DimPlot(scRNA1Azi,
        # split.by = 'tissue_type',
        group.by = "predicted.celltype.l2",
        # group.by = "celltype",
        reduction = "ref.umap",
        label = T, repel = T,
        cols = col_vector)+theme_bw()+NoLegend()
# print(p_top)
dev.off()



## 保存数据
saveRDS(scRNA1Azi,file = 'scRNA_anno注释azimuth.RDS')
scRNA1Azi

# table(scRNA$orig.ident)
# table(scRNA$celltype)
# scRNA
# scRNA1Azi <- scRNA
# scRNA1Azi <- readRDS('scRNA_anno.RDS')
scRNA1Azi <- readRDS('scRNA_anno注释azimuth.RDS')
scRNA1Azi


# scRNA1Azi$tissue_type3 <- scRNA1Azi$tissue_type
scRNA1Azi$tissue_type
table(scRNA1Azi$tissue_type)
# scRNA1Azi$tissue_type <- gsub("old","ctrl",)
# scRNA1Azi$tissue_type

table(scRNA1Azi$orig.ident)
table(scRNA1Azi$predicted.celltype.l2)
# View(scRNA1Azi@meta.data)



# scRNA1Azi$tissue_type <- factor(scRNA1Azi$tissue_type, levels = c("young","old","sepsis"))
scRNA1Azi$celltype <- NULL
scRNA1Azi$celltype <- scRNA1Azi$predicted.celltype.l2
scRNA1Azi$celltype


# DimPlot(scRNA1Azi, group.by = 'predicted.celltype.l1',reduction = "ref.umap",
DimPlot(scRNA1Azi, group.by = 'predicted.celltype.l2',reduction = "ref.umap",
        label = T, repel = T,
        split.by = 'tissue_type') + NoLegend()

# DimPlot(scRNA1Azi, group.by = 'predicted.celltype.l3',reduction = "ref.umap",
#         label = T, repel = T,
#         split.by = 'tissue_type') + NoLegend()

scRNA <- scRNA1Azi
table(scRNA$orig.ident)
table(scRNA$predicted.celltype.l1)


scRNA$celltype <- scRNA$predicted.celltype.l1
scRNA$celltype <- scRNA$predicted.celltype.l2

scRNA <- readRDS('scRNA_anno注释azimuth.RDS')
table(scRNA$celltype, scRNA$tissue_type)
table(scRNA$celltype, scRNA$orig.ident) 
table(scRNA$orig.ident)
table(scRNA$celltype)
 

scRNA
scRNAsub <- scRNA

table(scRNA$celltype)
## 细胞比例图
library(reshape2)
library(ggplot2)
library(dplyr)
prop_df <- table(scRNAsub@meta.data$celltype,
                 scRNAsub@meta.data$orig.ident) %>% melt()

# prop_df <- table(scRNAsub@meta.data$celltype,
#                  scRNAsub@meta.data$tissue_type) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 


# library(RColorBrewer)
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))
# mycolors <- getPalette(25)
# 
# col_vector <- sample(getPalette(30))
# 
# sample_color <- col_vector[1:10] 
# 
# saveRDS(col_vector, file = "col_vector用于eqtl.RDS")

# 作图
prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:30]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
dev.off()
prop


pdf(file="./figures/day1分布比例图.pdf", width=15, height=4)
print(prop)
dev.off()






##########################################   亚群分析            ############################               
############################ ##########################################

unique(scRNA$celltype)

table(scRNA$orig.ident, scRNA$celltype)
### 根据原文 我们直接提取T亚群
scRNA_T=subset(scRNA,celltype %in% c("NK"))
scRNA_T

pdf(file="./figures/day1_亚群.pdf", width=15, height=4)
# print(p_top)

DimPlot(scRNA_T, group.by = 'predicted.celltype.l2',reduction = "ref.umap",
        label = T, repel = T,
        split.by = 'tissue_type') + NoLegend()
dev.off()


scRNA_T
scRNA_T <- CreateSeuratObject(scRNA_T@assays$RNA@counts,
                              meta.data = scRNA_T@meta.data)
scRNA_T

library(Seurat)

# 提T细胞亚群,重新降维聚类
scRNAsub<- FindVariableFeatures(scRNA_T, selection.method = "vst", nfeatures = 3000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub<- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNAsub)


## 重新harmony
library(harmony)
set.seed(1000)
scRNAsub

scRNAsub <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")
DimPlot(scRNAsub, reduction = "harmony", group.by = "orig.ident")

ElbowPlot(scRNAsub,reduction = 'harmony')


scRNAsub <- FindNeighbors(scRNAsub, reduction = 'harmony',dims = 1:10)
scRNAsub <- FindClusters(scRNAsub)
scRNAsub <- RunUMAP(scRNAsub, reduction = 'harmony',dims = 1:10)

table(scRNAsub$predicted.celltype.l2)
table(scRNAsub$predicted.celltype.l3)
scRNAsub
table(scRNAsub$seurat_clusters)
table(scRNAsub$tissue_type)
table(scRNAsub$orig.ident)

# scRNAsub$tissue_type3 <- scRNAsub$tissue_type
# scRNAsub$tissue_type
# scRNAsub$tissue_type <- gsub("old","ctrl",scRNAsub$tissue_type)
table(scRNAsub$tissue_type)

# Idents(scRNAsub) <- "seurat_clusters"
# Idents(scRNAsub) <- "celltype"
dev.off()
table(scRNAsub$predicted.celltype.l2)
DimPlot(scRNAsub, label=TRUE,split.by = 'tissue_type',reduction = "umap") 
DimPlot(scRNAsub, label=TRUE,repel = T,
        split.by = 'tissue_type',
        group.by = "predicted.celltype.l2",
        reduction = "umap") 
DimPlot(scRNAsub, label=TRUE,repel = T,
        split.by = 'tissue_type',
        group.by = "seurat_clusters",
        reduction = "umap") 

dev.off()





colnames(scRNAsub@meta.data)
VlnPlot(scRNAsub, features = "mapping.score")
table(scRNAsub$seurat_clusters, scRNAsub$celltype)
scRNAsub <- subset(scRNAsub, seurat_clusters %in% c('12','13'), invert=T)
table(scRNAsub$tissue_type)
table(scRNAsub$orig.ident)

删掉细胞的话，再来一次harmony


# 提T细胞亚群,重新降维聚类
# scRNAsub <- subset(scRNAsub, subset = mapping.score > 0.5)
scRNAsub

scRNAsub<- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 3000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub<- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNAsub)

## 重新harmony
library(harmony)
set.seed(1000)
scRNAsub

scRNAsub <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")
DimPlot(scRNAsub, reduction = "harmony", group.by = "orig.ident")

ElbowPlot(scRNAsub,reduction = 'harmony')


scRNAsub <- FindNeighbors(scRNAsub, reduction = 'harmony',dims = 1:10)
scRNAsub <- FindClusters(scRNAsub)
scRNAsub <- RunUMAP(scRNAsub, reduction = 'harmony',dims = 1:10)

table(scRNAsub$predicted.celltype.l2)
table(scRNAsub$predicted.celltype.l3)
scRNAsub
table(scRNAsub$seurat_clusters)






pdf(file="./figures/day1分布-亚群.pdf", width=7, height=4)
DimPlot(scRNAsub, label=TRUE,repel = T,
        split.by = 'tissue_type',
        cols = col_vector,
        group.by = "seurat_clusters",
        reduction = "umap") 
dev.off()


table(scRNA$celltype)
table(scRNA$celltype, scRNA$tissue_type)
table(scRNA$celltype)



## 细胞比例图
library(reshape2)
library(ggplot2)
prop_df <- table(scRNAsub@meta.data$seurat_clusters,
                 scRNAsub@meta.data$tissue_type) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 


# library(RColorBrewer)
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))
# mycolors <- getPalette(25)
# 
# col_vector <- sample(getPalette(30))
# 
# sample_color <- col_vector[1:10] 
# 
# saveRDS(col_vector, file = "col_vector用于eqtl.RDS")

# 作图
prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:30]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
dev.off()
prop


pdf(file="./figures/day1分布比例图-亚群.pdf", width=15, height=4)
print(prop)
dev.off()



scRNAsub$celltype=ifelse(scRNAsub$seurat_clusters %in% c(0,3,5,6,9,10),
                         'HF_related_CD4_Tcm',
                         'Other_CD4_Tcm')




# 
# 
# 
# #
# # T细胞亚群注释参考
# # https://blog.csdn.net/weixin_52505487/article/details/126687526
# Idents(scRNAsub)=scRNAsub$seurat_clusters
# T_marker <- c("CCR7","LEF1", "TCF7",'SELL','KLF2', #CD4_Naive
#                 "ANXA1", "CXCR4", "IL2", #CD4_EM
#                 "BCL6", "CXCR5","ICA1", #CD4_FH
#                 'IL23R',"CCR6",'CAPG','RORC','IL17A', #TH17
#                 'FOXP3','IL2RA','IL1R2',#CD4_REG
#               'CD8A','GZMK',  # CD8_EM
#               'GZMA','CCL5',  #CD8_CM
#                'HAVCR2','PDCD1','LAG3', # CD8_exhau
#                 'EPCAM','CD19','CD3E')
# DotPlot(scRNAsub,
#         features = T_marker,
#         group.by = "seurat_clusters") + coord_flip()
# 
# 
# ## 不要求非常准确
# T_celltype=c('CD4_Naive',
#              'CD8_CM',
#              'CD8_EM',
#              'CD4_Naive',
#              'Th17',
#              'CD4_Naive',
#              'CD4_CM',
#              'CD4_Naive',
#              'CD8_EM',
#              'CD4_REG',
#              'CD8_Naive',
#              'CD8_Exhau',
#              'CD8_Exhau')
# 
# 
# 
# Idents(scRNAsub) <- scRNAsub@meta.data$seurat_clusters
# names(T_celltype) <- levels(scRNAsub)
# scRNAsub<- RenameIdents(scRNAsub, T_celltype)
# 
# scRNAsub@meta.data$T_celltype <- Idents(scRNAsub)
# 
# 
# #设置idents主要识别标签
# Idents(scRNAsub)=scRNAsub@meta.data$T_celltype
# 
# colors=c('#313c63','#b42e20','#ebc03e','#377b4c',
#          '#7bc7cd','#5d84a4','#bc3c29')
# DimPlot(scRNAsub, group.by="T_celltype", label=T, label.size=5,cols = col_vector[30:50],
#         pt.size = 1,split.by = 'tissue_type')
# 
# # ## 细胞比例图再次
# #
# # library(reshape2)
# # library(ggplot2)
# # prop_df <- table(scRNAsub@meta.data$T_celltype,scRNAsub@meta.data$tissue_type) %>% melt()
# # colnames(prop_df) <- c("Cluster","Sample","Number")
# # prop_df$Cluster <- factor(prop_df$Cluster)
# # library(RColorBrewer)
# # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# # #处理后有73种差异还比较明显的颜色，基本够用
# # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# #
# # sample_color <- col_vector[1:10]
# #
# # prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
# #   geom_bar(stat = "identity", width=0.8,position="fill")+
# #   scale_fill_manual(values=col_vector[1:20]) +
# #   theme_bw()+
# #   theme(panel.grid =element_blank()) +
# #   labs(x="",y="Ratio")+
# #   ####用来将y轴移动位置
# #   theme(axis.text.y = element_text(size=12, colour = "black"))+
# #   theme(axis.text.x = element_text(size=12, colour = "black"))+
# #   theme(
# #     axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
# #   )
# # prop


saveRDS(scRNAsub,file = 'scRNA_NKcells.RDS')
ent_blank()) +
# #   labs(x="",y="Ratio")+
# #   ####用来将y轴移动位置
# #   theme(axis.text.y = element_text(size=12, colour = "black"))+
# #   theme(axis.text.x = element_text(size=12, colour = "black"))+
# #   theme(
# #     axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
# #   )
# # prop


saveRDS(scRNAsub,file = 'scRNA_CD14_Mono.RDS')
