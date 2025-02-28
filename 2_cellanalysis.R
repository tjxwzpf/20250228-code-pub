# 关键细胞亚群的深入分析


## 读取数据
scRNA_T=readRDS('./scRNA_NKcells.RDS')
# BiocManager::install("slingshot")


library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)

## 细胞轨迹分析----------------------------

## 载入示例数据
# scRNA_T$T_celltype <- scRNA_T$celltype
scRNA_T$T_celltype <- scRNA_T$celltype
table(scRNA_T$T_celltype)

cellinfo <- scRNA_T@meta.data


## 构建SingleCellExperiment对象
scRNA_T

# scRNA_T1 <- CreateSeuratObject(counts = scRNA_T@assays$RNA@counts,
#                                meta.data = scRNA_T@meta.data)
sce <- as.SingleCellExperiment(scRNA_T,assay = "RNA")

ncol(scRNA_T)
nrow(scRNA_T)
scRNA_T
table(scrn)
## run
sce_slingshot <- slingshot(sce , clusterLabels = 'T_celltype', 
                           reducedDim = 'UMAP', 
                           # start.clus = c(3,5), shrink = 0.2)
                           start.clus = c("Other_CD4_Tcm"), shrink = 0.2)



#lin1 <- getLineages(sce_slingshot, 
#                    clusterLabels = "seurat_clusters", 
#                   start.clus ="0",#可指定起始细胞簇，用处不大
#                    end.clus="5",#可指定终点细胞簇,用处不大
#                    reducedDim = "UMAP")


sce_slingshot <- sim
如果上面失败了则跑下面：

# 在这里为了减少计算压力，我们把高变基因提取出来
# 你可以跟我一样从 scale.data 中提取，也可以直接在 seurat 对象中找出来
scRNA_T <- readRDS('./scRNA_CD14_Mono.RDS')
scobj <- scRNA_T
scale.data <- scobj@assays$RNA@scale.data
scale.gene <- rownames(scale.data)
length(scale.gene)

counts <- scobj@assays$RNA@counts
counts <- counts[scale.gene,]
dim(counts)
# 将表达矩阵转换为SingleCellExperiment对象
# 输入需要counts矩阵，否则影响下游分析
sim <- SingleCellExperiment(assays = List(counts = counts)) 


# 为了与前面画图一致，我们把降维坐标和细胞标签也提取并加入 sce 对象。
# umap reduction
umap = scobj@reductions$umap@cell.embeddings
colnames(umap) = c('UMAP-1', 'UMAP-2')
# 将降维的结果添加到SingleCellExperiment对象中
reducedDims(sim) = SimpleList(UMAP = umap)

# metadata
meta = scobj@meta.data
# colData(sim)相当于meta.data，但他不是data.frame格式
# 所以需要分步赋予
colData(sim)$sampleId = meta$sampleId
# colData(sim)$celltype = meta$celltype #这样的话，只有一条线
colData(sim)$celltype = meta$seurat_clusters

rd = umap
plot(rd, col = rgb(0,0,0,.5), pch=16, asp = 1)


table(sim$celltype)
# 一行命令就可生成，这里计算是比较快的

如果重新跑需要
rm(sim)

sim <- slingshot(sim, 
                 clusterLabels = 'celltype',  # 选择colData中细胞注释的列名
                 reducedDim = 'UMAP',  
                 # start.clus= "other_CD14_Mono",  # 选择起点
                 # start.clus= "0",  # 选择起点 seurat_clusters
                 start.clus= "4",  # 选择起点 seurat_clusters
                 end.clus = NULL     # 这里我们不指定终点
)     
colnames(colData(sim))


# 我们可以看到，sim 的 colData 中添加了一些时序信息，
# 每一个 lineage 添加了一列信息，每列信息以时序信息（拟时间）组成，
# 如果细胞对应不在这个 lineage 中，则以 NA 表示。

head(colData(sim)[,4:5])

colnames(colData(sim))
# 首先我们看看我们各个 lineage 是什么样的，我们使用官方提供的代码
pdf(file = "figures/day2_slingshot曲线.pdf",width = 4,height = 4)
plot(reducedDims(sim)$UMAP, 
     cex=.5, #cex就是点的大小
     pch=16, #pch就是点的形状
     # pch=16, 
     asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
legend("right",
       legend = paste0("lineage",1:6),
       col = unique(brewer.pal(6,"Set1")),
       # inset=0.8,
       pch = 16)
dev.off()



cellinfo <- meta
cl1 <- cellinfo$seurat_clusters
cl1

dev.off()


dev.off()
pdf(file = "figures/day2_slingshot折线.pdf",width = 5,height = 4)
plot(reducedDims(sce_slingshot)$UMAP,
     col = brewer.pal(12,"Paired")[cl1],
     cex=.5, #cex就是点的大小
     pch=16, #pch就是点的形状
     asp=1)

## 下面这行关键，否则容易报错！！（特别是Matrx>1.5-0的同学）
igraph::igraph.options(sparsematrices = FALSE)
## 曲线折线仍选一个即可
# lines(SlingshotDataSet(sce_slingshot), lwd=2,col = 'black')#,type = 'lineages'
lines(SlingshotDataSet(sce_slingshot), lwd=1, type = 'lineages', col = 'black')

# legend("right",legend = unique(sce$T_celltype),
legend("right",legend = levels(scRNA_T$seurat_clusters),
       col = unique(brewer.pal(12,"Paired")[cl1]),
       # inset=c(3,2,4), 
       pch = 16)

dev.off()




head(colData(sim)[,4:5])
colnames(colData(sim))
我们现在对 lineage3 非常感兴趣

# 这些 lineage 大致都符合我们先前注释的亚群背后的生物学知识，
# 也有轨迹出现部分重叠，可以设置 slingshot 参数 reweight = FALSE ，
# 让细胞不会在多个轨迹间重复评估。
# 接下来，假设我们现在对 lineage3 非常感兴趣，我们希望把 lineage3 进行可视化。
# 首先我们为细胞赋予颜色：
summary(sim$slingPseudotime_3)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   2.207   4.050   6.339  10.041  18.065    9446 

# 我们做一个绚丽的渐变色彩虹色
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) # 我们把这些颜色变成100个梯度，模拟渐变色
plotcol <- colors[cut(sim$slingPseudotime_3, breaks=100)] # 这里我们用cut函数把 lineage2 分割成100个区间，同一区间的细胞给与同一个颜色
plotcol[is.na(plotcol)] <- "lightgrey" # 不属于 lineage3 的点为NA，我们把他们变成灰色
plotcol

来看看效果：
（但是这些 lineage 数据怎么用 ggplot2 系统绘画，笔者还没找到方法）
# 保存时，把下面的注释运行
pdf(file = "figures/day2_slingshot_lineage3.pdf",width = 7,height = 4)
plot(reducedDims(sim)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
legend("right",
       legend = paste0("lineage",1:6),
       col = unique(brewer.pal(6,"Set1")),
       inset=0.8,
       pch = 16)
dev.off()










# 细胞通讯----------------------
# 在SLE和衰老中分别做，看看功能是不是一样的

## cell-cell chat----
scRNA=readRDS('./scRNA_anno注释azimuth.RDS')
scRNA
table(scRNA$predicted.celltype.l2)
table(scRNA$celltype)

scRNA_other=subset(scRNA,celltype %in% c( 'CD14 Mono'),invert=T)

rm(scRNA)

# 加一列！！！
scRNA_other$T_celltype =scRNA_other$celltype
table(scRNA_other$T_celltype)
table(scRNA_other$T_celltype)
table(scRNA_T$T_celltype)

DimPlot(scRNA_T, split.by = "tissue_type",group.by = "T_celltype")
scRNA_chat=merge(scRNA_other,scRNA_T)
rm(scRNA_T)
rm(scRNA_other)
gc()



# scRNA_chat <- scRNA
table(scRNA_chat$tissue_type)
table(scRNA_chat$T_celltype, scRNA_chat$tissue_type)
scRNA_chat_old=subset(scRNA_chat,tissue_type=='sepsis')



scRNA
table(scRNA_chat_old$orig.ident)


scRNA_chat_old
## 抽2000细胞做
set.seed(123)
a=sample(1:ncol(scRNA_chat_old),4000)
scRNA_chat_old=scRNA_chat_old[,a]

meta =scRNA_chat_old@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_old@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

table(meta$T_celltype)
library(CellChat)

cellchat <- createCellChat(object = data_input, meta = meta, 
                           group.by = "T_celltype")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)




cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
# cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- filterCommunication(cellchat, min.cells = 40)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))


saveRDS(cellchat, file = './figures/day2_cellchat.RDS')



library(CellChat)


cellchat <- readRDS('./figures/day2_cellchat.RDS')
##时常deff.off!!!!
dev.off()

p1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                       weight.scale = T, label.edge= T, 
                       # sources.use = 'HF_related_CD4_Tcm',
                       title.name = "Number of interactions")
p1




dev.off()

table(cellchat@idents)
p_bubble= netVisual_bubble(cellchat,
                           sources.use = c("sepsis_related_CD14_Mono",
                                           "other_CD14_Mono"),
                           remove.isolate = FALSE)+coord_flip()
p_bubble



table(scRNA_chat$tissue_type)
pdf(file="./figures/day2_cellchat_CD14Mono.pdf", width=7, height=5)
print(p_bubble)
dev.off()


saveRDS(cellchat, file = './figures/day2_cellchat.RDS')

# ## SLE
# scRNA_chat_SLE=subset(scRNA_chat,tissue_type=='HF')
# 
# set.seed(123)
# a=sample(1:ncol(scRNA_chat_SLE),2000)
# scRNA_chat_SLE=scRNA_chat_SLE[,a]
# 
# meta =scRNA_chat_SLE@meta.data # a dataframe with rownames containing cell mata data
# gc()
# data_input <- as.matrix(scRNA_chat_SLE@assays$RNA@data)
# #data_input=data_input[,rownames(meta)]
# identical(colnames(data_input),rownames(meta))
# 
# library(CellChat)
# cellchat <- createCellChat(object = data_input, meta = meta, 
#                            group.by = "celltype")
# 
# CellChatDB <- CellChatDB.human 
# groupSize <- as.numeric(table(cellchat@idents))
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# cellchat@DB <- CellChatDB.use 
# 
# dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# # 提取数据库支持的数据子集
# cellchat <- subsetData(cellchat)
# # 识别过表达基因
# cellchat <- identifyOverExpressedGenes(cellchat)
# # 识别配体-受体对
# cellchat <- identifyOverExpressedInteractions(cellchat)
# # 将配体、受体投射到PPI网络
# cellchat <- projectData(cellchat, PPI.human)
# unique(cellchat@idents)
# 
# cellchat <- computeCommunProb(cellchat)
# 
# # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
# cellchat <- filterCommunication(cellchat, min.cells = 10)
# cellchat <- computeCommunProbPathway(cellchat)
# 
# df.net<- subsetCommunication(cellchat)
# 
# cellchat <- aggregateNet(cellchat)
# groupSize <- as.numeric(table(cellchat@idents))
# 
# ##时常deff.off!!!!
# dev.off()
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
#                  weight.scale = T, label.edge= T, sources.use = 'NK',
#                  title.name = "Number of interactions")
# dev.off()
# 
# pp_bubble= netVisual_bubble(cellchat,
#                            sources.use = 'NK',
#                            remove.isolate = FALSE)+coord_flip()
# pp_bubble

# pdf(file="./figures/day2_cellchat_HF.pdf", width=15, height=4)
# print(pp_bubble)
# dev.off()






