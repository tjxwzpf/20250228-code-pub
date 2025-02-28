# https://gwas.mrcieu.ac.uk/datasets/?trait__icontains=pqtl

## 关键基因的孟德尔随机化分析

library(Seurat)
# 读取数据
scRNA=readRDS('./scRNA_anno注释azimuth.RDS')
table(Idents(scRNA))
scRNA
table(scRNA$tissue_type)
table(scRNA$orig.ident)
table(scRNA$predicted.celltype.l2)
# scRNA$celltype <- scRNA$predicted.celltype.l2
# Idents(scRNA) <- "celltype"
# table(Idents(scRNA))
# saveRDS(scRNA, file = "scRNA_anno注释azimuth.RDS")

scRNA=readRDS('./scRNA_anno注释azimuth.RDS')
library(Seurat)
genes <- readRDS('./hd_day26_单细胞多组差异基因logFC热图展示/GO_COAGULATION.RDS')
genes <- list(genes)
names(genes) <- "GO_COAGULATION"
scRNA <- AddModuleScore(scRNA,features = genes,name =names(genes) 
)
colnames(scRNA@meta.data)[ncol(scRNA@meta.data)] <- names(genes)


dev.off()

table(scRNA$orig.ident)
scRNA$orig.ident <- factor(scRNA$orig.ident, levels = c("before_1",
                                                        "before_2",
                                                        "after_1",
                                                        "after_2"))
library(ggpubr)
col_vector <- readRDS('col_vector用于eqtl.RDS')
p1 <- VlnPlot(scRNA, features = names(genes),
              # group.by = "tissue_type",
              group.by = "orig.ident",
              cols = col_vector,
              pt.size = 0)+
  geom_boxplot(width=.1,col="black",fill="white", outlier.alpha = 0)+
  stat_compare_means(label.x.npc = 0.1)+NoLegend()+labs(subtitle = "All cells")
p2 <- VlnPlot(scRNA, features = names(genes),
              group.by = "tissue_type",
              # group.by = "orig.ident",
              cols = col_vector,
              pt.size = 0)+
  geom_boxplot(width=.1,col="black",fill="white", outlier.alpha = 0)+
  stat_compare_means(label.x.npc = 0.1)+NoLegend()+labs(subtitle = "All cells")



p3 <- VlnPlot(scRNA, features = names(genes),
              group.by = "celltype",
              split.by = "tissue_type",
              # group.by = "orig.ident",
              cols = col_vector,
              pt.size = 0)+
  # geom_boxplot(width=.1,col="black",fill="white", outlier.alpha = 0)+
  stat_compare_means(label.x.npc = 0.1)+NoLegend()

p3

p1|p2


scRNAsub=readRDS('./scRNA_NKcells.RDS')
# 
# library(Seurat)
# genes <- readRDS('./hd_day26_单细胞多组差异基因logFC热图展示/GO_COAGULATION.RDS')
# genes <- list(genes)
# names(genes) <- "GO_COAGULATION"
# scRNAsub <- AddModuleScore(scRNAsub,features = genes,name =names(genes) 
#                             )
# colnames(scRNAsub@meta.data)[ncol(scRNAsub@meta.data)] <- names(genes)
# 

dev.off()
scRNAsub <- subset(scRNA, celltype %in% c("NK"))
table(scRNAsub$orig.ident)
scRNAsub$orig.ident <- factor(scRNAsub$orig.ident, levels = c("before_1",
                                                              "before_2",
                                                              "after_1",
                                                              "after_2"))
library(ggpubr)
col_vector <- readRDS('col_vector用于eqtl.RDS')
p4 <- VlnPlot(scRNAsub, features = names(genes),
        # group.by = "tissue_type",
        group.by = "orig.ident",
        cols = col_vector[7:10],
        pt.size = 0)+
  geom_boxplot(width=.1,col="black",fill="white", outlier.alpha = 0)+
  stat_compare_means(label.x.npc = 0.1)+NoLegend()+labs(subtitle = "NK cells")
p5 <- VlnPlot(scRNAsub, features = names(genes),
              group.by = "tissue_type",
              # group.by = "orig.ident",
              cols = col_vector[7:10],
              pt.size = 0)+
  geom_boxplot(width=.1,col="black",fill="white", outlier.alpha = 0)+
  stat_compare_means(label.x.npc = 0.1)+NoLegend()+labs(subtitle = "NK cells")


(p1|p2)/(p4|p5)

library(patchwork) 
pdf(file="./figures/拼图.pdf", width=7, height=7)
print((p1|p2)/(p4|p5))
dev.off()







# scRNA_other
# scRNA_other=subset(scRNA,celltype !='NK')
# # scRNA_other=subset(scRNA,celltype %in% c( 'CD14 Mono'),invert=T)
# 
# scRNAsub
# table(scRNAsub$celltype)
# scRNAsub$T_celltype <- scRNAsub$celltype
# table(scRNAsub$T_celltype)
# 
# scRNA_CM=subset(scRNAsub,T_celltype == 'sepsis_related_CD14_Mono')
# scRNA_CM
# gc()
# scRNA_other$T_celltype=scRNA_other$celltype
# table(scRNA_other$T_celltype)
# 
# scRNA_compare=merge(scRNA_other,scRNA_CM)
# table(scRNA_compare$T_celltype)

gc()
# rm()

scRNAsub



Idents(scRNAsub) <- "tissue_type"
table(scRNAsub$tissue_type)
## 1. CM相对其他亚群
df_CM=FindMarkers(scRNAsub,ident.1 = 'after',
                  # only.pos = T,
                  logfc.threshold = 0.1,
                  only.pos = T)


saveRDS(df_CM, file = "day3_细胞内部比的差异基因.RDS")

# Idents(scRNA_compare) <- "celltype"
# 
# table(scRNA_compare$T_celltype)
# table(scRNA_compare$celltype)
# ## 2. CM相对其他细胞
# df_T=FindMarkers(scRNA_compare,ident.1 = 'sepsis_related_CD14_Mono',
#                  only.pos = T,logfc.threshold = 0.25)
# 
# ss=intersect(rownames(df_CM),rownames(df_T))
# 
# ss
# saveRDS(df_T, file = "day3_细胞比其他所有类型的差异基因.RDS")

ss <- rownames(df_CM)
ss <- unique(ss)
ss

##########################################   上面失败了，继续用COSG            ############################              
############################ ##########################################

##########################################   COSG and 富集分析            ############################
############################ ##########################################
library(COSG)

pbmc <- readRDS("scRNA_NKcells.RDS")
# pbmc <- scRNA_compare
# pbmc <- scRNAsub
DimPlot(pbmc, split.by = "tissue_type")
# table(pbmc$T_celltype)
# 
# pbmc$T_celltype=ifelse(pbmc$T_celltype %in% c("HF_related_CD4_Tcm"),
#                        'HF_related_CD4_Tcm',
#                        'Other_celltypes')
table(pbmc$tissue_type)
Idents(pbmc) <- "tissue_type"
# # Run COSG:
pro = 'NKcells'

if(T){
  
  library(COSG)
  marker_cosg <- cosg(
    pbmc,
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=100)
  dir.create("./figures")
  save(marker_cosg,file = paste0('./figures/',pro,'_marker_cosg.Rdata'))
  head(marker_cosg)
  
  ## Top10 genes
  library(dplyr)
  top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10)))
  ## Top3 genes
  top_3 <- unique(as.character(apply(marker_cosg$names,2,head,3)))
  
}


library(dplyr)
symbols_list <-  as.list(as.data.frame(apply(marker_cosg$names,2,head,100)))
symbols_list_1 <- symbols_list$after
symbols_list_1
# symbols_list1 <- as.data.frame(unlist(symbols_list))
#
# # symbol <- data.frame()
# names(symbols_list)
# symbol <- c()
# for (i in names(symbols_list)) {
#   symbol1 <- rep(i,50)
#   symbol <- c(symbol,symbol1)
# }
# symbol
# symbols_list1$cluster <- symbol
# sig.markers <- symbols_list1


##########################################  接在day3mendelian后面             ############################
############################ ##########################################

pbmc <- scRNAsub

table(pbmc$T_celltype)

# pbmc$T_celltype=ifelse(pbmc$T_celltype %in% c("HF_related_CD4_Tcm"),
#                        'HF_related_CD4_Tcm',
#                        'Other_celltypes')
table(pbmc$T_celltype)
Idents(pbmc) <- "T_celltype"
# # Run COSG:
pro = 'CD4Tcm_比Other_CD4_Tcm_cosg_seurat_clusters'

if(T){
  
  library(COSG)
  marker_cosg <- cosg(
    pbmc,
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=200)
  dir.create("./figures")
  save(marker_cosg,file = paste0('./figures/',pro,'_marker_cosg.Rdata'))
  head(marker_cosg)
  
  ## Top10 genes
  library(dplyr)
  top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10)))
  ## Top3 genes
  top_3 <- unique(as.character(apply(marker_cosg$names,2,head,3)))
  
}


library(dplyr)
symbols_list <-  as.list(as.data.frame(apply(marker_cosg$names,2,head,200)))
symbols_list_2 <- symbols_list$after
symbols_list_2


# ss=intersect((symbols_list_1),(symbols_list_2))
# 
# ss
# 
# 
# ss <- rownames(df_CM)
# ss <- unique(ss)
# ss

ss <- symbols_list_2
ss

genes
gene <- genes$GO_COAGULATION

ss <- intersect(ss, gene)
ss



## 保存，后面要用
save(ss,file ='key_marker_gene.Rdata')
save(ss,file ='./figures/key_marker_gene.Rdata')

load('key_marker_gene.Rdata')
write.csv(ss,file ='./figures/key_marker_gene.csv',quote = F)


gene=ss
gene=genes$GO_COAGULATION
gene=gene[!duplicated(gene)]
gene

## 转ENSEMBL
library(clusterProfiler)
library(org.Hs.eg.db)
df=bitr(geneID = gene,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
id=df$ENSEMBL
id=id[!duplicated(id)]
id



## 载入包
# devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
# https://gwas.mrcieu.ac.uk/datasets/?trait__icontains=pqtl


# 暴露数据的SNP
exposure_id=paste0('eqtl-a-',id)
exposure_dat <- extract_instruments(exposure_id, p1=5e-08, 
                                    clump=TRUE)

exposure_id
# exposure_dat <- extract_instruments(exposure_id, p1=5e-06, 
#                                     clump=TRUE)

#报错可以试试这个
# try(tmp <- extract_instruments(ID,p1 = 5e-08,clump = F,access_token=NULL))

save(exposure_dat,file ='exposure_dat.Rdata')
load('exposure_dat.Rdata')
R2a=2*exposure_dat$beta.exposure*exposure_dat$beta.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2b=2*exposure_dat$se.exposure*exposure_dat$se.exposure*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2=R2a/(R2a+R2b)

exposure_dat$F_statistics=R2*(exposure_dat$samplesize.exposure-2)/(1-R2)

exposure_dat=exposure_dat[exposure_dat$F_statistics>10,]

exposure_dat$exposure
exposure_dat$exposure=stringr::str_sub(exposure_dat$exposure,1,15)
exposure_dat$exposure <- gsub(" || id:eqtl-a-","",exposure_dat$exposure)
exposure_dat$exposure <- gsub("\\||","",exposure_dat$exposure)
exposure_dat$exposure

exposure_dat$ENSEMBL =exposure_dat$exposure
exposure_dat$ENSEMBL
df
exposure_dat= left_join(exposure_dat,df, by = "ENSEMBL")

exposure_dat$exposure=NULL
exposure_dat$exposure=exposure_dat$SYMBOL
exposure_dat$SYMBOL=NULL
save(exposure_dat,file ='exposure_dat.Rdata')
save(exposure_dat,file ='./figures/exposure_dat.Rdata')

load('exposure_dat.Rdata')

unique(exposure_dat$SNP)
#这些基因的SNP就是我们的工具变量



#提取结局数据
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, 
                                    # outcomes="ebi-a-GCST90011866")
                                    outcomes="ebi-a-GCST009541") #这个可能还行
# outcomes="ukb-d-HEARTFAIL")
# outcomes="bbj-a-109")
# outcomes="uk <- b-d-I9_HEARTFAIL")

# outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, 
#                                     # outcomes="ebi-a-GCST90011866")
#                                     outcomes=c("ebi-a-GCST009541",
#                                                "ukb-d-I50",
#                                                "bbj-a-109",
#                                                "ukb-e-I50_CSA",
#                                                "ukb-e-428_CSA")) #这个可能还行

out_id='ukb-d-I50'#结局的IEU id     Diagnoses - main ICD10: I50 Heart failure
out_id='ebi-a-GCST009541'#结局的IEU id     Diagnoses - main ICD10: I50 Heart failure
out_id='ukb-d-HEARTFAIL'#结局的IEU id     Diagnoses - main ICD10: I50 Heart failure
out_id='ukb-d-I9_HEARTFAIL'#结局的IEU id     Diagnoses - main ICD10: I50 Heart failure
out_id='bbj-a-109'#结局的IEU id     Diagnoses - main ICD10: I50 Heart failure

# setwd('../5个GWAS数据集/')
dir.create("./figures")
saveRDS(outcome_dat, file = "./figures/outcome_dat.RDS")
outcome_dat <- readRDS('./figures/outcome_dat.RDS')
# 取交集, 筛选共有的SNP
exposure_dat=exposure_dat[exposure_dat$SNP %in% outcome_dat$SNP,]

#其实就是合并矩阵
harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)



## MR分析

mr_modified <- function(dat = harmonised_dat, prop_var_explained = T)
{
  mr_res <- mr(dat) #mr函数是核心
  #pve就是可解释的比例，可解释比例越高，暴露越重
  pve <- dat %>% 
    dplyr::select(id.exposure, beta.exposure, 
                  se.exposure, samplesize.exposure) %>% 
    dplyr::group_by(id.exposure) %>% 
    dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  
  if(prop_var_explained)
  {
    mr_res <- mr_res %>% 
      dplyr::left_join(pve, by = "id.exposure")
  }
  
  return(mr_res)
}



#运行mr
# Analysing 'eqtl-a-ENSG00000126353' on 'ukb-d-I50'即：分析这个暴露在这个结局上的情况，
# 一个个暴露的分析

mr_res <- mr_modified(harmonised_dat, prop_var_explained = T)


这个表格是最重要的，先看看pvalue，如果没有小于0.05的，放弃重来


# save(exposure_dat,outcome_dat,mr_res,harmonised_dat,file ='mr_input_res.Rdata')

# GWAS数据其实就是分析SNP
save(exposure_dat,outcome_dat,mr_res,harmonised_dat,file ='./figures/mr_input_res.Rdata')

load('./figures/mr_input_res.Rdata')

colnames(mr_res)

#一般有五种计算方法
unique(mr_res$method)

# 从这个数据可以看出来，其实mr的结果是不体现有多少个snp的，也就是不体现工具变量的，
# 只展示暴露和结局的关系


#  异质性 ，用处不大
# 异质性分析：看看工具变量之间差别大不大，如果差别太大说明结果可能不太可靠；
# 可以不汇报这个结果
# 对疾病的gwas数据可能有意义，对eqtl数据可能比较没用
# 因为snp比较少


heter_dat=mr_heterogeneity(harmonised_dat)
write.csv(heter_dat, file="table.heterogeneity.csv", row.names=F)
write.csv(heter_dat, file="./figures/table.heterogeneity.csv", row.names=F)

colnames(heter_dat)
unique(heter_dat$method)


# 水平多效性检验，用处不大
#水平多效性检验：这个snp可能不是先导致基因表达改变，才导致心衰的发生，可能本身就跟疾病有关系；
# 所以要检验看看
pleio_dat=mr_pleiotropy_test(harmonised_dat)
write.csv(pleio_dat, file="table.pleiotropy.csv", row.names=F)
write.csv(pleio_dat, file="./figures/table.pleiotropy.csv", row.names=F)

# 前两个都可以不看


最重要的是mr_res这个数据，主要是pvalue和beta这两列
p<0.05认为有价值
beta认为大于0是危险因素，是和心衰发生有直接关系的基因



### table
table1 <- mr_res %>% 
  filter(pval < 0.05,
         method %in% c("Wald ratio","Inverse variance weighted")) %>% 
  left_join(exposure_dat, by = "exposure")
saveRDS(table1,file ='./figures/mr.res.sig.RDS')

table1 <- table1 %>% 
  generate_odds_ratios()%>% 
  mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",or,or_lci95, or_uci95),
         `P value` = scales::scientific(pval),
         `PVE` = paste0(sprintf("%.2f",100*pve),"%"),
         `F statistics` = sprintf("%.2f",F_statistics)) %>% 
  dplyr::select(Gene = exposure, `ENSEMBL ID` = ENSEMBL,
                SNP, `Effect allele` = effect_allele.exposure, 
                `OR (95% CI)`, `P value`, 
                PVE, `F statistics`)

## 保存，后面常用
save(table1,file ='table1.Rdata')
save(table1,file ='./figures/table1.Rdata')



# ### 可视化
# volcano_plot <- function(.data, 
#                          number_comparasion = 1,
#                          title = "eQTL",
#                          col_beta = "b",
#                          col_size = "pve",
#                          col_label = "exposure",
#                          legend.position = "none")
# {
#   p_thershold <- 0.05/number_comparasion
#   
#   p <- .data %>% 
#     rename(beta := !!col_beta,
#            size := !!col_size,
#            label := !!col_label) %>% 
#     mutate(x = beta,
#            y = -log10(pval),
#            label = ifelse(pval < p_thershold, label, NA)) %>% 
#     ggplot(aes(x = x, y = y)) +
#     geom_point(aes(size = size), alpha = 0.5, color = "#0072b5") +
#     geom_vline(xintercept = 0, linetype = 2)+
#     geom_hline(yintercept = -log10(p_thershold), linetype = 2) +
#     theme_classic() +
#     theme(panel.grid = element_blank(),
#           legend.title = element_text(size = 6.5),
#           legend.text = element_text(size = 6.5),
#           legend.position = legend.position)+
#     labs(x = "ln(OR)", 
#          y = parse(text = "-log[10]*(italic(P)-value)"),
#          title = title) +
#     scale_size(name = "PVE",
#                breaks = c(0.2*1:3)) +
#     ggrepel::geom_label_repel(aes(label = label),size = 3)
#   plot(p)
# }
# 
# library(ggplot2)
# 
#   mr_res %>% 
#     filter(method %in% c("Wald ratio","Inverse variance weighted")) %>% 
#     volcano_plot(number_comparasion = 1)


### 如果包有冲突可以用下面
library(dplyr)
library(ggplot2)
library(ggrepel)

volcano_plot <- function(.data,
                         number_comparasion = 1,
                         title = "eQTL",
                         legend.position = "none") {
  
  p_thershold <- 0.05/number_comparasion
  
  p <- .data %>%
    mutate(y = -log10(pval),
           label = ifelse(pval < p_thershold, exposure, NA)) %>%
    ggplot(aes(x = b, y = y)) +
    # geom_point(aes(size = pve), alpha = 0.5, color = "#0072b5") +
    geom_point(aes(size = pve), alpha = 0.5, color = "#ba9839") +
    # geom_vline(xintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 3) +
    # geom_hline(yintercept = -log10(p_thershold), linetype = 2) +
    geom_hline(yintercept = -log10(p_thershold), linetype = 3) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6.5),
          legend.position = legend.position) +
    labs(x = "ln(OR)",
         y = parse(text = "-log[10]*(italic(P)-value)"),
         title = title) +
    scale_size(name = "PVE",
               breaks = c(0.2*1:3)) +
    ggrepel::geom_label_repel(aes(label = label), size = 1.5)  #字体大小
  plot(p)
}

# 画图
mr_res %>%
  filter(method %in% c("Wald ratio", "Inverse variance weighted")) %>%
  volcano_plot(number_comparasion = 1)



p_vol <- mr_res %>%
  filter(method %in% c("Wald ratio", "Inverse variance weighted")) %>%
  volcano_plot(number_comparasion = 1)
pdf(file="./figures/day3火山图.pdf", width=3.5, height=3.5)

p_vol
dev.off()


VlnPlot(scRNAsub, features = "AP3B1", group.by = "tissue_type")
VlnPlot(scRNAsub, features = "HLA-B", group.by = "tissue_type")


#将统计结果绘制森林图
mr_res2=mr_res[mr_res$exposure %in% table1$Gene,]
result_or=generate_odds_ratios(mr_res2)
write.table(result_or[,4:ncol(result_or)],"OR.txt",row.names = F,sep = "\t",quote = F)
write.table(result_or[,4:ncol(result_or)],"./figures/OR.txt",row.names = F,sep = "\t",quote = F)


library(grid)
library(forestploter)


mydata=read.table("OR.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
forest(mydata[,c(1:3,6,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)

pdf(file="./figures/day3_森林图.pdf", width=10, height=10,onefile = T)
forest(mydata[,c(1:3,6,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)
dev.off()



##########################################   单个基因的漏斗图森林图            ############################               
############################ ##########################################
## 其他可视化（一般不太好用，因为一个基因只有很有限的eqtl！！
# 取一个有多个工具变量的基因
mr_res1=mr_res[mr_res$exposure=='SNHG5',]
harmonised_dat1=harmonised_dat[harmonised_dat$exposure=='SNHG5',]

#绘制散点图
p_scatter <- mr_scatter_plot(mr_res1, harmonised_dat1)
p_scatter

#森林图
res_single=mr_singlesnp(harmonised_dat1)
p_forest <- mr_forest_plot(res_single)
p_forest

#漏斗图
p_funnel <- mr_funnel_plot(singlesnp_results = res_single)
p_funnel

#留一法敏感性分析
p_leave <- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harmonised_dat1))
p_leave



library(patchwork) 
pdf(file="./figures/拼图.pdf", width=10, height=8)
print(p_vol+p_scatter+p_forest+p_funnel+p_leave)
dev.off()


library(patchwork) 
pdf(file="./figures/拼图2.pdf", width=8, height=10)
print((p_vol+p_vol+p_forest+p_leave+p_scatter+p_funnel)+plot_layout(ncol = 2))
dev.off()

genes <- unique(table1$Gene)
genes

saveRDS(genes, file = "genes_mr_sig.RDS")

for (i in genes) {
  mr_res1=mr_res[mr_res$exposure==i,]
  harmonised_dat1=harmonised_dat[harmonised_dat$exposure==i,]
  
  #绘制散点图
  p_scatter <- mr_scatter_plot(mr_res1, harmonised_dat1)
  p_scatter
  
  #森林图
  res_single=mr_singlesnp(harmonised_dat1)
  p_forest <- mr_forest_plot(res_single)
  p_forest
  
  #漏斗图
  p_funnel <- mr_funnel_plot(singlesnp_results = res_single)
  p_funnel
  
  #留一法敏感性分析
  p_leave <- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harmonised_dat1))
  p_leave
  
  
  #必须要有p_vol, 不然绘图失败
  library(patchwork) 
  pdf(file=paste0("./figures/单基因孟德尔拼图_",i,".pdf"), width=8, height=10)
  print((p_vol+p_vol+p_forest+p_leave+p_scatter+p_funnel)+plot_layout(ncol = 2))
  dev.off()
}


pbmc <- readRDS('./scRNA_CD4_Tcm.RDS')

Idents(pbmc) <- "celltype"
DotPlot(pbmc, features = genes)
library(ggpubr)
VlnPlot(pbmc, features = genes, ncol = 4)&stat_compare_means(label.x.npc = 0.1)

pdf(file="./figures/day3_scRNA_genes.pdf", width=10, height=4)
VlnPlot(pbmc, features = genes, ncol = 4)&stat_compare_means(label.x.npc = 0.1)
dev.off()

pdf(file="./figures/day3_scRNA_genes放大.pdf", width=13, height=5)
VlnPlot(pbmc, features = genes, ncol = 4)&stat_compare_means(label.x.npc = 0.1)
dev.off()
