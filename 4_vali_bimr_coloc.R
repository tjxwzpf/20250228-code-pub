## 孟德尔后：验证集、双向孟德尔、共定位分析

library(tidyverse)

load('./day3_sc_mr_代码_xc所有差异表达基因做mr/figures/exposure_dat.Rdata')

gene=read.table('./day3_sc_mr_代码_xc所有差异表达基因做mr/figures/OR.txt',sep = '\t',header = T)
# gene=gene$exposure
gene=unique(gene$exposure)
gene

exposure_dat=exposure_dat[exposure_dat$exposure %in% gene,]

### 验证集----
library(TwoSampleMR)
#提取结局数据,换一个队列
exposure_dat$SNP

#训练用的是ukb-d-I50
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, 
                                    # outcomes="ukb-d-I50") #换一个队列
                                    # outcomes=c("bbj-a-109",
                                               # 'ukb-d-HEARTFAIL')) #换一个队列
                                    # outcomes="bbj-a-109") #换一个队列
                                    outcomes=c('ukb-d-HEARTFAIL')) #这个是成功的
                                    # outcomes=c("ukb-d-I9_HEARTFAIL_NS","finn-b-I9_HEARTFAIL_EXNONE")) #换一个队列PRKCQ-AS1
                                    # outcomes=c("ebi-a-GCST009541")) #换一个队列PRKCQ-AS1
                                    # outcomes="bbj-a-109") #换一个队列 失败
                                    # outcomes="ukb-e-428_CSA") #换一个队列
                                    # outcomes="finn-b-I9_HEARTFAIL_EXNONE") #换一个队列NSUN5P1
                                    # outcomes="finn-b-I9_HEARTFAIL_ALLCAUSE") #换一个队列
                                    # outcomes="finn-b-I9_HEARTFAIL_AND_HYPERTCARDIOM") #换一个队列
outcome_dat

out_id='ukb-d-I50'#结局的IEU id     Diagnoses - main ICD10: I50 Heart failure
out_id='ebi-a-GCST009541'#结局的IEU id     Diagnoses - main ICD10: I50 Heart failure
out_id='ukb-d-HEARTFAIL'#结局的IEU id     Diagnoses - main ICD10: I50 Heart failure
out_id='ukb-d-I9_HEARTFAIL'#结局的IEU id     Diagnoses - main ICD10: I50 Heart failure
out_id='bbj-a-109'#结局的IEU id     Diagnoses - main ICD10: I50 Heart failure


{
# 取交集
exposure_dat=exposure_dat[exposure_dat$SNP %in% outcome_dat$SNP,]
exposure_dat$exposure

harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)



## MR

mr_modified <- function(dat = harmonised_dat, prop_var_explained = T)
{
  mr_res <- mr(dat)
  
  pve <- dat %>% 
    dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>% 
    dplyr::group_by(id.exposure) %>% 
    dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  
  if(prop_var_explained)
  {
    mr_res <- mr_res %>% 
      dplyr::left_join(pve, by = "id.exposure")
  }
  
  return(mr_res)
}

}


mr_res_vali <- mr_modified(harmonised_dat, prop_var_explained = T)
View(mr_res_vali)
# save(mr_res_vali,harmonised_dat,file ='mr_input_res_vali.Rdata')
save(mr_res_vali,harmonised_dat,file ='./figures/mr_input_res_vali.Rdata')

load('./figures/mr_input_res_vali.Rdata')
result_or=generate_odds_ratios(mr_res_vali)
write.table(result_or[,4:ncol(result_or)],"./figures/OR_vali.txt",row.names = F,sep = "\t",quote = F)

#将统计结果绘制森林图


library(grid)
library(forestploter)


mydata=read.table("./figures/OR_vali.txt",header = T,sep = "\t")
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



p_forest <- forest(mydata[,c(1:3,6,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)
pdf(file="./figures/day4_验证-森林图.pdf", width=10, height=9)
p_forest
dev.off()


library(TwoSampleMR)
unique(mr_res_vali$exposure)
## 其他可视化（一般不太好用，因为一个基因只有很有限的eqtl
mr_res1=mr_res_vali[mr_res_vali$exposure=='BAX',]
harmonised_dat1=harmonised_dat[harmonised_dat$exposure=='BAX',]

dev.off()
#绘制散点图
mr_scatter_plot(mr_res1, harmonised_dat1)

#森林图
res_single=mr_singlesnp(harmonised_dat1)
mr_forest_plot(res_single)

#漏斗图
mr_funnel_plot(singlesnp_results = res_single)

#留一法敏感性分析
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harmonised_dat1))


##########################################    反向孟德尔           ############################               
############################ ##########################################
### 反向孟德尔
# 暴露
bimr_disea <- extract_instruments('ebi-a-GCST009541',
# bimr_disea <- extract_instruments('ukb-d-I50',
# bimr_disea <- extract_instruments('ukb-d-HEARTFAIL',
                                  p1=5e-08, clump=TRUE) #结果太少了
                                  # p1=5e-06, clump=TRUE)
saveRDS(bimr_disea,file = "./figures/day4_bimr_disea.RDS")

bimr_disea <- readRDS('./figures/day4_bimr_disea.RDS')

bimr_disea$SNP
bimr_SLE <- bimr_disea
#提取结局数据
colnames(exposure_dat) #在这里看看基因的ensg
unique(exposure_dat$exposure)
bimr_SLE$SNP
rm(outcome_gene)
outcome_gene<- extract_outcome_data(snps=bimr_SLE$SNP, 
                                    # outcomes="eqtl-a-ENSG00000203875") #SNHG5
                                    # outcomes="eqtl-a-ENSG00000120738") #EGR1
                                    outcomes="eqtl-a-ENSG00000087088") #BAX
                                    # outcomes="eqtl-a-ENSG00000170345") #FOS
outcome_gene
# saveRDS(outcome_gene, file = "./figures/outcome_gene_EGR1.RDS")
saveRDS(outcome_gene, file = "./figures/outcome_gene_BAX.RDS")

# 取交集
outcome_gene$SNP
bimr_SLE =bimr_SLE [bimr_SLE $SNP %in% outcome_gene$SNP,]


harmonised_SLE_gene <- harmonise_data(bimr_SLE, outcome_gene)

bimr_mr_SLE_gene <- mr(harmonised_SLE_gene)

result_or=generate_odds_ratios(bimr_mr_SLE_gene)
write.table(result_or[,4:ncol(result_or)],"./figures/bi_OR.txt",row.names = F,sep = "\t",quote = F)

#将统计结果绘制森林图


library(grid)
library(forestploter)


mydata=read.table("./figures/bi_OR.txt",header = T,sep = "\t")
## !!
mydata$outcome
# mydata$outcome='SNHG5'
# mydata$outcome='EGR1'
mydata$outcome='BAX'
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
forest(mydata[,c(1:3,6,12,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =6 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)

dev.off()

p_forest <- forest(mydata[,c(1:3,6,12,13,14)],
                   est = mydata$or,
                   lower =mydata$or_lci95, 
                   upper = mydata$or_uci95,
                   sizes =0.3,
                   ci_column =6 ,
                   ref_line = 1,
                   xlim = c(0.05, 3),
)
pdf(file="./figures/day4_双向孟德尔_森林图.pdf", width=12, height=4)
# print(p_top)
p_forest
dev.off()




########################################## 共定位分析 两种方法         ############################
############################ ##########################################
### 共定位分析----
# rm(list=ls())
# options(stringsAsFactors = F)



# load('mr_input_res.Rdata')
load('./day3_sc_mr_代码_xc所有差异表达基因做mr/figures/mr_input_res.Rdata')

#如果表型是二分类变量（case和control），输入文件二选一：
#1）rs编号`rs_id`、P值`pval_nominal`、SNP的效应值`beta`、效应值方差`varbeta`；（推荐）
#2）rs编号`rs_id`、P值`pval_nominal`、case在所有样本中的比例`s`，MAF也要，写在list最后

#如果表型是连续型变量，输入文件三选一：
#1）rs编号`rs_id`、P值`pval_nominal`、表型的标准差`sdY`；
#2）rs编号`rs_id`、P值`pval_nominal`、效应值`beta`,效应值方差 `varbeta`, 样本量`N`,次等位基因频率 `MAF`；
#3）rs编号`rs_id`、P值`pval_nominal`、次等位基因频率 `MAF`；(推荐)


## 下载
# data <- vcfR::read.vcfR("./eqtl-a-ENSG00000239713.vcf/eqtl-a-ENSG00000239713.vcf")


如果是Safari下载，一定要先另存为，不然文件大小有问题
https://gwas.mrcieu.ac.uk/datasets/?trait__icontains=pqtl

gse <- list.files('./figures/',  'vcf')
gse

# install.packages("vcfR")
paste0('./figures/',gse)
data <- vcfR::read.vcfR(paste0('./figures/',gse[1]))
# data <- vcfR::read.vcfR("./figures/eqtl-a-ENSG00000173221.vcf")

ensg <- "ENSG00000120738"
eqtl_ensg <- paste0("eqtl-a-",ensg)
#NSUN5P1
# data <- vcfR::read.vcfR("./figures/eqtl-a-ENSG00000223705.vcf")

# 1.SNP ID rs开头
# 2.effect allele，此处相当于ALT
# 3.other allele，此处相当于REF
# 4.beta
# 5.se
# 6.pval


#整理数据
## 处理gt
gt=data@gt
gt=as.data.frame(gt)

colnames(gt)
gt$FORMAT[1]

library(tidyverse)

# gt$`eqtl-a-ENSG00000239713`[1]
colnames(gt)
gt[1,1]
# gt=separate(gt,col='eqtl-a-ENSG00000239713',into = c('ES', 'SE',
gt=separate(gt,col=colnames(gt)[2],into = c('ES', 'SE',
                                            'LP','AF','SS',
                                            'ID'),sep = '\\:')

gc()
gt[1,]


gt=na.omit(gt)
colnames(gt)=c('format','beta','se','logpvalue','maf','samplesize','snp')
gt$beta=as.numeric(gt$beta)
gt$se=as.numeric(gt$se)
gt$logpvalue=as.numeric(gt$logpvalue)
gt$maf=as.numeric(gt$maf)
gt$samplesize=as.numeric(gt$samplesize)

gc()
gt$format=NULL
gt[1,]

fix=data@fix
fix=as.data.frame(fix)
colnames(fix)
colnames(fix)=c('chr','pos','snp','ref','alt')
fix=fix[,1:5]


## 合并gt fix
eqtl=left_join(fix,gt,by='snp')
eqtl=na.omit(eqtl)
dim(eqtl)
# https://gwas.mrcieu.ac.uk/datasets/eqtl-a-ENSG00000203875/  SNGH5基因

# load('./figures/exposure_dat.Rdata')
colnames(exposure_dat) #在这里看看基因的ensg
unique(exposure_dat$exposure)

# exposure_dat1 <- exposure_dat[exposure_dat$exposure %in% c("SNHG5"),]
exposure_dat1 <- exposure_dat[exposure_dat$exposure %in% c("EGR1"),]
View(exposure_dat1)

## 查找染色体和位置
#22
#39605007

#SNHG5
#6
# 86149161

#NSUN5P1
# 7
# 75147934

#MGST1
# 12
# 16502260

#EGR1
15
80260014

或
3
188135783

unique(eqtl$chr)  #必须要有22条染色体，不然就是下载错了

eqtl=eqtl[eqtl$chr==15,]
# eqtl=eqtl[eqtl$chr==3,]
eqtl$logpvalue=as.numeric(eqtl$logpvalue)
eqtl$p_value=10^(-eqtl$logpvalue)

eqtl$pos=as.numeric(eqtl$pos)

## 上下1mkb
eqtl=eqtl[eqtl$pos > 80260014-1000000 ,]
eqtl=eqtl[eqtl$pos < 80260014+1000000 ,]

# eqtl=eqtl[eqtl$pos > 188135783-1000000 ,]
# eqtl=eqtl[eqtl$pos < 188135783+1000000 ,]

my_eqtl=eqtl[,c('snp','p_value','maf')]

colnames(my_eqtl)=c('snp','pvalues','MAF')
my_eqtl=na.omit(my_eqtl)

my_eqtl=my_eqtl[my_eqtl$MAF>0 ,]
dim(my_eqtl)
#大概要有3000-8000个snp
saveRDS(my_eqtl,'./figures/day5_my_eqtl.RDS')


##-----接下来gwas疾病


library(TwoSampleMR)
my_eqtl$snp  #一般5000-8000个

outcomes_gwas
## ！！！！！

outcomes_gwas <- "ukb-d-I50"
# outcomes_gwas <- "ebi-a-GCST009541"
coloc_SLE_dat <- extract_outcome_data(snps = c(my_eqtl$snp),
                                      # outcomes = "ebi-a-GCST90011866",
                                      outcomes = outcomes_gwas,
                                      proxies = F) %>%
  mutate(chr.outcome = as.numeric(chr),
         pos.outcome = as.numeric(pos),
         outcome = "Heart failure",
         id.outcome = outcomes_gwas)

saveRDS(coloc_SLE_dat, file = "figures/coloc_SLE_dat.RDS")
coloc_SLE_dat <- readRDS('./figures/coloc_SLE_dat.RDS')
#gwas指的是暴露、疾病
gwas=coloc_SLE_dat
gwas$beta=as.numeric(gwas$beta.outcome)
gwas$se=as.numeric(gwas$se.outcome)
#X=gwas$beta/gwas$se
#P=2*pnorm(q=abs(X), lower.tail=FALSE) #Z=1.96  P=0.05
#gwas$pvalue=P
## check Lp=log10(P)

#计算varbeta用于而分类变量表型（心衰和正常）
#如果表型是二分类变量（case和control），输入文件二选一：
#1）rs编号`rs_id`、P值`pval_nominal`、SNP的效应值`beta`、效应值方差`varbeta`；（推荐）
#2）rs编号`rs_id`、P值`pval_nominal`、case在所有样本中的比例`s`，MAF也要，写在list最后

#如果表型是连续型变量，输入文件三选一：
#1）rs编号`rs_id`、P值`pval_nominal`、表型的标准差`sdY`；
#2）rs编号`rs_id`、P值`pval_nominal`、效应值`beta`,效应值方差 `varbeta`, 样本量`N`,次等位基因频率 `MAF`；
#3）rs编号`rs_id`、P值`pval_nominal`、次等位基因频率 `MAF`；(推荐)

gwas$varbeta=(gwas$se)^2

gwas=gwas[,c('SNP','pval.outcome',"beta",'varbeta')]

colnames(gwas)=c('snp','pvalues','beta','varbeta')

gwas=na.omit(gwas)
colnames(gwas)
# gwas#即暴露（心衰）的所有snp和beta值、pvalues值和varbeta值
saveRDS(gwas, file = "./figures/gwas.RDS")

gwas <- readRDS('./figures/gwas.RDS')
gwas <- gwas[!duplicated(gwas$snp),]

##---开始coloc
library(coloc)
# my_eqtl#基因的gwas
input <- merge(my_eqtl, gwas, by="snp", all=FALSE, suffixes=c("_eqtl","_gwas"))
head(input)
library(coloc)


print(outcomes_gwas)

outcomes_gwas
# https://gwas.mrcieu.ac.uk/datasets/?trait__icontains=pqtl  #上官网查看

# ebi-a-GCST009541
# ncase	47,309
# ncontrol	930,014
# Sample size	977,323

#cc指二分类表型  outcomes_hf
#s指病例数/所有样本的比例
47309/977323
#N指control样本

print(outcomes_gwas)

# ukb-d-I50
ncase	1088
ncontrol	360106
Sample size	361194

ncase <- 1088
s_pro <- ncase/361194
s_pro
ncontrol <- 360106

# ebi-a-GCST009541
# # ncase	47309
# # ncontrol	930014
# # Sample size	977323
# # 
# ncase <- 47309
# s_pro <- ncase/977323
# s_pro
# ncontrol <- 930014


# quant指连续型变量
View(exposure_dat1)
unique(exposure_dat1$samplesize.exposure) #基因的影响不大
unique(exposure_dat1$samplesize.exposure)[1]
unique(exposure_dat1$samplesize.exposure)[2]

coloc.abf(dataset1=list(snp=input$snp,pvalues=input$pvalues_gwas,
                        # type="cc", s=0.003897353, N=ncontrol), #疾病
                        type="cc", s=s_pro, N=ncontrol), #疾病
          # dataset2=list(snp=input$snp,pvalues=input$pvalues_eqtl, type="quant", N=18235), #基因
          dataset2=list(snp=input$snp,pvalues=input$pvalues_eqtl,
                        # type="quant", N=9173), #基因
                        type="quant", N=unique(exposure_dat1$samplesize.exposure)[1]), #基因
          MAF=input$MAF)

## 或者用下面的方法（第二种方法，一般推荐上面这种）

#my_eqtl=as.list(my_eqtl)
#my_eqtl[['type']]='quant'
#my_eqtl[['N']]=31430

#gwas=as.list(gwas)
#gwas[['type']]='cc'
#gwas[['N']]=12653
#gwas$pvalues=NULL
#coloc.abf(dataset1 = gwas,dataset2 = my_eqtl)


PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf
4.46e-304 2.57e-305  9.42e-01  5.43e-02  3.43e-03
## 解读------
#H0：该区域的两个性状都没有遗传关联
#H1 / H2：只有表型1或表型2在该区域具有遗传关联
#H3：两个特征都相关，但因果变量不同
#H4：两个特征都相关并且共享一个因果变量(共享同一个snp的可能性)

H2:有90%的可能这些snp只跟一个表型有关
H3:有5%的可能这些snp与两个表型都相关
H4:共享同一个snp的可能性是0.3%






## 或者用下面的方法（第二种方法，一般推荐上面这种）

my_eqtl=as.list(my_eqtl)
my_eqtl[['type']]='quant'
# my_eqtl[['N']]=18235
my_eqtl[['N']]=unique(exposure_dat1$samplesize.exposure)[1]

gwas=as.list(gwas)
gwas[['type']]='cc'
gwas[['N']]=ncontrol
gwas$pvalues=NULL
coloc.abf(dataset1 = gwas,dataset2 = my_eqtl)
