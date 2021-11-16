
data=read.delim("/home/anoop/Desktop/CCR6/Data.txt",row.names = 1,header = TRUE)
data_F=na.omit(data)
dim(data_F)

library(NormalyzerDE)
?normalyzer
normalyzer(jobName="Prot",designPath = "/home/anoop/Desktop/CCR6/Metadata.txt",
           dataPath = "/home/anoop/Desktop/CCR6/ImputedData.txt",
           outputDir = "/home/anoop/Desktop/CCR6/Normed"
)


library(PCAtools)
count=read.delim("/home/anoop/Desktop/CCR6/Normed/Quantile.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/CCR6/Normed/Metadata.txt",row.names = 1)
p <- PCAtools::pca(count, metadata = meta, removeVar = 0.1)
write.table(p$rotated,file = "/home/anoop/Desktop/CCR6/Normed/Quantile_PCA.txt",sep = "\t",col.names = NA,quote = FALSE)

head(p$variance)

dat=read.delim("/home/anoop/Desktop/CCR6/Normed/Quantile_PCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/CCR6/Normed/Quantile_PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group2)) + geom_point(size=5,aes(fill=Group2,shape=Status))+theme_gray()+
  scale_shape_manual(values=c(21, 8, 22))+
  scale_color_manual(values=c(HC="#457979",ART="#cc8400",EC="#0000e5",Pool="#666666"))+
  scale_fill_manual(values=c(HC="#5ca2a2",ART="#ffa500",EC="#1919ff",Pool="#808080"))+#geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
   labs(x="PC1, 46.19% variance",y="PC2, 25.19% variance")+
  theme(axis.title = element_text(size=10),legend.position = c(0.85, 0.88),plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 2,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()


library(sva)

prot=read.delim("/home/anoop/Desktop/CCR6/Normed/Quantile.txt",row.names = 1,check.names = FALSE)
mat <- as.matrix(prot)
des=read.table("/home/anoop/Desktop/CCR6/Normed/Metadata.txt",sep="\t",header=TRUE)
designCombat = model.matrix(~ des$group)
rnaseqCombat = ComBat(mat, batch = des$batch, mod = designCombat, par.prior = TRUE, prior.plots = TRUE)
write.table(rnaseqCombat,file="/home/anoop/Desktop/CCR6/Normed/Combat_Result.txt", sep="\t",quote=FALSE,col.names = NA)


count=read.delim("/home/anoop/Desktop/CCR6/Normed/Combat_Result.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/CCR6/Normed/Metadata.txt",row.names = 1)
p <- PCAtools::pca(count, metadata = meta, removeVar = 0.1)
write.table(p$rotated,file = "/home/anoop/Desktop/CCR6/Normed/Combat_PCA.txt",sep = "\t",col.names = NA,quote = FALSE)
head(p$variance)
library(ggrepel)
dat=read.delim("/home/anoop/Desktop/CCR6/Normed/Combat_PCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/CCR6/Normed/Combat_PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group2)) + geom_point(size=5,aes(fill=Group2,shape=Status))+theme_gray()+
  scale_shape_manual(values=c(21, 8, 22))+
  scale_color_manual(values=c(HC="#457979",ART="#cc8400",EC="#0000e5",Pool="#666666"))+
  scale_fill_manual(values=c(HC="#5ca2a2",ART="#ffa500",EC="#1919ff",Pool="#808080"))+geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  labs(x="PC1, 16.24% variance",y="PC2, 10.92% variance")+
  theme(axis.title = element_text(size=10),legend.position = c(0.18, 0.86),plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 2,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()


library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/CCR6/Normed/CCR6_Neg/Info.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.table("/home/anoop/Desktop/CCR6/Normed/CCR6_Neg/NegData.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(ARTvsEC=EC - ART,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ARTvsEC",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/Normed/CCR6_Neg/Result.txt",sep="\t",quote=FALSE,col.names = NA)


sampleinfo <- read.delim("/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/Info.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.table("/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/PosData.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(ARTvsEC=EC - ART,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ARTvsEC",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/Result.txt",sep="\t",quote=FALSE,col.names = NA)

sampleinfo <- read.delim("/home/anoop/Desktop/CCR6/Normed/ART/Info.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
Pair <- factor(sampleinfo$Pair)
design <- model.matrix(~Pair+group)


group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.table("/home/anoop/Desktop/CCR6/Normed/ART/ART_Data.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(NeqvsPos=groupPos,levels=design)
cont.matrix <- makeContrasts(NeqvsPos=Pos - Neg,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="NeqvsPos",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/Normed/ART/Result_1.txt",sep="\t",quote=FALSE,col.names = NA)




sampleinfo <- read.delim("/home/anoop/Desktop/CCR6/Normed/EC/Info.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
Pair <- factor(sampleinfo$Pair)
design <- model.matrix(~Pair+group)


group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.table("/home/anoop/Desktop/CCR6/Normed/EC/EC_Data.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(NeqvsPos=groupPos,levels=design)
cont.matrix <- makeContrasts(NeqvsPos=Pos - Neg,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="NeqvsPos",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/Normed/EC/Result_1.txt",sep="\t",quote=FALSE,col.names = NA)


xx=read.delim("/home/anoop/Desktop/CCR6/Normed/ART/MW/ART_Data.txt",row.names = 1,check.names = FALSE)
head(xx)
yy=t(xx)
write.table(yy,file="/home/anoop/Desktop/CCR6/Normed/ART/MW/Data.txt",sep = "\t",col.names = NA,quote = FALSE)

########

count=read.delim("/home/anoop/Desktop/CCR6/Normed/EC/EC_Data.txt",row.names = 1,check.names = FALSE)
countX=as.matrix(count)
calc_coef_var <- function(x) sd(x) / mean(x)
coef_var <- apply(countX, 1, calc_coef_var)
HVG_5 <- countX[rank(coef_var) / length(coef_var) > 0.99, ]
write.table(HVG_5,file="/home/anoop/Desktop/CCR6/Normed/EC/EC_Data_Top5.txt",sep="\t",col.names = NA,quote = FALSE)


library(gplots)
Dat=read.delim("/home/anoop/Desktop/CCR6/Normed/EC/EC_Data_Top5.txt",header = TRUE,row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/CCR6/Normed/EC/EC_heatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table((X$carpet),file="/home/anoop/Desktop/CCR6/Normed/EC/EC_Data_Top5_Z.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/CCR6/Normed/EC/EC_Data_Top5_Z.txt",header = TRUE, row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/CCR6/Normed/EC/Info.txt",row.names = 1,header = TRUE)


#col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#990000","#e50000" ,"#ff9999","white","#b5cde1","#315b7d","#23415a"))
col_fun1 = colorRamp2(c(2, 1, 0, -1, -2), c("#2a4e6c","#4682b4" ,"#fffcc9","#ffa700", "#ff5a00"))

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

col_fun1 = colorRamp2(c(2, 0, -3), c("#2a4e6c","yellow", "#ff5a00"))
library(ComplexHeatmap)
library(circlize)

colours <- list("Group"=c("Neg"="#48690E","Pos"="#ffd700"))
LFC=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HeatMap/FC_trans.txt",header = T,check.names = FALSE,row.names = 1)

head(LFC)
ncol(Zscore)
H1=Heatmap(as.matrix(t(Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(2, "cm"),column_title_gp =gpar(fontsize = 0),
           top_annotation  =columnAnnotation(Group = sampleinfo$Group,col=colours,show_legend=TRUE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(30, "cm"),width  = unit(18, "cm"))


pdf("/home/anoop/Desktop/CCR6/Normed/EC/EC_heatmap.pdf",height = 20,width =15)
draw(H1, merge_legend = TRUE)
dev.off()
?rowAnnotation


count=read.delim("/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/PosData.txt",row.names = 1,check.names = FALSE)
countX=as.matrix(count)
calc_coef_var <- function(x) sd(x) / mean(x)
coef_var <- apply(countX, 1, calc_coef_var)
HVG_5 <- countX[rank(coef_var) / length(coef_var) > 0.99, ]
write.table(HVG_5,file="/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/PosDataTop1.txt",sep="\t",col.names = NA,quote = FALSE)


library(gplots)
Dat=read.delim("/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/PosDataTop1.txt",header = TRUE,row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/PosDataTop1.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table((X$carpet),file="/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/PosDataTop1_Z.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/PosDataTop1_Z.txt",header = TRUE, row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/Info.txt",row.names = 1,header = TRUE)


#col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#990000","#e50000" ,"#ff9999","white","#b5cde1","#315b7d","#23415a"))
col_fun1 = colorRamp2(c(2, 1, 0, -1, -2), c("#2a4e6c","#4682b4" ,"#fffcc9","#ffa700", "#ff5a00"))

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

col_fun1 = colorRamp2(c(2, 0, -3), c("#2a4e6c","yellow", "#ff5a00"))
library(ComplexHeatmap)
library(circlize)

colours <- list("Group"=c("ART"="#48690E","EC"="#ffd700",HC="grey"))
LFC=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HeatMap/FC_trans.txt",header = T,check.names = FALSE,row.names = 1)

head(LFC)
ncol(Zscore)
H1=Heatmap(as.matrix(t(Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(2, "cm"),column_title_gp =gpar(fontsize = 0),
           top_annotation  =columnAnnotation(Group = sampleinfo$Group,col=colours,show_legend=TRUE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 10),height  = unit(30, "cm"),width  = unit(18, "cm"))


pdf("/home/anoop/Desktop/CCR6/Normed/CCR6_Pos/PosDataTop1.pdf",height = 20,width =15)
draw(H1, merge_legend = TRUE)
dev.off()



pval=read.delim("/home/anoop/Desktop/CCR6/Normed/ART/MW/MW_Result.txt",header = TRUE)
head(pval)
pval$AdjPval =p.adjust(pval$P_value,method = "BH")
write.table(pval,file="/home/anoop/Desktop/CCR6/Normed/ART/MW/MW_Result_ADj.txt",sep="\t",col.names = NA, quote = FALSE)


xx=read.delim("/home/anoop/Desktop/CCR6/Normed/EC/MW/EC_Data.txt",row.names = 1,check.names = FALSE)
head(xx)
yy=t(xx)
write.table(yy,file="/home/anoop/Desktop/CCR6/Normed/EC/MW/Data.txt",sep = "\t",col.names = NA,quote = FALSE)


pval=read.delim("/home/anoop/Desktop/CCR6/Normed/EC/MW/MW_Result.txt",header = TRUE)
head(pval)
pval$AdjPval =p.adjust(pval$P_value,method = "BH")
write.table(pval,file="/home/anoop/Desktop/CCR6/Normed/EC/MW/MW_Result_Adj.txt",sep="\t",col.names = NA, quote = FALSE)


#############################################

BiocManager::install("DEqMS")
n

library(DEqMS)
url <- "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/06/PXD004163/Yan_miR_Protein_table.flatprottable.txt"
download.file(url, destfile = "./miR_Proteintable.txt",method = "auto")
df.prot = read.table("miR_Proteintable.txt",stringsAsFactors = FALSE,header = TRUE, quote = "", comment.char = "",sep = "\t")

head(dat.log)
TMT_columns = seq(15,33,2)
dat = df.prot[df.prot$miR.FASP_q.value<0.01,TMT_columns]
rownames(dat) = df.prot[df.prot$miR.FASP_q.value<0.01,]$Protein.accession
dat.log = log2(dat)
dat.log = na.omit(dat.log)

boxplot(dat.log,las=2,main="TMT10plex data PXD004163")
cond = as.factor(c("ctrl","miR191","miR372","miR519","ctrl",
                   "miR372","miR519","ctrl","miR191","miR372"))

design = model.matrix(~0+cond)
colnames(design) = gsub("cond","",colnames(design))


x <- c("miR372-ctrl","miR519-ctrl","miR191-ctrl",
       "miR372-miR519","miR372-miR191","miR519-miR191")
contrast =  makeContrasts(contrasts=x,levels=design)
fit1 <- lmFit(dat.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)

library(matrixStats)
count_columns = seq(16,34,2)
psm.count.table = data.frame(count = rowMins(as.matrix(df.prot[,count_columns])), row.names =  df.prot$Protein.accession)

fit3$count = psm.count.table[rownames(fit3$coefficients),"count"]

head(psm.count.table,20)
fit4 = spectraCounteBayes(fit3)

write.table(df.prot,file="/home/anoop/Desktop/CCR6/Test_Data.txt",sep="\t",col.names = NA,quote = FALSE)


#########################

BiocManager::install("ROTS")
n
library(ROTS)
data(upsSpikeIn)
input = upsSpikeIn

groups = c(rep(1,3), rep(2,3))
groups=as.vector(groups)
class(groups)
results = ROTS(data = myData, groups = groups , B = 1000 , K = 2000 , seed = 1234,paired = TRUE)
?ROTS
myData=read.delim("/home/anoop/Desktop/CCR6/Normed/EC/EC_Data3.txt",row.names = 1)

library(reshape2)
tst=as.data.frame(do.call(cbind, results))

summary(results, fdr = 0.1)


pval=read.delim("/home/anoop/Desktop/CCR6/Normed/EC/Result_Names.txt",header = TRUE)
head(pval)
pval$AdjPval =p.adjust(pval$P.Value,method = "BH")
write.table(pval,file="/home/anoop/Desktop/CCR6/Normed/EC/Result_Names2.txt",sep="\t",col.names = NA, quote = FALSE)
?read.delim


tt=read.delim("/home/anoop/Desktop/CCR6/MW/EC_ART.txt",row.names = 1,check.names = FALSE)
YY=t(tt)
write.table(YY,file="/home/anoop/Desktop/CCR6/MW/Trans.txt",sep="\t",col.names = NA, quote = FALSE)

################################### NON Fraction ################################


library(NormalyzerDE)
?normalyzer
normalyzer(jobName="Prot",designPath = "/home/anoop/Desktop/CCR6/NonFraction/Analysis/Metadata.txt",
           dataPath = "/home/anoop/Desktop/CCR6/NonFraction/Analysis/Merged.txt",
           outputDir = "/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized"
)


library(impute)
?impute.knn
BiocManager::install("impute")
data=read.delim("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Quantile.txt",row.names = 1)
X=as.matrix(data)
KNN=impute.knn(data=X ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
write.table(KNN$data,file="/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Imputed.txt",sep="\t",col.names = NA,quote = FALSE)


prot=read.delim("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Imputed.txt",row.names = 1,check.names = FALSE)
mat <- as.matrix(prot)
des=read.table("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Metadata.txt",sep="\t",header=TRUE)
designCombat = model.matrix(~ des$group)
rnaseqCombat = ComBat(mat, batch = des$batch, mod = designCombat, par.prior = TRUE, prior.plots = TRUE)
write.table(rnaseqCombat,file="/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Batch_Corrected.txt", sep="\t",quote=FALSE,col.names = NA)


library(PCAtools)
count=read.delim("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Batch_Corrected.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Metadata.txt",row.names = 1)
p <- PCAtools::pca(count, metadata = meta, removeVar = 0.1)
write.table(p$rotated,file = "/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Batch_Corrected_PCA.txt",sep = "\t",col.names = NA,quote = FALSE)

head(p$variance)

dat=read.delim("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Batch_Corrected_PCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Batch_Corrected_PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=5,aes(fill=Group,shape=Status))+theme_gray()+
  scale_shape_manual(values=c(21, 8, 22))+
  scale_color_manual(values=c(HC="#457979",ART="#cc8400",EC="#0000e5",Pool="#666666"))+
  scale_fill_manual(values=c(HC="#5ca2a2",ART="#ffa500",EC="#1919ff",Pool="#808080"))+geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  labs(x="PC1, 23.95% variance",y="PC2, 11.68% variance")+
  theme(axis.title = element_text(size=10),legend.position = c(0.15, 0.12),plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 2,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()


sampleinfo <- read.delim("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Limma/EC/des.txt")
group <- paste(sampleinfo$group)
group <- factor(group)
Pair <- factor(sampleinfo$pair)
design <- model.matrix(~Pair+group)

data=read.table("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Limma/EC/Data.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(NeqvsPos=groupEC_pos,levels=design)

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="NeqvsPos",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Limma/EC/Result.txt",sep="\t",quote=FALSE,col.names = NA)

sampleinfo <- read.delim("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Limma/ART/design.txt")
group <- paste(sampleinfo$group)
group <- factor(group)
Pair <- factor(sampleinfo$pair)
design <- model.matrix(~Pair+group)

data=read.table("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Limma/ART/Data.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(NeqvsPos=groupEC_pos,levels=design)

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="NeqvsPos",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Limma/ART/Result.txt",sep="\t",quote=FALSE,col.names = NA)

tt=read.delim("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Limma/MW/Input.txt",row.names = 1,check.names = FALSE)
YY=t(tt)
write.table(YY,file="/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Limma/MW/Trans.txt",sep="\t",col.names = NA, quote = FALSE)



myData=read.delim("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Limma/EC/Data.txt",row.names = 1)
groups = c(rep(1,4), rep(2,4))
results = ROTS(data = myData, groups = groups , B = 1000 , K = 20 , seed = 1234,paired = TRUE)
summary(results, fdr = 0.05)

pval=as.data.frame(results$pvalue)
fdr=as.data.frame(results$FDR)
lfc=as.data.frame(results$logfc)
ff=plot(results, fdr = 0.05, type = "volcano")
Res=merge(lfc,pval,fdr,by.y="results$logfc")

myData=read.delim("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Limma/ART/Data.txt",row.names = 1)
groups = c(rep(1,6), rep(2,6))
results = ROTS(data = myData, groups = groups , B = 1000 , K = 20 , seed = 1234,paired = TRUE)
summary(results, fdr = 0.1)

#############

pval=read.delim("/home/anoop/Desktop/CCR6/MW/MW_Result.txt",header = TRUE)
head(pval)
pval$AdjPval =p.adjust(pval$P_value,method = "BH")
write.table(pval,file="/home/anoop/Desktop/CCR6/MW/MW_Result_Adj.txt",sep="\t",col.names = NA, quote = FALSE)


########################### NEW Analysis
sampleinfo=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/EC/Meta.txt")
pair <- factor(sampleinfo$pair)
group <- factor(sampleinfo$group, levels=c("Neg","Pos"))
design <- model.matrix(~pair+group)
data=read.table("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/EC/Pos_Neg.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat, design)
fit <- eBayes(fit)
limma.res=topTable(fit, coef="groupPos",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/EC/Pos_Neg_Results.txt",sep="\t",quote=FALSE,col.names = NA)


sampleinfo=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/EC/Info.txt")
f <- factor(sampleinfo$Group, levels=c("EC_Pos","EC_Neg","HC_Pos","HC_Neg"))
design <- model.matrix(~0+f)
colnames(design) <- c("EC_Pos","EC_Neg","HC_Pos","HC_Neg")
data=read.table("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/EC/EC.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat, design)
contrast.matrix <- makeContrasts(Pos=EC_Pos-HC_Pos, Neg=EC_Neg-HC_Neg,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
limma.res=topTable(fit2, coef="Pos",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/EC/Pos_HC_EC_Results.txt",sep="\t",quote=FALSE,col.names = NA)

limma.res=topTable(fit2, coef="Neg",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/EC/Neg_HC_EC_Results.txt",sep="\t",quote=FALSE,col.names = NA)


sampleinfo=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/Meta.txt")
pair <- factor(sampleinfo$pair)
group <- factor(sampleinfo$group, levels=c("Neg","Pos"))
design <- model.matrix(~pair+group)
data=read.table("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/ART_Neg_Pos.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat, design)
fit <- eBayes(fit)
limma.res=topTable(fit, coef="groupPos",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/ART_Neg_PosResults.txt",sep="\t",quote=FALSE,col.names = NA)


sampleinfo=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/Info.txt")
f <- factor(sampleinfo$Group, levels=c("ART_Pos","ART_Neg","HC_Pos","HC_Neg"))
design <- model.matrix(~0+f)
colnames(design) <- c("ART_Pos","ART_Neg","HC_Pos","HC_Neg")
data=read.table("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/ART.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat, design)
contrast.matrix <- makeContrasts(Pos=ART_Pos-HC_Pos, Neg=ART_Neg-HC_Neg,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
limma.res=topTable(fit2, coef="Pos",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/Pos_HC_ART.txt",sep="\t",quote=FALSE,col.names = NA)

limma.res=topTable(fit2, coef="Neg",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/Neg_HC_ART.txt",sep="\t",quote=FALSE,col.names = NA)

#############
sampleinfo=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/Info.txt")
f <- factor(sampleinfo$Group, levels=c("ART_pos","ART_neg","EC_pos","EC_neg"))
design <- model.matrix(~0+f)
colnames(design) <- c("ART_pos","ART_neg","EC_pos","EC_neg")
data=read.table("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/Data.txt",sep="\t",row.names=1,header=TRUE)
mat <- as.matrix(data)
fit <- lmFit(mat, design)
contrast.matrix <- makeContrasts(Pos=EC_pos-ART_pos, Neg=EC_neg-ART_neg,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
limma.res=topTable(fit2, coef="Pos",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/Pos_ART_EC.txt",sep="\t",quote=FALSE,col.names = NA)

limma.res=topTable(fit2, coef="Neg",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/Neg_ART_EC.txt",sep="\t",quote=FALSE,col.names = NA)

##########################

library(doParallel)
library(MUVR)

nCore=10
nRep=25
nOuter=4
varRatio=0.85
method='RF'

data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/MUVR/EC/Pos_Neg.txt",row.names = 1,check.names = FALSE)
TR=t(data)
rownames(TR) <- c()
write.table(TR,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/MUVR/EC/Pos_Neg_Tr.txt",sep="\t",quote=FALSE,col.names = NA)

XX=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/MUVR/EC/Pos_Neg_Tr.txt",check.names = FALSE)
head(XX)
Y=c(rep("Pos",4),rep("Neg",4))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, ID=XX$ID,nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar
?MUVR
vip=getVIP(classModel, model='min')
vip
write.table(vip,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/MUVR/EC/Variables.txt",sep="\t",col.names = NA,quote = FALSE)


###########

nCore=10
nRep=25
nOuter=4
varRatio=0.85
method='RF'

data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/MUVR/ART/ART_Neg_Pos.txt",row.names = 1,check.names = FALSE)
TR=t(data)
rownames(TR) <- c()
write.table(TR,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/MUVR/ART/ART_Neg_PosTr.txt",sep="\t",quote=FALSE,col.names = NA)

XX=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/MUVR/ART/ART_Neg_PosTr.txt",check.names = FALSE)
head(XX)
Y=c(rep("Pos",6),rep("Neg",6))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, ID=XX$ID,nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar
?MUVR
vip=getVIP(classModel, model='min')
vip
write.table(vip,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/MUVR/ART/Variables.txt",sep="\t",col.names = NA,quote = FALSE)



##############################################


BiocManager::install("piano")
n
library(piano)

data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/Pos.txt",header = TRUE,row.names = 9)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/KEGG_Pos_EC_ART.txt",sep="\t",col.names = NA)

#################

data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/Neg.txt",header = TRUE,row.names = 9)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/KEGG_Neg.txt",sep="\t",col.names = NA)

#####################

data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/ART.txt",header = TRUE,row.names = 9)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/ART_kegg.txt",sep="\t",col.names = NA)

######################

data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/EC/KEGG/EC.txt",header = TRUE,row.names = 9)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/EC/KEGG/EC_kegg.txt",sep="\t",col.names = NA)

################################################# Figures


data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/ART_Pos_Neg.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/ART_Pos_NegVolcano.pdf")
ggplot(data, aes(x=logFC, y=-log10(P.Value))) + 
  geom_point(data=subset(data, P.Value>=0.05),aes(x=logFC, y=-log10(P.Value),size=-log10(P.Value)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA)+
  geom_point(data=subset(data, P.Value<0.05 & logFC < 0),aes(x=logFC, y=-log10(P.Value),color=abs(logFC),size=-log10(P.Value)))+
  scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, P.Value<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(P.Value),fill=abs(logFC),size=-log10(P.Value)),color="transparent")+
  scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/EC_ART_neg.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/EC_ART_negVolcano.pdf")
ggplot(data, aes(x=logFC, y=-log10(P.Value))) + 
  geom_point(data=subset(data, P.Value>=0.05),aes(x=logFC, y=-log10(P.Value),size=-log10(P.Value)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA)+
  geom_point(data=subset(data, P.Value<0.05 & logFC < 0),aes(x=logFC, y=-log10(P.Value),color=abs(logFC),size=-log10(P.Value)))+
  scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, P.Value<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(P.Value),fill=abs(logFC),size=-log10(P.Value)),color="transparent")+
  scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/EC_ART_pos.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/EC_ART_posVocano.pdf")
ggplot(data, aes(x=logFC, y=-log10(P.Value))) + 
  geom_point(data=subset(data, P.Value>=0.05),aes(x=logFC, y=-log10(P.Value),size=-log10(P.Value)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA)+
  geom_point(data=subset(data, P.Value<0.05 & logFC < 0),aes(x=logFC, y=-log10(P.Value),color=abs(logFC),size=-log10(P.Value)))+
  scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, P.Value<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(P.Value),fill=abs(logFC),size=-log10(P.Value)),color="transparent")+
  scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/EC_Pos_Neg.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/EC_Pos_NegVolcano.pdf")
ggplot(data, aes(x=logFC, y=-log10(P.Value))) + 
  geom_point(data=subset(data, P.Value>=0.05),aes(x=logFC, y=-log10(P.Value),size=-log10(P.Value)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA)+
  geom_point(data=subset(data, P.Value<0.05 & logFC < 0),aes(x=logFC, y=-log10(P.Value),color=abs(logFC),size=-log10(P.Value)))+
  scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, P.Value<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(P.Value),fill=abs(logFC),size=-log10(P.Value)),color="transparent")+
  scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (Pvalue)")
dev.off()


###################################################

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/GSEA/HallMark.txt",row.names = 1,sep = "\t")
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/GSEA/HallMark1.pdf",family = "Arial")
ggplot(data, aes(x=Order, y=NES)) + theme_bw()+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=1,
                   box.padding=0.95,nudge_x = 1.8,nudge_y = -0.01)+
  geom_point(data=subset(data, NES>0),aes(x=Order, y=NES, size=Log),color="#90001c",fill="#ba0024")+geom_hline(yintercept=0,size=0.3)+
  geom_point(data=subset(data, NES<0),aes(x=Order, y=NES, size=Log),color="#004000",fill="#008000")+
  scale_size(range = c(3, 10))+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        axis.title.y=element_text(size=15),legend.position = "bottom",axis.ticks.x = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),axis.text.x =element_text(size=10,color="black",vjust =61),plot.margin = margin(3,0.5,3,0.5, "cm"))+
  guides(size=guide_legend(override.aes=list(color="grey",fill="grey")))+
  labs(y="NES(EC-vs-ART)",size="-Log10(FDR)")
dev.off()


#########################

data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/ART.txt",header = TRUE,row.names = 9)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/test.txt",sep="\t",col.names = NA)


##############
library(piano)
data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/Pos.txt",header = TRUE,row.names = 9)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/hallMarkGeneset.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/HallMark_Neg.txt",sep="\t",col.names = NA)


pdf("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/Neg_test.pdf")
nw <- networkPlot2(gsares,class="distinct",direction = "both",adjusted=TRUE,significance=0.1)
dev.off()

head(nw$x$edges)

write.table(as.data.frame(nw$x$edges),file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/Pos_edge.txt",sep="\t",col.names = NA)
write.table(as.data.frame(nw$x$nodes),file="/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/CCR6/Pos_node.txt",sep="\t",col.names = NA)
#######################################################

ip=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/Hall_pos.txt",header = TRUE)
ip$Name <- factor(ip$Name, levels = ip$Name)
pdf("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/Hall_pos.pdf")
ggplot(ip, aes(y=Name)) + 
  geom_point(data=ip,aes(x=1,y=Name,size=Lpval,color=Direction))+scale_x_discrete(limits=c("1"))+scale_color_manual(values=c(Up="#f58564",Down="#09b9b7"))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(4,5,5.2,4.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "-log10(Padj)"),color=guide_legend(title = "Direction"))
dev.off()


ip=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/Hall_neg.txt",header = TRUE)
ip$Name <- factor(ip$Name, levels = ip$Name)
pdf("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/Hall_neg.pdf")
ggplot(ip, aes(y=Name)) + 
  geom_point(data=ip,aes(x=1,y=Name,size=Lpval,color=Direction))+scale_x_discrete(limits=c("1"))+scale_color_manual(values=c(Up="#f58564",Down="#09b9b7"))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(2.7,5,5.2,4.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "-log10(Padj)"),color=guide_legend(title = "Direction"))
dev.off()

#########################


data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/EC_ART_pos.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/EC_ART_posVolcano_2.pdf")
ggplot(data, aes(x=logFC, y=-log10(P.Value ))) + geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA)+
  geom_point(data=subset(data, P.Value <.05 & logFC <= -1),aes(x=logFC,y=-log10(P.Value )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, P.Value <.05 & logFC >= 1),aes(x=logFC,y=-log10(P.Value )),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, P.Value <.05 & logFC > 0 & logFC < 1),aes(x=logFC,y=-log10(P.Value )),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, P.Value <.05 & logFC < 0 & logFC > -1),aes(x=logFC,y=-log10(P.Value )),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, P.Value >=.05),aes(x=logFC,y=-log10(P.Value )),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),plot.margin = margin(0.9,1,0.9,1, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/EC_ART_neg.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/Figures/EC_ART_negVolcano_2.pdf")
ggplot(data, aes(x=logFC, y=-log10(P.Value ))) + geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA)+
  geom_point(data=subset(data, P.Value <.05 & logFC <= -1),aes(x=logFC,y=-log10(P.Value )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, P.Value <.05 & logFC >= 1),aes(x=logFC,y=-log10(P.Value )),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, P.Value <.05 & logFC > 0 & logFC < 1),aes(x=logFC,y=-log10(P.Value )),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, P.Value <.05 & logFC < 0 & logFC > -1),aes(x=logFC,y=-log10(P.Value )),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, P.Value >=.05),aes(x=logFC,y=-log10(P.Value )),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),plot.margin = margin(0.9,1,0.9,1, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Pvalue)")
dev.off()
packageVersion("ggplot2")
?impute.knn::impute
