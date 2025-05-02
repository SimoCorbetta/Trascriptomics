## DATA PREPROCESSING
library(recount3)
library(ggplot2)
## checking samples on brain
rse_brain <- readRDS("rse_brain.RDS")
assays(rse_brain)$counts <- transform_counts(rse_brain) # traform counts as we know
rin_b<-c(colData(rse_brain)$gtex.smrin[13],colData(rse_brain)$gtex.smrin[14],colData(rse_brain)$gtex.smrin[15])# rin
df <- data.frame(x = 1:length(rin_b), y = rin_b,names = c("b_13", "b_14", "b_15"))
#barplot del RIN
x11()
ggplot(df, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "steelblue",width=0.35) +
  labs(title = "RIN in brain samples",y = "rin") +
  ylim(0, 10)+geom_hline(yintercept = 6, linetype = "dashed", color = "red",linewidth=1.5)+
  geom_text(aes(label = names), vjust = 1.5, size = 4,y = 0)+
  scale_x_discrete(labels = df$names)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))
#barplot rRNA %
rrna_b<-c(colData(rse_brain)$gtex.smrrnart[13],colData(rse_brain)$gtex.smrrnart[14],colData(rse_brain)$gtex.smrrnart[15])#rrna
df <- data.frame(x = 1:length(rin_b), y = rrna_b,names = c("b_13", "b_14", "b_15"))
x11()
ggplot(df, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "steelblue",width=0.35) +
  labs(title = "rRNA percentage in brain samples",y = "percentage of rRNA") +
  ylim(0, 0.15)+geom_hline(yintercept = 0.1, linetype = "dashed", color = "red",linewidth=1.5)+
  geom_text(aes(label = names), vjust = 1.5, size = 4,y = 0)+
  scale_x_discrete(labels = df$names)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))
#barplot mapping quality
mapping_b<-c(colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[13],
           colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[14],
           colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[15])#checking mapping quality
df <- data.frame(x = 1:length(rin_b), y = mapping_b,names = c("b_13", "b_14", "b_15"))
x11()
ggplot(df, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "steelblue",width=0.35) +
  labs(title = "percentage of uniquely mapping reads in brain samples",y = "% of uniquely mapping reads") +
  ylim(0, 100)+geom_hline(yintercept = 85, linetype = "dashed", color = "red",linewidth=1.5)+
  geom_text(aes(label = names), vjust = 1.5, size = 4,y = 0)+
  scale_x_discrete(labels = df$names)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))
rse_brain_selected <- rse_brain[,c(13,14,15)] # select my 3 replicates from brain since they passed my test
counts_brain_selected <- assays(rse_brain_selected)$counts # counts for my selected columns
colnames(counts_brain_selected)<-c("b_13","b_14","b_15")
location_b<-table(colData(rse_brain_selected)$gtex.smtsd) # which part of the organ my samples come from
barplot(location_b)
sex_b<-table(colData(rse_brain_selected)$gtex.sex) # sex of the patients
age_b<-table(colData(rse_brain_selected)$gtex.age) # age of my patients


#checking samples on colon
rse_colon <- readRDS("rse_colon.RDS")
assays(rse_colon)$counts <- transform_counts(rse_colon)
rin_c<-c(colData(rse_colon)$gtex.smrin[14],colData(rse_colon)$gtex.smrin[15],colData(rse_colon)$gtex.smrin[16])# rin
df_one <- data.frame(x = 1:length(rin_c), y = rin_c,names = c("c_14", "c_15", "c_16"))

#barplot del RIN
x11()
ggplot(df_one, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "steelblue",width=0.35) +
  labs(title = "RIN in colon samples",y = "rin") +
  ylim(0, 10)+geom_hline(yintercept = 6, linetype = "dashed", color = "red",linewidth=1.5)+
  geom_text(aes(label = names), vjust = 1.5, size = 4,y = 0)+
  scale_x_discrete(labels = df_one$names)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))

#barplot rRNA %
rrna_c<-c(colData(rse_colon)$gtex.smrrnart[14],colData(rse_colon)$gtex.smrrnart[15],colData(rse_colon)$gtex.smrrnart[16])#rrna
df_one <- data.frame(x = 1:length(rin_c), y = rrna_c,names = c("c_14", "c_15", "c_16"))
x11()
ggplot(df_one, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "steelblue",width=0.35) +
  labs(title = "rRNA percentage in colon samples",y = "percentage of rRNA") +
  ylim(0, 0.15)+geom_hline(yintercept = 0.1, linetype = "dashed", color = "red",linewidth=1.5)+
  geom_text(aes(label = names), vjust = 1.5, size = 4,y = 0)+
  scale_x_discrete(labels = df_one$names)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))
#barplot % of uniquely mapping reads
mapping_c<-c(colData(rse_colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[14],
           colData(rse_colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[15],
           colData(rse_colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[16])
df_one <- data.frame(x = 1:length(rin_b), y = mapping_c,names = c("c_14", "c_15", "c_16"))
x11()
ggplot(df_one, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "steelblue",width=0.35) +
  labs(title = "percentage of uniquely mapping reads in colon samples",y = "% of uniquely mapping reads") +
  ylim(0, 100)+geom_hline(yintercept = 85, linetype = "dashed", color = "red",linewidth=1.5)+
  geom_text(aes(label = names), vjust = 1.5, size = 4,y = 0)+
  scale_x_discrete(labels = df_one$names)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))
rse_colon_selected <- rse_colon[,c(14,15,16)] # select my 3 replicates from brain since they passed my test
counts_colon_selected <- assays(rse_colon_selected)$counts # counts for my selected columns
colnames(counts_colon_selected)<-c("c_14","c_15","c_16")
location_c<-table(colData(rse_colon_selected)$gtex.smtsd)  # which part of the organ my samples come from
sex_c<-table(colData(rse_colon_selected)$gtex.sex) # sex of the patients
age_c<-table(colData(rse_colon_selected)$gtex.age) # age of my patients

#### check samples on kidney
rse_kidney <- readRDS("rse_kidney.RDS")
assays(rse_kidney)$counts <- transform_counts(rse_kidney)
rin_k<-c(colData(rse_kidney)$gtex.smrin[17],colData(rse_kidney)$gtex.smrin[18],colData(rse_kidney)$gtex.smrin[20])# rin
rin_ks<-c(colData(rse_kidney)$gtex.smrin[13],colData(rse_kidney)$gtex.smrin[17],colData(rse_kidney)$gtex.smrin[18],colData(rse_kidney)$gtex.smrin[20])
#barplot rin

df_two <- data.frame(x = 1:length(rin_k), y = rin_k,names = c("k_17", "k_18", "k_20"))
x11()
ggplot(df_two, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "steelblue",width=0.35) +
  labs(title = "RIN in kidney samples",y = "rin") +
  ylim(0, 10)+geom_hline(yintercept = 6, linetype = "dashed", color = "red",linewidth=1.5)+
  geom_text(aes(label = names), vjust = 1.5, size = 4,y = 0)+
  scale_x_discrete(labels = df_two$names)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))
#barplot rin con errore
df_two <- data.frame(x = 1:length(rin_ks), y = rin_ks,names = c("k_13","k_17", "k_18", "k_20"))
x11()
p<-ggplot(df_two, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "blue",width=0.35) +
  labs(title = "RIN in kidney samples",y = "rin") +
  ylim(0, 10)+geom_hline(yintercept = 6, linetype = "dashed", color = "black",linewidth=1.5)+
  geom_text(aes(label = names), vjust = 1.5, size = 4,y = 0)+
  scale_x_discrete(labels = df_two$names)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))
ggsave("rin_conerrore.png",p,dpi=300)
#barplot rrna %

rrna_k<-c(colData(rse_kidney)$gtex.smrrnart[17],colData(rse_kidney)$gtex.smrrnart[18],colData(rse_kidney)$gtex.smrrnart[20])#rrna
df_two <- data.frame(x = 1:length(rin_c), y = rrna_k,names = c("k_17", "k_18", "k_20"))
x11()
ggplot(df_two, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "steelblue",width=0.35) +
  labs(title = "rRNA percentage in kidney samples",y = "percentage of rRNA") +
  ylim(0, 0.15)+geom_hline(yintercept = 0.1, linetype = "dashed", color = "red",linewidth=1.5)+
  geom_text(aes(label = names), vjust = 1.5, size = 4,y = 0)+
  scale_x_discrete(labels = df_two$names)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))

#barplot % of uniquely mapping reads
mapping_k<-c(colData(rse_kidney)$"recount_qc.star.uniquely_mapped_reads_%_both"[17],
             colData(rse_kidney)$"recount_qc.star.uniquely_mapped_reads_%_both"[18],
             colData(rse_kidney)$"recount_qc.star.uniquely_mapped_reads_%_both"[20])
df_two <- data.frame(x = 1:length(rin_b), y = mapping_k,names = c("k_17", "k_18", "k_20"))
x11()
ggplot(df_two, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "steelblue",width=0.35) +
  labs(title = "percentage of uniquely mapping reads in kidney samples",y = "% of uniquely mapping reads") +
  ylim(0, 100)+geom_hline(yintercept = 85, linetype = "dashed", color = "red",linewidth=1.5)+
  geom_text(aes(label = names), vjust = 1.5, size = 4,y = 0)+
  scale_x_discrete(labels = df_two$names)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5))
rse_kidney_selected <- rse_kidney[,c(17,18,20)] # select my 3 replicates from brain since they passed my test
counts_kidney_selected <- assays(rse_kidney_selected)$counts # counts for my selected columns
colnames(counts_kidney_selected)<-c("k_17","k_18","k_20")
location_k<-table(colData(rse_kidney_selected)$gtex.smtsd)  # which part of the organ my samples come from
sex_k<-table(colData(rse_kidney_selected)$gtex.sex) # sex of the patients
age_k<-table(colData(rse_kidney_selected)$gtex.age) # age of my patients
##### function for filtering genes
filter_rse <- function(rse) {
  rRNA_or_pseudogene <- rowData(rse)$gbkey %in% c('rRNA', 'Gene')
  mitochondrial <- rowRanges(rse)@seqnames == 'chrM'
  too_short <- rowData(rse)[, 'bp_length'] < 200
  # Get a vector with TRUE for any row for which at least one
  # of the filters is TRUE.
  filter_out <- rRNA_or_pseudogene | mitochondrial | too_short
  return(rse[!filter_out])
}

filtered_brain<-filter_rse(rse_brain_selected) # 13139 genes removed dal df con i miei 3 campioni
filtered_colon<-filter_rse(rse_colon_selected)
filtered_kidney<-filter_rse(rse_kidney_selected) 

counts_brain_filter<-assays(filtered_brain)$counts # count table dei 3 campioni avendo rimosso i geni non voluti
rownames(counts_brain_filter)<-rowRanges(filtered_brain)$gene_name # assegno alle righe il nome comune del gene
tot_reads_brain_filtered<-apply(counts_brain_filter,2,sum) # total number of reads for each sample in brain
colnames(counts_brain_filter)<-c("b_13","b_14","b_15")


counts_colon_filter<-assays(filtered_colon)$counts
rownames(counts_colon_filter)<-rowRanges(filtered_colon)$gene_name
tot_reads_colon_filtered<-apply(counts_colon_filter,2,sum) # total number of reads for each sample in colon
colnames(counts_colon_filter)<-c("c_14","c_15","c_16")

counts_kidney_filter<-assays(filtered_kidney)$counts # assegno nome ufficiale del gene come nome della riga
rownames(counts_kidney_filter)<-rowRanges(filtered_kidney)$gene_name
tot_reads_kidney_filtered<-apply(counts_kidney_filter,2,sum) # total number of reads for each sample in colon
colnames(counts_kidney_filter)<-c("k_17","k_18","k_20")


## rimozione di tutti i geni non annotati sui cromosomi tradizionali in kidney
can_chr<-levels(rowRanges(filtered_brain)@seqnames)[1:24] # solo i cromosomi canonici
#rowRanges(filtered_kidney)@seqnames[1]  %in% can_chr ## boolean operation
vk <- rep(F,40903)


for (i in 1: 40903) {
  vk[i] = (as.vector(rowRanges(filtered_kidney)@seqnames[i]  %in% can_chr)) # maschera che filtra i geni
}
counts_removed_kidney<-counts_kidney_filter[vk,] # count table con solo i geni annotati su cromosomi canonici
length(unique(rownames(counts_removed_kidney)))# there are 36141 uniquely named genes
dim(counts_removed_kidney) # while there are 36179 rows in the count table (ci sono 38 geni in piu comunque)
dupl_genes_k<-counts_removed_kidney[duplicated(rownames(counts_removed_kidney)),] # duplicated genes have counts = to 0
not_dupl_kidney<- counts_removed_kidney[!duplicated(rownames(counts_removed_kidney)),] # count table avendo rimossi i geni che compaiono due volte

## in brain
vb <- rep(F,40903)
for (i in 1: 40903) {
  vb[i] = (as.vector(rowRanges(filtered_brain)@seqnames[i]  %in% can_chr)) # maschera che filtra i geni
}
counts_removed_brain<-counts_brain_filter[vb,] # count table con solo i geni annotati su cromosomi canonici
length(unique(rownames(counts_removed_brain))) # there are 36141 uniquely named genes
dim(counts_removed_brain) # while there are 36179 rows in the count table (ci sono 38 geni in piu comunque)
dupl_genes_b<-counts_removed_brain[duplicated(rownames(counts_removed_brain)),] # duplicated genes have counts = to 0
not_dupl_brain<- counts_removed_brain[!duplicated(rownames(counts_removed_brain)),] # count table avendo rimossi i geni che compaiono due volte

### in colon
vc <- rep(F,40903)
for (i in 1: 40903) {
  vc[i] = (as.vector(rowRanges(filtered_colon)@seqnames[i]  %in% can_chr)) # maschera che filtra i geni
}
counts_removed_colon<-counts_colon_filter[vc,] # count table con solo i geni annotati su cromosomi canonici
length(unique(rownames(counts_removed_colon))) # there are 36141 uniquely named genes
dim(counts_removed_colon) # while there are 36179 rows in the count table (ci sono 38 geni in piu comunque)
dupl_genes_c<-counts_removed_colon[duplicated(rownames(counts_removed_colon)),] # duplicated genes have counts = to 0
not_dupl_colon<- counts_removed_colon[!duplicated(rownames(counts_removed_colon)),] # count table avendo rimossi i geni che compaiono due volte
##################################################################################
#ANALISI DATI
library(edgeR)
all_counts<-cbind(not_dupl_brain,not_dupl_colon,not_dupl_kidney) # colonne gia nominate con acronimo samples e geni chiamati con nomenclatura tradizionale
####### analisi colon specific
all_counts_colon<-not_dupl_colon
y_all_colon<-DGEList(counts=all_counts_colon)
group_c <- as.factor(c("Sigmoideo","traverso","traverso")) 
y_all_colon$samples$group<-group_c
keep.exprs <- filterByExpr(y_all_colon, group=group_c)
y_all_colon <- y_all_colon[keep.exprs,, keep.lib.sizes=FALSE]
y_all_colon<- calcNormFactors(y_all_colon, method = "TMM")
design <- model.matrix(~0+group, data=y_all_colon$samples)
y_col <- estimateDisp(y_all_colon, design)
colnames(design) <- levels(y_all_colon$samples$group)
fit_col <- glmQLFit(y_col, design)
#traverso (top) vs sigmoideo (bottom)
qlfTS <- glmQLFTest(fit_col, contrast=c(-1,1))
#sigmoideo (top) vs traverso (bottom)
qlfST<- glmQLFTest(fit_col, contrast=c(1,-1))
resultsTS <- topTags(qlfTS, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
upregulatedTS<-resultsTS[resultsTS$table$logFC>0,]
upregulatedsigmoideo<-resultsTS[resultsTS$table$logFC<=0,]
head(resultsTS)
head(resultsST)
write.table(resultsTS,"resultsTS.csv",sep=",",col.names=NA,row.names=T)
write.table(upregulatedsigmoideo,"upregulatedsigmoideo.csv",sep=",",col.names=NA,row.names=T)
#Brain
size_all <- colSums(all_counts)# library size
y_all<-DGEList(counts=all_counts)
group <- as.factor(c("B","B","B","C","C","C","K","K","K")) # assegno come gruppo l'organo
y_all$samples$group<-group
# it is advisable to remove altogether all genes with low or zero counts. edgeR does it by considering also library size,
#removing also genes with “near zero” counts.
y_all$samples$rin <- as.factor(c(colData(rse_brain_selected)$gtex.smrin,colData(rse_colon_selected)$gtex.smrin,colData(rse_kidney_selected)$gtex.smrin))
y_all$samples$slice <- as.factor(c(colData(rse_brain_selected)$gtex.smtsd,colData(rse_colon_selected)$gtex.smtsd,colData(rse_kidney_selected)$gtex.smtsd))
y_all$samples$sex <- as.factor(c(colData(rse_brain_selected)$gtex.sex,colData(rse_colon_selected)$gtex.sex,colData(rse_kidney_selected)$gtex.sex))
y_all$samples$age <- as.factor(c(colData(rse_brain_selected)$gtex.age,colData(rse_colon_selected)$gtex.age,colData(rse_kidney_selected)$gtex.age))

y_all$samples$rRNA <- as.factor(c(colData(rse_brain_selected)$gtex.smrrnart,colData(rse_colon_selected)$gtex.smrrnart,colData(rse_kidney_selected)$gtex.smrrnart))

y_all$samples$mapped <- as.factor(c(colData(rse_brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both", colData(rse_colon_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",colData(rse_kidney_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))

y_all$samples$chrm <- as.factor(c(colData(rse_brain_selected)$"recount_qc.aligned_reads%.chrm", colData(rse_colon_selected)$"recount_qc.aligned_reads%.chrm",colData(rse_kidney_selected)$"recount_qc.aligned_reads%.chrm"))

keep.exprs <- filterByExpr(y_all, group=group)
y_all <- y_all[keep.exprs,, keep.lib.sizes=FALSE] # rimossi 14314 geni 
dim(y_all) # rimango con 21827 geni
# Let us extract and store in a vector the log of the counts per million before normalization with the “cpm”
logcpm_before <- cpm(y_all, log=TRUE)
boxplot(logcpm_before,main="log(CPM) before normalization",ylab = "log(CPM)",col=c("red","red","red","green","green","green","blue","blue","blue")) # graphical inspection
#text(x = 1:length(logcpm_before), y = par("usr")[3] - 0.2, labels = labels, xpd = TRUE, srt = 45, adj = c(1, 0.5), pos = 2)
#We can now normalize the counts (TMM normalization)
y_all <- calcNormFactors(y_all, method = "TMM")
norm_factors<-y_all$samples$norm.factors
names<-c("b_13", "b_14", "b_15", "c_14", "c_15", "c_16", "k_17", "k_18", "k_20")
mat<-as.data.frame(cbind(names,round(norm_factors,digits=2)))
colnames(mat)<-c("samples","norm_factors")

#Distribution of (normalized) log-cpm values across samples
logcpm <- cpm(y_all, log=TRUE)
x11()
boxplot(logcpm,main="log(CPM) after normalization",ylab = "log(CPM)",col=c("red","red","red","green","green","green","blue","blue","blue"))
## design of linear model(Brain as base_class) (chiedi al prof come modificarla e qual è il modello piu adatto ai miei tessuti)
design <- model.matrix(~0+group, data=y_all$samples)
colnames(design) <- levels(y_all$samples$group)

x11()
plotMDS(logcpm, labels=group) # pcoA per i miei campioni, curioso che clusterizzino tutti molto bene tranne un campione di colon
plotMDS(logcpm,labels=y_all$samples$slice,cex=1.25,main="")
title("Clustering Based on Organ Portion", col.main = "blue", cex.main = 2)
plotMDS(logcpm,labels=y_all$samples$age,cex=1.5,main="")
title("Clustering Based on Age", col.main = "blue", cex.main = 2)
plotMDS(logcpm,labels=y_all$samples$sex,cex=1.5,main="") # 2 sono femmine
title("Clustering Based on Sex", col.main = "blue", cex.main = 2)
y_all$samples$sex<-ifelse(y_all$samples$sex==2,"F","M")



#We now estimate the NB dispersion and plot the BCV - common and gene-specific.
y <- estimateDisp(y_all, design)
plotBCV(y)
y$common.dispersion # red line

## fitting the model
fit <- glmQLFit(y, design)
#colon (top) vs brain (bottom)
qlfCB <- glmQLFTest(fit, contrast=c(-1,1,0))
#kidney (top) vs brain (bottom)
qlfKB <- glmQLFTest(fit, contrast=c(-1,0,1))
#kidney (top) vs colon (bottom)
qlfKC <- glmQLFTest(fit, contrast=c(0,-1,1))

resultsCB <- topTags(qlfCB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
upregulatedCB<-resultsCB[resultsCB$table$logFC>0,]
upregulated_genesCB<-rownames(upregulatedCB[upregulatedCB$table$FDR<0.01,])
downregulatedCB<-resultsCB[resultsCB$table$logFC<=0,]
downregulated_genesCB<-rownames(downregulatedCB[upregulatedCB$table$FDR<0.01,])
write.table(resultsCB,"resultsCB.csv",sep=",",col.names=NA,row.names=T)
write.table(upregulatedCB,"upregulatedCB.csv",sep=",",col.names=NA,row.names=T)
write.table(downregulatedCB,"downregulatedCB.csv",sep=",",col.names=NA,row.names=T)


summary(decideTests(qlfCB, p.value=0.05, lfc=0)) # p.value intende quello del FDR
summary(decideTests(qlfCB, p.value=0.05, lfc=1))
summary(decideTests(qlfCB, p.value=0.01, lfc=0))# choosing different will cause having more or less genes in the table
summary(decideTests(qlfCB, p.value=0.01, lfc=1))# having FDR threshold at 0.01 result in 669 genes downregulated in colon and 265 upregolati
##### kidney(top) vs brain(bottom)
resultsKB <- topTags(qlfKB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
upregulatedKB<-resultsKB[resultsKB$table$logFC>0,]
upregulated_genesKB<-rownames(upregulatedKB[upregulatedKB$table$FDR<0.0075,])
downregulatedKB<-resultsKB[resultsKB$table$logFC<=0,]
downregulated_genesKB<-rownames(downregulatedKB[downregulatedKB$table$FDR<0.0075,])
write.table(resultsKB,"resultsKB.csv",sep=",",col.names=NA,row.names=T)
write.table(upregulatedKB,"upregulatedKB.csv",sep=",",col.names=NA,row.names=T)
write.table(downregulatedKB,"downregulatedKB.csv",sep=",",col.names=NA,row.names=T)
summary(decideTests(qlfKB, p.value=0.05, lfc=0))
summary(decideTests(qlfKB, p.value=0.05, lfc=1))
summary(decideTests(qlfKB, p.value=0.01, lfc=0))# choosing different will cause having more or less genes in the table
summary(decideTests(qlfKB, p.value=0.0075, lfc=0)) # with FDR < 0.0075  i have 763 genes downregulated in kidney and 562 upre
#result for kidney vs colon
resultsKC <- topTags(qlfKC, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
upregulatedKC<-resultsKC[resultsKC$table$logFC>0,]
upregulated_genesKC<-rownames(upregulatedKC[upregulatedKC$table$FDR<0.05,])
downregulatedKC<-resultsKC[resultsKC$table$logFC<=0,]
downregulated_genesKC<-rownames(downregulatedKC[downregulatedKC$table$FDR<0.05,])
write.table(resultsKC,"resultsKC.csv",sep=",",col.names=NA,row.names=T)
write.table(upregulatedKC,"upregulatedKC.csv",sep=",",col.names=NA,row.names=T)
write.table(downregulatedKC,"downregulatedKC.csv",sep=",",col.names=NA,row.names=T)
summary(decideTests(qlfKC, p.value=0.05, lfc=0)) # scelgo questa threshold perche con FDR di 0.01 ho solo 29 geni upregolati
summary(decideTests(qlfKC, p.value=0.05, lfc=1))
summary(decideTests(qlfKC, p.value=0.01, lfc=0))# choosing different will cause having more or less genes in the table
summary(decideTests(qlfKC, p.value=0.01, lfc=1))
# DE genes specific of a class
up_kidney<-intersect(upregulated_genesKB,upregulated_genesKC)
up_brain<-intersect(downregulated_genesCB,downregulated_genesKB)
up_colon<-intersect(upregulated_genesCB,downregulated_genesKC)
which(rowData(rse_brain)$gene_name == "RHBG") #upregulated gene in kidney
# 10 most DE genes 
# table with all DE genes
#deg.kvsb <- topTags(qlf.kvsb, n=20000, adjust.method = "BH", sort.by = "PValue", p.value = a)$table
#up.genes.kvsb <- row.names(deg.kvsb[deg.kvsb$logFC > fc,]) # 773 genes higher expressed in kidney with respect to brain
#down.genes.kvsb <- row.names(deg.kvsb[deg.kvsb$logFC < fc,]) # 1015 genes lower expression in kidney with respect to brain

#or example, gene “PCDH10” is on the top of the DE lists, over-expressed in brain with respect to the other two tissues.
#Let us plot the distribution of expression across the three complete datasets, as TPM:
assays(rse_brain)$TPM <- recount::getTPM(rse_brain)
assays(rse_colon)$TPM <- recount::getTPM(rse_colon)
assays(rse_kidney)$TPM <- recount::getTPM(rse_kidney)
which(rowData(rse_brain)$gene_name == "RHBG")
#boxplot(assays(rse_brain_selected)$TPM[9091,]) #
#TPM of gene RHBG in brain
brain_count<-recount::getTPM(rse_brain)[9091,]
brain_count<-as.numeric(brain_count)
#TPM of gene RHBG in colon
colon_count<-recount::getTPM(rse_colon)[9091,]
colon_count<-as.numeric(colon_count)
#TPM of gene RHBG in kidney
kidney_count<-recount::getTPM(rse_kidney)[9091,]
kidney_count<-as.numeric(kidney_count)
x11()
boxplot(brain_count, colon_count, kidney_count, 
        names = c("Brain", "Colon", "Kidney"),
        ylab = "TPM",
        main = "TPM of gene RHBG in my samples",outline=F,
        col = c("red", "blue", "green"))
stars_x <- 3  
stars_y <- max(kidney_count)- 28 
par(xpd = TRUE)  # Enable drawing outside the plot area
text(stars_x, stars_y, labels = "***", cex = 2)

df<-data.frame(group=c(rep("B",2931),rep("C",822),rep("K",98)),y=c(brain_count,colon_count,kidney_count)) 

x<-aov(y~factor(group),data=df) # anova test to check stastical difference between TPM of geneRHBG in kidney,brain,colon
qq
