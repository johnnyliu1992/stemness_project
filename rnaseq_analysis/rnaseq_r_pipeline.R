library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(NMF)


setwd('/Users/johnny/Documents/17p_stemness/rnaseq_analysis')

# Read the data into R
seqdata <- read.delim("rnaseq_result_combined_symbol.csv", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("sample_info.csv", stringsAsFactors = TRUE)

# Remove first two columns from seqdata
countdata <- seqdata[,-(1:1)]

# Store Symbol as rownames
rownames(countdata) <- seqdata[,1]

y <- DGEList(countdata)


#ann <- select(org.Mm.eg.db,keys=rownames(y$counts),columns=c("SYMBOL","ENTREZID","GENENAME"))
#ann <- select(org.Mm.eg.db,keys=c("19888","27395"),columns=c("SYMBOL","ENTREZID","GENENAME"))
y$genes <- data.frame("symbol"=seqdata[,1])





group <- paste(sampleinfo$Cell,sampleinfo$Treatment,sep=".")
# Convert to factor
group <- factor(group)

y$samples$group <- group
y$samples


# Obtain CPMs
#myCPM <- cpm(countdata)
# Have a look at the output
#head(myCPM)

#QC filter out genes
keep <- filterByExpr(y)
y <- y[keep,]


# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")


########  MDS plot  ########
levels(sampleinfo$Status)
col.status <- c("black","dodgerblue","deeppink4","deeppink")[y$samples$group]
col.status

png("MDS_plot.png", width=1200, height=1200, res=200)
plotMDS(y,col=col.status)
legend("topleft",fill=c("black","dodgerblue","deeppink4","deeppink"),legend=levels(y$samples$group),cex=0.8)
title("Status")
dev.off()

########  top variance gene heatmap plot  ########

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("black","dodgerblue","deeppink4","deeppink")[y$samples$group]

# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")


########  top DE gene heatmap plot  ########

#before treat
select_var <- c('RTN1','GLUL','GPNMB','ETV5','TLL1','PSG4',
                'NFIB','CLCN4','RPS6KA2','SIM2','CABLES2',
                'ADAMTSL1','WWC3','CACNA2D2','HIP1','CYFIP2',
                'LINC01234','EFEMP1','GOLGA2','FEZ1','DHRS2','NT5E',
                'GALNT5','DMTN','ADRB2','CCN2','ICAM1','CYP2U1','RRAGD','YPEL5',
                'ADARB1','SPRY4','PXDC1','TFDP2','SIRPA','SUMO3','ARFGEF3','MIR4435-2HG',
                'AL391422.4','LAMA1','SFTA1P','NOVA2','CD177','BACE1','SCRN1','TMEM200A',
                'B4GALT4','SHANK1','SERPINB7','HUNK','LANCL2','BCAM','FAM120C','LTBP2',
                'SHISA3','OSBP2','SLC25A16','NCK2','FAXC','MPDZ','SPTAN1','H2AZ2','PCDH10',
                'AFF3','CCDC144NL-AS1','DCLK1','YIPF4','STXBP1','RYR2','MMP3','CCDC69','DNM1',
                'GLE1','INSYN2B','MAP3K7CL','FBLN5','LARP4B','INPP5D','BNC2','IGFBP7','TAGLN',
                'LMX1B','TTBK2','RALGPS1','RAI2','TTC39B','DRGX','NPC2','AJAP1','TSPYL4',
                'AL158206.1','ATP6V1D','SSTR1','ESM1','BEND7','FN3K','PLIN2','ACER2','LONRF2',
                'MIR31HG')[0:50]
#diff
select_var <- c('ATF3','BTG2','SULF2','CXCL8','TP53I11','RASSF4','MR1','TRAF1','TNFAIP3','DEPP1','RNF144B','P3H2','COL27A1','GDF15','IL1A','GPR87','XDH','CSF1','LRATD2','MDM2','TP53I3','CRISPLD2','BHLHE40','CFAP251','NECTIN4','CMBL','SAT1','BCL7A','GDNF','PLCXD2','ACTA2','CDKN1A','DGKA','KRT15','DUSP10','SRGAP3','TNFRSF10A','GREB1','TMEM40','PLAT','S1PR3','FAM13C','MCC','ZNF385A','AC245041.2','JUN','UNC13A','DRAM1','PAPLN','ABHD4','GADD45A','HAS3','NLRP1','MYOSLID','CDH10','ZFP36L2','CES2','OAS2','FBLIM1','MX2','C2orf88','PTPRU','IL1B','CYFIP2','ZNF79','PADI3','ANXA8','PLCL2','ELAPOR1','PDE4A','IRF1','FAS','NHLH2','RAP2B','DDB2','GABBR2','POLH','STX6','PRICKLE2','NTN1','MAST4','RIPK4','IFFO2','PPP1R14C','BTBD11','TRIM8','AL158206.1','HES2','AREG','ACER2','TGFA','KANK3','LAMC3','ACHE','JAG1','L3MBTL2-AS1','PRKAB1','PURPL','PVT1',
                'TRIM26')[0:30]
#stemness markers:: 15 markers
select_var <- c('TWIST1', 'NOTCH1', 'KLF4', 'CD44', 'ABCG2', 'CD34', 'EZH2', 'HIF1A', 'PROM1', 'MYC', 'EPAS1', 'BMI1', 'KDM5B', 'NES', 'CTNNB1')


#stemness markers:: 41 markers ==> sort by logFC
select_var <- c('WNT4','MDM2','PROM1','HHAT','CAMK2B','KLF4','NOTCH1','EPAS1','GJA1','PLCB1','ITGA6','ABCG2','GPR161','MYC','PRKACA','MAPK8','BMI1','DVL1','CTNNB1','SOX9','HIF1A','CD24','CD44','EZH2','RUNX2','TWIST1','PPP3CA','ROCK2','HHIP','LRP6','BCL2','EPCAM','KDM5B','FSTL1','NFATC1','NES','CDON','CD34','DHH','PTCH1','PRKCA')

#stemness markers:: 30 markers ==> sort by logFC
select_var <- c('WNT4','MDM2','PROM1','HHAT','CAMK2B','KLF4','NOTCH1','EPAS1','PLCB1','ITGA6','GPR161','MYC','PRKACA','MAPK8','HIF1A','CD44','EZH2','ROCK2','HHIP','LRP6','BCL2','EPCAM','KDM5B','FSTL1','NFATC1','NES','CDON','DHH','PTCH1','PRKCA')


# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("black","dodgerblue","deeppink4","deeppink")[y$samples$group]

# Plot the heatmap
heatmap.2(highly_variable_lcpm,
          col=rev(morecols(50)),
          trace="none", 
          main="Top100",
          ColSideColors=col.cell,
          scale="row")

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Specify colors
Cell = c("#EC2025", "#3953a3")
names(Cell) = c("cas9", "loss17p")
Treatment = c("#ff8084", "#a6a6a6")
names(Treatment) = c("Dox", "NoTreat")
ann_colors = list(Cell = Cell, Treatment = Treatment)

#reorder cols
highly_variable_lcpm <- highly_variable_lcpm[,c('TH_10_S101_L002','TH_11_S102_L002','TH_12_S103_L002','TH_4_S95_L002','TH_5_S96_L002','TH_6_S97_L002','TH_7_S98_L002','TH_8_S99_L002','TH_9_S100_L002','TH_1_S92_L002','TH_2_S93_L002','TH_3_S94_L002')]
#update sample info
sampleinfo_heatmap <- sampleinfo[colnames(highly_variable_lcpm),]

aheatmap(highly_variable_lcpm,
         col=rev(morecols(50)),
         main="Stemness Markers",
         annCol=sampleinfo_heatmap[, 2:3],
         annColors=ann_colors,
         labCol=c('10_S101','11_S102','12_S103','4_S95','5_S96','6_S97','7_S98','8_S99','9_S100','1_S92','2_S93','3_S94'),
         Colv=NA,
         Rowv=NA,
         scale="row",
         filename='stemness_markers(30)_in_diff.pdf'
         )

aheatmap(x = highly_variable_lcpm, scale = 'row', legend = F
         , main = 'Heatmap', filename = 'plot.pdf'
         , width = 10, height = 10) 



########  normalization  ########

y <- calcNormFactors(y)

########  DE  ########
# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)

par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

fit <- lmFit(v)
names(fit)


###########. BeforeTreat ############
cont.matrix <- makeContrasts(
  BeforeTreat=loss17p.NoTreat - cas9.NoTreat,
  levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

top.table <- topTable(fit.cont, sort.by = "P", n = Inf)
write.table(top.table, file = "DE_results_BeforeTreat.txt", row.names = F, sep = "\t", quote = F)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))




###########. Loss17pTreat ############
cont.matrix <- makeContrasts(
  Loss17pTreat=loss17p.Dox - loss17p.NoTreat,
  levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

top.table <- topTable(fit.cont, sort.by = "P", n = Inf)
write.table(top.table, file = "DE_results_Loss17pTreat.txt", row.names = F, sep = "\t", quote = F)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))




###########. Diff ############
cont.matrix <- makeContrasts(
  Diff=(loss17p.Dox - loss17p.NoTreat) - (cas9.Dox - cas9.NoTreat),
  levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)


top.table <- topTable(fit.cont, sort.by = "P", n = Inf)
write.table(top.table, file = "DE_results_Diff.txt", row.names = F, sep = "\t", quote = F)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))



###########. ........ ############
cont.matrix <- makeContrasts(
  BeforeTreat=cas9.NoTreat - loss17p.NoTreat,
  Loss17pTreat=loss17p.Dox - loss17p.NoTreat,
  Diff=(loss17p.Dox - loss17p.NoTreat) - (cas9.Dox - cas9.NoTreat),
  levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

#vacanoplot
volcanoplot(fit.cont,coef=1,style = "p-value",highlight=100,hl.col="red",names=fit.cont$genes$symbol, main="")


#####==========>>>> volcano plot <<<<<====== #######

## stemness markers
plot_data <- read.csv('DE_results_BeforeTreat_volcano.txt', sep = '\t', row.names = 1)
library(EnhancedVolcano)
#stemness markers
stemness_markers <- c('DHH','CAMK2B','WNT4','BCL2','CDON','HHAT','PTCH1','PRKACA','HHIP','MDM2','EPAS1','GPR161','MAPK8','FSTL1','KDM5B','NOTCH1','LRP6','PROM1')
nan_markers=c('BST1','NT5E','NMRK1','QPRT','ENPP3','SIRT1')

stemness_markers_plot <- c('DHH','CAMK2B','WNT4','BCL2','CDON','HHAT','PTCH1','PRKACA','HHIP','MDM2','EPAS1','GPR161','MAPK8','FSTL1','KDM5B','NOTCH1','LRP6','PROM1')

plotlabels <- c('DHH','CAMK2B','WNT4','BCL2','CDON','HHAT','PTCH1','PRKACA','HHIP','MDM2','EPAS1','GPR161','MAPK8','FSTL1','KDM5B','NOTCH1','LRP6','PROM1','BST1','NT5E','NMRK1','QPRT','ENPP3','SIRT1')

keyvals <- ifelse(plot_data$adj.P.Val > 5*10e-3, 'gray60',
         ifelse(plot_data$logFC < 0, ifelse(rownames(plot_data) %in% stemness_markers_plot, 'black', 'royalblue1'), 
                ifelse(rownames(plot_data) %in% stemness_markers_plot, 'black', 'red'))
  )
keyvals[is.na(keyvals)] <- 'red'
names(keyvals)[keyvals == 'red'] <- 'High'
names(keyvals)[keyvals == 'gray60'] <- 'Not Sig'
names(keyvals)[keyvals == 'royalblue1'] <- 'Low'
names(keyvals)[keyvals == 'black'] <- 'Stemness Markers'
#names(keyvals)[keyvals == 'green'] <- 'NANM Markers'

pdf("volcano_plot_beofretreat_stemness_marker.pdf", width = 10, height = 10)

EnhancedVolcano(plot_data,
  lab = rownames(plot_data),
  x = 'logFC',
  y = 'adj.P.Val', 
  ylim = c(0,13),
  title = 'Before Treatment', 
  selectLab = stemness_markers_plot,
  pCutoff = 5*10e-3,
  FCcutoff = 0, 
  colCustom = keyvals, 
  pointSize = 2.0,
  labSize = 3,
  #labCol = 'black',
  #labFace = 'bold',
  #boxedLabels = TRUE, 
  colAlpha = 0.8,
  legendPosition = 'right',
  legendLabSize = 14,
  legendIconSize = 4.0,
  #drawConnectors = TRUE,
  #widthConnectors = 1.0,
  #colConnectors = 'black', 
  gridlines.major = FALSE,
  gridlines.minor = FALSE
  
)
dev.off()




## NANM markers
plot_data <- read.csv('DE_results_BeforeTreat_volcano.txt', sep = '\t', row.names = 1)
library(EnhancedVolcano)
#stemness markers
stemness_markers <- c('DHH','CAMK2B','WNT4','BCL2','CDON','HHAT','PTCH1','PRKACA','HHIP','MDM2','EPAS1','GPR161','MAPK8','FSTL1','KDM5B','NOTCH1','LRP6','PROM1')
nan_markers=c('BST1','NT5E','NMRK1','QPRT','ENPP3','SIRT1')

stemness_markers_plot <- c('DHH','CAMK2B','WNT4','BCL2','CDON','HHAT','PTCH1','PRKACA','HHIP','MDM2','EPAS1','GPR161','MAPK8','FSTL1','KDM5B','NOTCH1','LRP6','PROM1')

plotlabels <- c('DHH','CAMK2B','WNT4','BCL2','CDON','HHAT','PTCH1','PRKACA','HHIP','MDM2','EPAS1','GPR161','MAPK8','FSTL1','KDM5B','NOTCH1','LRP6','PROM1','BST1','NT5E','NMRK1','QPRT','ENPP3','SIRT1')

keyvals <- ifelse(plot_data$adj.P.Val > 5*10e-3, 'gray60',
                  ifelse(plot_data$logFC < 0, ifelse(rownames(plot_data) %in% nan_markers, 'green', 'royalblue1'), 
                         ifelse(rownames(plot_data) %in% nan_markers, 'green', 'red'))
)
keyvals[is.na(keyvals)] <- 'red'
names(keyvals)[keyvals == 'red'] <- 'High'
names(keyvals)[keyvals == 'gray60'] <- 'Not Sig'
names(keyvals)[keyvals == 'royalblue1'] <- 'Low'
#names(keyvals)[keyvals == 'black'] <- 'Stemness Markers'
names(keyvals)[keyvals == 'green'] <- 'NANM Markers'

pdf("volcano_plot_beofretreat_nanm_marker.pdf", width = 10, height = 10)

EnhancedVolcano(plot_data,
                lab = rownames(plot_data),
                x = 'logFC',
                y = 'adj.P.Val', 
                ylim = c(0,13),
                title = 'Before Treatment', 
                selectLab = nan_markers,
                pCutoff = 5*10e-3,
                FCcutoff = 0, 
                colCustom = keyvals, 
                pointSize = 2.0,
                labSize = 3,
                #labCol = 'black',
                #labFace = 'bold',
                #boxedLabels = TRUE, 
                colAlpha = 0.8,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                #drawConnectors = TRUE,
                #widthConnectors = 1.0,
                #colConnectors = 'black', 
                gridlines.major = FALSE,
                gridlines.minor = FALSE
                
)
dev.off()















