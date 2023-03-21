setwd("/Users/johnny/Documents/17p_stemness")

#install packages
deps <- c("gelnet","dplyr","gdata","DT")
for(pkg in deps)  if (!pkg %in% installed.packages()) install.packages(pkg, dependencies = TRUE)
install.packages("xml2")

#install synapser
install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))

#install biomaRt for R 3.6
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library(gelnet)
library(dplyr)
library(biomaRt)
library(synapser)
synLogin(email='johnnyliu1992', password='Jiannan66095518')
library(data.table)

# Maps ENSEMBL IDs to HUGO
# Use srcType = "ensembl_gene_id" for Ensembl IDs
# Use srcType = "entrezgene" for Entrez IDs
genes2hugo <- function( v, srcType = "ensembl_gene_id" )
{
  ## Retrieve the EMSEMBL -> HUGO mapping
  ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl" )
  ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
  
  ## Make sure there was at least one mapping
  if( nrow(ID) < 1 ) top( "No IDs mapped successfully" )
  
  ## Drop empty duds
  j <- which( ID[,2] == "" )
  if( length(j) > 0 ) ID <- ID[-j,]
  stopifnot( all( ID[,1] %in% v ) )
  
  ID
}


# Load RNAseq data
synRNA <- synGet( "syn2701943", downloadLocation = NULL )
X <- read.delim( "c:/pancanstem_data/rnaseq_norm.tsv" ) %>%
  tibble::column_to_rownames( "tracking_id" ) %>% as.matrix
X[1:3,1:3]


synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
Y_file <- fread(synMeta[["filepath"]],select=c(3,4))
Y <- Y_file %>%
  mutate( UID = gsub("-", ".", UID) ) %>%
  tibble::column_to_rownames( "UID" )
Y[1:4,]

# Retrieve the labels from the metadata
y <- Y[colnames(X),]

names(y) <- colnames(X)
# Fix the missing labels by hand
y["SC11.014BEB.133.5.6.11"] <- "EB"
y["SC12.039ECTO.420.436.92.16"] <- "ECTO"

## Drop the splice form ID from the gene names
v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()

rownames(X) <- v
head(y)

# Map Ensembl IDs to HUGO
V <- genes2hugo( rownames(X) )
head(V)

X <- X[V[,1],]
rownames(X) <- V[,2]
X[1:3,1:3]

#optional
if(!is.null(fnGenes)){
  vGenes <- read.delim( fnGenes, header=FALSE ) %>% as.matrix() %>% drop()
  VE <- genes2hugo( vGenes, "entrezgene" )
  X <- X[intersect( rownames(X), VE[,2] ),]
}

m <- apply( X, 1, mean )
X <- X - m
X[1:3,1:3]

j <- which( y == "SC" )
X.tr <- X[,j]
X.tr[1:3,1:3]

X.bk <- X[,-j]
X.bk[1:3,1:3]

mm <- gelnet( t(X.tr), NULL, 0, 1 )

write.table(mm$w, file = 'output_signature.txt', sep = "\t", quote = FALSE, col.names = FALSE)
















#predict=================================>>>>>=================>>>>>>>>>>>>

setwd("/Users/johnny/Documents/17p_stemness")
w <- read.delim('output_signature.txt', header = FALSE, row.names = 1 ) %>% as.matrix() %>% drop()

#s <- synGet( "syn4976369", downloadLocation = "/data/pancan" )

# Auxiliary function: Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )

file_list=list.files("/Users/johnny/Documents/17p_stemness/exp_matrix_divided_symbol")

for (curr_file in file_list) {
X <- read.delim( paste("exp_matrix_divided_symbol/", curr_file, sep=""), as.is=TRUE, check.names=FALSE ) %>%    ## Read the raw values
  #filter( !grepl( "\\?", id ) ) %>%     ## Drop genes with no mapping to HUGO
  #mutate( id = genes2hugo( id ) ) %>%       ## Clip gene ids to HUGO
  filter( id %in% names(w) )   %>% 
  filter( id != "MIA2" )

j <- grep( "SLC35E2", X[,1] )
if( length(j) > 1 ) X <- X[-j[-1],]

rownames(X) <- NULL
X <- X %>% tibble::column_to_rownames( "id" ) %>% as.matrix()
X[1:3,1:3]

stopifnot( all( rownames(X) %in% names(w) ) )
w <- w[ rownames(X) ]
w[1:5]

s <- apply( X, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
s[1:5]

s <- s - min(s)
s <- s / max(s)
s[1:5]

#file_name=unlist(strsplit( curr_file, "." ))[1]
write.table( cbind(s), file=paste("final_stemness/", unlist(strsplit( curr_file, ".g" ))[1], ".txt", sep=""), sep="\t", quote=FALSE, col.names=FALSE )

}



#singel cell test
setwd("/Users/johnny/Documents/singlecell")
library(magrittr)
w <- read.delim('output_signature.txt', header = FALSE, row.names = 1 ) %>% as.matrix() %>% drop()

#s <- synGet( "syn4976369", downloadLocation = "/data/pancan" )

# Auxiliary function: Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )



X <- read.delim( paste("", "inputed_exp_sciv_10k_final.csv", sep=""), as.is=TRUE, check.names=FALSE ) %>%    ## Read the raw values
  #filter( !grepl( "\\?", id ) ) %>%     ## Drop genes with no mapping to HUGO
  #mutate( id = genes2hugo( id ) ) %>%       ## Clip gene ids to HUGO
  dplyr::filter( id %in% names(w) )   %>% 
  dplyr::filter( id != "MIA2" )

j <- grep( "SLC35E2", X[,1] )
if( length(j) > 1 ) X <- X[-j[-1],]

rownames(X) <- NULL
X <- X %>% tibble::column_to_rownames( "id" ) %>% as.matrix()
X[1:3,1:3]

stopifnot( all( rownames(X) %in% names(w) ) )
w <- w[ rownames(X) ]
w[1:5]

s <- apply( X, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
s[1:5]

s <- s - min(s)
s <- s / max(s)
s[1:5]

#file_name=unlist(strsplit( curr_file, "." ))[1]
write.table( cbind(s), file=paste("", unlist(strsplit( "inputed_exp_sciv_10k_final_stemness", ".g" ))[1], ".txt", sep=""), sep="\t", quote=FALSE, col.names=FALSE )











