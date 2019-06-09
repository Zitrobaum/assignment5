# install standard R packages:
install.packages("ggplot2") # for plotting
install.packages("reshape") # for Venn diagrams
install.packages("Vennerable", repos="http://R-Forge.R-project.org") # for Venn diagrams

# vector of Bioconductor package names:
pkgs <- c(
  "rtracklayer",
  "TFBSTools",
  "JASPAR2014",
  "BSgenome.Hsapiens.UCSC.hg18",
  "TxDb.Hsapiens.UCSC.hg18.knownGene"
)

# install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# install packages
source("http://bioconductor.org/biocLite.R")
biocLite(pkgs, suppressUpdates=TRUE)

# set two lists with the names of the Biocodnuctor and the CRAN packages
bioc.libs <- c("Biobase", "affy", "annotate", "hgu95av2.db", "hgu95av2cdf",
               "estrogen", "limma", "genefilter", "vsn", "GEOquery", "annaffy",
               "GO.db","KEGG.db")
cran.libs <- c("RColorBrewer", "R.utils", "gplots")

library(BiocInstaller)
biocLite(bioc.libs)
biocLite(cran.libs)

# install the package "ChIPpeakAnno" 
biocLite("ChIPpeakAnno")

# read the input data directly from the URL
FOXA1.df <- read.table("http://www.carroll-lab.org.uk/FreshFiles/Data/Data_Sheet_3/MCF7_FOXA1%20binding.bed",
                       header=TRUE)
# show the first six lines of the data.frame
head(FOXA1.df)

# load the GenomicRanges package
require(GenomicRanges) 
FOXA1 <- GRanges(
  FOXA1.df$chr,
  IRanges(FOXA1.df$star, FOXA1.df$end),
  strand="*"
)
# we can add more data to each peak subsequently
names(FOXA1) <- paste("FOXA1_peak", 1:nrow(FOXA1.df), sep="_")
score(FOXA1) <- FOXA1.df[,"X.10.log10.pvalue."]

# show the first and last lines of the GRanges object
FOXA1

# number of peaks
length(FOXA1)

# maximum length of peak
max(width(FOXA1))

# average value
mean(width(FOXA1))

# median value
median(width(FOXA1))

# removing all peaks larger 2000
FOXA1 <- FOXA1[width(FOXA1) <= 2000]

# create histogramm of peak size
hist(width(FOXA1), xlab="Peak size")

# get the -log_10 transformed p-values from the score column
mlPvals <- score(FOXA1)
hist(mlPvals, xlab="-log_10(p-value)", col="gray")

# load the rtracklayer package
require(rtracklayer) 
ER <- import.bed("http://www.carroll-lab.org.uk/FreshFiles/Data/Data_Sheet_3/MCF7_ER_binding.bed")

# assign scores by converting the 'name' feald to type numeric
score(ER) <- as.numeric(ER$name)
# overwrite the name column
ER$name <- paste("ER_peaks", 1:length(ER), sep="_")
# use the names() function
names(ER) <- ER$name
ER

# barplot with number of peaks
bp <- barplot(c(length(ER), length(FOXA1)), names=c("ER", "FOXA1"))
# add actual values as text lables to the plot
text(bp, c(length(ER), length(FOXA1)), labels=c(length(ER), length(FOXA1)), pos=1)

# overlaps
findOverlaps(ER,FOXA1)

# get subsets of binding sites
ER.uniq <- setdiff(ER, FOXA1)
FOXA1.uniq <- setdiff(FOXA1, ER)

# plot overlap of binding sites as Venn diagram
require(Vennerable)
# build objects with the numbers of sites in the subsets
venn <- Venn(SetNames=c("ER", "FOXA1"),
             Weight=c(
               '10'=length(ER.uniq),
               '11'=length(ovl),
               '01'=length(FOXA1.uniq)
             )
)
# plot Venn Diagram
plot(venn)

# load gene annotation from UCSC for human genome build hg18
require(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene # just a shortcut

require(ChIPpeakAnno) # laod package to annotate peaks with genes
# calculate the overlap with features
ER.features <- assignChromosomeRegion(ER, TxDb=txdb, nucleotideLevel=FALSE)
# show the results
ER.features

# make pie-chart
pie(ER.features$percentage)

# make barplot
bp <- barplot(ER.features$percentage, ylab="%")
text(bp, ER.features$percentage, signif(ER.features$percentage, 4),
     pos=1)

# get all genes from the data base as GRanges object
genes <- genes(txdb)
# take the region around the gene start as promoter
prom <- promoters(genes, upstream=2000, downstream=200)
prom

# get only those ER peaks that overlap a promoter region
ERatProm <- subsetByOverlaps(ER, prom)
# subset size
length(ERatProm)

# percent of all ER peaks
length(ERatProm) / length(ER) * 100

# We search for overlap between ER peaks and promoters
ERatProm.Hits = findOverlaps(ER, prom)

ERprom = genes[subjectHits(ERatProm.Hits)]
# take only unique ids
gene.ids <- unique(names(ERprom))
# write names to an output file
write.table(gene.ids, file="ER_regulated_genes.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)

# get TSS as only the start coordinate of genes.
tss <- resize(genes, width=1, fix="start")
# calculate the distances from peaks to tss
d <- distanceToNearest(ER, tss)
# show object d
d

# show the metacolumns as DataFrame object
mcols(d)
## DataFrame with 12999 rows and 1 column
dist <- mcols(d)[,1]
# get the average distance in kb
mean(dist) * 10^-3

# subset hits object by distance
close.Hits <- d[dist <= 10000,]
# show the subset
close.Hits

# get the indeces of genes
ER.genes <- genes[subjectHits(close.Hits)]

# extract the vector of names from the GRanges object
gene.ids <- names(ER.genes)
# take only unique ids
gene.ids <- unique(gene.ids)
# write names to an output file
write.table(gene.ids, file="ER_regulated_genes.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)

export.bed(ER, "ER_peaks.bed")
export.bed(FOXA1, "FOXA1_peaks.bed")

