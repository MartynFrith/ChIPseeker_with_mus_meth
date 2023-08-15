library(tidyverse) #load tidyverse
library(ChIPseeker) #load ChIPseeker
library(TxDb.Mmusculus.UCSC.mm10.knownGene) #txdb data
library(EnsDb.Mmusculus.v79) #annotation data 
library(GenomicRanges) #load GenomicRanges 

setwd("User/ChIPseeker_directory") #set working directory to ChIPseeker folder


mouse_manifest <- read_tsv("mouse_manifest_data.tsv") #read in your mouse manifest data
mouse_manifest_cleaned <- drop_na(mouse_manifest, c("CpG_beg", "CpG_end"))  #drop rows with NAs in start coordinate or end coordinate columns


mouse_manifest_cleaned_GRanges <- makeGRangesFromDataFrame(mouse_manifest_cleaned, keep.extra.columns = F, ignore.strand = F,                #create GRanges object using loaded data, drop extra columns, don't specify strand
                                       seqnames.field = "CpG_chrm", start.field = "CpG_beg", end.field = "CpG_end") #specify the columns with the chromosome names, start coordinates, and end coordinates


mouse_manifest_annotated <- annotatePeak(peak = mouse_manifest_cleaned_GRanges, tssRegion = c(-1500, 1500),                                       
                         TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db", verbose = T)      
#use annotatPeak to run annotation of peaks in manifest data. Must give a TxDB and annotation database. Optional arguments include tssRegions (default is -3000, 3000 for 3000 base pairs 
#down- and upstream), level for specifiying "transcript" or "gene" level annotation (default is transcript), genomicAnnotationPriority for specifying order of prioritization for annotations
#(default is "Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"), addFlankGeneInfo to add flankning gene information to eack peak, and various ignore arguments to allow 
# conditions like overlap of TSS regions with peaks to be ignored


mouse_manifest_annotated <- as_tibble(mouse_manifest_annotated) #transforming the CSanno object created by annotatePeak to a tibble for ease of use in future
mouse_manifest_annotated <- arrange(mouse_manifest_annotated, seqnames) #arranging by seq names to get annotations in order of chromosome number   
