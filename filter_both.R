library(dartRverse)


nt <- readRDS("./data/NTCrocs.rds")
ql <- readRDS("./data/QlandCrocs.rds")



ntscaf <- nt@other$loc.metrics$Chrom_Crocodile_ncbi_v2
qldcaf <- ql@other$loc.metrics$Chrom_Crocodile_ncbi_v3


ntchrsnp <- paste(nt@other$loc.metrics$Chrom_Crocodile_ncbi_v2, nt@other$loc.metrics$Pos_Crocodile_ncbi_v2, sep = ":")
qldcrsnp <- paste(ql@other$loc.metrics$Chrom_Crocodile_ncbi_v3, ql@other$loc.metrics$ChromPos_Crocodile_ncbi_v3, sep = ":")


table(ntscaf %in% qldcaf, useNA = "always")
table(qldcaf %in% ntscaf)



nt2 <- gl.filter.callrate(nt, method = "loc",0.9)
nt3 <- gl.filter.callrate(nt2, method = "ind",0.9)


ql2 <- gl.filter.callrate(ql, method = "loc",0.9)
ql3 <- gl.filter.callrate(ql2, method = "ind",0.9)


#filter for length of trimmed sequence==69

nt4 <- nt3[, nchar(as.character(nt3@other$loc.metrics$TrimmedSequence))==69]
ql4 <- ql3[, nchar(as.character(ql3@other$loc.metrics$TrimmedSequence))==69]

#saveRDS(nt4, "d:/bernd/Projects/Elise_crocs/data/NTCrocs_filtered.rds")
#saveRDS(ql4, "d:/bernd/Projects/Elise_crocs/data/QlandCrocs_filtered.rds")

nt4 <- readRDS("./data/NTCrocs_filtered.rds")
ql4 <- readRDS("./data/QlandCrocs_filtered.rds")

#run the blast against one genome


ntchrsnp <- paste(nt@other$loc.metrics$Chrom_Crocodile_ncbi_v2, nt@other$loc.metrics$Pos_Crocodile_ncbi_v2, sep = ":")
qldcrsnp <- paste(ql@other$loc.metrics$Chrom_Crocodile_ncbi_v3, ql@other$loc.metrics$ChromPos_Crocodile_ncbi_v3, sep = ":")





#ql4b <- gl.blast(ql4[,], ref_genome = "d:/bernd/Projects/Elise_crocs/crocs2/crocs.fasta", task = "megablast", number_of_threads = 8 )
#saveRDS(ql4b, "d:/bernd/Projects/Elise_crocs/crocs2/QlandCrocs_blasted.rds")
ql4b <- readRDS("d:/bernd/Projects/Elise_crocs/crocs2/QlandCrocs_blasted.rds")

#nt4b <- gl.blast(nt4[,], ref_genome = "d:/bernd/Projects/Elise_crocs/crocs2/crocs.fasta", task = "megablast", number_of_threads = 8 )
#saveRDS(nt4b, "d:/bernd/Projects/Elise_crocs/crocs2/NTCrocs_blasted.rds")
nt4b <- readRDS("d:/bernd/Projects/Elise_crocs/crocs2/NTCrocs_blasted.rds")


nt4b$snpp <- paste(nt4b$sacc, nt4b$sstart, sep = "_")
ql4b$snpp <- paste(ql4b$sacc, ql4b$sstart, sep = "_")

#take only the unique ones [though not 100%]
nt4b <- nt4b[!duplicated(nt4b$snpp),]
ql4b <- ql4b[!duplicated(ql4b$snpp),]



#check how often and if sequences are matched
yy <- sapply(1:nrow(nt4b),  function(x) which(nt4b$qseq[x]==nt4@other$loc.metrics$TrimmedSequence)[1])



###Todo
#1. check the duplicates in the blasts (bs) and nt4, ql4. Are they all genuine, snps are different positiion? What about triple hits?
#2. find out if the ones with indels are safeable
#3. filter snsp in nt4 that are good hits
#4. same for ql4
#5. make sure loc.all is the same in both (or potentially reverse e.g. AT->TA one of them)
#6. rbind snps which are in both datasets

#6




