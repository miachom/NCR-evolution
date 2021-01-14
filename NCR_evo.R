## NCR evolution analysis in comparison to two control gene families
setwd("/Users/ma2292/Documents/warwick/NCR_dist_matrix/analysis")
library("seqinr")
AtPPR <- read.fasta(file = "AtPPR.fasta")
nt_length <- getLength(AtPPR)
hist(nt_length, xlab = "CDS length", ylab = "Frequency", breaks = 100, main = "AtPPR gene family- expanded and stopped")
abline(col = "red", v=2700)
abline(col = "red", v=1200)
f1 <- AtPPR[which(getLength(AtPPR) >= 1200)]
f2 <- AtPPR[which(getLength(AtPPR) <= 2700)]
length(intersect(f1,f2)) # 
f <- intersect(f1,f2)
length(f)
f_id <- getAnnot(f)
f_seq <- getSequence(f)
write.fasta(f_seq, f_id, "AtPPR_opt.fasta")

#####
AtPG <- read.fasta(file = "AtPGs.fasta")
nt_length <- getLength(AtPG)
hist(nt_length, xlab = "CDS length", ylab = "Frequency", breaks = 50, main = "AtPG gene family- continuously expanding")
#####

Gm_PGs <- read.fasta(file = "Gmax_PGs.fasta")
nt_length <- getLength(Gm_PGs)
abline(col = "red", v=1000)
f <- Gm_PGs[which(getLength(Gm_PGs) >= 1000)]
length(f)
f_id <- getAnnot(f)
f_seq <- getSequence(f)
write.fasta(f_seq, f_id, "Gm_PGs_opt.fasta")

#### pairwise alignment score for G. max polygalacturonase gene family with cds length ranging from 1200-2700 bp using Needleman-Wunsch algorithm with Bioconductor package Biostrings
library(Biostrings)
gm_nt <- readDNAStringSet("Gm_PGs_opt.fasta")
id <- names(gm_nt)
id_vec <- paste(gm_nt)
DNA_df <- data.frame(id, id_vec)
length(DNA_df)
combinations <- combn(100,2)
sigma <- nucleotideSubstitutionMatrix(match = 0, mismatch = -1, baseOnly = T)
choose(n=100,k=2) #4950

# looping for 4950 possible combinations
scores <- c()
gm_seq1 <- c()
gm_seq2 <- c()
for (i in 1:4950) {
  pwa <- pairwiseAlignment(toString(DNA_df[combinations[1,i],2]), toString(DNA_df[combinations[2,i],2]), substitutionMatrix = sigma, gapOpening = 0, gapExtension = -1, scoreOnly = F, type = "global")
  scores[i] <- c(pwa@score)
  gm_seq1[i] <- c(as.character(DNA_df[combinations[1,i],1]))
  gm_seq2[i] <- c(as.character(DNA_df[combinations[2,i],1]))
}
scoreDF <- data.frame(gm_seq1, gm_seq2, scores)
colnames(scoreDF) <- c("gm_seq1", "gm_seq2", "scores")
write.csv(scoreDF, file = "gm_PG_pwa_score.csv")

hist(scores, xlab = "Similarity", ylab = "Frequency", breaks = 200, probability = T, main = "NCR gene family PWA score distribution")
lines(density(scores), col = "Blue", lwd = 2)

#########
#### pairwise alignment score for At PPR gene family with cds length ranging from > 1000 bp using Needleman-Wunsch algorithm with Bioconductor package Biostrings
library(Biostrings)
AtPPR_nt <- readDNAStringSet("AtPPR_opt2.fasta")
id <- names(AtPPR_nt)
id_vec <- paste(AtPPR_nt)
DNA_df <- data.frame(id, id_vec)
length(DNA_df)
combinations <- combn(265,2)
sigma <- nucleotideSubstitutionMatrix(match = 0, mismatch = -1, baseOnly = T)
choose(n=265,k=2)

# looping for 34980 possible combinations
scores <- c()
AtPPR_seq1 <- c()
AtPPR_seq2 <- c()
for (i in 1:34980) {
  pwa <- pairwiseAlignment(toString(DNA_df[combinations[1,i],2]), toString(DNA_df[combinations[2,i],2]), substitutionMatrix = sigma, gapOpening = 0, gapExtension = -1, scoreOnly = F, type = "global")
  scores[i] <- c(pwa@score)
  AtPPR_seq1[i] <- c(as.character(DNA_df[combinations[1,i],1]))
  AtPPR_seq2[i] <- c(as.character(DNA_df[combinations[2,i],1]))
}
scoreDF <- data.frame(AtPPR_seq1, AtPPR_seq2, scores)
colnames(scoreDF) <- c("AtPPR_seq1", "AtPPR_seq2", "scores")
write.csv(scoreDF, file = "AtPPR_pwa_score.csv")

### for NCRs
f1<- ncrs[which(getLength(ncrs) >= 170)] 
f2<- ncrs[which(getLength(ncrs) <= 190)]
length(intersect(f1,f2)) # 240 NCRs
f <- (intersect(f1,f2))
length(f)
geneID <- getAnnot(f) # store geneIDs
f_seq <- getSequence(f)
write.fasta(f_seq, geneID, "ncrs_opt.fasta")

ncrs_nt <- readDNAStringSet("ncrs_opt.fasta")
id <- names(ncrs_nt) # holding all names as vector
id_vec <- paste(ncrs_nt)
DNA_df <- data.frame(id, id_vec) # 240 NCRs
combinations <- combn(240,2)
sigma <- nucleotideSubstitutionMatrix(match = 0, mismatch = -1, baseOnly = T)
choose(n=240,k=2) # 28680 combinations possible for 41 NCRs

# looping for 28680 possible combinations
scores <- c()
ncr_seq1 <- c()
ncr_seq2 <- c()
for (i in 1:28680) {
  pwa <- pairwiseAlignment(toString(DNA_df[combinations[1,i],2]), toString(DNA_df[combinations[2,i],2]), substitutionMatrix = sigma, gapOpening = 0, gapExtension = -1, scoreOnly = F, type = "global")
  scores[i] <- c(pwa@score)
  ncr_seq1[i] <- c(as.character(DNA_df[combinations[1,i],1]))
  ncr_seq2[i] <- c(as.character(DNA_df[combinations[2,i],1]))
}
scoreDF <- data.frame(ncr_seq1, ncr_seq2, scores)
colnames(scoreDF) <- c("ncr_seq1", "ncr_seq2", "scores")
write.csv(scoreDF, file = "ncr_pwa2_score.csv")
hist(scores, xlab = "Similarity", ylab = "Frequency", breaks = 100, probability = T, main = "NCR gene family PWA score distribution")
lines(density(scores), col = "Blue", lwd = 2)

### AtRLK

AtRLK_nt <- readDNAStringSet("At_rlk_opt2.fasta")
id <- names(AtRLK_nt) # holding all names as vector
id_vec <- paste(AtRLK_nt)
DNA_df <- data.frame(id, id_vec)
combinations <- combn(115,2)
sigma <- nucleotideSubstitutionMatrix(match = 0, mismatch = -1, baseOnly = T)
choose(n=115,k=2) ## 6555 possibilities

## looping
scores <- c()
rlk_seq1 <- c()
rlk_seq2 <- c()
for (i in 1:6555) {
  pwa <- pairwiseAlignment(toString(DNA_df[combinations[1,i],2]), toString(DNA_df[combinations[2,i],2]), substitutionMatrix = sigma, gapOpening = 0, gapExtension = -1, scoreOnly = F, type = "global")
  scores[i] <- c(pwa@score)
  rlk_seq1[i] <- c(as.character(DNA_df[combinations[1,i],1]))
  rlk_seq2[i] <- c(as.character(DNA_df[combinations[2,i],1]))
}
scoreDF <- data.frame(rlk_seq1, rlk_seq2, scores)
colnames(scoreDF) <- c("rlk_seq1", "rlk_seq2", "scores")
write.csv(scoreDF, file = "AtRLK_pwa2_score.csv")

hist(scores, xlab = "Similarity", ylab = "Frequency", breaks = 100, probability = T, main = "AtRLK gene family PWA score distribution")
lines(density(scores), col = "Blue", lwd = 2)

#### overlay histogram of pwa scores-
ncr_scores <- read.table("ncr_pwa2_scores.txt", sep="\t", header=T)
AtPPR_scores <- read.table("AtPPR_scores.txt", sep="\t", header=T)
Gm_scores <- read.table("Gm_pg_score.txt", sep="\t", header=T)
ncr_scores$pwa_score <- 'NCRs' 
AtPPR_scores$pwa_score <- 'AtPPR' 
Gm_scores$pwa_score <- 'Gm_PGs' 
pwa_scores <- rbind(ncr_scores, AtPPR_scores, Gm_scores)
library(ggplot2)
ggplot(pwa_scores, aes(scores, fill = pwa_score)) + geom_density(alpha = 0.2)

###### Replacing Gm with AtRLKs
AtRLK_scores <- read.table("AtRLK_scores.txt", sep="\t", header=T)
AtRLK_scores$pwa_score <- 'AtRLKs'
pwa_scores2 <- rbind(ncr_scores, AtPPR_scores, AtRLK_scores)
ggplot(pwa_scores2, aes(scores, fill = pwa_score)) + geom_density(alpha = 0.2)

######
ggplot(pwa_scores, aes(pwa_score, fill = pwa_score)) + geom_density(alpha = 0.2)
```
