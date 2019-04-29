###
###
###

# set working directory to where file is saved
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))

# get working directory
getwd()

# load in libraries
library(dplyr)

# load files
dataset <- read.table(file.path(getwd(), "Datasets", "merged_data_final.txt"), header = T, sep = "\t")
inflGenes <- read.table("P:/FSE_MACSBIO/2018-Laurent-Winckers/Final_files/merged_infl_genes.txt", header = T, sep = "\t")
pwGenes <- read.table("P:/FSE_MACSBIO/2018-Laurent-Winckers/Final_files/nodes_final_entrezgene.txt", header = T, sep = "\t")
pwGenes <- pwGenes[c(-1:-10),]

dataset1 <- dataset[dataset$entrezgene %in% pwGenes$entrezgene,]
dataset2 <- dataset[dataset$entrezgene %in% inflGenes$entrezgene,]

# breast cancer
datBC <- dataset[,-5:-16]
datBC <- na.omit(datBC)

# lung cancer
datLC <- dataset[,c(-3,-4,-7:-16)]
datLC <- na.omit(datLC)

# metabolically unhealthy obese
datMUO <- dataset[,c(-3:-6,-9:-16)]
datMUO <- na.omit(datMUO)

# rheumatoid arthritis
datRA <- dataset[,c(-3:-8,-11:-16)]
datRA <- na.omit(datRA)

# DCM
datDCM <- dataset[,c(-3:-10,-13:-16)]
datDCM <- na.omit(datDCM)

# ICM
datICM <- dataset[,c(-3:-12,-15:-16)]
datICM <- na.omit(datICM)

# systemic lupus erythematosus
datSLE <- dataset[,c(-3:-14)]
datSLE <- na.omit(datSLE)

###
###
###

diseases <- c("BC","LC","MUO","RA","DCM","ICM","SLE")
frames <- list(datBC, datLC, datMUO, datRA, datDCM, datICM, datSLE)

y <- 1
repeat {

x <- 1
repeat {
assign(paste0("dat", diseases[x], "1"), subset(frames[[x]], frames[[x]][[paste0("logFC_", diseases[x])]] <= -0.26))
  x = x + 1
  if (x == length(diseases)+1){
    break}
}
frames1 <- list(datBC1, datLC1, datMUO1, datRA1, datDCM1, datICM1, datSLE1)

x <- 1
repeat {
assign(paste0("dat", diseases[x], "2"), subset(frames[[x]], frames[[x]][[paste0("logFC_", diseases[x])]] >= 0.26))
  x = x + 1
  if (x == length(diseases)+1){
    break}
}
frames2 <- list(datBC2, datLC2, datMUO2, datRA2, datDCM2, datICM2, datSLE2)

x <- 1
repeat {
assign(paste0("dat", diseases[x], "3"), rbind(frames1[[x]], frames2[[x]]))
  x = x + 1
  if (x == length(diseases)+1){
    break}
}
frames3 <- list(datBC3, datLC3, datMUO3, datRA3, datDCM3, datICM3, datSLE3)

x <- 1
repeat {
assign(paste0("dat", diseases[x], "4"), subset(frames3[[x]], frames3[[x]][[paste0("PValue_", diseases[x])]] <= 0.05)) # data frame with significant differentially expressed genes
  x = x + 1
  if (x == length(diseases)+1){
    break}
}
frames4 <- list(datBC4, datLC4, datMUO4, datRA4, datDCM4, datICM4, datSLE4)

x <- 1
repeat {
assign(paste0("n", diseases[x]), nrow(frames4[[x]]))
assign(paste0("dat", diseases[x], "5"), frames4[[x]][frames4[[x]][,"entrezgene"] %in% pwGenes$entrezgene,])
x = x + 1
if (x == length(diseases)+1){
  break}
}
frames5 <- list(datBC5, datLC5, datMUO5, datRA5, datDCM5, datICM5, datSLE5)
nframes <- c(nBC, nLC, nMUO, nRA, nDCM, nICM, nSLE)

x <- 1
repeat {
assign(paste0("dat", diseases[x], "5"), unique(frames5[[x]]))
assign(paste0("n", diseases[x], "pw"), nrow(frames5[[x]]))
x = x + 1
if (x == length(diseases)+1){
  break}
}
framesPW <- c(nBCpw, nLCpw, nMUOpw, nRApw, nDCMpw, nICMpw, nSLEpw)

x <- 1
repeat {
assign(paste0("dat", diseases[x], "6"), frames4[[x]][frames4[[x]][,"entrezgene"] %in% inflGenes$entrezgene,])
  x = x + 1
  if (x == length(diseases)+1){
    break}
}
frames6 <- list(datBC6, datLC6, datMUO6, datRA6, datDCM6, datICM6, datSLE6)

x <- 1
repeat {
assign(paste0("dat", diseases[x], "6"), unique(frames6[[x]]))
assign(paste0("n", diseases[x], "infl"), nrow(frames6[[x]]))
  x = x + 1
  if (x == length(diseases)+1){
    break}
}
framesInfl <- c(nBCinfl, nLCinfl, nMUOinfl, nRAinfl, nDCMinfl, nICMinfl, nSLEinfl)

x <- 1
repeat {
assign(paste0("perc", diseases[x]), (as.numeric(nframes[x]) / as.numeric(nrow(frames[[x]]) * 100)))
  x = x + 1
  if (x == length(diseases)+1){
    break}
}
percframes <- c(percBC, percLC, percMUO, percRA, percDCM, percICM, percSLE)

x <- 1
repeat {
assign(paste0("perc", diseases[x], "pw"), (as.numeric(framesPW[x]) / as.numeric(nrow(frames[[x]]) * 100)))
  x = x + 1
  if (x == length(diseases)+1){
    break}
}
percframesPW <- c(percBCpw, percLCpw, percMUOpw, percRApw, percDCMpw, percICMpw, percSLEpw)

x <- 1
repeat {
assign(paste0("perc", diseases[x], "infl"), (as.numeric(framesInfl[x]) / as.numeric(nrow(frames[[x]]) * 100)))
  x = x + 1
  if (x == length(diseases)+1){
    break}
}
percframesInfl <- c(percBCinfl, percLCinfl, percMUOinfl, percRAinfl, percDCMinfl, percICMinfl, percSLEinfl)

y = y + 1
if (y == length(diseases)+1){
  break}
}

