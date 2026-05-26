#This script performs PCAIR on the European subset (can be applied to any population) of the MESA WGS data after QC filtering. It uses the KING matrix to account for relatedness in the PCA calculation. The resulting eigenvectors and eigenvalues are saved for downstream use in adjusting for population structure in association analyses.

library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(dplyr)
library(tibble)

#Convert PLINK files to GDS
snpgdsBED2GDS(bed.fn = "~/path/to/europop_filtered.bed", 
              bim.fn = "~/path/to/europop_filtered.bim", 
              fam.fn = "~/path/to/europop_filtered.fam", 
              out.gdsfn = "~/path/to/genotype.gds")

#Close all previous GDS files
showfile.gds(closeall=TRUE)

#Open GDS file
gdsfile <- "~/path/to/genotype.gds"
gds <- snpgdsOpen(gdsfile)

#Calculate the kinship coefficients 
king <- snpgdsIBDKING(gds)

kingMat <- king$kinship 
colnames(kingMat)<-king$sample.id
row.names(kingMat)<-king$sample.id

#Save king RDS file for future use (if needed)
saveRDS(king, file = "~/path/to/King_matrix.RDS")

#Close all previous GDS files
showfile.gds(closeall=TRUE)

#Create genotype object needed for pcair function
geno <- GdsGenotypeReader(filename = "~/path/to/genotype.gds") 
genoData <- GenotypeData(geno)

#Run PCAIR 
mypcair <- pcair(genoData, kinobj = kingMat, divobj = kingMat)

# plot top 2 PCs
plot(mypcair)
# plot PCs 2 and 3
plot(mypcair, vx = 2, vy = 3)
# plot PCs 3 and 4
plot(mypcair, vx = 3, vy = 4)

#Convert to PCA vec and val files
eigenvec<-mypcair$vectors %>% as.data.frame() %>% rownames_to_column(var="sample_id")
str(eigenvec)
val<-mypcair$values %>% as.data.frame()
fwrite(eigenvec,"~/path/to/PCAIR.eigenvec",col.names = T,row.names = F,sep='\t')
fwrite(val,"~/path/to/PCAIR.eigenval")