
# Script to visualize aspects of data in VCF files.
# Dependencies: 'variantannotation' package in 'bioconductor'
# Written by Jamie Walters (jrw47 [at] stanford __dot__ edu)

# Run as:
#% Rscript analyze_vcf.R  myvariants.vcf output_file

#Version 2 is amended to deal with missing data which crops up when running on large data files.



args <- commandArgs(TRUE)
library(VariantAnnotation)

#setwd("~/Documents/PostDoc_Research/MK-tests_Hmelpomene/vcf_R_scripts")
#file = "~/Documents/PostDoc_Research/MK-tests_Hmelpomene/vcf_R_scripts/test_amaagl.vcf"
#basename = "out"
file <- args[1]
basename <- args[2]


vcf <- readVcf(file, "Hmel")
hdr <- exptData(vcf)[["header"]]



####################################
# Plot distribution of QUAL score  #
####################################

quals <- vcf@fixed$QUAL

mean.qual   <- mean(quals)
median.qual <- median(quals)
sd.qual     <- sd(quals)

bin <-100
png(paste(basename,"_QualHist.png", sep=""), height=5,width=5,units="in", res=200)
hist(quals, breaks=seq(0,max(quals)+bin,bin), col="dark green", border="dark green", xlim=c(0, mean.qual+5*sd.qual),
	main="Distribution of SNP Quality Scores", xlab="SNP Quality Score", cex.axis=0.8)
dev.off()

###############################
# Plot distribution of depth  #
###############################

depth <- vcf@info$DP

depth.mean  <- mean(depth)
depth.median <- median(depth)
depth.sd     <- sd(depth)


bin <-5
png(paste(basename,"_DepthHist.png", sep=""), height=5,width=5,units="in", res=200)
hist(depth, breaks=seq(0,max(depth)+bin,bin), col="dark green", border="dark green", xlim=c(0, depth.mean+5*depth.sd),
	main="Distribution of Total Depth per SNP", xlab="Summed Coverage per SNP", cex.axis=0.8)
dev.off()

#col2rgb("black")
png(paste(basename,"_QualByDepth.png", sep=""), height=5,width=5,units="in", res=200)
plot (quals ~ depth, xlim=c(0, depth.mean+5*depth.sd), ylim=c(0, mean.qual+5*sd.qual), pch=19, cex=.3, col=rgb(0,0,0,.05),
xlab="Summed Coverage per SNP", ylab="SNP Quality Score", main="Coverage Depth ~ SNP Quality", cex.axis=0.8 ) 
dev.off()


##################################
# Count & Plot numbers of empty  #
# genotypes per sample and site  #
##################################
count.empty.genos <- function (x) {
	length(grep('./.',x, fixed=T))
}
empty.genos.sample <- apply(geno(vcf)$GT, 2, count.empty.genos )  # apply to rows
empty.genos.sites  <- apply(geno(vcf)$GT, 1, count.empty.genos )  # apply to cols

#plot frequency of sites missing N samples
png(paste(basename,"_EmptyGeno-Sites.png", sep=""), height=5,width=5,units="in", res=200)
barplot(table(empty.genos.sites[empty.genos.sites > 0]), main = "Missing genotypes per site", ylab="Number of sites", xlab="Count of missing genotypes", sub=paste(dim(vcf)[1], "sites counted"), cex.lab=1.2, col="light blue" )
dev.off()

#Plot counts sites missing data per sample
png(paste(basename,"_EmptyGeno-Samples.png", sep=""), height=5,width=5,units="in", res=200)
par(oma=c(4,1,0.5,1))
barplot(empty.genos.sample, las=2, main = "Missing genotypes per sample", ylab="Number of sites", cex.lab=1.2, col="light blue", cex.names=.7, )
mtext( paste(dim(vcf)[1], "sites counted"), outer=T, padj=4, 1)
dev.off()

write.table(table(empty.genos.sites),file=paste(basename,"_EmptyGeno-Sites.txt", sep=""), quote=F, sep="\t", row.names=F ) 
write.table(empty.genos.sample, file=paste(basename,"_EmptyGeno-Samples.txt", sep=""), quote=F, sep="\t", row.names=T, col.names=F) 

##################################
#Plot allele frequencies & counts#
##################################
af <- unlist(vcf@info$AF)  # get allele frequncies
ac <- unlist(vcf@info$AC)  # get allele counts


#Allele freqs
af.table.all <- table(af)
af.table.complete <- table(af[empty.genos.sites == 0])

png(paste(basename,"_SFS_completeGenos.png", sep=""), height=5,width=5,units="in", res=200)
barplot(af.table.complete/sum(af.table.complete), main="Allele Frequency Spectrum", ylab= "Proportion of SNPs", xlab="Allele Frequency", col="light blue", las=2, cex.lab=1.2, cex.names=.5, sub=paste(sum(af.table.complete), "sites with complete genotypes counted") )
dev.off()

out.af.complete <- cbind(names(af.table.complete),af.table.complete, af.table.complete/sum(af.table.complete))  # create object for writing to file
out.af.all <- cbind(names(af.table.all),af.table.all, af.table.all/sum(af.table.all))  # create object for writing to file

write.table(out.af.complete, file=paste(basename,"_SFS_completeGenos.txt", sep=""), quote=F, sep="\t", 
			row.names=F, col.names=c("Allele_Frequency","Count","Proportion"))
write.table(out.af.all, file=paste(basename,"_SFS_allsites.txt", sep=""), quote=F, sep="\t", 
			row.names=F, col.names=c("Allele_Frequency","Count","Proportion"))


#Allele counts
ac.table.all <- table(ac)
ac.table.complete <- table(ac[empty.genos.sites == 0])

png(paste(basename,"_AlleleCountSpectrum.png", sep=""), height=5,width=5,units="in", res=200)
barplot(ac.table.complete/sum(ac.table.complete), main="Allele Frequency Spectrum", ylab= "Proportion of SNPs", xlab="Allele Count", col="light blue", las=1, cex.lab=1.2, cex.names=.9, sub=paste(sum(ac.table.complete), "sites with complete genotypes counted") )
dev.off()

out.ac.complete <- cbind(names(ac.table.complete),ac.table.complete, ac.table.complete/sum(ac.table.complete))  # create object for writing to file
out.ac.all <- cbind(names(ac.table.all),ac.table.all, ac.table.all/sum(ac.table.all))  # create object for writing to file
write.table(out.ac.complete, file=paste(basename,"_AlleleCounts_completeGenos.txt", sep=""), quote=F, sep="\t", 
			row.names=F, col.names=c("Allele_Count","Site_Count","Proportion"))
write.table(out.af.all, file=paste(basename,"_AlleleCounts_allsites.txt", sep=""), quote=F, sep="\t", 
			row.names=F, col.names=c("Allele_Count","Site_Count","Proportion"))



