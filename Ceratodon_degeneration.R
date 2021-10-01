#### R script for the manuscript by S.F.McDaniel and S.B.Carey "On the degeneration of UV sex chromosomes"
### Script by Sarah B. Carey
## Analyses below run using R v4.1.1

######### Install packages #########
## install all packages used in this script (versions used are next to package name)

# BiocManager v1.30.16
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#karyoploteR v1.18.0
#BiocManager::install("karyoploteR")

# ggplot2 v3.3.5
#install.packages("ggplot2")

# ggtree v3.0.4
#BiocManager::install("ggtree")


######### Figures #########
##### Fig. 1A, Gene vs TE density ######

library(karyoploteR)

### read in density data
gene.gtf <- read.csv("plotData/gene_density_filtered.csv",header=T)
gtf.GR <- toGRanges(data.frame(chr=gene.gtf$Chromosome, start=gene.gtf$Start, end=gene.gtf$Stop))

TE.gtf <- read.csv("plotData/TE_all_density_filtered.csv",header=T)
TE.gtf.GR <- toGRanges(data.frame(chr=TE.gtf$Chromosome, start=TE.gtf$Start, end=TE.gtf$Stop))

### read in genome lengths
r40 <- read.table("plotData/R40_genome_lengths_sex.txt", header=T)
custom.genome <- toGRanges(data.frame(chr=r40$Chromosome, start=r40$ChromStart, end=r40$ChromEnd))


png("Figures/Gene_TE_densities_ceratodon.png", width = 10, height = 5, units = 'in', res = 300)

pp <- getDefaultPlotParams(plot.type = 4)
pp$ideogramlateralmargin <- 0.015
pp$leftmargin <- 0.1

kp <- plotKaryotype(genome=custom.genome,pin=8, plot.type = 4,labels.plotter = NULL, plot.params=pp,
                    chromosomes=c("Chr01","U","V"))
kp <- kpDataBackground(kp, data.panel = 1, color="gray95")
kpAddChromosomeNames(kp, chr.names=c("Chr01","U","V"), yoffset=-5, cex=1.25, font=2)
kpAddLabels(kp, "Genes", r0=0.65, srt=90, label.margin=0.07,font=2, cex=1.5)
kpAddLabels(kp, "TEs", r0=-0.5, srt=90, label.margin=0.07,font=2, cex=1.5)

kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 5, tick.col="black", cex=0.75,
                 minor.tick.dist = 5000000, minor.tick.len = 5, minor.tick.col = "black", add.units = F)

# gene density
kp <- kpLines(kp, data=gtf.GR, data.panel = 1, col=("gray50"),y=gene.gtf$Density,r0=0.6,r1=1,lwd=2, ymax=0.75)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.6, r1=1, cex=1, numticks=5, labels=c("0","","","","0.75"))
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.6, r1=1, cex=1, numticks=5, labels=c("0","","","","0.75"), side=3)
kpAbline(kp, h=0.5, col="black", r0=0.6, r1=1, lwd=3)


# all TE density
kp <- kpLines(kp, data=TE.gtf.GR, data.panel = 1, col=("gray50"),y=TE.gtf$Density,r0=0,r1=0.4,lwd=2)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.4, cex=1, numticks=5, labels=c("0","","","","1"))
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.4, cex=1, numticks=5,  labels=c("0","","","","1"), side=3)
kpAbline(kp, h=0.5, col="black", r0=0, r1=0.4, lwd=3)


kpAbline(kp, h=0.5, col="white",lwd=5)

dev.off()



####### Fig. 1B, Scatterplot of frequency of optimal codon vs gene expression ########

library(ggplot2)

codonbias_cerat <- read.csv("plotData/codon_analyses_nonZero.csv",header=TRUE)

auto <- subset(codonbias_cerat, chromosome=="Autosomal")
chrU <- subset(codonbias_cerat, chromosome=="U-linked")
chrV <- subset(codonbias_cerat, chromosome=="V-linked")


png("Figures/fopVsExpr.png", width = 6, height = 6, units = 'in', res = 300)

p <- ggplot(auto, aes(x=log2(auto$mean_count), y=auto$fop)) + 
  geom_point(col="gray50", alpha=0.1, size=2) +
  geom_point(aes(x =log2(chrU$mean_count), y= chrU$fop), col="orange", data=chrU, alpha=0.1, pch=17, size=2) +
  geom_point(aes(x =log2(chrV$mean_count), y =chrV$fop), col="blue", data=chrV,alpha=0.1, pch=15, size=2) +
  geom_line(method = "lm", size = 4, col="midnightblue", aes(x =log2(chrV$mean_count), y =chrV$fop), data=chrV, stat="smooth", alpha=0.75) + 
  geom_line(method = "lm", size = 4, col="peru", aes(x =log2(chrU$mean_count), y= chrU$fop), data=chrU, stat="smooth", alpha=0.75) +
  geom_line(method = "lm", size = 4, col="gray10", stat="smooth", alpha=0.75) +
  xlab("Gene expression") +
  ylab("fop") +
  theme(axis.title.x = element_text(size=25)) +
  theme(axis.title.y = element_text(size=25))+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20)) 

p

dev.off()

#### get slope and intercept of lm 

Auto.lm <- lm(formula=fop~log2(mean_count),data=auto)
summary(Auto.lm)

U.lm <- lm(formula=fop~log2(mean_count),data=chrU)
summary(U.lm)

V.lm <- lm(formula=fop~log2(mean_count),data=chrV)
summary(V.lm)



############ Fig. 2A, Moss sex chromosome evolution tree ############

library(ggtree)

tree_file <- read.tree("plotData/moss_sex_chrom_tree.tre")

genus <- c("Brachythecium", "Hylocomium", "Aulacomnium", "Syntrichia","Ceratodon", "Scouleria", 
           "Physcomitrium", "Buxbaumia","Sphagnum", "Liverworts")
species <- c("rivulare", "splendens", "palustre", "princeps","purpureus", "aquatica", 
             "patens","aphylla", "palustre","")
d <- data.frame(label = tree_file$tip.label, genus = genus,
                species = species)


tree_plot <- ggtree(tree_file, branch.length = "none", color="black", size=1.25) %<+% d +
  geom_tiplab(size=6, color="black", aes(label=paste0('bolditalic(', genus, ')~bolditalic(', species, ')')), parse=T) +
  xlim(NA,15) +
  theme_tree("white") 
tree_plot

ggsave("Figures/moss_sex_chrom_tree.png", tree_plot, units="in", width=8, height=4, dpi=300,
       device="png")




############ Fig. 2B, Ks of Ancestral Element D genes ############

library(ggplot2)

ks_data <- read.csv("plotData/ceratodon_ks.csv", header=TRUE)
ks_D <- subset(ks_data, ks_data$Ancestral_element =="D")


png("Figures/Ks_AncEleD.png", width = 6, height = 6, units = 'in', res = 300)

p <- ggplot(ks_D, aes(reorder(Gene.tree.clusterID, Ks), Ks)) + 
  geom_point(col="chocolate4", alpha=0.75, size=3) +
  xlab("Gene") +
  ylab("Ks") +
  theme(axis.title.x = element_text(size=25)) +
  theme(axis.title.y = element_text(size=25))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=20)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) + 
  theme(panel.background = element_rect(color="black")) + 
  scale_y_continuous(minor_breaks = seq(1, 10, 1)) +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x = element_blank())
p

dev.off()




####### Fig. 2C, Example tree showing coalescence in pleurocarps #######

library(ggtree)
library(treeio)

tree_file <- read.tree("plotData/sex_chrom_tree_exampleConvergence.tre")

genus <- c("Brachythecium", "Hylocomium", "Brachythecium", "Hylocomium","Aulacomnium", "Syntrichia","Ceratodon", "Scouleria", "Buxbaumia", "Aulacomnium", "Syntrichia","Ceratodon","Scouleria", "Buxbaumia","Sphagnum")
species <- c("rivulare", "splendens","rivulare", "splendens", "palustre", "princeps","purpureus", "aquatica", "aphylla", "palustre", "princeps","purpureus", "aquatica", "aphylla", "palustre")
sex <- c("male", "male", "female", "female", "male", "male","male", "male", "male", "female", "female","female", "female", "female"," ")

d <- data.frame(label = tree_file$tip.label, genus = genus,
                species = species, sex = sex)


tree_plot <- ggtree(tree_file, branch.length = "none", color="black", size=1.25) %<+% d +
  theme_tree("white") +
  geom_hilight(node=26, fill="#D55E00", alpha=0.75, extend=22)+
  geom_hilight(node=29, fill="#D55E00", alpha=0.25, extend=22) +  
  geom_hilight(node=18, fill="#0072B2", alpha=0.25, extend=22) +  
  geom_hilight(node=24, fill="#0072B2", alpha=0.75, extend=22)+
  geom_tiplab(size=6, color="black", aes(label=paste0('bolditalic(', genus, ')~bolditalic(', species, ')~bold(', sex, ')')), parse=T) +
  xlim(NA,40)
tree_plot

ggsave("Figures/sex_chrom_tree_exampleConvergence.png", tree_plot, units="in", width=8, height=4, dpi=300,
       device="png")




####### Fig. 2D, Gene tree showing coalescence in pleurocarps #######

library(ggtree)
library(treeio)


tree_file <- read.raxml("plotData/RAxML_bipartitionsBranchLabels.cluster350.tree")

##get node numbers
tree_plot <- ggtree(tree_file, color="black", size=1) + 
  geom_tiplab(size=2, color="black") +
  geom_nodelab(aes(label=node),size=2, color="red")+
  theme_tree("white") 
tree_plot


png("Figures/ancient_convergence.png", width = 6, height = 10, units = 'in', res = 300)

tree_plot <- ggtree(tree_file, color="black", size=1) 

tree_plot <- collapse(tree_plot,node=232) + 
  geom_point2(aes(subset=(node==232)), shape=23, size=5, fill='black')

tree_plot <- collapse(tree_plot,node=207) + 
  geom_point2(aes(subset=(node==207)), shape=23, size=5, fill='black')

tree_plot <- collapse(tree_plot,node=405) + 
  geom_point2(aes(subset=(node==405)), shape=23, size=5, fill='black')


tree_plot  +  
  geom_hilight(node=279, fill="#D55E00", alpha=0.25, extend=1.2) +
  geom_hilight(node=380, fill="#D55E00", alpha=0.75, extend=1.4) +
  geom_hilight(node=336, fill="#0072B2", alpha=0.25, extend=1.2) +
  geom_hilight(node=358, fill="#0072B2", alpha=0.75, extend=1.3) +
  geom_treescale(color="gray50", fontsize=5,x=0, y=60) + 
  geom_tiplab(size=2, color="black", hjust=-0.05) +
  geom_nodelab(aes(label=bootstrap),size=2, color="red") +
  theme_tree("white") +
  xlim(0,3)

dev.off()

