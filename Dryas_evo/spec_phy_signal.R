#Determine phylogenetic signal between species and populations
################################################################################
#Set up
################################################################################

library(spectrolab)
library(ape)
library(phytools)
library(geiger)

#load spectra, normalize, and take average spectra per individual
spec = readRDS("Data/clean_all.rds")
spec = spec[!meta(spec)$Species_ID == "DX",]
spec = normalize(spec)
spec = aggregate(spec, by = meta(spec)$Name, mean, try_keep_txt(mean))
spec = resample(spec, seq(400, 2400, by = 10))
spec.m = as.matrix(spec)

meta_new = read.csv("Data/mean_meta_new.csv", stringsAsFactors = F)
rownames(spec.m) <- meta_new[, "genetic_name"]

#load tree
tree = read.tree(file = "Data/dryas_phylogeny.treefile")

################################################################################
#Match tree to data
################################################################################

matched_treedata = geiger::treedata(phy = tree, data = spec.m)

spec_data = matched_treedata$data
phylo = matched_treedata$phy

################################################################################
#Calculate phylogenetic signal of spectra 
################################################################################

k = c()
p = c()

for (i in 1:201){
sig = phylosig(phylo, spec_data[,i], method = 'K', test = T, nsim = 999)

sig.k = assign(paste0("sig.k",i), sig[[1]])
k <- append(k, get('sig.k'))

sig.p = assign(paste0('sig.p',i), sig[[2]])
p <- append(p, get('sig.p'))
}

################################################################################
#Plot phylogenetic signal
################################################################################
k1 = matrix(nrow = 201)
rownames(k1) <- seq(400,2400, by = 10)
k1 = cbind(k1, k)
k.m = k1[,-1]
k2 = as.matrix(k.m)
k2$num = seq(1,201, 1)
plot(k2)


plot(k, xaxt = "n", xlab = 'Wavelength (nm)')
axis(1, at = c(0,50,100, 150,200) )















