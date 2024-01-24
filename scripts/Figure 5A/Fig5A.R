# loading the Vegan library
library(vegan)
library(pracma)

BATS <- read.table("BATS_short_reads_2031_coverage_rm2v.csv",check.names=FALSE, header=T,sep=",",row.names=1)
BATS <- as.data.frame(t(BATS))

VPs <- data.matrix(BATS)
d.clr = nthroot(VPs, 3)
Bray=vegdist(d.clr,"bray")

VPs_bray=capscale(Bray~-1)
VPs_bray$CA$eig
VPs_bray_eig = eigenvals(VPs_bray)
percentage_variance_explained <- VPs_bray_eig / sum(VPs_bray_eig)
sum_percentage_variance_explained <- cumsum(VPs_bray_eig / sum(VPs_bray_eig))

xlabel= as.numeric(format(round((percentage_variance_explained[1]*100), 2), nsmall = 2))
xlabel= sprintf("%.2f %%", xlabel)
xlabel= paste ("PCo1 (", xlabel, ")")
ylabel= as.numeric(format(round((percentage_variance_explained[2]*100), 2), nsmall = 2))
ylabel= sprintf("%.2f %%", ylabel)
ylabel= paste ("PCo2 (", ylabel, ")")

env=read.table("BATS_env.csv", sep=",",check.names=FALSE, header=T,row.names=1)

# unique(env$depth)
depth_shape=rep(17,nrow(env))
depth_shape[env$depth=="200m"]=16

# unique(env$type)
type_col=rep('#f7a520',nrow(env))
type_col[env$type=="c"]='#7568ad'

env$depth_type <- paste(env$depth, '_', env$type)

plot(VPs_bray, main= "PCoA (Bray-Curtis)",type="n", xlab=xlabel, ylab=ylabel,xlim=c(-2, 2),cex=0)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
points(VPs_bray,cex= 1.0, pch=depth_shape, col=type_col)
# text(VPs_bray, cex = 0.2)

ordiellipse(VPs_bray, env$depth_type, conf=0.95, lty=5,col="#d1cec9")
ordiellipse(VPs_bray, env$depth_type, conf=0.98, lty=5,col="#ebe8e4")

legend('topleft',c("95%","98%"),pch=3,col=c('#d1cec9', '#ebe8e4'),box.lty=0)

# https://chrischizinski.github.io/rstats/adonis/ 
adonis2(d.clr ~ depth_type, data = env, permutations = 999,method = "bray") # *** 0.001


