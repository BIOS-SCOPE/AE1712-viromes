#### Figure 3 ####
library(vegan)
library(pracma)
library(tibble)
library(plyr)

BATS_GOV2 <- read.table("GOV2.0_BATS_coverage.csv",
                        check.names=FALSE, header=T,sep=",",row.names=1)
BATS_GOV2 <- as.data.frame(t(BATS_GOV2))

VPs <- data.matrix(BATS_GOV2)
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

env=read.table("BATS_GOV2.0_env.csv",
               check.names=FALSE, header=T,sep=",")

BATS_GOV2 <- tibble::rownames_to_column(BATS_GOV2, "Site")
env=join(BATS_GOV2["Site"],env, by='Site')

col=rep('#000000',nrow(env))
col[env$Zone=="TT-EPI"]="#f7e092"
col[env$Zone=="TT-MES"]="#f98ce4"
col[env$Zone=="ANT"]="#b3b3b3"
col[env$Zone=="ARC"]="#b4d7f0"
col[env$Zone=="v"]="#f7a520"
col[env$Zone=="c"]="#7568ad"

shape=rep(15,nrow(env))
shape[env$Shape=="200m"]=16
shape[env$Shape=="80m"]=17

type=col
type[env$Zone=="v"]="#000000"
type[env$Zone=="c"]="#000000"

plot(VPs_bray, main= "PCoA (Bray-Curtis)",type="n", xlab=xlabel, ylab=ylabel,xlim=c(-2, 2),cex=0)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
points(VPs_bray,cex= 1.0, pch=shape, col=col)

ordiellipse(VPs_bray, type, conf=0.95, lty=5,col="#d1cec9")
ordiellipse(VPs_bray, type, conf=0.98, lty=5,col="#ebe8e4")

legend('topleft',c("95%","98%"),pch=3,col=c('#d1cec9', '#ebe8e4'),box.lty=0)

adonis2(d.clr ~ type, permutations = 999,method = "bray") # *** 0.001

