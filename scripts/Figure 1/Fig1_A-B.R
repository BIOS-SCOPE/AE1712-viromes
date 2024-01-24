#### Fig 1A ####
library(ggplot2)

rank_abundance<-read.table("rank_abundance_2301.csv",check.names=FALSE, header=T,sep=",")


ggplot(rank_abundance,aes(x=rank,y=log10(total_rm2v),
                          fill=type_use),show.legend = TRUE)+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#83c77c","#ae8dc1"),
                    labels = c("VirION-read Assembled Viral Populations", 
                               "Short-read Assembled Viral Populations"))+
  theme_classic()

#### Fig 1B ####
library(dplyr)

virus_length=read.table("length_2301.txt", sep="\t", 
                        row.names =1, check.names=FALSE, header=T)

ggplot(virus_length,aes(x=type))+
  geom_boxplot(aes(y=log(length),colour = type))+
  scale_colour_manual(values=c("#83c77c","#ae8dc1"),
                      labels = c("VirION-read Assembled Viral Populations", 
                                 "Short-read Assembled Viral Populations"))+
  theme_classic()
wilcox.test(length~type, data=virus_length)$p.value

group_by(virus_length, type) %>%
  summarise(
    count = n(),
    median = median(length, na.rm = TRUE),
    IQR = IQR(length, na.rm = TRUE),
    mean = mean(length, na.rm = TRUE)
  )

