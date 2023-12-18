library(tidyverse)
library(cowplot)
library(colorspace)


okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999")

adjust_data<-function(df, cutoff){
  return(df %>% mutate(reads_per_kb_of_genome=ifelse(pct_covered>=cutoff, reads_per_kb_of_genome, NA), cutoff=cutoff))
}



cube_root<-function(x){
  x^(1/3)
}

cubed<-function(x){
  x^3
}


make_plot<-function(df){
  p = ggplot(df, aes(x=sample_label, y=reads_per_kb_of_genome)) +
    geom_point(aes(color=ecological_zone, fill=ecological_zone), size=4, alpha=0.7, shape=21) +
    scale_y_continuous("Reads per kB of genome", trans=scales::trans_new('cube root', cube_root, cubed)) +
    scale_x_discrete("Sample") +
    scale_color_manual(values = darken(okabe_ito, 0.3)) +
    scale_fill_manual(values = okabe_ito) +
    facet_wrap(~contig, ncol=1) +
    theme_bw(16) +
    theme(axis.text.x = element_text(size=6, angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
  return(p)
  
}


make_viral_plot<-function(df, cutoff=40){
  sub_df = adjust_data(df, cutoff)
  p = ggplot(sub_df, aes(x=contig, y=reads_per_kb_of_genome)) +
    geom_point(aes(color=phage_type, fill=phage_type), size=3, alpha=0.7, shape=21) +
    scale_y_continuous("Reads per kB of genome",
                       trans=scales::trans_new('cube root', cube_root, cubed),
                       breaks=c(100, 500, 2000, 4000, 8000, 12000)) +
    scale_x_discrete("Phage Genome") +
    scale_color_manual('Phage type', values = darken(okabe_ito, 0.3)) +
    scale_fill_manual('Phage type', values = okabe_ito) +
    facet_grid(phage_type~depth, scales='free_y') +
    theme_bw(16) +
    theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size=11),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none") +
    coord_flip()
  return(p)
}

## we want to remove any genomes that have 0 coverage at 40% in ALL samples
df = read_tsv('data/biller_merged_dataframe.tsv.gz')
keep_these = df %>% mutate(reads_per_kb_of_genome=ifelse(pct_covered>=40, reads_per_kb_of_genome, 0),
                             read_count=ifelse(pct_covered>=40, read_count, 0)) %>% 
  group_by(contig) %>%
  summarise(total_reads=sum(read_count), count=sum(read_count>0)) %>%
  filter(total_reads>0) %>%
  pull(contig)

df = df %>% filter(contig %in% keep_these)


#now figure out what the most abundantly recruited samples are
contig_order = df %>% group_by(contig) %>%
  summarise(total_abundance=sum(reads_per_kb_of_genome)) %>%
  arrange(total_abundance) %>%
  pull(contig)

#now figure out what the most abundantly recruited samples are



df = df %>% mutate(contig=factor(contig, levels=contig_order))
df = df %>% mutate(depth=factor(depth,levels=c(1,10, 40, 60, 80, 100, 120, 160))) %>% filter(depth!='NA')



p = make_viral_plot(df, cutoff=40)
ggsave('plots/Biller-combined-40pct.pdf', plot = p, device=pdf, width=210*2, height=297, units='mm', limitsize = FALSE)
p2 = make_viral_plot(df, cutoff=70)
ggsave('plots/Biller-combined-70pct.pdf', plot = p2, device=pdf, width=210*2, height=297, units='mm', limitsize = FALSE)


