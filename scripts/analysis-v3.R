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



## we want to remove any genomes that have 0 coverage at 40% in ALL samples
df = read_tsv('data/merged_dataframe.tsv.gz') %>% filter(sample_type=='GOV2')
keep_these = df %>% mutate(reads_per_kb_of_genome=ifelse(pct_covered>=40, reads_per_kb_of_genome, 0),
                             read_count=ifelse(pct_covered>=40, read_count, 0)) %>% 
  group_by(contig) %>%
  summarise(total_reads=sum(read_count), count=sum(read_count>0)) %>%
  filter(total_reads>0) %>%
  pull(contig)

df = df %>% filter(contig %in% keep_these)


#now figure out what the most abundantly recruited samples are
sample_order = df %>% group_by(sample_label, ecological_zone) %>%
  summarise(total_abundance=sum(reads_per_kb_of_genome)) %>%
  group_by(ecological_zone) %>%
  arrange(desc(total_abundance), .by_group = TRUE) %>%
  pull(sample_label)

df = df %>% mutate(sample_label=factor(sample_label, levels=sample_order))

sub_df = adjust_data(df, 40) %>% filter(phage_type=='pelagiphage')
make_plot(sub_df)
ggsave('plots/GOV2-Pelagiphage-40pct.pdf', device=pdf, width=3*210, height=8*297, units='mm', limitsize = FALSE)

make_plot(sub_df)
sub_df = adjust_data(df, 70) %>% filter(phage_type=='pelagiphage')
ggsave('plots/GOV2-Pelagiphage-70pct.pdf', device=pdf, width=3*210, height=8*297, units='mm', limitsize = FALSE)


sub_df = adjust_data(df, 40) %>% filter(phage_type=='cyanophage')

make_plot(sub_df)
ggsave('plots/GOV2-Cyanophage-40pct.pdf', device=pdf, width=3*210, height=8*297, units='mm', limitsize = FALSE)

sub_df = adjust_data(df, 70) %>% filter(phage_type=='cyanophage')
make_plot(sub_df)
ggsave('plots/GOV2-Cyanophage-70pct.pdf', device=pdf, width=3*210, height=8*297, units='mm', limitsize = FALSE)

