# Is there a statistical signal that BATS viral contigs are enriched in samples with low phosphate and nitrate in GOV2?

library(tidyverse)
map_df = read_tsv('./data/Figure 2/BATS-contigs-only-GOV2-mapping-tpmean.txt') %>%
  select(c(1:9, 11:23, 25:100, 102:118, 120:136))

gov2_metadata = read_tsv('./data/Figure 2/gov2.metadata.tsv') %>%
  group_by(sample) %>% slice(1)

viral_contigs_found_elsewhere = map_df %>% 
  pivot_longer(!Contig, names_to='sample', values_to='rpkm') %>%
  group_by(Contig) %>%
  summarise(total_in_all_samples=sum(rpkm)) %>%
  filter(total_in_all_samples > 0) %>% pull(Contig)

filtered_map_df = map_df %>%
  filter(Contig %in% viral_contigs_found_elsewhere) %>%
  pivot_longer(!Contig, names_to='sample', values_to='rpkm') %>%
  mutate(presence = rpkm > 1) %>%
  left_join(gov2_metadata) %>%
  filter(ecological_zone=='TT_EPI')

ggplot(filtered_map_df, aes(ecological_zone, fill=presence)) +
  geom_bar() +
  coord_flip()

logistic = glm(presence ~ phosphate + nitrate_nitrite, data=filtered_map_df, family='binomial')
print(summary(logistic),show.residuals=TRUE)