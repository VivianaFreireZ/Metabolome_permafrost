library(tidyverse)
library(ggpubr)
library(rstatix)
library(patchwork)



# Beta NTI feature 

# Barplot elemental composition
# 
palsa <- read_csv("../Assembly_feature_level_metabolome/output/Palsa_feature_TWCD_mol_8-2-22.csv")

bog <- read_csv("../Assembly_feature_level_metabolome/output/bog_feature_TWCD_mol_8-2-22.csv")

fen <- read_csv("../Assembly_feature_level_metabolome/output/fen_feature_TWCD_mol_8-2-22.csv")


palsa_elcomp <- palsa %>% 
  group_by(El_comp) %>% 
  count(Direction) %>%  
  mutate(percent_count = n/sum(n)*100) %>% 
  mutate(habitat = 'Palsa')

bog_elcomp <- bog %>% 
  group_by(El_comp) %>% 
  count(Direction) %>%  
  mutate(percent_count = n/sum(n)*100) %>% 
  mutate(habitat = 'Bog')


fen_elcomp <- fen %>% 
  group_by(El_comp) %>% 
  count(Direction) %>% 
  mutate(percent_count = n/sum(n)*100) %>% 
  mutate(habitat = 'Fen')

color_pal <- get_palette(palette = 'RdBu', 5)
color_pal[3] <- '#F9FAE0'
names(color_pal) <- c('Sig. Divergence', 'Divergence', 'Insignificant', 'Convergence', 'Sig. Convergence')

bar_table <- rbind(bog_elcomp, palsa_elcomp, fen_elcomp) %>% 
  filter(!El_comp == 'CH') %>% 
  filter(!El_comp == 'CHS') %>% 
  filter(!El_comp == 'CHP') %>% 
  filter(!El_comp == 'CHSP')


bar_table_t <- bar_table %>% 
  mutate(Direction = str_replace_all(Direction, "Sig. High", "Sig. Divergence")) %>% 
  mutate(Direction = str_replace_all(Direction, "High", "Divergence")) %>% 
  mutate(Direction = str_replace_all(Direction, "Low", "Convergence")) %>% 
  mutate(Direction = str_replace_all(Direction, "Sig. Low", "Sig. Convergence"))

bar_plot <- bar_table_t %>% 
  mutate(habitat = factor(habitat, levels = c('Palsa', 'Bog', 'Fen')),
         Direction = factor(Direction, levels = c('Sig. Divergence', 'Divergence', 'Insignificant', 'Convergence', 'Sig. Convergence'))) %>% 
  ggplot(aes(x = habitat, y = percent_count, fill = Direction)) +
  geom_col() +
  scale_fill_manual(values = color_pal)+
  labs(title = "Beta NTI feature - elemental composition",
       x = "",
       y = "Percentage (%)",
       fill = "Contribution")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ El_comp)+
  theme(legend.position = "bottom")


bar_plot



## Boxplot per elemental composition

report <- read_csv("../Assembly_Metabolome/output/Processed_Abisko_Data.csv") %>% 
  pivot_longer(!Mass, names_to = 'SampleID', values_to = 'abundance')

mol <- read_csv("../Assembly_Metabolome/output/Processed_Abisko_Mol.csv")

metadata <- read_csv('../Assembly_Metabolome/input/metadata_updated.csv')

metadata$Habitat <- factor(metadata$Habitat, levels = c("Palsa", "Bog", "Fen"))

my_colors <- c('Palsa' = '#703C1B', 'Bog' = '#058000', 'Fen' = '#0001FF')

df <- left_join(report, mol, by = "Mass") %>% 
  filter(!is.na(El_comp),
         abundance > 0) %>% 
  select(El_comp, SampleID) %>% 
  group_by(SampleID) %>% 
  count(El_comp) %>% 
  mutate(perc_counts = n/sum(n) *100) %>% 
  left_join(metadata, by = 'SampleID')

stat_table <- df %>% 
  filter(El_comp != 'CH',
         El_comp != 'CHNSP',
         El_comp != 'CHP',
         El_comp != 'CHSP',
         El_comp != 'CHS') %>% 
  group_by(El_comp) %>% 
  wilcox_test(perc_counts ~ Habitat) %>% 
  add_y_position()

stat_table[1, 'y.position'] <- 30
stat_table[10, 'y.position'] <- 1
stat_table[16, 'y.position'] <- 3.5
stat_table[17, 'y.position'] <- 4
stat_table[19, 'y.position'] <- 15
stat_table[20, 'y.position'] <- 20

box_plot <- df %>% 
  filter(El_comp != "CHP", 
         El_comp !="CHSP",
         El_comp != "CHP", 
         El_comp !="CHNSP",
         El_comp !="CH",
         El_comp !="CHS") %>% 
  ggplot()+
  geom_boxplot(aes(x = Habitat, y = perc_counts, fill = Habitat)) +
  facet_wrap(~El_comp, scales = 'free_y')+
  scale_fill_manual(values = my_colors) +
  stat_pvalue_manual(data = stat_table, hide.ns = TRUE)+
  labs(title = "Relative richness of elemental compositions", 
       y = "Relative richness (%)") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom")

box_plot


# Clusters - elemental composition


habitat <- c('Bog', 'Fen', 'Palsa')

bnti_counts_files <- list.files(path = '../Metabolite_clusters/output/',
                                pattern = 'rep_features',
                                recursive = TRUE,
                                full.names = TRUE)

bnti_counts_files <- bnti_counts_files[str_detect(bnti_counts_files, 'bnti', negate = TRUE)]

names(bnti_counts_files) <- habitat

bnti_counts <- imap(bnti_counts_files, function(x, y){
  df <- read_csv(x) %>% 
    group_by(cluster, Direction) %>% 
    count() %>% 
    group_by(cluster) %>% 
    mutate(perc = n/sum(n) * 100,
           Habitat = y)
})

rep_df <- imap(bnti_counts_files, function(x, y){
  df <- read_csv(x) %>% 
    mutate(Habitat = y)
})

classes <-  sort(unique(rep_df$Palsa$bs2_class))

class_colors <- set_names(get_palette(palette = 'Paired', 10),
                          nm = classes)


el_plots <- imap(rep_df, function(x, y){
  plot <- x %>% 
    group_by(cluster, Habitat) %>% 
    count(bs2_class) %>% 
    mutate(bs2_class = factor(bs2_class, levels = classes),
                              perc = n / sum(n)) %>% 
    ggplot() +
    geom_col(aes(x = cluster,
                 y = perc,
                 fill = bs2_class),
             color = 'black') +
    labs(fill = 'Elemental Composition',
         title = y) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = class_colors, drop = FALSE) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = 'bold', hjust = 0.5))
})

# Ordination plots















final <- (bar_plot + box_plot) / (el_plots$Palsa + el_plots$Bog + el_plots$Fen + 
                                    plot_layout(guides = 'collect')) +
  
  plot_annotation(tag_levels = list(c('A', 'B', 'C'))) &
  theme(legend.position = 'bottom')

final

ggsave("output/Supplementary_fig_1.png", final, dpi = 300, height = 12, width = 15)

