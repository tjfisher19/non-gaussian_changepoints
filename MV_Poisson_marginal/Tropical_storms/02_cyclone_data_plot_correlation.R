library(tidyverse)
library(patchwork)

load("./data/tropicalStormCounts.RData")

storm_counts %>%
  group_by(BASIN) %>%
  summarize(Mean = mean(Counts),
            Var = var(Counts))
global_counts <- storm_counts %>%
  group_by(SEASON) %>%
  summarize(Counts = sum(Counts) ) %>%
  ungroup() 

global_counts %>%
  summarize(Mean = mean(Counts),
            Var = var(Counts))

##################
## Counts in time

p_counts <- ggplot(storm_counts, aes(x=SEASON, y=Counts, group=BASIN)) +
  geom_line(color="gray50", linewidth=0.5) + 
  geom_point(color="gray40", size=1) +
  facet_wrap(~BASIN) +
  labs(title="Tropical Cyclones stratified by basin",
      # subtitle=paste("All storms achieving at least 34 knots from 1980 -", max(storm_counts$SEASON)),
       caption="Source: International Best Track Archive for Climate Stewardship\nhttps://www.ncei.noaa.gov/products/international-best-track-archive") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        plot.caption = element_text(family="mono"),
        plot.title.position = "plot",
        panel.spacing = unit(1.5, "lines"))

p_counts
# ggsave(filename="cyclone_counts.png", plot=p_counts,
#        width=8, height=6, bg="white")


p_global_counts <- ggplot(global_counts, aes(x=SEASON, y=Counts)) +
  geom_line(color="gray50", linewidth=0.5) + 
  geom_point(color="gray40", size=1) +
  labs(title="Tropical Cyclones per season",
       subtitle=paste("All storms achieving at least 34 knots from 1980 -", max(storm_counts$SEASON)) ) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        plot.caption = element_text(family="mono"),
        plot.title.position = "plot")

p_global_counts

# ggsave(filename="cyclone_counts_global.png", plot=p_global_counts,
#        width=8, height=4.5, bg="white")


p_count_paper <- p_global_counts / p_counts +
  plot_layout(heights=c(1,2.5))

ggsave(filename="cyclone_counts_plot.png", p_count_paper,
       width=9, height=7, bg="white")


#######################
## Correlations of storms counts
##   by basin
cor_mat <- cor(storm_counts_wide[,-1])

cor_mat

#######################
## Heatmap plot of correlations
cor_mat[upper.tri(cor_mat, diag=TRUE)] <- NA
P_val_mat <- Hmisc::rcorr(as.matrix(storm_counts_wide[,-1]))$P
P_val_mat[upper.tri(P_val_mat, diag=TRUE)] <- NA

cor_pval_df <- as.data.frame(P_val_mat) %>%
  mutate(Basin1 = rownames(P_val_mat)) %>%
  pivot_longer(-Basin1, values_to = "p_value")

cor_heatmap_df <- as.data.frame(cor_mat) %>%
  mutate(Basin1 = rownames(cor_mat)) %>%
  pivot_longer(-Basin1) %>%
  left_join(cor_pval_df, by=c("Basin1", "name")) %>%
  mutate(value = round(value, 2),
         Basin1 = factor(Basin1, levels=levels(storm_counts$BASIN)),
         name = factor(name, levels=levels(storm_counts$BASIN)),
         pval = round(p_value, 3),
         sig = case_when(pval <= 0.01 ~ "**",
                         pval <= 0.05 ~ "*",
                         .default = ""), 
         bf_sig = ifelse(pval <=0.05, 2, 1),
         p_val_label = paste0("(", sprintf("%0.3f", pval), sig, ")")) %>%
  drop_na()

p_cor <- ggplot(cor_heatmap_df, aes(x=Basin1, y=name)) +
  geom_tile(aes(fill=value), color="white") +
  scale_fill_gradient2(high="firebrick", low="royalblue4", mid="gray99", 
                       midpoint = 0, limits=c(-0.6,0.6))+
  geom_text(aes(label=sprintf("%0.2f", round(value, digits = 2)), fontface=bf_sig),
            color="gray10", nudge_y=0.15) +
  geom_text(aes(label=p_val_label, fontface=bf_sig), 
            color="gray10", nudge_y=-0.15, size=2) +
  theme_minimal() +
  theme(legend.position="none",
        axis.title= element_blank(),
        axis.text.x=element_text(angle=0, hjust = 0.5, vjust=0.85))

p_cor

ggsave(filename = "cyclone_count_correlations.png", plot=p_cor,
       width=8, height=5, bg="white")



