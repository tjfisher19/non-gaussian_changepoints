#####################################################
##
## Here we explore the results of the MDL search
##    based on the joint multivariate poisson. 


## Load the data
load("./data/tropicalStormCounts.RData")
## Load the MDL Search results
load("tropicalStorm_gaFindings.RData")

source("joint_distro_code.R")

## Refit with more precision (MLEs)





nochange <- mle_mv_pois_chpt(storm_counts_wide[,-1], approx=FALSE)
mdl_model <- mle_mv_pois_chpt(storm_counts_wide[,-1], approx=FALSE,
                              chpts=(which(joint.pois.ga.mdl@solution==1)+1) )
bic_model <- mle_mv_pois_chpt(storm_counts_wide[,-1], approx=FALSE,
                              chpts=(which(joint.pois.ga.bic@solution==1)+1) )

nochange$MDL
bic_model$MDL
mdl_model$MDL
joint.pois.ga.mdl@fitnessValue

nochange$BIC
bic_model$BIC
mdl_model$BIC
joint.pois.ga.bic@fitnessValue

save(nochange, mdl_model, bic_model,
     file="best_fit_joint_models.RData")

load("best_fit_joint_models.RData")

######
## Lambda Estimates?

round(matrix(bic_model$Lambda, ncol=6, byrow=TRUE), 3)

round(matrix(bic_model$R, ncol=6, byrow=TRUE), 2)


################################################################################
##
## Plot of fitted model
##
################################################################################

bic_mod_df <- matrix(bic_model$Lambda, nrow=length(bic_model$chpts)+1, byrow=TRUE) %>%
  as.data.frame(make.names=TRUE, col.names=names(storm_counts_wide[,-1]) ) %>%
  mutate(Regime = 0:length(bic_model$chpts) )

names(bic_mod_df)[1:6] <- names(storm_counts_wide[,-1])

bic_mod_df <- bic_mod_df %>%
  pivot_longer(-Regime, names_to="BASIN", values_to="Mean_value") %>%
  mutate(Method = "BIC")

mdl_mod_df <- matrix(mdl_model$Lambda, nrow=length(mdl_model$chpts)+1, byrow=TRUE) %>%
  as.data.frame(make.names=TRUE, col.names=names(storm_counts_wide[,-1]) ) %>%
  mutate(Regime = 0:length(mdl_model$chpts) )

names(mdl_mod_df)[1:6] <- names(storm_counts_wide[,-1])

mdl_mod_df <- mdl_mod_df %>%
  pivot_longer(-Regime, names_to="BASIN", values_to="Mean_value") %>%
  mutate(Method = "MDL")

storm_basin_regimes <- bind_rows(
  storm_counts %>%
    ungroup() %>%
    arrange(SEASON) %>%
    mutate(Regime = rep(findInterval(1:nrow(storm_counts_wide), 
                                     (bic_model$chpts) ), each=6 ) ) %>%
    group_by(BASIN, Regime) %>%
    summarize(Min_year = min(SEASON)-0.5,
              Max_year = max(SEASON)+0.5) %>%
    inner_join(bic_mod_df, by=c("BASIN", "Regime")) %>%
    mutate(BASIN = factor(BASIN, levels=levels(storm_counts$BASIN))),
  storm_counts %>%
    ungroup() %>%
    arrange(SEASON) %>%
    mutate(Regime = rep(findInterval(1:nrow(storm_counts_wide), 
                                     (mdl_model$chpts) ), each=6 ) ) %>%
    group_by(BASIN, Regime) %>%
    summarize(Min_year = min(SEASON)-0.5,
              Max_year = max(SEASON)+0.5) %>%
    inner_join(mdl_mod_df, by=c("BASIN", "Regime")) %>%
    mutate(BASIN = factor(BASIN, levels=levels(storm_counts$BASIN)))
)

p_marg_segments <- ggplot(storm_counts, aes(x=SEASON, y=Counts)) + 
  geom_line(color="gray70", linewidth=0.5) + 
  geom_point(color="gray60", size=1) +
  geom_segment(data=storm_basin_regimes, 
               aes(x=Min_year, xend=Max_year,
                   y=Mean_value, yend=Mean_value, group=Regime, color=Method),
              linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Tropical Cyclones per season stratified by basin",
       subtitle=paste("Segmentation based on the multivariate Poisson model"),
       caption="Source: International Best Track Archive for Climate Stewardship\nhttps://www.ncei.noaa.gov/products/international-best-track-archive") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        plot.caption = element_text(family="mono"),
        plot.title.position = "plot",
        panel.spacing = unit(1.5, "lines"),
        legend.position="inside",
        legend.position.inside = c(0.93, 0.32),
        legend.box.background = element_rect(color="gray20") ) +
  scale_color_manual(values=c("gray35", "gray5"),
                     name="Search Criteria")
p_marg_segments

ggsave(filename="cyclone_counts_with_segmentation.png", plot=p_marg_segments,
       width=9, height=6, bg="white")



p_bic_segments <- ggplot(storm_counts, aes(x=SEASON, y=Counts)) + 
  geom_line(color="gray70", linewidth=0.5) + 
  geom_point(color="gray60", size=1) +
  geom_segment(data=dplyr::filter(storm_basin_regimes, Method=="BIC"), 
               aes(x=Min_year, xend=Max_year,
                   y=Mean_value, yend=Mean_value, group=Regime, color=Method),
               linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Tropical Cyclones per season stratified by basin",
       subtitle=paste("BIC segmentation based on the multivariate Poisson model"),
       caption="Source: International Best Track Archive for Climate Stewardship\nhttps://www.ncei.noaa.gov/products/international-best-track-archive") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        plot.caption = element_text(family="mono"),
        plot.title.position = "plot",
        panel.spacing = unit(1.5, "lines"),
        legend.position="inside",
        legend.position.inside = c(0.93, 0.32),
        legend.box.background = element_rect(color="gray20") ) +
  scale_color_manual(values=c("gray25", "gray5"),
                     name="Search Criteria", guide="none")
p_bic_segments

ggsave(filename="cyclone_counts_with_bic_segmentation.png", plot=p_bic_segments,
       width=9, height=6, bg="white")



p_mdl_segments <- ggplot(storm_counts, aes(x=SEASON, y=Counts)) + 
  geom_line(color="gray70", linewidth=0.5) + 
  geom_point(color="gray60", size=1) +
  geom_segment(data=dplyr::filter(storm_basin_regimes, Method=="MDL"), 
               aes(x=Min_year, xend=Max_year,
                   y=Mean_value, yend=Mean_value, group=Regime, color=Method),
               linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Tropical Cyclones per season stratified by basin",
       subtitle=paste("MDL segmentation based on the multivariate Poisson model"),
       caption="Source: International Best Track Archive for Climate Stewardship\nhttps://www.ncei.noaa.gov/products/international-best-track-archive") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        plot.caption = element_text(family="mono"),
        plot.title.position = "plot",
        panel.spacing = unit(1.5, "lines"),
        legend.position="inside",
        legend.position.inside = c(0.93, 0.32),
        legend.box.background = element_rect(color="gray20") ) +
  scale_color_manual(values=c("gray5"),
                     name="Search Criteria", guide="none")
p_mdl_segments

ggsave(filename="cyclone_counts_with_mdl_segmentation.png", plot=p_mdl_segments,
       width=9, height=6, bg="white")
