###########################################################
###########################################################
## 
## This code performs the analysis of the baseball
##    data provided in the paper.
## We effectively do the following
##   - Read in the data and do some processing
##   - a little EDA and plots for the paper
##   - the optimal segmentations are found via a 
##     genetic algorithm search
##   - We analyze our fitted models
##   - Perform a residual analysis/goodness-of-fit
##

library(tidyverse)
library(GA)
library(patchwork)

library(moments)
library(fitdistrplus)
library(betareg)

library(WeightedPortTest)
library(ggtext)
library(qqplotr)
library(nortest)

## The function that fits our model
setwd("../Code_implementation/")
source("beta_ar1_model_ga_search_functions.R")
setwd("../baseball_analysis/")

###########################################################
###########################################################
## Data processing
###########################################################
###########################################################

batting <- read_csv("./lahman_1871-2023_csv/Batting.csv")

yearly_players <- batting %>%
  group_by(yearID, playerID) %>%
  summarize(AB = sum(AB) ) %>%
  group_by(playerID) %>%
  mutate(PrevYear = lag(yearID) ) %>%
  group_by(yearID) %>%
  summarize(Pct_return = mean(!is.na(PrevYear))*100) %>%
  filter(yearID >= 1920)

batting %>%
  filter(yearID >= 1920) %>%
  group_by(yearID, playerID) %>%
  summarize(AB = sum(AB) ) %>%
  filter(AB > 0 ) %>%
  group_by(playerID) %>%
  summarize(Season = n() ) %>%
  ungroup() %>%
  summarize(AvgSeasons = mean(Season),
            SD_Seasons = sd(Season),
            MedianSeasons = median(Season))


yearly_batting <- batting %>%
 # filter(teamID=="BAL") %>%
  group_by(yearID) %>%
  summarize(Most_HR = max(HR),
            Most_H = max(H),
            Most_3B = max(`3B`),
            H = sum(H),
            AB = sum(AB),
            HR = sum(HR),
            `3B` = sum(`3B`)
  ) %>%
  mutate(Avg = H/AB,
         HR_avg = HR/AB) %>%
  filter(yearID >= 1920)

###########################################################
###########################################################
##  Plots of data, and player continuation
###########################################################
###########################################################

p_players <- ggplot(yearly_players, aes(x=yearID, y=Pct_return) ) +
  geom_line(color="gray50", linewidth=0.5) + 
  geom_point(color="gray40", size=1) +
  scale_x_continuous(breaks=seq(1920, 2025, 15), 
                     minor_breaks = seq(1915, 2025, 5)) +
  scale_y_continuous(breaks=seq(65, 95, 10),
                     minor_breaks = seq(67.4, 92.5, 2.5), limits=c(65,95)) + 
  labs(x="Season", y="Percent (%)",
       title="Batters returning from previous season",
       caption="Source: Lahman's Baseball Database, 1871-2023") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        plot.caption = element_text(family="mono"),
       # plot.title = element_text(size=12),
        plot.title.position = "plot" )


p_hr <- ggplot(yearly_batting, aes(x=yearID, y=HR_avg) ) +
  geom_line(color="gray50", linewidth=0.5) + 
  geom_point(color="gray40", size=1) +
  scale_x_continuous(breaks=seq(1920, 2025, 15), 
                     minor_breaks = seq(1915, 2025, 5)) +
  scale_y_continuous(breaks=seq(0.01, 0.04, 0.01),
                     minor_breaks = seq(0.0075, 0.0425, 0.0025),
                     labels=scales::number_format(accuracy = 0.001)) + 
  labs(x="Season", y="Proportion of HR (HR/AB)",
       title="Homerun \"Average\" during the \'Live Ball\' era",
       caption="") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        plot.caption = element_text(family="mono"),
       # plot.title = element_text(size=12),
        plot.title.position = "plot" )


p_hr + plot_spacer() + p_players + plot_layout(widths = c(4, 0.5 ,4))

ggsave("baseball_player_hr_props.png", width=8, height=3, bg="white" )

ggplot(yearly_batting, aes(x=HR_avg) ) + 
  geom_density(fill="gray70", alpha=0.3, linewidth=1) +
  scale_x_continuous(breaks=seq(0.01, 0.04, 0.01),
                     minor_breaks = seq(0.0075, 0.0425, 0.0025),
                     labels=scales::number_format(accuracy = 0.001)) +
  labs(x="Proportion of AB resulting in a HR",
       y="Empirical Density") +
  theme_minimal()

###########################################################
###########################################################
## Baseline calculations, checking things
###########################################################
###########################################################

## Baseline comparisons (Quick sanity check as our method
##    should have a larger log-likelihood)

fit_mle <- fitdist(yearly_batting$HR_avg, "beta")
fit_reg <- betareg(yearly_batting$HR_avg ~ 1)

fit_mle$loglik
fit_reg$loglik
fit_beta_ar1_chpt(yearly_batting$HR_avg, type="shape")$loglik
fit_beta_ar1_chpt(yearly_batting$HR_avg, type="mean")$loglik
## Less accurate approximation (uses MMEs)
fit_beta_ar1_chpt(yearly_batting$HR_avg, type="shape", approx=TRUE)$loglik
fit_beta_ar1_chpt(yearly_batting$HR_avg, type="mean", approx=TRUE)$loglik


TIME <- 1:dim(yearly_batting)[1]
fit_trend <- betareg(yearly_batting$HR_avg ~ TIME )
fit_trend$loglik
fit_beta_ar1_chpt(yearly_batting$HR_avg, type="trend")$loglik

###########################################################
###########################################################
## Change point search for changes in mean proportion
###########################################################
###########################################################

min.regime <- 8

N <- nrow(yearly_batting)
N*3
PopSize <- 350
suggest_matrix <- matrix(0,
                         ncol= N-1,
                         nrow=PopSize)
for(i in 2:(NCOL(suggest_matrix)-2*min.regime+3) ) {
  suggest_matrix[i, min.regime + i-2] <- 1
}
for(i in (NCOL(suggest_matrix)-min.regime+2):PopSize) {
  suggest_matrix[i, sample(min.regime:((NCOL(suggest_matrix)-min.regime+1)/2), size=1)] <- 1
  suggest_matrix[i, sample(((NCOL(suggest_matrix)-min.regime+1)/2+1):(NCOL(suggest_matrix)-min.regime+1), size=1)] <- 1
}

ga.search.hr_avg.mean.min8.bic <- ga(type="binary", 
                                     fitness=ga_search_function_beta_ar1_chpt_bic,
                                     nBits=N-1,
                                     series= yearly_batting$HR_avg,
                                     model_type="mean",
                                     approx=FALSE,
                                     min.regime=min.regime,
                                     max.chpts=floor(N/min.regime),
                                     popSize=PopSize,
                                     suggestions = suggest_matrix,
                                     parallel = 8,
                                     maxiter=5000,
                                     run=500)

ga.search.hr_avg.mean.min8.mdl <- ga(type="binary", 
                                     fitness=ga_search_function_beta_ar1_chpt_mdl,
                                     nBits=N-1,
                                     series= yearly_batting$HR_avg,
                                     model_type="mean",
                                     approx=FALSE,
                                     min.regime=min.regime,
                                     max.chpts=floor(N/min.regime),
                                     popSize=PopSize,
                                     suggestions = suggest_matrix,
                                     parallel = 8,
                                     maxiter=5000,
                                     run=500)

save(ga.search.hr_avg.mean.min8.bic, ga.search.hr_avg.mean.min8.mdl,
     file="baseball_ga_search_results.RData")

load("baseball_ga_search_results.RData")
(which(ga.search.hr_avg.mean.min8.bic@solution==1)+1)
(which(ga.search.hr_avg.mean.min8.mdl@solution==1)+1)

###########################################################
###########################################################
## Compare model fits
###########################################################
###########################################################

ga.search.hr_avg.mean.min8.bic@fitnessValue
ga.search.hr_avg.mean.min8.mdl@fitnessValue

hr_avg_mean_fit_bic <- fit_beta_ar1_chpt(yearly_batting$HR_avg, chpts=(which(ga.search.hr_avg.mean.min8.bic@solution==1)+1), type="mean")
hr_avg_mean_fit_mdl <- fit_beta_ar1_chpt(yearly_batting$HR_avg, chpts=(which(ga.search.hr_avg.mean.min8.mdl@solution==1)+1), type="mean")

1920 + hr_avg_mean_fit_bic$chpts
1920 + hr_avg_mean_fit_mdl$chpts

hr_avg_mean_fit_bic$loglik
hr_avg_mean_fit_mdl$loglik

hr_avg_mean_fit_bic$BIC
hr_avg_mean_fit_bic$MDL
hr_avg_mean_fit_mdl$BIC
hr_avg_mean_fit_mdl$MDL

fit_beta_ar1_chpt(yearly_batting$HR_avg, type="trend")$BIC
fit_beta_ar1_chpt(yearly_batting$HR_avg, type="trend")$MDL

fit_beta_ar1_chpt(yearly_batting$HR_avg, type="trend")
fit_beta_ar1_chpt(yearly_batting$HR_avg, type="mean")
hr_avg_mean_fit_bic
hr_avg_mean_fit_mdl

###########################################################
###########################################################
## Plot of BIC & MDL fitted models
###########################################################
###########################################################

yearly_batting_fits <- bind_rows(
  yearly_batting %>%
    arrange(yearID) %>%
    mutate(RegimeMean = hr_avg_mean_fit_bic$Regimes,
           AvgHR_mean = hr_avg_mean_fit_bic$fitted.values) %>%
    group_by(RegimeMean) %>%
    slice(c(1,n())) %>%
    mutate(yearID = ifelse(row_number()==1, yearID-0.5, yearID+0.5) ) %>%
    dplyr::select(yearID, AvgHR_mean, RegimeMean) %>%
    mutate(Method = "BIC"),
  yearly_batting %>%
    arrange(yearID) %>%
    mutate(RegimeMean = hr_avg_mean_fit_mdl$Regimes,
           AvgHR_mean = hr_avg_mean_fit_mdl$fitted.values) %>%
    group_by(RegimeMean) %>%
    slice(c(1,n())) %>%
    mutate(yearID = ifelse(row_number()==1, yearID-0.5, yearID+0.5) ) %>%
    dplyr::select(yearID, AvgHR_mean, RegimeMean) %>%
    mutate(Method = "MDL") 
) %>%
  mutate(Grouping=paste0(Method, RegimeMean) )

ggplot(yearly_batting, aes(x=yearID, y=HR_avg) ) +
  geom_line(color="gray70", linewidth=0.5) + 
  geom_point(color="gray60", size=1) +
  scale_x_continuous(breaks=seq(1920, 2025, 15), 
                     minor_breaks = seq(1915, 2025, 5)) +
  scale_y_continuous(breaks=seq(0.01, 0.04, 0.01),
                     minor_breaks = seq(0.0075, 0.0425, 0.0025),
                     labels=scales::number_format(accuracy = 0.001)) + 
  labs(x="Season", y="Proportion of HR (HR/AB)",
       title="Homerun \"Average\" during the \'Live Ball\' era",
       subtitle="Segmentations with 8-year minimum regime length",
       caption="Source: Lahman's Baseball Database, 1871-2023") +
  theme_minimal() + 
  geom_line(data=yearly_batting_fits, 
               aes(x=yearID, y=AvgHR_mean, group=Grouping, color=Method), 
               linewidth=1.25) +
  theme(axis.title.x = element_blank(),
        plot.caption = element_text(family="mono"),
        plot.title.position = "plot",
        legend.position="inside",
        legend.position.inside = c(0.15, 0.7),
        legend.box.background = element_rect(color="gray20")) +
  scale_color_manual(values=c("gray35", "gray5"),
                     name="Search Criteria")

ggsave("baseball_hr_props.png",
       width=8, height=3, bg="white")

###########################################################
###########################################################
## Exploring the BIC selected model
##     Model Parameters & Residuals analysis
###########################################################
###########################################################

hr_avg_mean_fit_bic$Parameters


mu_vec <- rep(hr_avg_mean_fit_bic$Parameters$Means,
              hr_avg_mean_fit_bic$Regime_lengths)

alpha_vec <- mu_vec*hr_avg_mean_fit_bic$Parameters$Precision
beta_vec <- (1-mu_vec)*hr_avg_mean_fit_bic$Parameters$Precision

Z_series <- qnorm(pbeta(yearly_batting$HR_avg, shape1=alpha_vec, shape2=beta_vec) )

###########################################################
## ACF and PACF of Z_series
###########################################################

correlation_z_series <- data.frame(
  lag = 1:10,
  ACF = acf(Z_series, plot=FALSE, lag.max = 10)$acf[-1],
  PACF = pacf(Z_series, plot=FALSE, lag.max=10)$acf
)

cor_z_series_long <- correlation_z_series |>
  pivot_longer(-lag, names_to="Measure") |>
  mutate(Measure = ifelse(Measure=="ACF", 
                          "Autocorrelation Function", 
                          "Partial Autocorrelation"))

df_ci <- as.data.frame(Z_series) |> 
  summarise(ci = qnorm((1 + 0.95)/2)/sqrt(n()))

df_portmanteau <- bind_rows(
  data.frame(value = Z_series) %>%
    summarize(Q_stat = Weighted.Box.test(value, lag=10, type="Ljung")$statistic,
              Q_pval = Weighted.Box.test(value, lag=10, type="Ljung")$p.value) %>%
    mutate(Measure = "Autocorrelation Function",
           Label = paste0('FG~Test \n', 'italic(p)-value', '==', sprintf('%0.3f', Q_pval)),
           Label = paste0("FG Test<br/><i>p</i>-value = ", sprintf("%0.3f", Q_pval)),
           x=6, y=0.4),
  data.frame(value = Z_series) %>%
    summarize(Q_stat = Weighted.Box.test(value, lag=10, type="Monti")$statistic,
              Q_pval = Weighted.Box.test(value, lag=10, type="Monti")$p.value) %>%
    mutate(Measure = "Partial Autocorrelation",
           Label = paste0('FG~Test \n', 'italic(p)-value', '==', sprintf('%0.3f', Q_pval)),
           Label = paste0("FG Test<br/><i>p</i>-value = ", sprintf("%0.3f", Q_pval)),
           x=6, y=0.4) 
)


ggplot(cor_z_series_long, aes(x=lag, y=value) ) +
  geom_hline(yintercept=0, color="gray10") +
  geom_hline(data = df_ci, aes(yintercept = -ci), color="gray5", linetype="dotted") +
  geom_hline(data = df_ci, aes(yintercept = ci), color="gray5", linetype="dotted") +
  geom_col(width=.5, color="gray40", fill="gray40") +
  geom_richtext(data=df_portmanteau, aes(x=x, y=y, label=paste(Label)), 
                hjust=0, color="gray20", size=3 ) +
  facet_wrap(~Measure) +
  scale_x_continuous(breaks=1:10, minor_breaks = NULL) +
  labs(x="Lag",
       y="Correlation Value",
       title="Correlograms of latent normal series from the baseball marginal Beta fit",
       subtitle="Results of Weighted Portmanteau Tests reported") +
  theme_minimal() +
  theme(panel.spacing = unit(1.5, "lines"),
        plot.subtitle = element_markdown(),
        plot.title.position = "plot") 
  

ggsave("baseball_z_series_acf.png", width=8, height=3, bg="white" )

###########################################################
## Residuals after modeling the AR(1)
###########################################################

phi_est <- hr_avg_mean_fit_bic$Parameters$AR_coef
resid <- Z_series[2:N] - phi_est*Z_series[1:(N-1)]

correlation_resid <- data.frame(
  lag = 1:10,
  ACF = acf(resid, plot=FALSE, lag.max = 10)$acf[-1],
  Method = "Autocorrelation Function"
)


df_ci <- as.data.frame(resid) |> 
  summarise(ci = qnorm((1 + 0.95)/2)/sqrt(n()))

df_portmanteau <- 
  data.frame(value = Z_series) %>%
    summarize(Q_stat = Weighted.Box.test(resid, lag=10, type="Ljung", fitdf=1)$statistic,
              Q_pval = Weighted.Box.test(resid, lag=10, type="Ljung", fitdf=1)$p.value) %>%
    mutate(Measure = "Autocorrelation Function",
           Label = paste0('FG~Test \n', 'italic(p)-value', '==', sprintf('%0.3f', Q_pval)),
           Label = paste0("FG Test<br/><i>p</i>-value = ", sprintf("%0.3f", Q_pval)),
           x=6, y=0.15)


p_resid_acf <- ggplot(correlation_resid, aes(x=lag, y=ACF) ) +
  geom_hline(yintercept=0, color="gray10") +
  geom_hline(data = df_ci, aes(yintercept = -ci), color="gray5", linetype="dotted") +
  geom_hline(data = df_ci, aes(yintercept = ci), color="gray5", linetype="dotted") +
  geom_col(width=.5, color="gray40", fill="gray40") +
  geom_richtext(data=df_portmanteau, aes(x=x, y=y, label=paste(Label)), 
                hjust=0, color="gray20", size=3 ) +
  facet_wrap(~Measure) +
  theme_minimal() +
  theme(panel.spacing = unit(0.5, "in") ) + 
  scale_x_continuous(breaks=1:10, minor_breaks = NULL) +
  labs(y="Correlation Value",
       x="Lag")

tmp <- data.frame(value=resid) |>
  mutate(Measure = "Normal Q-Q Plot")

ad_test_df <- data.frame(
  x=0.90*(min(tmp$value)+0), 
  y=0.75*(max(tmp$value)),
  Label = paste0("AD Test<br/><i>p</i>-value = ", sprintf("%0.3f", 
                                                          ad.test(resid)$p.value))
)

p_qq <- ggplot() + 
  stat_qq_line(data=tmp, aes(sample=value),
               distribution = "norm", color="gray10") +
  stat_qq_point(data=tmp, aes(sample=value),
                distribution = "norm", color="gray40") +
  geom_richtext(data=ad_test_df, aes(x=x, y=y, label=Label), 
                hjust=0, color="gray20", size=3) +
  facet_wrap(~Measure) + 
  theme_minimal() +
  labs(x="Theoretical Quantiles", 
       y="Empirical Quantiles")
  


p_resid_acf + plot_spacer() + p_qq + plot_layout(widths = c(4, 0.15 ,4)) +
  plot_annotation(title="Correlogram and Q-Q plot of residuals from the baseball home run fit",
                  subtitle="Results of Weighted Portmanteau and Anderson-Darling tests reported") +
  theme(plot.subtitle = element_markdown())

ggsave("baseball_resid_plots.png", width=8, height=3, bg="white" )


