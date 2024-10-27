#########################################
## This code performs the residual analysis
##    for our fitted multivariate Poisson model



## Load the results and necessary packages

load("best_fit_joint_models.RData")
load("./data/tropicalStormCounts.RData")

library(tidyverse)
library(patchwork)
library(qqplotr)
library(ggtext)
library(MVN)
library(forecast)
library(WeightedPortTest)
library(tmvtnorm)


## This code has the Multivarite Hosking test (from Fisher & Robbins, 2018)
source("mvTestStats.R")


###########################################
## Get residuals...
##
## First we need a matrix with the lambda values
##   based on our segmentation

N <- dim(storm_counts_wide)[1]
mean_vals <- as.data.frame(matrix(bic_model$Lambda, ncol=6, byrow=TRUE) )
names(mean_vals) <- names(bic_model$Lambda)[1:6]
mean_vals$Counts <- c(20, 24)  ## Length of regimes
Mean_vec_df <- mean_vals %>%
  uncount(Counts)

##################################################
## This function calculates the univariate conditional truncated Normal expectation
##    It works with column vectors, so it effectively goes down the rows of our 44x6 data
##    matrix and matrix of lambda/mean values. 
## So we call this function with sapply() traversing across the 6 columns.
get_res_uni_cond_norm <- function(col=1) {
  (exp(-qnorm(ppois((as.numeric(as.data.frame(storm_counts_wide)[,(col+1)])-1), lambda=Mean_vec_df[,col]) )^2/2) -
     exp(-qnorm(ppois(as.numeric(as.data.frame(storm_counts_wide)[,(col+1)]), lambda=Mean_vec_df[,col]))^2/2) ) /
    (sqrt(2*pi)*(ppois(as.numeric(as.data.frame(storm_counts_wide)[,(col+1)]), lambda=Mean_vec_df[,col]) -
                   ppois(as.numeric(as.data.frame(storm_counts_wide)[,(col+1)])-1, lambda=Mean_vec_df[,col])))
}

##################################################
## This function calculates the multivariate conditional truncated Normal expectation
##   It works across each row vector, so effectively acrosse each column of our 44x6 data matrix
##   and matrix of lambda/mean values
## So we call this function with sapply() traversing across all 44 rows.
get_res_mv_cond_norm <- function(row=1) {
  mtmvnorm(sigma=bic_model$R, 
           lower=qnorm(ppois(as.numeric(storm_counts_wide[row,-1])-1, lambda=as.numeric(Mean_vec_df[row,1:6]) )),
           upper=qnorm(ppois(as.numeric(storm_counts_wide[row,-1]), lambda=as.numeric(Mean_vec_df[row,1:6] ) ) ) )$tmean
}

## Get the residuals
resids <- as.data.frame(t(sapply(1:N, get_res_mv_cond_norm)))
names(resids) <- names(storm_counts_wide)[-1]

resids_tall <- resids %>%
  pivot_longer(everything(),
               names_to="BASIN") %>%
  mutate(BASIN = factor(BASIN, levels=levels(storm_counts$BASIN)))


#########################
## Residual Normal QQ Plots
##  We also provide the Anderson Darling test p-value
mvn_resids <- mvn(resids)
get_qq_plot <- function(basin) {
  
  tmp <- resids_tall %>%
    filter(BASIN == basin)
  
  ad_test_df <- data.frame(
    x=0.90*(min(tmp$value)+0), 
    y=0.75*(max(tmp$value)),
    Label = paste0("AD Test<br/><i>p</i>-value = ", sprintf("%0.3f", 
                                                            as.numeric(mvn_resids$univariateNormality[str_trim(mvn(resids)$univariateNormality$Variable)==basin,4]))),
    BASIN = basin
  )
  
  ggplot() + 
    stat_qq_line(data=tmp, aes(sample=value),
                 distribution = "norm", color="gray10") +
    stat_qq_point(data=tmp, aes(sample=value),
                  distribution = "norm", color="gray40") +
    facet_wrap(~BASIN, scales="free") +
    geom_richtext(data=ad_test_df, aes(x=x, y=y, label=Label), 
                  hjust=0, color="gray20", size=3) +
    theme_minimal() +
    theme(axis.title = element_blank() ) 
}

qq_plots <- lapply(unique(storm_counts$BASIN), get_qq_plot)

wrap_plots(qq_plots)

#########
## Trying to label -- had to piece this together based on some examples online

dummy_plot <- ggplot() + labs(x="Theoretical Quantiles", y="Empirical Quantiles")

qq_plots$x_axis <- cowplot::get_plot_component(dummy_plot, "xlab-b")
qq_plots$y_axis <- cowplot::get_plot_component(dummy_plot, "ylab-l")

designs <- "
HABC
HDEF
#GGG
"

p_resid_qqplot <- wrap_plots(qq_plots) +
  plot_annotation(title="Normal Q-Q plots for residuals from each basin",
                  subtitle="Anderson Darling Test *p*-value reported for each basin",
                  theme=theme(plot.subtitle = element_markdown() ) ) +
  plot_layout(design=designs, heights = c(40, 40, 1), widths = c(1, 50, 50, 50))

p_resid_qqplot


ggsave(filename="cyclone_residuals_joint_model_qqplot.png", plot=p_resid_qqplot,
       width=8, height=5, bg="white")

############################################
## ACF Plot of marginal residuals

ggAcf(resids)

df_acf <- resids_tall %>% 
  group_by(BASIN) %>% 
  summarise(list_acf=list(acf(value, plot=FALSE, lag.max=6))) %>%
  mutate(acf_vals=purrr::map(list_acf, ~as.numeric(.x$acf))) %>% 
  dplyr::select(-list_acf) %>% 
  unnest(acf_vals) %>% 
  group_by(BASIN) %>% 
  mutate(lag=row_number() - 1) %>%
  filter(lag>0) %>%
  mutate(BASIN = factor(BASIN, levels=levels(storm_counts$BASIN)))

df_ci <- resids_tall %>% 
  group_by(BASIN) %>% 
  summarise(ci = qnorm((1 + 0.95)/2)/sqrt(n()))

df_portmanteau <- resids_tall %>%
  group_by(BASIN) %>%
  summarize(Q_stat = Weighted.Box.test(value, lag=6, type="Ljung")$statistic,
            Q_pval = Weighted.Box.test(value, lag=6, type="Ljung")$p.value) %>%
  mutate(Label = paste0('FG~Test \n', 'italic(p)-value', '==', sprintf('%0.3f', Q_pval)),
         Label = paste0("FG Test<br/><i>p</i>-value = ", sprintf("%0.3f", Q_pval)),
         x=0.95, y=0.24) 


p_resid_acf <- ggplot(df_acf, aes(x=lag, y=acf_vals)) +
  geom_hline(yintercept = 0, color="gray10") +
  geom_hline(data = df_ci, aes(yintercept = -ci), color="gray5", linetype="dotted") +
  geom_hline(data = df_ci, aes(yintercept = ci), color="gray5", linetype="dotted") +
  geom_bar(stat="identity", width=.15, color="gray40", fill="gray40") +
  geom_richtext(data=df_portmanteau, aes(x=x, y=y, label=paste(Label)), 
                hjust=0, color="gray20", size=3 ) +
  labs(x="Lag", y="Autocorrelation Function",
       title="Autocorrelation Function for residuals from each Basin",
       subtitle="Weighted Portmanteau Test *p*-value reported for each basin") +
  scale_x_continuous(breaks=1:6) + 
  facet_wrap(~BASIN) +
  theme_minimal() +
  theme(plot.subtitle = element_markdown(),
        panel.spacing = unit(1.5, "lines") )

p_resid_acf

LBHoskingTest(resids, lag=6, fitdf=0, weighted=TRUE)


ggsave(filename="cyclone_residuals_joint_model_acf.png", plot=p_resid_acf,
       width=8, height=5, bg="white")


d_vals <- sapply(1:N, function(i) {as.numeric(resids[i,] - colMeans(resids))%*%solve(cov(resids))%*%as.numeric(resids[i,] - colMeans(resids)) } )
d_vals_df <- data.frame(Sample=d_vals)
hz_test_df <- data.frame(x=3, y=13, label=paste0("HZ Test<br/><i>p</i>-value = ", sprintf("%0.3f", as.numeric(mvn_resids$multivariateNormality[3]))))

p_chi_sq <- ggplot( ) + 
  geom_abline(slope=1, intercept=0, color="gray10") +
  stat_qq(data=d_vals_df, aes(sample=Sample),
          distribution = qchisq, dparams=list(df=6), color="gray40") +
  geom_richtext(data=hz_test_df, aes(x=x, y=y, label=label), 
                hjust=0, color="gray20", size=3) +
  labs(title=expression(italic(chi)**2~" Q-Q Plot for multivariate mormality"),
       x="Theoretical Quantiles", y="Observed Mahalanobis Distance",
       subtitle = "Henze-Zirkler Test for multivariate normality reported") +
  theme_minimal() +
  theme(plot.subtitle = element_markdown(),
        plot.title.position = "plot") +
  scale_x_continuous(limits=c(0,17)) +
  scale_y_continuous(limits=c(0,17))

p_chi_sq



df_acf <- d_vals_df %>% 
  summarise(list_acf=list(acf(Sample, plot=FALSE, lag.max=6))) %>%
  mutate(acf_vals=purrr::map(list_acf, ~as.numeric(.x$acf))) %>% 
  dplyr::select(-list_acf) %>% 
  unnest(acf_vals) %>% 
  mutate(lag=row_number() - 1) %>%
  filter(lag>0)

df_ci <- d_vals_df %>% 
  summarise(ci = qnorm((1 + 0.95)/2)/sqrt(n()))

df_portmanteau <- resids_tall %>%
  summarize(Q_stat = LBHoskingTest(resids, lag=6, fitdf=0, weighted=TRUE)[1],
            Q_pval = LBHoskingTest(resids, lag=6, fitdf=0, weighted=TRUE)[2]) %>%
  mutate(Label = paste0('FR~Test \n', 'italic(p)-value', '==', sprintf('%0.3f', Q_pval)),
         Label = paste0("FR Test<br/><i>p</i>-value = ", sprintf("%0.3f", Q_pval)),
         x=0.95, y=0.24) 

p_chi_acf <- ggplot(df_acf, aes(x=lag, y=acf_vals)) +
  geom_hline(yintercept = 0, color="gray10") +
  geom_hline(data = df_ci, aes(yintercept = -ci), color="gray5", linetype="dotted") +
  geom_hline(data = df_ci, aes(yintercept = ci), color="gray5", linetype="dotted") +
  geom_bar(stat="identity", width=.15, color="gray40", fill="gray40") +
  geom_richtext(data=df_portmanteau, aes(x=x, y=y, label=paste(Label)), 
                hjust=0, color="gray20", size=3 ) +
  labs(x="Lag", y="Autocorrelation Function",
       title="ACF of Mahalanobis distance of residuals",
       subtitle="Weighted Portmanteau Test reported") +
  scale_x_continuous(breaks=1:6) + 
  theme_minimal() +
  theme(plot.subtitle = element_markdown(),
        plot.title.position = "plot")

p_chi_acf

p_chi_plots <- p_chi_sq + plot_spacer() + p_chi_acf +
  plot_layout(widths=c(4, 0.25, 4) )
p_chi_plots

ggsave(p_chi_plots, filename="cyclone_residuals_chi_sq_qqplot.png",
       width=8, height=3.5, bg="white")


