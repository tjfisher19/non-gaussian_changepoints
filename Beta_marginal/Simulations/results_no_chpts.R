library(tidyverse)
library(ggtext)
# devtools::install_github('stefano-meschiari/latex2exp')

#############################
## Parameter estimates

source("beta_ar1_model_fit_function.R")
load("results/beta_no_change_sims_baseballParams.RData")

fits_list <- lapply(fake_series, fit_beta_ar1_chpt, approx=FALSE)
param_df <- bind_rows(lapply(fits_list, function(x) {x$Parameters}))

param_df_tall <- param_df %>%
  pivot_longer(everything()) %>%
  mutate(name = factor(name, levels=c("Means", "Precision", "AR_coef"),
                       labels=c("Mean *<span style='font-family: serif; font-size:12pt'>\u03BC</span>*",
                                "Precision *<span style='font-family: serif; font-size:12pt'>\u03BA</span>*",
                                "AR Coefficient *<span style='font-family: serif; font-size:12pt'>&phi;</span>*") ) )
                                
true_vals <- data.frame(name = c("Means", "Precision", "AR_coef"),
                        value = c(parameters$Alpha/(parameters$Alpha+parameters$Beta),
                                  parameters$Alpha+parameters$Beta,
                                  parameters$phi)) %>%
  mutate(name = factor(name, levels=c("Means", "Precision", "AR_coef"),
                       labels=c("Mean *<span style='font-family: serif; font-size:12pt'>\u03BC</span>*",
                                "Precision *<span style='font-family: serif; font-size:12pt'>\u03BA</span>*",
                                "AR Coefficient *<span style='font-family: serif; font-size:12pt'>&phi;</span>*") ) )

                          


ggplot(param_df_tall, aes(x=value, y="") ) +
  geom_vline(data=true_vals, aes(xintercept=value), linetype=2) +
  geom_boxplot(width=0.35) +
  stat_summary(geom="point", fun="mean", shape=23, fill="gray70") + 
  facet_wrap(name~., scales="free_x", nrow=1) + # , labeller=label_parsed) +
  theme_minimal() + 
  labs(title="Distribution of estimated parameters",
       subtitle="Based on estimated values from 100 simulated datasets") +
  theme(axis.title = element_blank(),
        #text = element_text(family="serif"),
        panel.spacing = unit(1.5, "lines"),
        strip.text = ggtext::element_markdown(),
        plot.title.position = "plot")

ggsave(filename="beta_ar1_parameterEstimates.png",
       width=6.5, height=2, bg="white", units="in", device=png,
       symbolfamily="serif")


############################################
## Number of changepoints found
############################################

## phi=0

load("results/beta_no_change_sims_noCorrelation.RData")
chpt_found <- sapply(full_fit_fake_series, function(x) { c(sum(x$iid_mdl_mean_model@solution), sum(x$ar1_mdl_mean_model@solution),
                                                           sum(x$iid_bic_mean_model@solution), sum(x$ar1_bic_mean_model@solution)) } )
                                                        #   sum(x$iid_mdl_shape_model@solution), sum(x$ar1_mdl_shape_model@solution),
                                                        #   sum(x$iid_bic_shape_model@solution), sum(x$ar1_bic_shape_model@solution))})
df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")
               #"MDL_iid_Shape", "MDL_ar1_Shape", "BIC_iid_Shape", "BIC_ar1_Shape")
df_tall <- df %>% pivot_longer(everything() )

no_change_results <- df_tall %>%
  group_by(name) %>%
  summarize(Mean = mean(value),
            SD = sd(value),
            PctZero = mean(value==0)*100 ) %>%
  mutate(Phi = 0)


## phi = 0.3

load("results/beta_no_change_sims_baseballParams.RData")

chpt_found <- sapply(full_fit_fake_series, function(x) { c(sum(x$iid_mdl_mean_model@solution), sum(x$ar1_mdl_mean_model@solution),
                                                           sum(x$iid_bic_mean_model@solution), sum(x$ar1_bic_mean_model@solution) ) } ) #,
                                                          # sum(x$iid_mdl_shape_model@solution), sum(x$ar1_mdl_shape_model@solution),
                                                           #sum(x$iid_bic_shape_model@solution), sum(x$ar1_bic_shape_model@solution))})

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")
               #"MDL_iid_Shape", "MDL_ar1_Shape", "BIC_iid_Shape", "BIC_ar1_Shape")

df_tall <- df %>% pivot_longer(everything() )

no_change_results <- bind_rows(no_change_results,
                               df_tall %>%
                                 group_by(name) %>%
                                 summarize(Mean = mean(value),
                                           SD = sd(value),
                                           PctZero = mean(value==0)*100 ) %>%
                                 mutate(Phi = 0.3) )


## Phi = 0.6

load("results/beta_no_change_sims_strongCorrelation.RData")

chpt_found <- sapply(full_fit_fake_series, function(x) { c(sum(x$iid_mdl_mean_model@solution), sum(x$ar1_mdl_mean_model@solution),
                                                           sum(x$iid_bic_mean_model@solution), sum(x$ar1_bic_mean_model@solution) ) } )
                                                           #sum(x$iid_mdl_shape_model@solution), sum(x$ar1_mdl_shape_model@solution),
                                                           #sum(x$iid_bic_shape_model@solution), sum(x$ar1_bic_shape_model@solution))})

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")
              # "MDL_iid_Shape", "MDL_ar1_Shape", "BIC_iid_Shape", "BIC_ar1_Shape")

df_tall <- df %>% pivot_longer(everything() )

no_change_results <- bind_rows(no_change_results,
                               df_tall %>%
                                 group_by(name) %>%
                                 summarize(Mean = mean(value),
                                           SD = sd(value),
                                           PctZero = mean(value==0)*100 ) %>%
                                 mutate(Phi = 0.6) )



## Phi = -0.3

load("results/beta_no_change_sims_baseballNegativeCorrelation.RData")

chpt_found <- sapply(full_fit_fake_series, function(x) { c(sum(x$iid_mdl_mean_model@solution), sum(x$ar1_mdl_mean_model@solution),
                                                           sum(x$iid_bic_mean_model@solution), sum(x$ar1_bic_mean_model@solution) ) } ) 
                                                      #     sum(x$iid_mdl_shape_model@solution), sum(x$ar1_mdl_shape_model@solution),
                                                      #     sum(x$iid_bic_shape_model@solution), sum(x$ar1_bic_shape_model@solution))})

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean" )
        #       "MDL_iid_Shape", "MDL_ar1_Shape", "BIC_iid_Shape", "BIC_ar1_Shape")

df_tall <- df %>% pivot_longer(everything() )

no_change_results <- bind_rows(no_change_results,
                               df_tall %>%
                                 group_by(name) %>%
                                 summarize(Mean = mean(value),
                                           SD = sd(value),
                                           PctZero = mean(value==0)*100 ) %>%
                                 mutate(Phi = -0.3) )

## Phi = -0.6

load("results/beta_no_change_sims_strongNegativeCorrelation.RData")

chpt_found <- sapply(full_fit_fake_series, function(x) { c(sum(x$iid_mdl_mean_model@solution), sum(x$ar1_mdl_mean_model@solution),
                                                           sum(x$iid_bic_mean_model@solution), sum(x$ar1_bic_mean_model@solution) ) } )
                                                        #   sum(x$iid_mdl_shape_model@solution), sum(x$ar1_mdl_shape_model@solution),
                                                        #   sum(x$iid_bic_shape_model@solution), sum(x$ar1_bic_shape_model@solution))})

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")
  #               "MDL_iid_Shape", "MDL_ar1_Shape", "BIC_iid_Shape", "BIC_ar1_Shape")

df_tall <- df %>% pivot_longer(everything() )

no_change_results <- bind_rows(no_change_results,
                               df_tall %>%
                                 group_by(name) %>%
                                 summarize(Mean = mean(value),
                                           SD = sd(value),
                                           PctZero = mean(value==0)*100 ) %>%
                                 mutate(Phi = -0.6) )


ggplot(filter(no_change_results, name %in% c("MDL_iid_Mean", "MDL_ar1_Mean") ) ) +
  geom_point(aes(x = Phi, y=Mean, color=name)) +
  geom_errorbar(aes(x=Phi, ymin=(Mean-SD), ymax=(Mean+SD)))

ggplot(filter(no_change_results, name %in% c("MDL_iid_Mean", "MDL_ar1_Mean") ) ) +
  geom_point(aes(x = Phi, y=PctZero, color=name))

no_change_results %>% filter(Phi==-0.6)

no_change_results
