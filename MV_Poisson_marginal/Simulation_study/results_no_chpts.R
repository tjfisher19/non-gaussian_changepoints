library(tidyverse)
library(ggtext)


############################################
## Number of changepoints found
############################################

## phi=0

load("results/mvpois_no_change_sims_pos_pos.RData")
chpt_found <- sapply(full_fit_fake_series, function(x) { c(sum(x$iid_mdl_mvpois_model@solution), sum(x$ar1_mdl_mvpois_model@solution),
                                                           sum(x$iid_bic_mvpois_model@solution), sum(x$ar1_bic_mvpois_model@solution)) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")
               
df_tall <- df %>% pivot_longer(everything() )

no_change_results <- df_tall %>%
  group_by(name) %>%
  summarize(Mean = mean(value),
            SD = sd(value),
            PctZero = mean(value==0)*100 ) %>%
  mutate(Phi = "PosPos")


## phi = 0.3

load("results/mvpois_no_change_sims_pos_neg.RData")

chpt_found <- sapply(full_fit_fake_series, function(x) { c(sum(x$iid_mdl_mvpois_model@solution), sum(x$ar1_mdl_mvpois_model@solution),
                                                           sum(x$iid_bic_mvpois_model@solution), sum(x$ar1_bic_mvpois_model@solution)) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

no_change_results <- bind_rows(no_change_results,
                               df_tall %>%
                                 group_by(name) %>%
                                 summarize(Mean = mean(value),
                                           SD = sd(value),
                                           PctZero = mean(value==0)*100 ) %>%
                                 mutate(Phi = "PosNeg") )


## Phi = 0.6

load("results/mvpois_no_change_sims_neg_pos.RData")

chpt_found <- sapply(full_fit_fake_series, function(x) { c(sum(x$iid_mdl_mvpois_model@solution), sum(x$ar1_mdl_mvpois_model@solution),
                                                           sum(x$iid_bic_mvpois_model@solution), sum(x$ar1_bic_mvpois_model@solution)) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")
              

df_tall <- df %>% pivot_longer(everything() )

no_change_results <- bind_rows(no_change_results,
                               df_tall %>%
                                 group_by(name) %>%
                                 summarize(Mean = mean(value),
                                           SD = sd(value),
                                           PctZero = mean(value==0)*100 ) %>%
                                 mutate(Phi = "NegPos") )



## Phi = -0.3

load("results/mvpois_no_change_sims_neg_neg.RData")

chpt_found <- sapply(full_fit_fake_series, function(x) { c(sum(x$iid_mdl_mvpois_model@solution), sum(x$ar1_mdl_mvpois_model@solution),
                                                           sum(x$iid_bic_mvpois_model@solution), sum(x$ar1_bic_mvpois_model@solution)) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean" )
        

df_tall <- df %>% pivot_longer(everything() )

no_change_results <- bind_rows(no_change_results,
                               df_tall %>%
                                 group_by(name) %>%
                                 summarize(Mean = mean(value),
                                           SD = sd(value),
                                           PctZero = mean(value==0)*100 ) %>%
                                 mutate(Phi = "NegNeg") )



ggplot(filter(no_change_results, name %in% c("MDL_iid_Mean", "MDL_ar1_Mean") ) ) +
  geom_point(aes(x = Phi, y=Mean, color=name)) +
  geom_errorbar(aes(x=Phi, ymin=(Mean-SD), ymax=(Mean+SD)))

ggplot(filter(no_change_results, name %in% c("MDL_iid_Mean", "MDL_ar1_Mean") ) ) +
  geom_point(aes(x = Phi, y=PctZero, color=name))

no_change_results <- no_change_results |>
  mutate(name=factor(name, levels=c("BIC_iid_Mean", "BIC_ar1_Mean", "MDL_iid_Mean", "MDL_ar1_Mean"))) |>
  dplyr::select(name, PctZero, Mean, Phi) |>
  arrange(name) 

no_change_results |>
  dplyr::filter(Phi=="PosPos")

no_change_results |>
  dplyr::filter(Phi=="PosNeg")

no_change_results |>
  dplyr::filter(Phi=="NegPos")

no_change_results |>
  dplyr::filter(Phi=="NegNeg")
