
library(tidyverse)
library(lpSolve)


################################
## Metric for how accurate our change point is
cpt.dist = function(C1, C2, N){
  m = length(C1)
  k = length(C2)
  
  ##Generate Cost Matrix via all paired distance
  pair = expand.grid(C1, C2)
  if(m==k){ 
    cost.mat = matrix(abs(pair[,1]-pair[,2]), 
                      nrow=m,ncol=k,byrow=T )
  }else if(m > k){ #C1 has more changepoints than C2
    cost.mat = cbind(matrix(abs(pair[,1]-pair[,2]), 
                            nrow=m,ncol=k,byrow=T ), 
                     matrix(0, nrow=m, ncol=(m-k), byrow=T))
  }else{ #C1 has less changepoints than C2, nrow < ncol
    cost.mat = rbind(matrix(abs(pair[,1]-pair[,2]), 
                            nrow=m,ncol=k,byrow=F ), 
                     matrix(0, nrow=(k-m), ncol=k, byrow= T ))
  }
  cpt.asgn = lp.assign(cost.mat, direction = "min")
  cpt.asgn$objval/N + abs(m-k)
}

############################################################
############################################################
## Single Change - In lambda vector parameter
############################################################

rep_series <- 100
N <- 50
real_chpt_config <- 21

############################################################
## Positive cross-correlation and positive correlation

load("results/mvpois_change_pts_sims_pos_pos.RData")
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100) ) } )
      
df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- df_tall %>%
  mutate(Case = "PosPos") 

############################################################
## Positive cross-correlation and negative correlation

load("results/mvpois_change_pts_sims_pos_neg.RData")
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100) ) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- bind_rows(change_results, 
                            df_tall %>%
  mutate(Case = "PosNeg") 
)

############################################################
## Negative cross-correlation and positive correlation


load("results/mvpois_change_pts_sims_neg_pos.RData")
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100) ) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- bind_rows(change_results, 
                            df_tall %>%
  mutate(Case = "NegPos") 
)

############################################################
## Negative cross-correlation and Negative correlation

load("results/mvpois_change_pts_sims_neg_neg.RData")
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100) ) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- bind_rows(change_results, 
                            df_tall %>%
  mutate(Case = "NegNeg") 
)


############################################################
############################################################
## Two Changes - In lambda vector parameter
############################################################

rep_series <- 100
N <- 50
real_chpt_config <- c(16, 26)

############################################################
## Positive cross-correlation and negative correlation

load("results/mvpois_mult_change_pts_sims_pos_neg.RData")
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100) ) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- bind_rows(change_results, 
                            df_tall %>%
                              mutate(Case = "MultPosNeg") 
)
  

############################################################
############################################################
## Three Changes - In lambda vector parameter
############################################################

rep_series <- 100
N <- 50
real_chpt_config <- c(16, 26, 41 ) 

############################################################
## Positive cross-correlation and negative correlation

load("results/mvpois_mult_change_pts_sims_pos_neg2.RData")
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mvpois_model@solution==1)+1, real_chpt_config, N=100) ) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- bind_rows(change_results, 
                            df_tall %>%
                              mutate(Case = "MultPosNeg2") 
)


change_results_plot <- change_results %>%
  filter(str_detect(name, "ar1")) %>%
  mutate(`Penalty Structure` = factor(name,
                                      levels=c("MDL_ar1_Mean", "BIC_ar1_Mean"),
                                      labels=c("MDL", "BIC") ) ) |>
  mutate(Case = factor(Case, levels=c("PosPos", "PosNeg", "NegPos", "NegNeg", "MultPosNeg", "MultPosNeg2"),
                       labels=paste("MV Poisson Configuration", 1:6)))

ggplot(change_results_plot, aes(y=`Penalty Structure`, x=value)) +
  geom_boxplot( ) + 
  stat_summary(fun="mean", geom="point", shape=23, fill="gray40") + 
  facet_wrap(~Case, dir="v", nrow=2) + 
  labs(x="Distance in estimated vs true changepoint configuration",
       title="Comparison of BIC and MDL penalty structures in change point detection") +
  theme_minimal() + 
  theme(panel.spacing = unit(1.5, "lines"),
        plot.title.position = "plot")

ggsave("mv_pois_ar1_chptFindings.png", 
       width=6.5, height=3, bg="white")
