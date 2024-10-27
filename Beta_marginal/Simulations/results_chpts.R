
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


#########################################
## Single Change - Mean only

load("results/beta_change_sims_singleChangePoint.RData")

rep_series <- 100
N <- 100

fake_series_df1 <- data.frame(Series=unlist(fake_series),
                              Num = rep(1:rep_series, each=N),
                              Time = rep(1:N, rep_series),
                              Config = "Beta Configuration 1")

mu_series1 <- data.frame(Mu=parameters$Alpha/(parameters$Alpha+parameters$Beta),
                         SD=sqrt(with(parameters, Alpha*Beta/((Alpha+Beta)^2*(Alpha+Beta+1)) ) ),
                         Time = 1:N) %>%
  mutate(Group=as.factor(Mu)) %>%
  group_by(Group) %>%
  slice(c(1,n())) %>%
  mutate(New_Time = ifelse(row_number()%%2==0, Time+0.4, Time-0.4),
         Config = "Beta Configuration 1")


real_chpt_config <- which(diff(parameters$Alpha)!=0) + 1

real_chpt_config
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mean_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mean_model@solution==1)+1, real_chpt_config, N=100) ) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- df_tall %>%
  mutate(Case = "Beta Configuration 1") 


##################################
## Single change - mean and precision

load("results/beta_change_sims_meanPrecisionChangePoint.RData")

fake_series_df2 <- data.frame(Series=unlist(fake_series),
                              Num = rep(1:rep_series, each=N),
                              Time = rep(1:N, rep_series),
                              Config = "Beta Configuration 2")

mu_series2 <- data.frame(Mu=parameters$Alpha/(parameters$Alpha+parameters$Beta),
                         SD=sqrt(with(parameters, Alpha*Beta/((Alpha+Beta)^2*(Alpha+Beta+1)) ) ),
                         Time = 1:N) %>%
  mutate(Group=as.factor(SD)) %>%
  group_by(Group) %>%
  slice(c(1,n())) %>%
  mutate(New_Time = ifelse(row_number()%%2==0, Time+0.4, Time-0.4),
         Config = "Beta Configuration 2")


real_chpt_config <- which(diff(parameters$Alpha)!=0) + 1
real_chpt_config
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mean_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mean_model@solution==1)+1, real_chpt_config, N=100) ) } )

df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- bind_rows(change_results,
                            df_tall %>%
                              mutate(Case = "Beta Configuration 2") )

#############################################
## Two changes - mean only

load("results/beta_change_sims_twoChangePoints.RData")

fake_series_df3 <- data.frame(Series=unlist(fake_series),
                              Num = rep(1:rep_series, each=N),
                              Time = rep(1:N, rep_series),
                              Config = "Beta Configuration 3")

mu_series3 <- data.frame(Mu=parameters$Alpha/(parameters$Alpha+parameters$Beta),
                         SD=sqrt(with(parameters, Alpha*Beta/((Alpha+Beta)^2*(Alpha+Beta+1)) ) ),
                         Time = 1:N) %>%
  mutate(Group=as.factor(SD)) %>%
  group_by(Group) %>%
  slice(c(1,n())) %>%
  mutate(New_Time = ifelse(row_number()%%2==0, Time+0.4, Time-0.4),
         Config = "Beta Configuration 3")


real_chpt_config <- which(diff(parameters$Alpha)!=0) + 1
real_chpt_config
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mean_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mean_model@solution==1)+1, real_chpt_config, N=100) ) } )
                                                        
df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- bind_rows(change_results,
                            df_tall %>%
                              mutate(Case = "Beta Configuration 3") )

#################################################
## Two changes, mean and precision

load("results/beta_change_sims_twoMeanPrecisionChangePoints.RData")

fake_series_df4 <- data.frame(Series=unlist(fake_series),
                                       Num = rep(1:rep_series, each=N),
                                       Time = rep(1:N, rep_series),
                                       Config = "Beta Configuration 4")

mu_series4 <- data.frame(Mu=parameters$Alpha/(parameters$Alpha+parameters$Beta),
                        SD=sqrt(with(parameters, Alpha*Beta/((Alpha+Beta)^2*(Alpha+Beta+1)) ) ),
                        Time = 1:N) %>%
  mutate(Group=as.factor(SD)) %>%
  group_by(Group) %>%
  slice(c(1,n())) %>%
  mutate(New_Time = ifelse(row_number()%%2==0, Time+0.4, Time-0.4),
         Config = "Beta Configuration 4")


real_chpt_config <- which(diff(parameters$Alpha)!=0) + 1
real_chpt_config
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mean_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mean_model@solution==1)+1, real_chpt_config, N=100) ) } )
                                                       
df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- bind_rows(change_results,
                            df_tall %>%
                              mutate(Case = "Beta Configuration 4") )


###############################################
## Four changes - mean only

load("results/beta_change_sims_fourChangePoints.RData")

fake_series_df5 <- data.frame(Series=unlist(fake_series),
                              Num = rep(1:rep_series, each=N),
                              Time = rep(1:N, rep_series),
                              Config = "Beta Configuration 5")

mu_series5 <- data.frame(Mu=parameters$Alpha/(parameters$Alpha+parameters$Beta),
                         SD=sqrt(with(parameters, Alpha*Beta/((Alpha+Beta)^2*(Alpha+Beta+1)) ) ),
                         Time = 1:N) %>%
  mutate(Diff = SD - lag(SD, default=0),
         Diff = ifelse(Diff==0, NA, Diff)) %>%
  fill(Diff) %>%
mutate(Group=as.factor(Diff)) %>%
  group_by(Group) %>%
  slice(c(1,n())) %>%
  mutate(New_Time = ifelse(row_number()%%2==0, Time+0.4, Time-0.4),
         Config = "Beta Configuration 5")


real_chpt_config <- which(diff(parameters$Alpha)!=0) + 1
real_chpt_config
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mean_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mean_model@solution==1)+1, real_chpt_config, N=100) ) } )


df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- bind_rows(change_results,
                            df_tall %>%
                              mutate(Case = "Beta Configuration 5") )

############################################
## Four changes, mean and precision

load("results/beta_change_sims_fourMeanPrecisionChangePoints.RData")

fake_series_df6 <- data.frame(Series=unlist(fake_series),
                              Num = rep(1:rep_series, each=N),
                              Time = rep(1:N, rep_series),
                              Config = "Beta Configuration 6")

mu_series6 <- data.frame(Mu=parameters$Alpha/(parameters$Alpha+parameters$Beta),
                         SD=sqrt(with(parameters, Alpha*Beta/((Alpha+Beta)^2*(Alpha+Beta+1)) ) ),
                         Time = 1:N) %>%
  mutate(Group=as.factor(SD)) %>%
  group_by(Group) %>%
  slice(c(1,n())) %>%
  mutate(New_Time = ifelse(row_number()%%2==0, Time+0.4, Time-0.4),
         Config = "Beta Configuration 6")


real_chpt_config <- which(diff(parameters$Alpha)!=0) + 1
real_chpt_config
chpt_found <- sapply(full_fit_fake_series, function(x) { c(cpt.dist(which(x$iid_mdl_mean_model@solution==1)+1, real_chpt_config, N=100), 
                                                           cpt.dist(which(x$ar1_mdl_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$iid_bic_mean_model@solution==1)+1, real_chpt_config, N=100),
                                                           cpt.dist(which(x$ar1_bic_mean_model@solution==1)+1, real_chpt_config, N=100) ) } ) 
                                                       
df <- data.frame(t(chpt_found) )
names(df) <- c("MDL_iid_Mean", "MDL_ar1_Mean", "BIC_iid_Mean", "BIC_ar1_Mean")

df_tall <- df %>% pivot_longer(everything() )

change_results <- bind_rows(change_results,
                            df_tall %>%
                              mutate(Case = "Beta Configuration 6") )


fake_series_df <- bind_rows(fake_series_df1, fake_series_df2,
                            fake_series_df3, fake_series_df4,
                            fake_series_df5, fake_series_df6)
mu_series <- bind_rows(mu_series1, mu_series2,
                       mu_series3, mu_series4,
                       mu_series5, mu_series6)

ggplot(fake_series_df) +
  geom_line(aes(x=Time, y=Series, group=Num), 
            color="gray50", linewidth=0.2, alpha=0.5) +
  geom_line(data=mu_series, aes(x=New_Time, y=Mu, group=Group),
            color="gray5", linewidth=1.0) +
  geom_line(data=mu_series, aes(x=New_Time, y=(Mu-2*SD), group=Group),
            color="gray5", linewidth=0.6, linetype=2) + 
  geom_line(data=mu_series, aes(x=New_Time, y=(Mu+2*SD), group=Group),
            color="gray5", linewidth=0.6, linetype=2) + 
  facet_wrap(~Config, dir="v", nrow=2) + 
  theme_minimal() +
  labs(title="Simulated Beta marginal time series",
       subtitle="Expected value with 2 standard deviation bands provided") +
  theme(axis.title = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        plot.title.position = "plot") 

ggsave("beta_ar1_chptConfigurations.png", 
       width=6.5, height=3, bg="white")

change_results_plot <- change_results %>%
  filter(str_detect(name, "ar1")) %>%
  mutate(`Penalty Structure` = factor(name,
                                      levels=c("MDL_ar1_Mean", "BIC_ar1_Mean"),
                                      labels=c("MDL", "BIC") ) )
ggplot(change_results_plot, aes(y=`Penalty Structure`, x=value)) +
  geom_boxplot( ) + 
  stat_summary(fun="mean", geom="point", shape=23, fill="gray40") + 
  facet_wrap(~Case, dir="v", nrow=2) + 
  labs(x="Distance in estimated vs true changepoint configuration",
       title="Comparison of BIC and MDL penalty structures in change point detection") +
  theme_minimal() + 
  theme(panel.spacing = unit(1.5, "lines"),
        plot.title.position = "plot")

ggsave("beta_ar1_chptFindings.png", 
       width=6.5, height=3, bg="white")
