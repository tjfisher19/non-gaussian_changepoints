

## Some function calls to build datasets for our simulations
source("mv_pois_ar1_generation.R")

set.seed(123)

##############################################################
## No Change point data

## Positve AR cross correlation, Positive at lag 0
X_series <- gen_mv_pois_series_no_change(n=50, n_series=100,
                                         ar_cor="pos", sigma_cor="pos",
                                         lambda=c(13.6, 17.8) )

save(X_series, file="data/nochange_pos_pos.RData")

## Positve AR cross correlation, Negative at lag 0
X_series <- gen_mv_pois_series_no_change(n=50, n_series=100,
                                         ar_cor="pos", sigma_cor="neg",
                                         lambda=c(13.6, 17.8) )

save(X_series, file="data/nochange_pos_neg.RData")

## Negative AR cross correlation, Positive at lag 0
X_series <- gen_mv_pois_series_no_change(n=50, n_series=100,
                                         ar_cor="neg", sigma_cor="pos",
                                         lambda=c(13.6, 17.8) )

save(X_series, file="data/nochange_neg_pos.RData")

## Negative AR cross correlation, Negative at lag 0
X_series <- gen_mv_pois_series_no_change(n=50, n_series=100,
                                         ar_cor="neg", sigma_cor="neg",
                                         lambda=c(13.6, 17.8) )

save(X_series, file="data/nochange_neg_neg.RData")


##############################################################
## Change point data

## Positve AR cross correlation, Positive at lag 0
set.seed(321)
X_series <- gen_mv_pois_series_changes(n=50, n_series=100,
                                       ar_cor="pos", sigma_cor="pos",
                                       lambda=rbind(c(10, 19), c(16, 17)),
                                       chpts=21)

save(X_series, file="data/change_pos_pos.RData")


## Positve AR cross correlation, Negative at lag 0
set.seed(321)
X_series <- gen_mv_pois_series_changes(n=50, n_series=100,
                                       ar_cor="pos", sigma_cor="neg",
                                       lambda=rbind(c(10, 19), c(16, 17)),
                                       chpts=21)

save(X_series, file="data/change_pos_neg.RData")

## Negative AR cross correlation, Negative at lag 0
set.seed(321)
X_series <- gen_mv_pois_series_changes(n=50, n_series=100,
                                       ar_cor="neg", sigma_cor="neg",
                                       lambda=rbind(c(10, 19), c(16, 17)),
                                       chpts=21)

save(X_series, file="data/change_neg_neg.RData")


## Negative AR cross correlation, Positive at lag 0
set.seed(321)
X_series <- gen_mv_pois_series_changes(n=50, n_series=100,
                                       ar_cor="neg", sigma_cor="pos",
                                       lambda=rbind(c(10, 19), c(16, 17)),
                                       chpts=21)

save(X_series, file="data/change_neg_pos.RData")



##############################################################
## Multiple Change point


set.seed(321)
X_series <- gen_mv_pois_series_changes(n=50, n_series=100,
                                       ar_cor="pos", sigma_cor="neg",
                                       lambda=rbind(c(10, 19), c(16, 15), c(20, 19)),
                                       chpts=c(16, 26) )

save(X_series, file="data/mult_change_pos_neg.RData")


set.seed(321)
X_series <- gen_mv_pois_series_changes(n=50, n_series=100,
                                       ar_cor="pos", sigma_cor="pos",
                                       lambda=rbind(c(10, 19), c(16, 15), c(20, 19), c(22, 16)),
                                       chpts=c(16, 26, 41 ) )

save(X_series, file="data/mult_change_pos_neg2.RData")
