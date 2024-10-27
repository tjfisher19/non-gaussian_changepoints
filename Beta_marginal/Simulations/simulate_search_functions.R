

####################################################
## We build a suggestion matrix with some initial
##   change point configurations. By default
##   we should be leaning towards no segmentations
##   so the suggestion matrix will include the 
##   following suggestions (rows)
##   - No change point - row of all zeroes
##   - Single changepoint - row with a single 1
##   - Two changepoints - rows with 2 change points
## For the case of no & one change point we look
##   at all viable configurations. For the case of
##   two changepoints we randomly select two
##   locations but in such a way that both
##   guarantee the minimum regime

suggest_matrix <- matrix(0,
                         ncol= N-1,
                         nrow= PopSize)
for(i in 2:(NCOL(suggest_matrix)-2*min.regime+3) ) {
  suggest_matrix[i, min.regime + i-2] <- 1
}
for(i in (NCOL(suggest_matrix)-min.regime+2):PopSize) {
  suggest_matrix[i, sample(min.regime:((NCOL(suggest_matrix)-min.regime+1)/2), size=1)] <- 1
  suggest_matrix[i, sample(((NCOL(suggest_matrix)-min.regime+1)/2+1):(NCOL(suggest_matrix)-min.regime+1), size=1)] <- 1
}



fit_fake_series <- function(x) {
  cnt <<- cnt + 1
  cat(cnt)
  ga.search.mdl.iid.mean <- ga(type="binary",
                               fitness=ga_search_function_beta_chpt_mdl,
                               nBits=N-1,
                               series= x,
                               model_type="mean",
                               quick_eval=TRUE,
                               min.regime=min.regime,
                               max.chpts=floor(N/min.regime),
                               popSize=PopSize,
                               suggestions = suggest_matrix,
                               parallel = cores,
                               monitor = FALSE,
                               maxiter=maxit,
                               run=max.run)
  cat(".")
  ga.search.mdl.ar1.mean <- ga(type="binary", 
                               fitness=ga_search_function_beta_ar1_chpt_mdl,
                               nBits=N-1,
                               series= x,
                               model_type="mean",
                               approx=TRUE,
                               min.regime=min.regime,
                               max.chpts=floor(N/min.regime),
                               popSize=PopSize,
                               suggestions = suggest_matrix,
                               parallel = cores,
                               monitor = FALSE,
                               maxiter=maxit,
                               run=max.run)
  cat(".")
  ga.search.bic.iid.mean <- ga(type="binary", 
                               fitness=ga_search_function_beta_chpt_bic,
                               nBits=N-1,
                               series= x,
                               model_type="mean",
                               quick_eval=TRUE,
                               min.regime=min.regime,
                               max.chpts=floor(N/min.regime),
                               popSize=PopSize,
                               suggestions = suggest_matrix,
                               parallel = cores,
                               monitor = FALSE,
                               maxiter=maxit,
                               run=max.run)
  cat(".")
  ga.search.bic.ar1.mean <- ga(type="binary", 
                               fitness=ga_search_function_beta_ar1_chpt_bic,
                               nBits=N-1,
                               series= x,
                               model_type="mean",
                               approx=TRUE,
                               min.regime=min.regime,
                               max.chpts=floor(N/min.regime),
                               popSize=PopSize,
                               suggestions = suggest_matrix,
                               parallel = cores,
                               monitor = FALSE,
                               maxiter=maxit,
                               run=max.run)
  # cat(".")
  # ga.search.mdl.iid.shape <- ga(type="binary",
  #                               fitness=ga_search_function_beta_chpt_mdl,
  #                               nBits=N-1,
  #                               series= x,
  #                               model_type="shape",
  #                               quick_eval=TRUE,
  #                               min.regime=min.regime,
  #                               max.chpts=floor(N/min.regime),
  #                               popSize=PopSize,
  #                               suggestions = suggest_matrix,
  #                               parallel = cores,
  #                               monitor = FALSE,
  #                               maxiter=maxit,
  #                               run=max.run)
  # cat(".")
  # ga.search.mdl.ar1.shape <- ga(type="binary", 
  #                               fitness=ga_search_function_beta_ar1_chpt_mdl,
  #                               nBits=N-1,
  #                               series= x,
  #                               model_type="shape",
  #                               approx=TRUE,
  #                               min.regime=min.regime,
  #                               max.chpts=floor(N/min.regime),
  #                               popSize=PopSize,
  #                               suggestions = suggest_matrix,
  #                               parallel = cores,
  #                               monitor = FALSE,
  #                               maxiter=maxit,
  #                               run=max.run)
  # cat(".")
  # ga.search.bic.iid.shape <- ga(type="binary", 
  #                               fitness=ga_search_function_beta_chpt_bic,
  #                               nBits=N-1,
  #                               series= x,
  #                               model_type="shape",
  #                               quick_eval=TRUE,
  #                               min.regime=min.regime,
  #                               max.chpts=floor(N/min.regime),
  #                               popSize=PopSize,
  #                               suggestions = suggest_matrix,
  #                               parallel = cores,
  #                               monitor = FALSE,
  #                               maxiter=maxit,
  #                               run=max.run)
  # cat(".")
  # ga.search.bic.ar1.shape <- ga(type="binary", 
  #                               fitness=ga_search_function_beta_ar1_chpt_bic,
  #                               nBits=N-1,
  #                               series= x,
  #                               model_type="shape",
  #                               approx=TRUE,
  #                               min.regime=min.regime,
  #                               max.chpts=floor(N/min.regime),
  #                               popSize=PopSize,
  #                               suggestions = suggest_matrix,
  #                               parallel = cores,
  #                               monitor = FALSE,
  #                               maxiter=maxit,
  #                               run=max.run)
  cat(".\n")
  list(iid_mdl_mean_model = ga.search.mdl.iid.mean,
       ar1_mdl_mean_model = ga.search.mdl.ar1.mean,
       iid_bic_mean_model = ga.search.bic.iid.mean,
       ar1_bic_mean_model = ga.search.bic.ar1.mean)
       # iid_mdl_shape_model = ga.search.mdl.iid.shape,
       # ar1_mdl_shape_model = ga.search.mdl.ar1.shape,
       # iid_bic_shape_model = ga.search.bic.iid.shape,
       # ar1_bic_shape_model = ga.search.bic.ar1.shape)
}
