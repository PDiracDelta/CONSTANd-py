# see http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html
library(limma)

moderated_ttest <- function(dat, design) {  # dat = datp[, c(tr,ct)]
  # This version makes MAXIMAL use of limma. It is an extension of eb.fit, but 
  # with terminology from moderated_ttest_extracted, which returns 
  # output for more than 1 sample comparison.
  ngenes <- dim(dat)[1]
  fit <- eBayes(lmFit(dat, design))
  logFC <- fit$coefficients # estimate of the log2-fold-change corresponding to the effect size
  df.r <- fit$df.residual # residual degrees of freedom assiciated with ordinary t-statistic and p-value
  df.0 <- rep(fit$df.prior, ngenes) # degrees of freedom associated with s2.0
  s2.0 <- rep(fit$s2.prior, ngenes) # estimated prior value for the variance
  s2 <- (fit$sigma)^2 # sample variance
  s2.post <- fit$s2.post # posterior value for the variance
  t.ord <- fit$coefficients / fit$stdev.unscaled / fit$sigma # vector of ordinary t-statistic: using sigma vs s2.post
  t.mod <- fit$t # moderated t-statistic
  p.ord <- 2*pt(-abs(t.ord), fit$df.residual) # ordinary p-value corresonding to the ordinary t-statistic
  p.mod <- fit$p.value # moderated p-value corresonding to the moderated t-statistic
  if(ngenes>1) q.ord <- apply(X = p.ord, MARGIN = 2, FUN = p.adjust, method='BH') # ordinary q-value corresponding to the ordinary t-statistic
  # if(ngenes>1) q.ord <- qvalue(p.ord)$q#, pi0=1)$q # avoid qvalue library when using BH correction
  else q.ord <- p.ord
  if(ngenes>1) q.mod <- apply(X = p.mod, MARGIN = 2, FUN = p.adjust, method='BH') # moderated q-value corresponding to the moderated t-statistic
  # if(ngenes>1) q.mod <- qvalue(p.mod)$q#, pi0=1)$q # avoid qvalue library when using BH correction
  else q.mod <- p.mod
  # incorporate entity type into colnames to overwrite identical factor names
  colnames(logFC) <- paste0('logFC_', colnames(logFC))
  colnames(t.ord) <- paste0('t.ord', '_', colnames(t.ord))
  colnames(t.mod) <- paste0('t.mod', '_', colnames(t.mod))
  colnames(p.ord) <- paste0('p.ord', '_', colnames(p.ord))
  colnames(p.mod) <- paste0('p.mod', '_', colnames(p.mod))
  colnames(q.ord) <- paste0('q.ord', '_', colnames(q.ord))
  colnames(q.mod) <- paste0('q.mod', '_', colnames(q.mod))
  results <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  # results <- results[order(results$p.mod), ]  # ordering doesn't work when you have mutiple groups and thus columns
  return(results)
}