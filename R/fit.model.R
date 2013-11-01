fit.model <-
	function(obs.fit, snm.obj, basisSplineFunction, lmer.max.iter=500)
{
  snm.obj$M <- snm.obj$dat - obs.fit$res1
  snm.obj$M[snm.obj$nulls,] <- obs.fit$fit0[snm.obj$nulls,]
	
# Split the data into nbins bins based on their mean intensities
  bins <- getSpanningSet(snm.obj)
  
# Build the matrix of weighted raw data and matrix of weighted fitted values for each bin.
  lnp <- length(bins)
  np <- 1:lnp
  Y.pooled <- 0*snm.obj$dat[np,]
  M.pooled <- 0*snm.obj$M[np,]
  for(i in 1:lnp) {
    Y.pooled[i,] = apply(matrix(snm.obj$r.dat[as.vector(bins[[i]]),], ncol=ncol(snm.obj$dat)),2,
        weighted.mean, w=snm.obj$weights[as.vector(bins[[i]])])
    M.pooled[i,] = apply(matrix(snm.obj$M[as.vector(bins[[i]]),], ncol=ncol(snm.obj$M)),2,
        weighted.mean, w=snm.obj$weights[as.vector(bins[[i]])])
  }
  
# Build the basis spline matrix for the pooled coefficients.
  bSM.model <- buildBasisSplineMatrix(M.pooled, basisSplineFunction)
  exp <- new.env()
# Build the data object and fit the mixed effects model
  expObj <- makeDataObject(Y.pooled, np, snm.obj, exp,bins)
  expObj$sp <- as.matrix(bSM.model)
  model.objects <- make.ref.model.matrices(snm.obj, exp)

  lmerCtrl <- lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun = lmer.max.iter))
  yourWarn <- getOption("warn")
  options(warn=-1)
  lf <- do.call("lmer", list(formula=model.objects$ZF, data=expObj, REML=TRUE, control = lmerCtrl, 
                             verbose=FALSE, subset=1:nrow(expObj),weights=expObj$weights))
  options(warn=yourWarn)
  
#  for(i in 1:ncol(snm.obj$int.var)) { 
#    lf$FL$trms[[i]]$ST <- matrix(0,nr=1,nc=1)
#    rownames(lf$FL$trms[[i]]$ST) <- colnames(lf$FL$trms[[i]]$ST) <- paste("spline",i,sep="")
#  }
#  rff <- do.call(lme4:::lmer_finalize,lf)

  rff <- ranef(lf)
	
# Add useful variables to snm.obj
  snm.obj$E.pooled <- matrix(residuals(lf), nrow=dim(Y.pooled)[1])
  snm.obj$Y.pooled <- Y.pooled
  snm.obj$M.pooled <- M.pooled
  snm.obj$rff <- rff
  snm.obj$bin.densities <- sapply(bins,length)
  
  snm.obj$array.fx <- calcArrayEffects(rff, basisSplineFunction, snm.obj, model.objects, snm.obj$M, lf)
  return(snm.obj)
}
