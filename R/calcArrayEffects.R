calcArrayEffects <-
  	function(rff, basisSplineFunction, snm.obj, model.objects, M.matrix,lf) 
{
#  offset <- 1;
#  rfx <- list();
#  for(i in 1:ncol(snm.obj$int.var)) { 
#    rfx[[i]] <- matrix(rff@ranef[offset:(-1 + offset + nrow(lf$FL$trms[[i]]$Zt))], nr=length(unique(snm.obj$int.var[,i])))
#    rfx[[i]] <- rff[[i]]
#    offset <- offset + nrow(rff[[i]])
#  }

  ars <- sapply(1:ncol(M.matrix), function(i) {
    		mREFs <- sapply(names(rff), function(j) {
      				model.objects$F.mats[[j]][i, ] %*% as.matrix(rff[[j]])
    				})
    		bSM <- predict(basisSplineFunction, M.matrix[, as.numeric(i)])
    		arsL <- bSM %*% mREFs
    		rowSums(arsL)
  		})
  ars
}
