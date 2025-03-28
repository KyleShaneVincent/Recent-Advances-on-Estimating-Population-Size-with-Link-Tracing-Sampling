#The "Gelman_Rubin" function will apply the Gelman-Rubin test for convergence. MCMC chains commence with "overdispersed" and "underdispersed" seed sample reorderings.

Gelman_Rubin <- function(s.overdispersed, s.underdispersed, num_mcmc, use_beta_est, var_est, var_est_new)
	{
	print("Working on first part (overdispersed)...")
	MCMC(s=s.overdispersed, num_mcmc, num_chains=1, search_iterations=0, use_beta_est=use_beta_est, var_est=var_est, var_est_new=var_est_new)
	GR_N.hat.strata_overdispersed = mcmc(N.hat.strata.mcmc) 
	GR_N.hat.strata.new_overdispersed = mcmc(N.hat.strata.new.mcmc)
	GR_N.hat.strata.1_overdispersed = mcmc(N.hat.strata.1.mcmc) 
	GR_N.hat.strata.1.new_overdispersed = mcmc(N.hat.strata.1.new.mcmc)
	print("Mean acceptance rate for first chain (overdispersed) is:")
	print(mean(accept_mcmc))	

	print("Working on second part (underdispersed)...")
	MCMC(s=s.underdispersed, num_mcmc, num_chains=1, search_iterations=0, use_beta_est=use_beta_est, var_est=var_est, var_est_new=var_est_new)
	GR_N.hat.strata_underdispersed = mcmc(N.hat.strata.mcmc) 
	GR_N.hat.strata.new_underdispersed = mcmc(N.hat.strata.new.mcmc)
	GR_N.hat.strata.1_underdispersed = mcmc(N.hat.strata.1.mcmc) 
	GR_N.hat.strata.1.new_underdispersed = mcmc(N.hat.strata.1.new.mcmc)
	print("Mean acceptance rate for second chain (underdispersed) is:")
	print(mean(accept_mcmc))

	GR_N.hat.strata = list()
	GR_N.hat.strata.new = list()
	GR_N.hat.strata.1 = list()
	GR_N.hat.strata.1.new = list()
	for(k in 1:num_strata)
		{
		GR_N.hat.strata[[k]] = mcmc.list(GR_N.hat.strata_overdispersed[,k], GR_N.hat.strata_underdispersed[,k])
		GR_N.hat.strata.new[[k]] = mcmc.list(GR_N.hat.strata.new_overdispersed[,k], GR_N.hat.strata.new_underdispersed[,k])
		GR_N.hat.strata.1[[k]] = mcmc.list(GR_N.hat.strata.1_overdispersed[,k], GR_N.hat.strata.1_underdispersed[,k])
		GR_N.hat.strata.1.new[[k]] = mcmc.list(GR_N.hat.strata.1.new_overdispersed[,k], GR_N.hat.strata.1.new_underdispersed[,k])
		}
		
	GR_N.hat.strata <<- GR_N.hat.strata
	GR_N.hat.strata.new <<- GR_N.hat.strata.new	
	GR_N.hat.strata.1 <<- GR_N.hat.strata.1
	GR_N.hat.strata.1.new <<- GR_N.hat.strata.1.new	
			
	for(k in 1:num_strata)
		{
		cat('\n')
		cat("Gelman-Rubin diagnostic for original estimators, stratum", k,'\n')
		cat("Estimator based on linkage counts", '\n')
		print(gelman.diag(GR_N.hat.strata[[k]]))
		cat("Estimator based on node counts", '\n')
		print(gelman.diag(GR_N.hat.strata.1[[k]]))
		}
	for(k in 1:num_strata)
		{
		cat('\n')
		cat("Gelman-Rubin diagnostic for new estimators, stratum", k,'\n')
		cat("Estimator based on linkage counts", '\n')
		print(gelman.diag(GR_N.hat.strata.new[[k]]))
		cat("Estimator based on node counts", '\n')
		print(gelman.diag(GR_N.hat.strata.1.new[[k]]))
		}
		
	X11()	
	par(mfrow=c(num_strata,1))	
	for(k in 1:num_strata)
		{	
		plot(as.vector(GR_N.hat.strata[[k]][[1]]), typ='l', main=paste("MCMC Chain Trace Plots for Strata", k, "Population Size Estimator Based on Linkage Counts"),
		xlab="MCMC Iteration (overdispersed seeded reording in blue and underdispersed seeded reordering in red)", 
		ylab="Strata Pop Size", col="blue", ylim=range(c(GR_N.hat.strata[[k]][[1]], GR_N.hat.strata[[k]][[2]])))
		points(as.vector(GR_N.hat.strata[[k]][[2]]), typ='l', col="red")
		}		
		
	X11()	
	par(mfrow=c(num_strata,1))	
	for(k in 1:num_strata)
		{	
		plot(as.vector(GR_N.hat.strata.1[[k]][[1]]), typ='l', main=paste("MCMC Chain Trace Plots for Strata", k, "Population Size Estimator Based on Node Counts"),
		xlab="MCMC Iteration (overdispersed seeded reording in blue and underdispersed seeded reordering in red)", 
		ylab="Strata Pop Size", col="blue", ylim=range(c(GR_N.hat.strata.1[[k]][[1]], GR_N.hat.strata.1[[k]][[2]])))
		points(as.vector(GR_N.hat.strata.1[[k]][[2]]), typ='l', col="red")
		}	
		
	X11()	
	par(mfrow=c(num_strata,1))	
	for(k in 1:num_strata)
		{	
		plot(as.vector(GR_N.hat.strata.new[[k]][[1]]), typ='l', main=paste("MCMC Chain Trace Plots for Strata", k, "New Population Size Estimator Based on Linkage Counts"),
		xlab="MCMC Iteration (overdispersed seeded reording in blue and underdispersed seeded reordering in red)", 
		ylab="Strata Pop Size", col="blue", ylim=range(c(GR_N.hat.strata.new[[k]][[1]], GR_N.hat.strata.new[[k]][[2]])))
		points(as.vector(GR_N.hat.strata.new[[k]][[2]]), typ='l', col="red")
		}			

	X11()	
	par(mfrow=c(num_strata,1))	
	for(k in 1:num_strata)
		{	
		plot(as.vector(GR_N.hat.strata.1.new[[k]][[1]]), typ='l', main=paste("MCMC Chain Trace Plots for Strata", k, "New Population Size Estimator Based on Node Counts"),
		xlab="MCMC Iteration (overdispersed seeded reording in blue and underdispersed seeded reordering in red)", 
		ylab="Strata Pop Size", col="blue", ylim=range(c(GR_N.hat.strata.1.new[[k]][[1]], GR_N.hat.strata.1.new[[k]][[2]])))
		points(as.vector(GR_N.hat.strata.1.new[[k]][[2]]), typ='l', col="red")
		}	

	X11()
	par(mfrow=c(num_strata,1))
	for(k in 1:num_strata)
		{
		hg1 = hist(GR_N.hat.strata[[k]][[1]], breaks=200, plot=FALSE); c1=rgb(0,0,1,1/4)
		hg2 = hist(GR_N.hat.strata[[k]][[2]], breaks=200, plot=FALSE); c2=rgb(1,0,0,1/4)
		xlim=range(unlist(GR_N.hat.strata[[k]]))
		ylim=range(0,max(c(hg1$counts,hg2$counts)))
		plot(hg1, col=c1, main=paste("Histogram of MCMC Strata", k, "Original Estimates Based on Linkage Counts with First Chain in Blue and Second Chain in Red", sep=" "), 
		xlab="Estimate with MCMC Seeds in Circles and Improved Estimate as Triangles",xlim=xlim,ylim=ylim)
		abline(v=mean(GR_N.hat.strata[[k]][[1]], col="blue", lwd=2))
		plot(hg2, col=c2, add=TRUE)
		abline(v=mean(GR_N.hat.strata[[k]][[2]], col="blue"))
		points(c(GR_N.hat.strata[[k]][[1]][1],GR_N.hat.strata[[k]][[2]][1]),c(0,0), pch=21, cex=2.5, bg=c("blue", "red"))	
		points(c(mean(GR_N.hat.strata[[k]][[1]]),mean(GR_N.hat.strata[[k]][[2]])[1]),c(0,0), pch=24, cex=2.5, bg=c("blue", "red"))		
		}
		
	X11()
	par(mfrow=c(num_strata,1))
	for(k in 1:num_strata)
		{
		hg1 = hist(GR_N.hat.strata.1[[k]][[1]], breaks=200, plot=FALSE); c1=rgb(0,0,1,1/4)
		hg2 = hist(GR_N.hat.strata.1[[k]][[2]], breaks=200, plot=FALSE); c2=rgb(1,0,0,1/4)
		xlim=range(unlist(GR_N.hat.strata.1[[k]]))
		ylim=range(0,max(c(hg1$counts,hg2$counts)))
		plot(hg1, col=c1, main=paste("Histogram of MCMC Strata", k, "Original Estimates Based on Node Counts with First Chain in Blue and Second Chain in Red", sep=" "), 
		xlab="Estimate with MCMC Seeds in Circles and Improved Estimate as Triangles",xlim=xlim,ylim=ylim)
		abline(v=mean(GR_N.hat.strata.1[[k]][[1]], col="blue", lwd=2))
		plot(hg2, col=c2, add=TRUE)
		abline(v=mean(GR_N.hat.strata.1[[k]][[2]], col="blue"))
		points(c(GR_N.hat.strata.1[[k]][[1]][1],GR_N.hat.strata.1[[k]][[2]][1]),c(0,0), pch=21, cex=2.5, bg=c("blue", "red"))	
		points(c(mean(GR_N.hat.strata.1[[k]][[1]]),mean(GR_N.hat.strata.1[[k]][[2]])[1]),c(0,0), pch=24, cex=2.5, bg=c("blue", "red"))		
		}			
		
	X11()
	par(mfrow=c(num_strata,1))
	for(k in 1:num_strata)
		{
		hg1 = hist(GR_N.hat.strata.new[[k]][[1]], breaks=200, plot=FALSE); c1=rgb(0,0,1,1/4)
		hg2 = hist(GR_N.hat.strata.new[[k]][[2]], breaks=200, plot=FALSE); c2=rgb(1,0,0,1/4)
		xlim=range(unlist(GR_N.hat.strata.new[[k]]))
		ylim=range(0,max(c(hg1$counts,hg2$counts)))
		plot(hg1, col=c1, main=paste("Histogram of MCMC Strata", k, "New Estimates Based on Linkage Counts with First Chain in Blue and Second Chain in Red", sep=" "), 
		xlab="Estimate with MCMC Seeds in Circles and Improved Estimate as Triangles",xlim=xlim,ylim=ylim)
		plot(hg2, col=c2, add=TRUE)
		points(c(GR_N.hat.strata.new[[k]][[1]][1],GR_N.hat.strata.new[[k]][[2]][1]),c(0,0), pch=21, cex=2.5, bg=c("blue", "red"))	
		points(c(mean(GR_N.hat.strata.new[[k]][[1]]),mean(GR_N.hat.strata.new[[k]][[2]])[1]),c(0,0), pch=24, cex=2.5, bg=c("blue", "red"))	
		}	
	
	X11()
	par(mfrow=c(num_strata,1))
	for(k in 1:num_strata)
		{
		hg1 = hist(GR_N.hat.strata.1.new[[k]][[1]], breaks=200, plot=FALSE); c1=rgb(0,0,1,1/4)
		hg2 = hist(GR_N.hat.strata.1.new[[k]][[2]], breaks=200, plot=FALSE); c2=rgb(1,0,0,1/4)
		xlim=range(unlist(GR_N.hat.strata.1.new[[k]]))
		ylim=range(0,max(c(hg1$counts,hg2$counts)))
		plot(hg1, col=c1, main=paste("Histogram of MCMC Strata", k, "New Estimates Based on Node Counts with First Chain in Blue and Second Chain in Red", sep=" "), 
		xlab="Estimate with MCMC Seeds in Circles and Improved Estimate as Triangles",,xlim=xlim,ylim=ylim)
		plot(hg2, col=c2, add=TRUE)
		points(c(GR_N.hat.strata.1.new[[k]][[1]][1],GR_N.hat.strata.1.new[[k]][[2]][1]),c(0,0), pch=21, cex=2.5, bg=c("blue", "red"))	
		points(c(mean(GR_N.hat.strata.1.new[[k]][[1]]),mean(GR_N.hat.strata.1.new[[k]][[2]])[1]),c(0,0), pch=24, cex=2.5, bg=c("blue", "red"))	
		}			
	}
	