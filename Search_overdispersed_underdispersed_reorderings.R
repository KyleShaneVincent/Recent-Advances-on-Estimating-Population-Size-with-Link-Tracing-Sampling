#The "Search_overdispersed_underdispersed_reorderings" function will conduct a search for an overdispersed sample reordering and an undispersed sample reordering
#to use in the Gelman_Rubin function.
#The first searches for one that is smaller in probability of selection relative to the original reordering.
#The second searches for one that is larger in probability of selection relative to the original reordering.

Search_overdispersed_underdispersed_reorderings <- function(s, search_iterations, use_beta_est, var_est, var_est_new)
	{
	#Use the appropriate values for the beta parameters, i.e., if they are the true or estimated values	
	if(use_beta_est == TRUE)
		{
		beta_strata_calc = beta_estimates_list$Strata_Parameter_Estimates
		if(is.na(beta_sex_worker))
			beta_sex_worker_calc = NA else
			beta_sex_worker_calc = beta_estimates_list$Beta_Sex_Worker_Parameter_Estimate		
		beta_net_size_calc = beta_estimates_list$Log_Net_Size_Parameter_Estimate
		} else
			{
			beta_strata_calc = beta_strata
			if(is.na(beta_sex_worker))
				beta_sex_worker_calc = NA else
				beta_sex_worker_calc = beta_estimates_list$Beta_Sex_Worker_Parameter_Estimate
			beta_net_size_calc = beta_net_size
			}

	#Determine the proportional empirical log probability of obtaining the original sample reordering
	P.prob.list.start = list()
	for(cluster in 1:Graph_s_numclusters)
		{
		P.prob.list.start[[cluster]] = log(1)
		s.cluster = as.numeric(names(which(Graph_s_clusters == cluster)))
		s.cluster.list = list()
		
		
		for(k in 1:(num_waves+1))
			s.cluster.list[[k]] = s.cluster[which(s.cluster %in% s[[k]])]
			
		
		for(k in 2:(num_waves+1))
			{
			candidates.cluster = setdiff(which(matrix(pop$x[s.cluster.list[[k-1]],],nrow=length(s.cluster.list[[k-1]]))==1,arr.ind=TRUE)[,2], 
									unlist(s.cluster.list[1:(k-1)]))
			candidates.cluster_selected = intersect(candidates.cluster,s.cluster.list[[k]])
			candidates.cluster_not_selected = setdiff(candidates.cluster,s.cluster.list[[k]])		
			for(i in candidates.cluster_selected)
				{
				strata_i = pop$strata[[i]]		
				ngbrs_i = s.cluster.list[[k-1]][which(matrix(pop$x[s.cluster.list[[k-1]],i],nrow=length(s.cluster.list[[k-1]]))==1)]
				strata_ngbrs_i = pop$strata[ngbrs_i]
				sex_worker_ngrs_i = pop$sex.worker[ngbrs_i]
				if(length(ngbrs_i) == 1)
					net_size_ngbrs_i = pop$net_size[ngbrs_i] else
						net_size_ngbrs_i = pop$net_size[ngbrs_i]
				P.prob.list.start[[cluster]] = P.prob.list.start[[cluster]] + log(1 - prod(1-inv.logit(beta_strata_calc[strata_ngbrs_i,strata_i] + 
												ifelse(is.na(beta_sex_worker),0,beta_sex_worker_calc*sex_worker_ngrs_i) +
												beta_net_size_calc*log(net_size_ngbrs_i))))
				}
				
			for(i in candidates.cluster_not_selected)
				{
				strata_i = pop$strata[[i]]		
				ngbrs_i = s.cluster.list[[k-1]][which(matrix(pop$x[s.cluster.list[[k-1]],i],nrow=length(s.cluster.list[[k-1]]))==1)]
				strata_ngbrs_i = pop$strata[ngbrs_i]
				sex_worker_ngrs_i = pop$sex.worker[ngbrs_i]
				if(length(ngbrs_i) == 1)
					net_size_ngbrs_i = pop$net_size[ngbrs_i] else
						net_size_ngbrs_i = pop$net_size[ngbrs_i]
				P.prob.list.start[[cluster]] = P.prob.list.start[[cluster]]  + log(prod(1-inv.logit(beta_strata_calc[strata_ngbrs_i,strata_i] + 
										ifelse(is.na(beta_sex_worker),0,beta_sex_worker_calc*sex_worker_ngrs_i) +
										beta_net_size_calc*log(net_size_ngbrs_i))))
				}	
			}
		}

	P.prob.list.overdispersed = P.prob.list.start			
	
	accept_search_overdispersed = accept_search_underdispersed = numeric(search_iterations)	
	s.resample = s	
	for(k in 1:search_iterations)
		{
		Resample_selection(s.resample=s.resample, beta_strata_calc=beta_strata_calc, beta_sex_worker_calc=beta_sex_worker_calc, 
							beta_net_size_calc=beta_net_size_calc)
		if(Q.prob > log(0))
			if(Q.prob < P.prob.list.overdispersed[[cluster_select]])
				{
				P.prob.list.overdispersed[[cluster_select]] = Q.prob
				s.resample = s.resample.proposed
				accept_search_overdispersed[k] = 1
				}
		}
	P.prob.overdispersed <<- sum(unlist(P.prob.list.overdispersed))
	s.overdispersed <<- s.resample
	accept_search_overdispersed <<- accept_search_overdispersed

	P.prob.list.underdispersed = P.prob.list.start
	s.resample = s	
	for(k in 1:search_iterations)
		{
		Resample_selection(s.resample=s.resample, beta_strata_calc=beta_strata_calc, beta_sex_worker_calc=beta_sex_worker_calc, 
							beta_net_size_calc=beta_net_size_calc)
		if(Q.prob > P.prob.list.underdispersed[[cluster_select]])
			{
			P.prob.list.underdispersed[[cluster_select]] = Q.prob
			s.resample = s.resample.proposed
			accept_search_underdispersed[k] = 1
			}
		}
	P.prob.underdispersed <<- sum(unlist(P.prob.list.underdispersed))
	s.underdispersed <<- s.resample
	accept_search_underdispersed <<- accept_search_underdispersed
	
	cat("The log probability of the original reordering is", sum(unlist(P.prob.list.start)), '\n')
	cat('\n')
	
	cat("The log probabilities of the overdispersed and underdispersed reorderings are, respecitvely", P.prob.overdispersed, P.prob.underdispersed, '\n')
	cat('\n')
	
	cat("The total number of \"accepted\" overdispersed and underdispersed reorderings in the search algorithms are, respectively,", 
			sum(accept_search_overdispersed), "and", sum(accept_search_underdispersed), '\n')
	cat('\n')			
			
	cat("The original estimates corresponding with these reorderings (and which will act as MCMC chain seeds) are, respectively,", '\n')
	cat("Overdispersed reordering:", '\n')
	FS_estimate(s0=s.overdispersed$Wave_0, var_est) 
	cat("Estimator based on linkage counts", '\n')
	print(N.hat.strata)
	cat("Estimator based on node counts", '\n')
	print(N.hat.strata.1)
	cat("Underdispersed reordering:", '\n')
	FS_estimate(s0=s.underdispersed$Wave_0, var_est_new) 
	cat("Estimator based on linkage counts", '\n')
	print(N.hat.strata)
	cat("Estimator based on node counts", '\n')
	print(N.hat.strata.1)
	
	cat('\n')
	cat("The new estimates corresponding with these reorderings (and which will act as MCMC chain seeds) are, respectively,", '\n')
	cat("Overdispersed reordering:", '\n')
	FS_new_estimate(s0=s.overdispersed$Wave_0, var_est) 
	cat("Estimator based on linkage counts", '\n')
	print(N.hat.strata.new)
	cat("Estimator based on node counts", '\n')
	print(N.hat.strata.1.new)
	cat("Underdispersed reordering:", '\n')
	FS_new_estimate(s0=s.underdispersed$Wave_0, var_est_new) 
	cat("Estimator based on linkage counts", '\n')
	print(N.hat.strata.new)
	cat("Estimator based on node counts", '\n')
	print(N.hat.strata.1.new)	
	}
	