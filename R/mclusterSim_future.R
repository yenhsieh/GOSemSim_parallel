mclusterSim_future = function(named_gslist, semData, measure = "Rel", drop = "IEA", combine = "BMA", n_job = 5){
	require(furrr)
	#options(future.globals.maxSize= 891289600)
	plan(multisession, workers = n_job)
	message('making gene2GO...')
	cluster_gos = future_map(1:length(named_gslist), function(i){
		map(named_gslist[[i]], function(g) gene2GO(g, semData, dropCodes = drop))
	})
	uniqueGO <- unique(unlist(cluster_gos))
	message(paste0('unique GO: ', length(uniqueGO) ))
	message('creating GO matrix...')
	go_matrix <- mgoSim(uniqueGO, uniqueGO, semData, measure = measure, combine = NULL)
	message('GO similarity scoring...')
	pairwise_df = 
		combn(1:length(named_gslist), m =2) %>%
		t() %>%
		as_tibble() %>%
		mutate(
			gs1 = names(named_gslist)[V1],
			gs2 = names(named_gslist)[V2],
			gos1 = future_map(V1, function(i){
				gs <- unlist(cluster_gos[[i]])
				gs <- gs[!is.na(gs)]
				return(gs)
			}),
			gos2 = future_map(V2, function(j){
				gs <- unlist(cluster_gos[[j]])
				gs <- gs[!is.na(gs)]
				return(gs)
			}),
			score = future_map2_dbl(gos1, gos2, function(go1, go2){
				combineScores(go_matrix[go1, go2], combine = combine)
			})
		) %>%
		dplyr::select(-V1, -V2)
    
	if (!inherits(plan(), "sequential")) plan(sequential)	
	rm(cluster_gos, go_matrix, uniqueGO)
	return(pairwise_df)
}
