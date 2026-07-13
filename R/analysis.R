

#phase things to record:
#time to N fully contracted particles
#rate of change of links (i.e., how fast are links being added?)




Network = function(links){
	links = links[links[,1] > 0,] #only existing links. still indexed by [lid]

	network = cbind(
		links[,3], # X
		links[,4], # Y

		links[,7], #parent lid (lid = row index in links[]/network[])
		links[,8], #child 1 lid
		links[,9], #child 2 lid

		links[,6], #generation
		0, # strahler TODO
		0, # leaves
		0  # descendants
	)

	#calculate leaves
	leaf_lids = which(network[,4] + network[,5] == 0) #childless nodes
	counted = rep(1,nrow(network))
	for (lid in leaf_lids){ #for each leaf, follow parents to progenitor
		node = lid
		d = 0
		while(network[node,3] != 0){ #the initial link has a parent of "0"
			#parents are strictly older than children, so this cannot hang
			if (counted[node]){
				d = d + counted[node]
				counted[node] = FALSE #mark node as counted for descendants
			}
			p = network[node,3] #parent lid
			network[p,7] = network[p,7] + 1 #add leaf count
			network[p,8] = network[p,8] + d #add descendants count
			node = p
		}
	}

	network
}


