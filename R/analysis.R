

#phase things to record:
#time to N fully contracted particles
#rate of change of links (i.e., how fast are links being added?)
#quantifying the banding: error term on `distance from initial link ~ generation + error`. "How radially/consistently does it expand?"


# --- key ruleset statistics: ---
# --- all calculated when the model hits 1k links
# tortuosity ~ (distance along structure to origin) / (distance to origin)
# 	> how winding is the structure? aka efficiency!
# expansion speed = ticks to get to 1k links
# 	> how fast does the structure grow?
# density = average distance from origin when it hits 1k links
# 	> how densely packed are the links?
# constriction = average distance between a particle and its nearest neighbor (or 2).
# 	> how clumped are the particles?
# [leafiness] = (number of linked vertices) / (number of links (1k)). >1 for dendrites, <1 for porous
#		> how branching is the structure?
# surface area = number of non-node links
# 	> how much of the structure is internal?

# --- pore-related (put these on a 2d heatmap)
# porosity = number of cycles at 1k links
# 	> how holey is the structure?
# openness = average length of pore cycles (exclude cycles of less than 5 or so nodes; HISTOGRAM this to determine cutoff)
# 	> how big are the pores?
# [roundness] = average ratio of (shortest diameter) / (longest diameter) for each pore
# 	> how circular are the pores?



# --- other stats
# for each of these, HAVE A HISTOGRAM of the distribution (instead of just the average/sd) for default output (+default porous)
# image methods: fractal dimension, lacunarity,


#questions:
# tortuosity I believe does work properly with center-out calculation, but we should have a catch to PICK THE SHORTEST OPTION
# how do we get around fast growth, slow contraction making it hard to know when to collect data


Network = function(links){
	links = links[links[,1] > 0,] #only existing links. still indexed by [lid]

	stats = list()

	network = cbind(
		links[,3], # X
		links[,4], # Y

		links[,7], #parent lid (lid = row index in links[]/network[])
		links[,8], #child 1 lid
		links[,9], #child 2 lid

		links[,6], #generation
		0, # leaves, i.e. order
		0 # descendants
	)

	terminal_lids = which(network[,4] + network[,5] == 0)
	network[terminal_lids,] = NA #remove links with no children; they are not new triangles in the network
	#it's important to keep these rows so that network still follows [lid] indexing
	leaf_lids = which(network[,4] %in% terminal_lids & network[,5] %in% terminal_lids) #childless nodes => leaves

	#calculate leaves
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


	#now center-out search to calculate distance along trunk to each node


	network
}






