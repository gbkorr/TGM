

#phase things to record:
#time to N fully contracted particles
#rate of change of links (i.e., how fast are links being added?)
#quantifying the banding: error term on `distance from initial link ~ generation + error`. "How radially/consistently does it expand?"


# --- key ruleset statistics: ---
# --- all calculated when the model hits 1k links
# tortuosity ~ (distance along structure to origin) / (distance to origin)
# 	> how winding is the structure? aka efficiency!
# expansion speed = age of the 1000th link
# 	> how fast does the structure grow?
# density = average distance from origin
# 	> how densely packed are the links?
# constriction = average link length
# 	> how clumped are the particles?
# [leafiness] = (number of linked vertices) / (number of links). >1 for dendrites, <1 for porous
#		> how branching is the structure?
# surface area = number of non-node links
# 	> how much of the structure is internal?

# --- pore-related (put these on a 2d heatmap)
# porosity = number of cycles
# 	> how holey is the structure?
# openness = average length of pore cycles (exclude cycles of less than 5 or so nodes; HISTOGRAM this to determine cutoff)
# 	> how big are the pores?
# [roundness] = average ratio of (shortest diameter) / (longest diameter) for each pore
# 	> how circular are the pores?



# --- other stats
#
# for each of these, HAVE A HISTOGRAM of the distribution (instead of just the average/sd) for default output (+default porous)
# image methods: fractal dimension, lacunarity,


#questions:
# tortuosity I believe does work properly with center-out calculation, but we should have a catch to PICK THE SHORTEST OPTION
# how do we get around fast growth, slow contraction making it hard to know when to collect data?
# 	> anything based on first 1k is going to basically ignore contraction, making the whole thing kinda dumb
# 	> solution: wait quite a bit longer, but only analyze the n (1k) lowest-generation links
# 	> this makes collection work at any time
# some of these methods use links, some use the network... it might be a bit of a mess to handle both nicely.


# maybe replace the current network with just a copy links, but adding some columns and a "is_node" bool (replace is_active!)
# ^^^^^ ^^^^^^ ^^^^^ ^^^^^^

Network = function(links){
	links = links[links[,1] > 0,] #only existing links. still indexed by [lid]

	stats = list()

	network = cbind(
		links[,3], # X
		links[,4], # Y

		links[,7], #parent lid (lid = row index in links[]/network[])
		links[,8], #child 1 lid
		links[,9], #child 2 lid

		links[,6], #unused
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

Stats = function(sim,n=1000){
	list2env(sim$state,environment()) #unpack state into local variables to edit
	rules = lapply(sim$rules,parameter)

	nodes = order(links[links[,1] > 0,11])[1:n] #[lid] of n oldest links
	links[nodes,6] = 1 #these get a flag in links[,6]

	links = cbind(links,0,0,0,0,0,0,0,0)

	# 1 endpoint1: [pid] of first particle in link
	# 2 endpoint2: [pid] of second particle
	# 3 avg. x: link midpoint. only ever calculated twice: on generation, and when the link attempts to grow
	# 4 avg. y:
	# 5 active?:
	# 6 top n?: is the link in nodes?
	# 7 parent: [lid] of the link that spawned this
	# 8 child1: [lid] of the first link spawned by this
	# 9 child2: [lid] of the second link spawned by this
	# 10 age: tick of creation
	# 11 generation: parent generation + 1

	# 12 distance: from origin
	# 13 path: from origin along structure (calc by spreading out from first)
	# 14 constriction: link length
	# 15 leaves
	# 16 descendants
	# 17 node?: does the node have any children?
	# 18 leaf?: does the node have


	# origin-first calculations
	origin = links[1,3:4]
	nodes_to_check = c(links[1,8:9])
	path_dists = c(0,0)
	while (length(nodes_to_check)) {
			lid = nodes_to_check[1]
			node = links[lid,]

			# calculate stats
			dist_to_origin = mag(node[3:4] - origin)
			dist_to_parent = mag(node[3:4] - links[node[7],3:4])
			link_length = mag(particles[node[1],1:2] - particles[node[2],1:2])

			# update links with stat
			links[lid,12] = dist_to_origin #distance to origin
			links[lid,13] = path_dists[1] + dist_to_parent #distance to origin along path
			links[lid,14] = link_length

			# pop node
			path_dists = path_dists[-1]
			nodes_to_check = nodes_to_check[-1]

			# push children
			valid_children = c(node[8] %in% nodes, node[9] %in% nodes)
			if (any(valid_children)){
				path_dists = c(path_dists,rep(links[lid,13],2)[valid_children])
				nodes_to_check = c(nodes_to_check,node[8:9][valid_children])
			}
	}

	# leaf-inward calculations
	terminal_lids = which(!(links[nodes,8] %in% nodes) | !(links[nodes,9] %in% nodes)) #links which do not have children in nodes.
	real_nodes = nodes[-terminal_lids] #remove them; they are not real triangles in the network
	leaf_nodes = which(links[nodes,8] %in% terminal_lids & links[nodes,9] %in% terminal_lids) #childless nodes => leaves

	#calculate descendants
	counted = rep(1,nrow(links))
	for (lid in leaf_nodes){ #for each leaf, follow parents to progenitor
		node = lid
		d = 0 #descendants
		while(links[node,7] != 0){ #the initial link has a parent of "0"
			#parents are strictly older than children, so this cannot hang (will always trace back to origin)
			if (counted[node]){
				d = d + counted[node]
				counted[node] = 0 #mark node as counted for descendants
			}
			node = links[node,7] #parent lid
			links[node,15] = links[node,15] + 1 #add leaf count to parent
			links[node,16] = links[node,16] + d #add descendants count
		}
	}

	#return total stats, not the raw
	links[nodes,]
}




