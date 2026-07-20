
# there seems to be an unexplored phase change in 9b:3rd from bottom right; the model acts dendritically on a larger scale.
# this is quite important to explore and measure similarly to the others
# the solution? visual methods like fractal dim and lacunarity still apply; it's only the network methods that don't work as well

# we need a native fractal dim + boxcounting algorithm :(
# ideally performed on a render, since checking for within tris is hard

#NEED TO FIGURE OUT: how much do we care about the structure itself, vs. the network?? I'd like to only focus on the network, and generally ignore the triangular output...


# --- key ruleset statistics: ---
# tortuosity = average (distance along structure to origin) / (distance to origin)
# 	> how winding is the structure? aka efficiency!
#		> appears to be mode-invariant (mode=dendritic vs. porous)
# [leafiness] = (number of leaves) / (number of nodes)
#		> how branching is the structure?
# density = nodes per unit^2
# 	> how densely packed are the nodes?


# --- porous-only ----
# surface area = number of non-node links / total
# 	> how much of the structure is internal?
#		> surface_area = length(terminal_lids) / length(nodes)
# porosity = number of cycles
# 	> how holey is the structure?
# openness = average length of pore cycles (exclude cycles of less than 5 or so nodes; HISTOGRAM this to determine cutoff)
# 	> how big are the pores?
# [roundness] = average ratio of (shortest diameter) / (longest diameter) for each pore
# 	> how circular are the pores?



# --- other stats
#
# for each of these, HAVE A HISTOGRAM of the distribution (instead of just the average/sd) for default output (+default porous)
# image methods: fractal dimension, lacunarity...
# expansion speed = age of the 1000th link
#		> expansion_speed = links[nodes[length(nodes)],10],
# 	> how fast does the structure grow?
# constriction = (average link length) / link_range
#		> constriction = mean(links[nodes,14])/rules$link_range(origin), #NOT accurate if link_range is a heatmap
# 	> how clumped are the particles?


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

Stats = function(sim,n=5000){
	list2env(sim$state,environment()) #unpack state into local variables to edit
	rules = lapply(sim$rules,parameter)

	links = cbind(links,0,0,0,0,0,0,0,0)

	# 1 endpoint1: [pid] of first particle in link
	# 2 endpoint2: [pid] of second particle
	# 3 avg. x: link midpoint. only ever calculated twice: on generation, and when the link attempts to grow
	# 4 avg. y:
	# 5 active?:
	# 6 unused: is the link in nodes?
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

	# ---- Get Real Network Nodes ----
	nodes = 1:n #[lid] of n oldest links
	#order(links[links[,1] > 0,11])[1:n]  n oldest by generation, but this doesn't make as much sense

	terminal_lids = which(!(links[nodes,8] %in% nodes) | !(links[nodes,9] %in% nodes)) #links which do not have children in nodes.
	nodes = nodes[-terminal_lids] #remove them; they are not real triangles in the network
	leaf_nodes = which(links[nodes,8] %in% terminal_lids & links[nodes,9] %in% terminal_lids) #childless nodes => leaves

	links[nodes,17] = 1 #record as real vertex
	links[leaf_nodes,18] = 1 #record as leaf vertex

	# ---- Leaf-inward calculations ----
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

	# ---- Origin-outward calculations ----
	origin = links[1,3:4]
	nodes_to_check = c(links[1,8:9])
	path_dists = c(0,0)
	best_generation = rep(Inf,nrow(links)) #we want to the shortest distance, ~ path with the fewest links
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
			counted[lid] = 0

			# push children
			valid_children = c(node[8] %in% nodes && node[11] < best_generation[node[8]], node[9] %in% nodes && node[11] < best_generation[node[9]]) #ignore children past the top n
			if (any(valid_children)){
				path_dists = c(path_dists,rep(links[lid,13],2)[valid_children])
				nodes_to_check = c(nodes_to_check,node[8:9][valid_children])
				best_generation[node[8:9][valid_children]] = node[11]
			}
	}

	# ---- Calculate Stats ----
	radius = median(links[nodes,12]) #radius of circle encompassing half the points
	area = pi * radius^2 #area of circle encompassing half the points, in units^2
	density = 0.5*length(nodes) / area #points per unit^2

	#stats about the NETWORK, not the actual structure. network is what we care about
	list(
			tortuosity = mean(links[nodes[-1],13]/links[nodes[-1],12]), #ignore origin link to avoid /0
			leafiness = length(leaf_nodes) / length(nodes),
			density = density

	)

	#area = pi * r^2. r = mean(dist)


	#links[nodes,]
}


#euler characteristic should be around two, with minimal-ish crossing

#use something like
#data=c();for (i in 1:10) data = rbind(data,unlist(Stats(sim,100*i)))
#to show that these are invariant to size over a certain amount of links

debug_draw = function(sim,n=5000){
	sim$state$links	= sim$state$links[1:n,]
	sim
}


#function to extract info from a simset
#function to multithread that




