

# ---- Parameter Initialization ----
# rules prototype
Rules = function(
	link_range = 0.18,
	mobility = 0.04,
	contraction = 1.0,
	cohesion = 0.0,
	branching = 1.0,

	seed = 1,
	growth_rate = 0.05 #usually constant
) list(as.list(environment()))[[1]]

#evaluates parameter as either constant or a function of position
parameter = function(rule) ifelse(is.function(rule),rule,function(pos) rule(pos))

# ---- Particle Initialization ----
#particle initialization prototype
P_rules = function(
	size = 4, #length and width of bounding region
	density = 200, #particles/unit. TODO: add ability for heatmap
	particle_seed = 1, #seed for particle distribution
	chunksize = 0.25 #initial chunksize; if smaller than link_range, particles are re-chunked
) list(as.list(environment()))[[1]]

#initialize particles
Init = function(p_rules=P_rules()){
	n_particles = p_rules$size^2 * p_rules$density

#particles: MATRIX indexed by [pid,], i.e. one row per particle
# x position
# y position
# status: -1 = unlinked, 0 = frozen, >0 = mobile
#--- looped through via chunking

	set.seed(p_rules$particle_seed)
	particles = matrix(c(runif(2*n_particles,0,p_rules$size),rep(1,n_particles)),ncol=3)

#chunks: MATRIX indexed by [chunkX,chunkY]
# entry: LIST of [pid] of particles in that chunk
#--- particles move sufficiently little such that we don't have to update this dynamically

	chunkrows = ceiling(p_rules$size/p_rules$chunksize)
	chunks = matrix(list(),chunkrows + 2,chunkrows + 2) #empty buffer chunks in the +x +y direction to prevent subscript OOB
	for (pid in 1:n_particles){
		chunk_pos = ceiling(particles[pid,1:2]/p_rules$chunksize) #get which chunk the particle is in
		chunks[chunk_pos[1],chunk_pos[2]][[1]] = c(chunks[chunk_pos[1],chunk_pos[2]][[1]], pid) #add its pid to that chunk
	}

#particles_neighbors: LIST indexed by [pid], i.e. one row per particle
# entry: list of particles [pid] connected to this particle
#--- never looped through, only accessed by [pid]

	particle_neighbors = rep(list(NULL),n_particles)

#particle_links: LIST indexed by [pid], i.e. one row per particle
# entry: list of links [lid] connected to this particle
#--- this is ONLY used to check and avoid creating duplicate links when cohesion > 0
#--- never looped through, only accessed by [pid]

	particle_links = rep(list(NULL),n_particles)

#links: MATRIX indexed by [lid,], i.e. one row per link
#endpoint1: [pid] of first particle in link
#endpoint2: [pid] of second particle
#avg. x: link midpoint. only ever calculated twice: on generation, and when the link attempts to grow
#avg. y:
#active?: bool if the link can grow
#parent: [lid] of the link that spawned this
#child1: [lid] of the first link spawned by this
#child2: [lid] of the second link spawned by this
#--- full loops through this matrix (e.g. links[,3] == TRUE) should ABSOLUTELY MINIMIZED, ideally no more than once per tick

	links = matrix(0,nrow=n_particles*2,ncol=8) #links gets doubled in size any time it's almost full

	list(particles=particles,chunks=chunks,particle_neighbors=particle_neighbors,particle_links=particle_links,links=links)
}
default_init = Init()

# ---- Sim Initialization ----
Sim = function(rules=Rules(),state=default_init){
	#check chunksize and rechunk if it's above TODO
	rules$chunksize = 0.25

	# ---- Spawn Initial Link ----


	set.seed(rules$seed) #set seed for model growth
	list(state=state,rules=rules,time=1,rng=.Random.seed)
}




# ---- Tick ----
#unfortunately very delicate and optimized; would not recommend editing.
tick = function(sim, n_ticks=1){
	list2env(state) #unpack state into local variables to edit
	rules = lapply(sim$rules,parameter) #unpack rules into local functions
	chunksize = sim$rules$chunksize

	.Random.seed <<- sim$rng #restore rng state

	#loop in here instead of a separate function to reduce overhead
	for (t in 1:n_ticks){
		# ---- Grow Links (in series) ----
		active_lids = which(links[,5] > 0)
		if (is.function(sim$rules$growth_rate)) growth_rates = rules$growth_rate(links[active_lids,3:4]) #growth rate heatmap
		else growth_rates = rules$growth_rate()
		lids_to_grow = active_lids[runif(length(active_lids)) < growth_rates]
		links[lids_to_grow,5] == 0 #deactivate growing links

		#update link midpoints (x,y)
		links[lids_to_grow,3] = (particles[links[lids_to_grow,1],1] + particles[links[lids_to_grow,2],1]) / 2
		links[lids_to_grow,4] = (particles[links[lids_to_grow,1],2] + particles[links[lids_to_grow,2],2]) / 2

		#grow in series
		for (lid in lids_to_grow){
			link = links[lid,]
			# ---- Get Closest Particle ----
			closest_particle = 0 #pid
			link_range = rules$link_range(link[3:4])
			best_dist = link_range #can't bond with particles further than bond range
			preexisting_links = c(FALSE,FALSE) #which links are already there?

			nearby_pids = unlist(chunks[get_adjacent_chunks(ceiling(link[3:4]/chunksize))])

			#search through all nearby particles
			for (pid in nearby_pids){

			}
		}

		# ---- Contract Particles (in parallel)----


		# ---- Cleanup ----
		#grow links matrix TODO
	}


	sim$rng = .Random.seed #save rng state
	sim$time = sim$time + n_ticks
	sim$state = list(particles=particles,particle_neighbors=particle_neighbors,particle_links=particle_links,links=links) #repack into sim object
	sim #return sim
}


get_adjacent_chunk_ids = function(cpos){
	chunks = rbind(
		cpos,
		cpos + c(1,0),
		cpos + c(1,1),
		cpos + c(1,-1),
		cpos + c(-1,0),
		cpos + c(-1,1),
		cpos + c(-1,-1),
		cpos + c(0,1),
		cpos + c(0,-1)
	)

	culllist = c() #remove invalid chunk ids
	if (cpos[1] <= 0) {
		if (cpos[2] <= 0) chunks = chunks[-c(4,5,6,7,9),]
		else chunks = chunks[-c(5,6,7),]
	}
	else if (cpos[2] <= 0) chunks = chunks[-c(4,7,9),]

	return(chunks)
}



