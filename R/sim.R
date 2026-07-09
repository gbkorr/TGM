

# ---- Parameter Initialization ----
# rules prototype
Rules = function(
	link_range = 0.18,
	mobility = 0.04,
	contraction = 1.0,
	cohesion = 0.0,
	branching = 1.0,

	seed = 1,
	seed_pos = NULL, #if defined, exact coordinates of starting link
	growth_rate = 0.05 #usually constant
) list(as.list(environment()))[[1]]

#evaluates parameter as either constant or a function of position
parameter = function(rule) ifelse(is.function(rule),rule,function(pos) rule)

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
	particles = matrix(c(runif(2*n_particles,0,p_rules$size),rep(-1,n_particles)),ncol=3)

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
#generation:
#parent: [lid] of the link that spawned this
#child1: [lid] of the first link spawned by this
#child2: [lid] of the second link spawned by this
#--- full loops through this matrix (e.g. links[,3] == TRUE) should ABSOLUTELY MINIMIZED, ideally no more than once per tick

	links = matrix(0,nrow=n_particles*2.5,ncol=9) #links gets doubled in size any time it's almost full

	list(particles=particles,chunks=chunks,particle_neighbors=particle_neighbors,particle_links=particle_links,links=links)
}
default_init = Init()

# ---- Sim Initialization ----
Sim = function(rules=Rules(),state=default_init){
	#check chunksize and rechunk if it's above TODO
	rules$chunksize = 0.25

	# ---- Spawn Initial Link ----
	if (!is.null(rules$seed_pos)) state$particles[1,1:2] = rules$seed_pos
	else state$particles[1,1:2] = c(max(state$particles[,1]),max(state$particles[,2]))/2
	state$links[1,] = c(1,1,state$particles[1,1],state$particles[1,2],1,1,0,0,0)

	set.seed(rules$seed) #set seed for model growth
	list(state=state,rules=rules,time=1,rng=.Random.seed)
}




# ---- Tick ----
#unfortunately very delicate and optimized; would not recommend editing.
tick = function(sim, n_ticks=1){
	list2env(sim$state,environment()) #unpack state into local variables to edit
	rules = lapply(sim$rules,parameter) #unpack rules into local functions
	chunksize = sim$rules$chunksize

	.Random.seed <<- sim$rng #restore rng state

	total_lids = sum(links[,1] != 0) + 1 #first empty slot in links

	#loop in here instead of a separate function to reduce overhead
	for (t in 1:n_ticks){
		# ---- Grow Links (in series) ----
		active_lids = which(links[,5] == 1)
		if (is.function(sim$rules$growth_rate)) growth_rates = rules$growth_rate(links[active_lids,3:4]) #growth rate heatmap
		else growth_rates = rules$growth_rate()
		lids_to_grow = active_lids[runif(length(active_lids)) < growth_rates]
		if (length(lids_to_grow)){ #if there are any new links to grow
			links[lids_to_grow,5] = 0 #deactivate growing links

			#update link midpoints (x,y)
			links[lids_to_grow,3] = (particles[links[lids_to_grow,1],1] + particles[links[lids_to_grow,2],1]) / 2
			links[lids_to_grow,4] = (particles[links[lids_to_grow,1],2] + particles[links[lids_to_grow,2],2]) / 2

			start_of_new_lids = total_lids

			#grow in series
			for (lid in lids_to_grow){
				link = links[lid,]

				#get parameters for current link pos
				pos = link[3:4]
				link_range = rules$link_range(pos)
				cohesion = rules$cohesion(pos)
				contraction_timer = rules$contraction(pos) / (rules$mobility(pos) * rules$growth_rate(pos))

				# ---- Get Closest Particle ----
				nearby_pids = unlist(chunks[get_adjacent_chunk_ids(ceiling(pos/chunksize))])

				#search through all nearby particles
				nearby_pids = nearby_pids[nearby_pids != link[1] & nearby_pids != link[2]] #ignore your constituents
				vecs = cbind(
					particles[nearby_pids,1] - link[3],
					particles[nearby_pids,2] - link[4]
				) #vectors from link midpoint to particles
				dists = apply(vecs,1,\(xy)sqrt(sum(xy^2)))
				already_linked = particles[nearby_pids,3] != -1
				if (cohesion > 0) dists[already_linked] = dists[already_linked]/cohesion #apply cohesion penalty
				else dists[already_linked] = Inf

				closest_particle = which.min(dists)
				preexisting_links = c(0,0) #which links are already there? lid of preexisting link w/ the closest particle connected to endpoint1, endpoint2. 0 = does not exist
				if (dists[closest_particle] > link_range) break #invalid; closest particle is out of link_range
				else {
					pid = nearby_pids[closest_particle]
					if (particles[pid,3] != -1) { #particle is already linked with something
						bp = particle_links[[pid]] #get all links connected to that particle
						#we want to make sure we don't duplicate a link
						preexisting_links[1] = link[1] %in% bp #shared links
						preexisting_links[2] = link[2] %in% bp

						sum_pre = sum(preexisting_links)
						if (sum_pre == 2) break #already bonded with both endpoints; invalid
						else if (sum_pre == 1){
							#expensive, unavoidable intersect() call to see what link is already there
							if (l1) preexisting_links = c(intersect(particle_links[[endpoints[1]]],particle_links[[pdx]]),0)
							else if (l2) preexisting_links = c(0,intersect(particle_links[[endpoints[2]]],particle_links[[pdx]]))
						}
					}
				}

				# ---- Link with Closest Particle
				if (particles[pid,3] == -1) particles[pid,3] = contraction_timer #assign contraction timer if particle was free. uses location of link rather than particle; this doesn't really matter

				active_sides = c(TRUE,TRUE) #for branching, TODO

				if (preexisting_links[1] == 0){ #if there's not already a link between link[1] and the closest particle
					links[lid,8] = total_lids #set child
					links[total_lids,] = c(link[1],pid,0,0,active_sides[1],link[6]+1,lid,0,0) #create child
					particle_neighbors[[link[1]]] = c(particle_neighbors[[link[1]]],pid) #record new particle neighbors
					particle_neighbors[[pid]] = c(particle_neighbors[[pid]],link[1])
					particle_links[[link[1]]] = c(particle_links[[link[1]]],total_lids) #record links for particles
					particle_links[[pid]] = c(particle_links[[pid]],total_lids)
					total_lids = total_lids + 1
				}
				else links[lid,8] = preexisting_links[1] #IMPORTANT: records this new connection in the network; the preexisting link becomes a child because it's included in the triangle.
				#without the above, the network would not accurately represent the interconnected triangles
				if (preexisting_links[2] == 0){ #if there's not already a link between link[1] and the closest particle
					links[lid,8] = total_lids #set child
					links[total_lids,] = c(link[2],pid,0,0,active_sides[2],link[6]+1,lid,0,0)
					particle_neighbors[[link[2]]] = c(particle_neighbors[[link[2]]],pid)
					particle_neighbors[[pid]] = c(particle_neighbors[[pid]],link[2])
					particle_links[[link[2]]] = c(particle_links[[link[2]]],total_lids)
					particle_links[[pid]] = c(particle_links[[pid]],total_lids)
					total_lids = total_lids + 1
				}
				else links[lid,8] = preexisting_links[2]
			}

			#update position of new links
			if (start_of_new_lids != total_lids){
				new_lids = start_of_new_lids:(total_lids - 1)
				links[new_lids,3] = (particles[links[new_lids,1],1] + particles[links[new_lids,2],1]) / 2
				links[new_lids,4] = (particles[links[new_lids,1],2] + particles[links[new_lids,2],2]) / 2
			}
		}

		# ---- Contract Particles (in parallel)----
		active_particles = which(particles[,3] > 0)
		n_active = length(active_particles)
		if (n_active){ #if there are any particles to contract
			new_particles = matrix(0,nrow=length(active_particles),ncol=3) #will overwrite particles later to apply changes in parallel

			for (p in 1:n_active){
				pid = active_particles[p]
				particle = particles[pid,]
				mobility = rules$mobility(particle[1:2]) * rules$growth_rate(particle[1:2])
				neighbors = particle_neighbors[[pid]]
				centroid = sum_coords(particles[neighbors,1:2,drop=FALSE])

				particle[1:2] = particle[1:2] * (1 - mobility) + mobility * centroid / length(neighbors)

				particle[3] = particle[3] - 1 #decrement contraction timer
				new_particles[p,] = particle
			}

			particles[active_particles,] = new_particles
		}
		# ---- Cleanup ----
		#if link matrix is almost full, double its size
		if (total_lids > (nrow(links) - length(active_lids) * 4)) links = rbind(links,matrix(0,nrow(links),ncol(links)))
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
} #chunkpos
sum_coords = function(M) c(sum(M[,1]),sum(M[,2]))


test = function() Sim() |> tick(1000) |> draw()
#TESTING ONLY draw function
draw = function(sim){
	particles=sim$state$particles
	links=sim$state$links
	plot(NULL,xlim=range(particles[,1]),ylim=range(particles[,2]))
	for (l in 1:sum(links[,1] != 0)){
			p1 = particles[links[l,1],1:2]
			p2 = particles[links[l,2],1:2]
			lines(c(p1[1],p2[1]),c(p1[2],p2[2]))
	}
}

