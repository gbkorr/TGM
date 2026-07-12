


draw = function(sim,type){
	region = c(0,sim$state$p_rules$size)
	plot(NULL,xlab='',ylab='',axes=FALSE,xlim = region, ylim = region)

	state = sim$state
	switch(type,
		l=,links = draw_links(state), #wow I can't believe switch works like this
		t=,tris = draw_tris(state),
		n=,network = draw_network(state),
		warning("Unknown draw type. Options: 'links', 'tris', 'network'") #"l", "t", "n" also work
	)
}

draw_links = function(state){
	particles=state$particles
	links=state$links
	for (l in 1:sum(links[,1] != 0)){
		pid1 = links[l,1]
		pid2 = links[l,2]
		lines(particles[c(pid1,pid2),1:2])
	}
}

draw_tris = function(state){
	particles=state$particles
	links=state$links
	for (l in 1:sum(links[,1] != 0)){
		link = links[l,]
		pid1 = link[1]
		pid2 = link[2]
		pid3 = links[link[7],2] #get the parent's endpoint that's not included in this link
		#this double-draws most links, which isn't great.

		polygon(particles[c(pid1,pid2,pid3),1:2],border=NA,col='black')
	}
}

draw_network = function(state){
	links=state$links
	for (l in 1:sum(links[,1] != 0)){
		link = links[l,]
		parent_link = links[link[7],]

		lines(rbind(link[3:4],parent_link[3:4]))
	}
}


