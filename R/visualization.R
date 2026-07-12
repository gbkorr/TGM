


draw = function(sim,type,fun=NULL){
	region = c(0,sim$state$p_rules$size)
	plot(NULL,xlab='',ylab='',axes=FALSE,xlim = region, ylim = region)

	state = sim$state
	switch(type,
		l=,links = draw_links(state,fun), #wow I can't believe switch works like this
		t=,tris = draw_tris(state,fun),
		n=,network = draw_network(state,fun),
		warning("Unknown draw type. Options: 'links', 'tris', 'network'") #"l", "t", "n" also work
	)
}

draw_links = function(state,fun=NULL){
	particles=state$particles
	links=state$links
	for (l in 1:sum(links[,1] != 0)){
		pid1 = links[l,1]
		pid2 = links[l,2]
		lines(particles[c(pid1,pid2),1:2])
	}
}

draw_tris = function(state,col=NULL){
	if (is.null(col)) col = function(link, state) "black"
	particles=state$particles
	links=state$links
	for (l in 1:sum(links[,1] != 0)){
		link = links[l,]
		pid1 = link[1]
		pid2 = link[2]
		pid3 = links[link[7],2] #get the parent's endpoint that's not included in this link
		#this double-draws most links, which isn't great.

		polygon(particles[c(pid1,pid2,pid3),1:2],border=NA,col=col(link,state))
	}
}

draw_network = function(state,lwd=NULL){
	if (is.null(lwd)) lwd = function(link, state) 1
	links=state$links
	for (l in 1:sum(links[,1] != 0)){
		link = links[l,]
		parent_link = links[link[7],]

		lines(rbind(link[3:4],parent_link[3:4]),lwd=lwd(link,state))
	}
}



#example: draw tris colored by generation
#draw(sim,'t',\(link,state)rainbow(max(state$links[,6]))[link[6]])



