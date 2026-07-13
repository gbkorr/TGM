

# ---- Basic Drawing ----
draw = function(sim,type='l'){
	region = c(0,sim$state$p_rules$size)
	plot(NULL,xlab='',ylab='',axes=FALSE,xlim = region, ylim = region)

	state = sim$state

	switch(type,
		l=,links = draw_links(state), #wow I can't believe switch works like this
		t=,tris = draw_tris(state),
		n=,network = draw_network(state),
		bands = {
			cols = color_bands(10)
			draw_tris(state,\(link)cols(link[6]))
		},
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

draw_tris = function(state,color=function(link)'black'){
	particles=state$particles
	links=state$links
	for (l in 1:sum(links[,1] != 0)){
		link = links[l,]
		pid1 = link[1]
		pid2 = link[2]
		pid3 = links[link[7],2] #get the parent's endpoint that's not included in this link
		#this double-draws most links, which isn't great.

		polygon(particles[c(pid1,pid2,pid3),1:2],border=NA,col=color(link))
	}
}

draw_network = function(state,weight=function(link)1){
	links=state$links
	for (l in 1:sum(links[,1] != 0)){
		link = links[l,]
		child1 = links[link[8],]
		child2 = links[link[9],]

		lines(rbind(link[3:4],child1[3:4]),lwd=weight(link))
		lines(rbind(link[3:4],child2[3:4]),lwd=weight(link))
	}
}



# ---- Color ----
#this is useful for analyzing the evenness of growth
color_bands = function(period,col1='black',col2='white')\(x)colorRampPalette(c(col1,col2,col1))(period)[1 + floor(x %% period)]

draw_color = function(state){
		color = color_bands(50)
		particles=state$particles
		links=state$links
		for (l in 1:sum(links[,1] != 0)){
			link = links[l,]
			pid1 = link[1]
			pid2 = link[2]
			pid3 = links[link[7],2] #get the parent's endpoint that's not included in this link
			#this double-draws most links, which isn't great.

			polygon(particles[c(pid1,pid2,pid3),1:2],border=NA,col=color(link[6]))
		}
}


