

# ---- Basic Drawing ----
draw = function(sim,type='l',args=NULL){
	region = c(0,sim$state$p_rules$size)
	plot(NULL,xlab='',ylab='',axes=FALSE,xlim = region, ylim = region)

	state = sim$state

	#get line weighting
	if (type %in% c('d','desc','o','order')){
		if(is.null(args)) args = list(NULL)
		args = modifyList(list(
			min = 1,
			max = 10,
			scale = 0.01
		),args)
		weight = root_weight(args$min,args$max,args$scale)
	}

	switch(type,
		l=,links = draw_links(state), #wow I can't believe switch works like this
		t=,tris = draw_tris(state),
		n=,network = draw_network(state),
		b=,bands = { #color by generation
			cols = color_bands(ifelse(is.null(args),100,args))
			draw_tris(state,\(link)cols(link[6]))
		},
		d=,desc = draw_network(state,\(node)weight(node[8])), #number of descendants
		o=,order = draw_network(state,\(node)weight(node[7])), #number of descendant leaves
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
		pid3 = links[links[l,7],1:2] #of the parent's points, choose the not-shared one
		pid3 = pid3[!(pid3 %in% c(pid1,pid2))][1]
		#this double-draws most links, which isn't great.

		polygon(particles[c(pid1,pid2,pid3),1:2],border=NA,col=color(link))
	}
}

draw_network = function(state,weight=function(node)1){
	network = Network(state$links)
	for (l in 1:nrow(network)){
		node = network[l,]
		child1 = network[node[4],]
		child2 = network[node[5],]

		lines(rbind(node[1:2],child1[1:2]),lwd=weight(node))
		lines(rbind(node[1:2],child2[1:2]),lwd=weight(node))
	}
}

# ---- Color ----
#this is useful for analyzing the evenness of growth
root_weight = function(min,max,s)\(d)(max-(max-min)*(1+s)^(-d))
color_bands = function(period,col1='black',col2='white')\(x)colorRampPalette(c(col1,col2,col1))(period)[1 + floor(x %% period)]



