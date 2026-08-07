


# ---- Basic Drawing ----
draw = function(sim,type='l',args=NULL,f=NULL){
	size = sim$state$p_rules$size
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

	if (is.null(f)) plot(NULL,xlab='',ylab='',axes=FALSE,xlim = c(0,size), ylim = c(0,size))
	else {
		resolution = 100
		draw_function(f,size,resolution)
		#normalize positions
		adj_size = size/resolution
		state$links[,3] = state$links[,3]/adj_size + 0.5
		state$links[,4] = state$links[,4]/adj_size + 0.5
		state$particles[,1] = state$particles[,1]/adj_size + 0.5
		state$particles[,2] = state$particles[,2]/adj_size + 0.5
	}

	switch(type,
		l=,links = draw_links(state), #wow I can't believe switch works like this
		t=,tris = draw_tris(state),
		n=,network = draw_network(state),
		b=,bands = { #color by generation
			cols = color_bands(ifelse(is.null(args),100,args))
			draw_tris(state,\(link)cols(link[11]))
		},
		d=,descendants = draw_network(state,\(node)weight(node[8])), #number of descendants
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


# ---- Utility ----
#draw a spatially varying function
draw_function = function(f,size,resolution=100,label=TRUE){
	pm = par('mar')
	on.exit(par(mar=pm)) #restore previous margins
	par(mar=c(pm + c(0,0,0,4)))

	M = matrix(0,resolution,resolution)
	for (Y in 1:resolution){
		for (X in 1:resolution){
			x = seq(0,size,size/resolution)[X]
			y = seq(0,size,size/resolution)[Y]
			M[Y,X] = f(c(x,y))
		}
	}

	M = t(M)[,nrow(M):1] #rotate because image() flips the axes

  image(1:resolution,
  			1:resolution,
  			M,
  			col = hcl.colors(resolution, "YlOrRd"),
  			axes = FALSE,xlab='',ylab='')


  #colorbar hack
  if (resolution == 100 && label){
	  res = 100
	  half_res = floor(res/2)
	  for (i in 1:res - 1) points(resolution + 4, 2 + 2 * i * 20/res,pch=15,cex=3,xpd=TRUE,col=rev(hcl.colors(res,'YlOrRd',rev=TRUE))[i])
	  points(resolution + 4,2+2*res,pch=15,cex=3,xpd=TRUE,col='#fff')
	  text(resolution + 8,2+2*c(0:4/5 * 25),labels=seq(min(M),max(M),length.out=5) |> round(4),xpd=TRUE,cex=1,adj=0)
  }
}
