require(igraph)

#' Generate a random graph with desired size, clustering coefficient, and degree distribution
#'
#' This uses a fast network growth algorithm. It is not an exact 
#' algorithm, so the resulting clustering coefficient may differ 
#' slightly from the requested value. Users may provide their own 
#' R functions to sample from the degree distribution.
#'
#' Citation: Volz, EM, Physical Review E, 2004
#'
#' @param n Integer size of random graph
#' @param CC Desired clustering coefficient. Value in (0,1). 
#' @param rdegdist This is either a character string specifiying the desired degree distribution or an R function which generates the desired distribution.
#' @param ... Additional arguments are passed to rdegdist (such as lambda if using Poisson)
#' @return A graph in igraph format 
#' @examples 
#' rclustnet( n=500, CC = .25, rdegdist='pois', lambda=5 )
rclustnet <- function( n, CC, rdegdist = 'pois', ...)
{
	if (class(rdegdist)=='character'){
		if (rdegdist == 'pois'){
			rdegdist <- function(n) rpois( n, ... )
		} else if (rdegdist == 'nbinom'){
			rdegdist <- function(n) rnbinom(n, ...)
		} else if( rdegdist=='geom'){
			rdegdist <- function(n) rgeom(n, ... )
		} else{
			stop( 'Given rdegdist is not supported. Provide function to generate degree distribution (pois, geom, nbinom).')
		}
	}
	if (class( rdegdist)!='function') 
	  stop( 'Given rdegdist is not supported. Provide function to generate degree distribution (pois, geom, nbinom).')
	
	CC <- min(1, max(0, CC))
	
	d <- rdegdist( n )
	nstubs <- d  # will decrease as links formed 

	nbrs <- lapply( 1:n, function(i) rep( NA, d[i] ) )

	md <- max( d ) 
	p_tauij_x <- sapply( 1:md, function(k) sapply( 0:(md-1) , function(x){
		dbinom( x, size=k-1, prob=CC )
	}))

	W <- sapply( 1:md, function(k) sapply( 1:md, function(kk){
		sum(   p_tauij_x[, k] * p_tauij_x[, kk]   )  
	}))

	nzi <- setdiff( 1:n, which( d == 0 ) ) # non zero indices 

	while( length( nzi ) > 1 )
	{
	wave <- c( nzi[1] )
	nextwave <- c()
		iwave <- 0
	repeat { 
		for (u in wave ){
			# form cluster links
			# two steps away
			tsa <- suppressWarnings( na.omit( ( unlist( nbrs[ nbrs[[u]] ]  )  ) ) )
			tsa <- tsa[ nstubs[tsa] > 0 ]
			tsa <- setdiff( setdiff( tsa, nbrs[[u]] ), u)
			ntsa <- length(tsa)
			nc <- min( nstubs[u], rbinom( 1, size = ntsa, prob = CC*2 ) )
			if (nc > 0 ){
				newtsa <- unlist( sample( as.list(tsa), size = nc , replace=FALSE)  )
				for (v in newtsa){
					if (nstubs[v] > 0)
					{
						nbrs[[u]][ (d[u] - nstubs[u] + 1) ] <- v
						nstubs[u] <- nstubs[u] - 1
						nbrs[[v]][ (d[v] - nstubs[v] + 1) ] <- u
						nstubs[v] <- nstubs[v] - 1
					} else{
						stop('Error encountered. Could not add 3-clique to graph.')
					}
				}
			}
			
			if ( nstubs[u] > 0 ){
				notconn <- setdiff( setdiff( nzi, u ) , na.omit( nbrs[[u]] ) )  
				w <- pmax(0, W[ d[notconn], d[u] ] * nstubs[notconn] )
				if (length(notconn)==1){
					newnbrs <- ifelse( w > 0, notconn, NULL)
				} else{
					newnbrs <- tryCatch( { sample( notconn, size=nstubs[u], prob = w , replace=FALSE )
					}, error = function(e) NULL)
				}
				if (length( newnbrs) > 0) {
					nbrs[[u]][ (d[u] - nstubs[u] + 1):(d[u]) ] <- newnbrs 
					nstubs[u] <- 0
					for ( v in newnbrs ){
						nbrs[[v]][ d[v] - nstubs[v] + 1 ] <- u
						nstubs[v] <- nstubs[v]  - 1
					}
					nextwave <- c( nextwave , setdiff( newnbrs, wave ))
				}
			}
		}
		
		wave <- unique( nextwave )
		nextwave <- c() 
		iwave <- iwave + 1
		
		if ( length( wave ) < 1 ) break 
	}
	nzi <- which( nstubs > 0 ) # non zero indices 
	}
	
	el <- matrix( NA, nrow=(length( unlist(nbrs))), ncol=2)
	k <- 1
	for (u in 1:length(nbrs)){
		for (v in na.omit(nbrs[[u]]) ){
			el[k,1] <- u
			el[k,2] <- v
			k <- k + 1
		}
	}
	el <- el[ !is.na(el[,1]) , ]
	simplify( graph_from_edgelist( el, directed=FALSE) )
}

