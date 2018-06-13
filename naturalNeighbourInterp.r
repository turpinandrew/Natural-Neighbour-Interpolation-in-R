#
# Natural neighbour interpolation
# Defines: 
#   nnInterp(xs, ys, zs, x, y, plotIt=FALSE)    # return estimate of (x,y)
#   getNN(xs, ys, i)                            # return indexes of NN of (x,y) in (xs,ys)
#   getNNWeights(xs, ys, x, y)            # return indexes and weights if (x,y) is added to (xs,ys)
#
# Andrew Turpin
# Tue 26 Mar 2013 09:12:37 EST
#
# Modified 6 Jun 2013: added check to nnInterp to see if (x,y) already in zs
# Modified 7 Jul 2014: added rw to first deldir in case x or y is outside xs or ys
# Modified 13 Jun 2018: altered algorithm to be correct!
#

require(deldir)

##############################################################
# Return value for (x,y) interp'ing on (xs,ys,zs)
##############################################################
nnInterp <- function(xs, ys, zs, x, y, plotIt=FALSE) {
    n <- length(xs)

    if (length(ys) != n || length(zs) != n) {
        warning("xs, ys and zs not same length in nnInterp")
        return(NA)
    }

        # check if (x,y) is not already in zs
    ixy <- intersect(which(xs == x), which(ys == y))
    if (length(ixy) != 0) {
        return(zs[ixy])
    }

    rw <- c(range(xs,x)[1], range(xs,x)[2], range(ys,y)[1], range(ys,y)[2])

        # Get Voronoi tiles without (x,y)  and plot
    d1<-deldir(xs,ys,plotit=plotIt, wlines="tess", pch=19, rw=rw, suppressMsge=TRUE)

        # Get Voronoi tiles with (x,y)
    d2<-deldir(c(xs,x),c(ys,y),rw=rw, plotit=FALSE, suppressMsge=TRUE)

        # Work out the changed tiles between d1 and d2 (ignoring new one),
        # compute the differences in their area, and then 
        # take the weighted sum of their zs values.
    changedTiles <- d2$dirsgs[d2$dirsgs$ind1 == n+1, "ind2"]

    deltaAreas <- d1$summary[,"dir.area"][changedTiles] - d2$summary[,"dir.area"][changedTiles]

    z <- sum(zs[changedTiles] * deltaAreas)/sum(deltaAreas)
    if (plotIt) {
        newSegs <- d2$dirsgs[d2$dirsgs$ind1 == n+1, 1:4]
        text(xs, ys, paste("(",1:n, ") ",zs,"dB",sep=""), pos=4)
        points(x,y,col="red",pch=19)
        #segments(newSegs[,1], newSegs[,2], newSegs[,3], newSegs[,4], col="red", lty=3)
        polygon(c(newSegs[,1], newSegs[,3]), c(newSegs[,2], newSegs[,4]), col=rgb(1,0,0,0.2), border=NA)

        text(x, y, paste("(",n+1, ") ",round(z,3),"dB",sep=""), pos=4, col="red")
    }

    return(z)
}

##############################################################
# Return indexes of NN of (xs[i],ys[i]) in (xs,ys)
##############################################################
getNN <- function(xs, ys, i=length(xs)) {

    if (length(ys) != length(xs)) {
        warning("xs and ys not same length in nnInterp")
        return(NA)
    }

    d <- deldir(xs,ys,plotit=FALSE, suppressMsge=TRUE)
    nn1 <- d$dirsgs[d$dirsgs$ind1 == i, "ind2"]
    nn2 <- d$dirsgs[d$dirsgs$ind2 == i, "ind1"]

    return(c(nn1,nn2))
}

##############################################################
# Return the weights for each neighbour if (x,y) is added to (xs,ys)
##############################################################
getNNWeights <- function(xs, ys, x, y) {

    if (length(ys) != length(xs)) {
        warning("xs and ys not same length in nnInterp")
        return(NA)
    }

    d1<-deldir(xs,ys,plotit=FALSE, suppressMsge=TRUE)# Get Voronoi tiles without (x,y) 

    d2<-deldir(c(xs,x),c(ys,y),plotit=FALSE, suppressMsge=TRUE) # Get Voronoi tiles with (x,y)

        # Work out the changed tiles between d1 and d2 (ignoring new one),
        # compute the differences in their area, and then 
        # take the weighted sum of their zs values.
    changedTiles <- d2$dirsgs[d2$dirsgs$ind1 == length(xs)+1, "ind2"]

    deltaAreas <- d1$summary[,"dir.area"][changedTiles] - d2$summary[,"dir.area"][changedTiles]

    return(list(neighbours=changedTiles, weights=deltaAreas/sum(deltaAreas)))
}
