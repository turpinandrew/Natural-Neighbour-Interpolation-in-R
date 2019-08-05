#
# Natural neighbour interpolation
# Defines: 
#   nnInterp(xs, ys, zs, x, y, plotIt=FALSE)    # return estimate of (x,y)
#   getNNWeights(xs, ys, x, y)            # return indexes and weights if (x,y) is added to (xs,ys)
#
# Andrew Turpin
# Tue 26 Mar 2013 09:12:37 EST
#
# Modified 6 Jun 2013: added check to nnInterp to see if (x,y) already in zs
# Modified 7 Jul 2014: added rw to first deldir in case x or y is outside xs or ys
# Modified 13 Jun 2018: altered algorithm to be correct!
# Modified 31 Jul 2019: complete rewrite to use areas, not segments
# Modified  6 Aug 2019: allow for points well outside of the area of xs and ys, and matrix argument
#

require(deldir)

##############################################################
# Return the weights for each neighbour of (x,y) 
# if it is added to (xs,ys).
##############################################################
getNNWeights <- function(xs, ys, x, y, plotIt=FALSE) {
    if (length(ys) != length(xs)) {
        warning("xs and ys not same length in nnInterp")
        return(list(neighbours=NA, weights=NA))
    }

    if (length(x) != 1 || length(y) != 1) {
        warning("getNNWeights not vectorised: only one x and y allowed")
        return(list(neighbours=NA, weights=NA))
    }

    rw <- c(range(xs), range(ys))
    if (x < rw[1] || x > rw[2])
        warning(sprintf("getNNWeights: x (%f) is outside the rectangle bounding xs.  Extending rectangle to include.",x))
    if (y < rw[3] || y > rw[4])
        warning(sprintf("getNNWeights: y (%f) is outside the rectangle bounding ys.  Extending rectangle to include.",y))

    rw <- c(range(xs, x), range(ys,y))
    d1<-deldir(xs,ys,rw=rw,plotit=plotIt, asp=1, suppressMsge=TRUE)# Get Voronoi tiles without (x,y) 

    d2<-deldir(c(xs,x),c(ys,y),rw=rw,plotit=FALSE, suppressMsge=TRUE) # Get tiles with (x,y)

        # Work out the changed tiles between d1 and d2 (ignoring new one),
        # compute the differences in their area, and then normalise.
    areaChange <- d1$summary$dir.area - head(d2$summary$dir.area, -1)

    changed <- which(areaChange > 0)

    if (plotIt) {
        n <- nrow(d1$summary)
        newSegs <- rbind(d2$dirsgs[d2$dirsgs$ind1 == n+1, 1:4],
                     d2$dirsgs[d2$dirsgs$ind2 == n+1, 1:4])
        #text(xs, ys, paste("(",1:n, ") ",zs,"dB",sep=""))
        points(x,y,col="red",pch=19)
        
        segments(newSegs[,1], newSegs[,2], newSegs[,3], newSegs[,4], col=rgb(1,0,0,0.2))
    }

    return(list(neighbours=changed, weights=areaChange[changed]/sum(areaChange)))
}

##############################################################
# xs can be a 3 column matrix of xs,ys,zs or
# xs,ys,zs are vectors of x y and z values supplied separately.
#
# x can be a length 2 vector of the point (x,y) to compute, or
# x and y can be provided separately.
#
# If plotIt is true, the original tiling of the space is plotted in black,
# and the new tile for the added point (x,y) is plotted in red.
#
# Return value for (x,y) interp'ing on (xs,ys,zs)
##############################################################
nnInterp <- function(xs, ys=NULL, zs=NULL, x, y=NULL, plotIt=FALSE) {

    if (is.matrix(xs) || is.data.frame(xs)) {
        ys <- xs[,2]
        zs <- xs[,3]
        xs <- xs[,1]
    }

    if (length(x) == 2) {
        y <- x[2]
        x <- x[1]
    }

    if (is.null(ys) || is.null(zs))
        stop("nnInterp cannot proceed: xs should be a 3 column matrix, or specify xs,ys,zs separately.")
    if (is.null(y))
        stop("nnInterp cannot proceed: x can be a length 2 vector of (x,y), or specify x and y separately.")

        # check if (x,y) is not already in zs
    ixy <- intersect(which(xs == x), which(ys == y))
    if (length(ixy) != 0) {
        return(zs[ixy])
    }

    nn <- getNNWeights(xs, ys, x, y, plotIt)

    if (length(nn$weights) == 1 && is.na(nn$weights)) {
        warning("xs, ys and zs not same length in nnInterp")
        return(NA)
    }

    z <- sum(zs[nn$neighbours] * nn$weights)
    return(z)
}

#########################################
# Tests
#########################################
###test <- function(xs, ys, zs, x, y, answer) {
###    #print(getNNWeights(xs,ys,x,y,plotIt=TRUE))
###    #print(nnInterp(xs,ys,zs,x,y,plotIt=TRUE))
###    print(abs(nnInterp(xs,ys,zs,x,y,plotIt=TRUE) - answer) < 1e-05)
###}
###
###xs <- c(-1,-1,1,1)
###ys <- c(-1,1,1,-1)
###zs <- rep(30, 4)
###
###test(xs,ys,zs,  0,  0, 30)
###test(xs,ys,zs,0.5,0.5, 30)
###
###zs <- 1:4
###test(xs,ys,zs,  0,  0, sum(0.25*zs))
###
###test(xs,ys,zs,0.5,0.5, 2.85)
###
###abline(h=c(-1.2, 0.5, -0.17, -0.5), v=c(-0.5, -0.17,0.5, -1.2), col="green")
#### Pt 1 - tri with area 0.5 * 0.5^2 = 1/8 = 0.125
#### Pt 2 - trap area 0.5 * (0.5 + 0.15) * 1 = 0.325
#### Pt 3 - rec - tri = 1 - 0.5^3 = 0.875
#### Pt 4 - Pt 2 = 0.325
#### sum(c(1/8,0.325,0.875,0.325)/sum(c(1/8,0.325,0.875,0.325)) * 1:4) = 2.85
###
###test(xs,ys,zs,20,20, 3)  # should warn
###
###test(c(-15,-9,-9),c(-9,-9,-15),c(30,30,30),-13,-15, 30)
###
###test(c(-21,21,-15,15), c(3,3,21,21), c(10,15,20,25), 1,1, 16.95234)
###test(c(-21,21,-15,15), c(3,3,21,21), c(10,15,20,25), -1,-1, 15.59268)
