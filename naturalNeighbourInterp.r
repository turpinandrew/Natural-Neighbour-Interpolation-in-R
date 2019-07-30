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

    d1<-deldir(xs,ys,plotit=plotIt, suppressMsge=TRUE)# Get Voronoi tiles without (x,y) 
    if (x < d1$rw[1] || x > d1$rw[2])
        warning(sprintf("getNNWeights: x (%f) is outside the rectangle bounding xs",x))
    if (y < d1$rw[3] || y > d1$rw[4])
        warning(sprintf("getNNWeights: y (%f) is outside the rectangle bounding ys",y))

    d2<-deldir(c(xs,x),c(ys,y),plotit=FALSE, suppressMsge=TRUE) # Get tiles with (x,y)

        # Work out the changed tiles between d1 and d2 (ignoring new one),
        # compute the differences in their area, and then normalise.
    areaChange <- d1$summary$dir.area - head(d2$summary$dir.area, -1)

    changed <- which(areaChange > 0)

    if (plotIt) {
        n <- nrow(d1$summary)
        newSegs <- rbind(d2$dirsgs[d2$dirsgs$ind1 == n+1, 1:4],
                     d2$dirsgs[d2$dirsgs$ind2 == n+1, 1:4])
        text(xs, ys, paste("(",1:n, ") ",zs,"dB",sep=""))
        points(x,y,col="red",pch=19)
        
        segments(newSegs[,1], newSegs[,2], newSegs[,3], newSegs[,4], col=rgb(1,0,0,0.2))
    }

    return(list(neighbours=changed, weights=areaChange[changed]/sum(areaChange)))
}

##############################################################
# Return value for (x,y) interp'ing on (xs,ys,zs)
##############################################################
nnInterp <- function(xs, ys, zs, x, y, plotIt=FALSE) {
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
###test(xs,ys,zs,0.5,0.5, 2.867021)
###
####abline(h=c(-1.2, 1.2, -0.5, 0.3), v=c(-0.5, 0.3, 1.2, -1.2), col="green")
#### Pt 1 - tri with area 0.5 * 0.5^2 = 1/8
#### Pt 2 - trap area 0.5 * (0.5 + 0.1) * 1.2 = 0.36
#### Pt 3 - rec + rec + tri = 0.3*1.2 + (1.2-0.3)*0.3 + 0.5*(1.2-0.3)^2 = 1.035
#### Pt 4 - trap 0.5 * (0.5 + 0.1) * 1 = 0.36
#### sum(c(1/8,0.36,1.035,0.36)/sum(c(1/8, 0.36,1.035,0.36)) * 1:4) = 2.867021
###
###test(xs,ys,zs,20,20, 0)  # should warn
###
###test(c(-15,-9,-9),c(-9,-9,-15),c(30,30,30),-13,-15, 30)
