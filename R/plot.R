
# -----------------------------------
#' mplot
#'
#' Convenient function for plotting multiple lines/points on the same graph. Multiple x, y, col etc. arguments can be specified using lists of values. Two legends can also be specified manually.
#'
#' @export

mplot <- function(x, y, type="l", col=1, lty=1, lwd=1, pch=1, cex=1, xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL, xlab=NULL, ylab=NULL, legend_margin=5, legend_inset=0.2, legend1_title="legend1_title", legend1_pos="topright", legend1_label=NULL, legend1_col=1, legend1_lty=1, legend1_lwd=1, legend1_pch=NULL, legend1_cex=1, legend1_boxOn=FALSE, legend2_title="legend2_title", legend2_pos="right", legend2_label=NULL, legend2_col=1, legend2_lty=1, legend2_lwd=1, legend2_pch=NULL, legend2_cex=1, legend2_boxOn=FALSE, ...) {
	
	# if no y variable then first argument becomes y, and x becomes sequential numbers of same length
    if (missing(y)) {
        y <- x
        if (is.list(y)) {
            x <- mapply(function(x){1:length(x)}, y, SIMPLIFY=FALSE)
        } else {
            x <- 1:length(y)
        }
    }
    
    # try to set axis names automatically from names of input arguments
    if (is.null(xlab) && !is.list(x)) {
        xlab <- substitute(x)
    }
    if (is.null(ylab) && !is.list(y)) {
        ylab <- substitute(y)
    }
    
    # force y to list
    if (!is.list(y)) {
        y <- list(y)
    }
    
    # replicate x over list
    if (!is.list(x)) {
        x <- replicate(length(y), x, simplify=FALSE)
    }
    
    # replicate other attributes over list
    if (!is.list(type)) {
        type <- replicate(length(y), type, simplify=FALSE)
    }
    if (!is.list(col)) {
        col <- replicate(length(y), col, simplify=FALSE)
    }
    if (!is.list(lty)) {
        lty <- replicate(length(y), lty, simplify=FALSE)
    }
    if (!is.list(lwd)) {
        lwd <- replicate(length(y), lwd, simplify=FALSE)
    }
    if (!is.list(pch)) {
        pch <- replicate(length(y), pch, simplify=FALSE)
    }
    if (!is.list(cex)) {
        cex <- replicate(length(y), cex, simplify=FALSE)
    }
    
    # check that all lists same length
    l <- c(length(x), length(y), length(type), length(col), length(lty), length(lwd), length(pch), length(cex))
    if(length(unique(l))!=1) {
        stop("input lists are of different lengths")
    }
    
    # check that x and y inputs same length
    if (!all(mapply(length,x)==mapply(length,y))) {
        stop("x and y inputs must be corresponding lengths")
    }
    
    # define min and max functions that return NA if all values are NA
    safeMin <- function(x) {
    	if (all(is.na(x))) {
    		return(NA)
    	} else {
    		return(min(x,na.rm=TRUE))
    	}
    }
    safeMax <- function(x) {
    	if (all(is.na(x))) {
    		return(NA)
    	} else {
    		return(max(x,na.rm=TRUE))
    	}
    }
    
    # set default plotting limits
    if (is.null(xmin)) {
    	xmin <- safeMin(mapply(safeMin, x))
    }
    if (is.null(xmax)) {
    	xmax <- safeMax(mapply(safeMax, x))
    }
    if (is.null(ymin)) {
    	ymin <- safeMin(mapply(safeMin, y))
    }
    if (is.null(ymax)) {
    	ymax <- safeMax(mapply(safeMax, y))
    }
    
    # change margins
    oldmar <- newmar <- par(mar=rep(0,4), xpd=TRUE)
    anyLegend <- (!is.null(legend1_label) | !is.null(legend2_label))
    if (anyLegend) {
    	newmar$mar[4] <- newmar$mar[4] + legend_margin
    	newmar$xpd <- TRUE
    }
    par(newmar)
    
    # produce empty plot and add lines
	plot(0, type="n", xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=xlab, ylab=ylab, ...)
    for (i in 1:length(x)) {
    	lines(x[[i]], y[[i]], type=type[[i]], col=col[[i]], lty=lty[[i]], lwd=lwd[[i]], pch=pch[[i]], cex=cex[[i]])
    }
	
	# add first legend
	if (!is.null(legend1_label)) {		
		legend(x=legend1_pos, legend=legend1_label, col=legend1_col, lty=legend1_lty, lwd=legend1_lwd, pch=legend1_pch, cex=legend1_cex, inset=c(-legend_inset,0), title=legend1_title, title.adj=1, bty=switch(legend1_boxOn+1,"n","o"))
	}
	
	# add second legend
	if (!is.null(legend2_label)) {
		legend(x=legend2_pos, legend=legend2_label, col=legend2_col, lty=legend2_lty, lwd=legend2_lwd, pch=legend2_pch, cex=legend2_cex, inset=c(-legend_inset,0), title=legend2_title, title.adj=1, bty=switch(legend2_boxOn+1,"n","o"))
	}
	
	# change margins back
    par(oldmar)
}
