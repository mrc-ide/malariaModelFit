
# -----------------------------------
#' plot_nice1
#'
#' Produces plot similar to base plot() function, but in ggplot2. Has the advantage of being able to specify legends outside of plotting area more easily. Multiple x, y and col arguments can be specified using lists of values.
#'
#' @export

plot_nice1 <- function(x, y, col=1, xlab=NULL, ylab=NULL, main=NULL, pch=1, cex=1, legend=FALSE, legend_position="right", legend_title="legend_title", legend_label=NULL) {
    
    # if no y variable then first argument becomes y
    if (missing(y)) {
        y <- x
        if (is.list(y)) {
            x <- mapply(function(x){1:length(x)}, y)
        } else {
            x <- 1:length(y)
        }
    }
    
    # try to set axis names automatically
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
    
    # replicate col over list
    if (!is.list(col)) {
        col <- replicate(length(y), col, simplify=FALSE)
    }
    
    # check that all lists same length
    if(!(length(x)==length(y) & length(x)==length(col))) {
        stop("input lists of different lengths")
    }
    
    # check that x and y inputs same length
    if (!all(mapply(length,x)==mapply(length,y))) {
        stop("x and y inputs must be corresponding lengths")
    }
    
    # add each list element to data frame
    df1 <- NULL
    for (i in 1:length(x)) {
        
        # if col is numeric then convert to default colours
        if (is.numeric(col[[i]])) {
            col[[i]] <- (col[[i]]-1)%%8+1	# recycle colours
            col[[i]] <- palette()[col[[i]]]
        }
        
        # create data frame
        col[[i]] <- rep(col[[i]], length(x[[i]]))
        df1 <- rbind(df1, data.frame(x=x[[i]], y=y[[i]], col=col[[i]]))
    }
    
    # get unique colour levels
    col_levels <- levels(df1$col)
    
    # create default legend labels
    if (is.null(legend_label)) {
        legend_label <- paste0("legend_label", 1:length(col_levels))
    }
    
    # --------------------
    
    # create plot and theme
    plot1 <- ggplot()
    plot1 <- plot1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    plot1 <- plot1 + theme(panel.background=element_blank())
    plot1 <- plot1 + theme(axis.line=element_line(colour="black", size=0.25))
    plot1 <- plot1 + theme(plot.title=element_text(hjust=0.5))
    
    # add legend
    if (!legend) {
        plot1 <- plot1 + theme(legend.position="none")
    } else {
        plot1 <- plot1 + theme(legend.position=legend_position)
    }
    
    # add points
    plot1 <- plot1 + geom_point(aes(x, y, colour=col), shape=pch, size=2*cex, data=df1)
    
    # add titles
    plot1 <- plot1 + labs(x=xlab, y=ylab)
    plot1 <- plot1 + ggtitle(main)
    
    # add legend
    plot1 <- plot1 + scale_colour_manual(labels=legend_label, values=col_levels, name=legend_title)
    
    # plot
    plot1
}
