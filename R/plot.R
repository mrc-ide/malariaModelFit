
# -----------------------------------
#' plot_nice1
#'
#' Produces plot similar to base plot() function, but in ggplot2. Has the advantage of being able to specify legends outside of plotting area more easily. Multiple x, y and col arguments can be specified using lists of values.
#'
#' @export

plot_nice1 <- function(x, y, col=1, lty=1, pch=1, cex=1, type="p", xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, main=NULL, legend_position="right", legend_col=FALSE, legend_col_title="legend_col_title", legend_col_label=NULL, legend_lty=FALSE, legend_lty_title="legend_lty_title", legend_lty_label=NULL, logx=FALSE, logy=FALSE) {
    
    # check input formats
    stopifnot(type%in%c("p","l"))
    
    # if no y variable then first argument becomes y
    if (missing(y)) {
        y <- x
        if (is.list(y)) {
            x <- mapply(function(x){1:length(x)}, y, SIMPLIFY=FALSE)
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
    
    # replicate lty over list
    if (!is.list(lty)) {
        lty <- replicate(length(y), lty, simplify=FALSE)
    }
    
    # check that all lists same length
    if(!(length(x)==length(y) & length(x)==length(col) & length(x)==length(lty))) {
        stop("input lists are of different lengths")
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
        
        # if lty is numeric then convert to string
        if (is.numeric(lty[[i]])) {
            lty_strings <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")
            lty[[i]] <- lty_strings[lty[[i]]]
        }
        
        # create data frame
        col[[i]] <- rep(col[[i]], length(x[[i]]))
        lty[[i]] <- as.factor(rep(lty[[i]], length(x[[i]])))
        df1 <- rbind(df1, data.frame(x=x[[i]], y=y[[i]], col=col[[i]], lty=lty[[i]], group=i))
    }
    
    # get unique colour levels
    col_levels <- levels(df1$col)
    
    # get unique linetype levels
    lty_levels <- levels(df1$lty)
    
    # set default plotting limits
    if (is.null(xlim)) {
        xlim <- range(df1$x, na.rm=TRUE)
    }
    if (is.null(ylim)) {
        ylim <- range(df1$y, na.rm=TRUE)
    }
    
    # create default legend labels (not necessarily used)
    if (is.null(legend_col_label)) {
        legend_col_label <- paste0("legend_col_label", 1:length(col_levels))
    }
    
    # --------------------
    
    # create plot and theme
    plot1 <- ggplot()
    plot1 <- plot1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    plot1 <- plot1 + theme(panel.background=element_blank())
    plot1 <- plot1 + theme(axis.line=element_line(colour="black", size=0.25))
    plot1 <- plot1 + theme(plot.title=element_text(hjust=0.5))
    
    # set plotting limits
    plot1 <- plot1 + coord_cartesian(xlim=xlim, ylim=ylim)
    
    # add legends
    plot1 <- plot1 + theme(legend.position=legend_position, legend.key = element_blank())
    if (!legend_col) {
        plot1 <- plot1 + guides(colour=FALSE)
    }
    if (!legend_lty) {
        plot1 <- plot1 + guides(lty=FALSE)
    }
    plot1 <- plot1 + scale_colour_manual(labels=legend_col_label, values=col_levels, name=legend_col_title)
    plot1 <- plot1 + scale_linetype_manual(labels=legend_lty_label, values=lty_levels, name=legend_lty_title)
    
    # add points or lines
    if (type=="p") {
        plot1 <- plot1 + geom_point(aes(x, y, colour=col, group=group), shape=pch, size=2*cex, data=df1)
    } else if (type=="l") {
        plot1 <- plot1 + geom_line(aes(x, y, colour=col, linetype=lty, group=group), data=df1)
    }
    
    # add titles
    plot1 <- plot1 + labs(x=xlab, y=ylab)
    plot1 <- plot1 + ggtitle(main)
    
    # log axes
    if (logx) {
        plot1 <- plot1 + scale_x_continuous(trans='log10')
    }
    if (logy) {
        plot1 <- plot1 + scale_y_continuous(trans='log10')
    }
    
    # return plot
    plot1
}

# -----------------------------------
#' points_nice1
#'
#' Add points to plot_nice1.
#'
#' @export

points_nice1 <- function(plot1, x, y, col=1, pch=1, cex=1) {
    
    # if no y variable then first argument becomes y
    if (missing(y)) {
        y <- x
        if (is.list(y)) {
            x <- mapply(function(x){1:length(x)}, y, SIMPLIFY=FALSE)
        } else {
            x <- 1:length(y)
        }
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
        stop("input lists are of different lengths")
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
        df1 <- rbind(df1, data.frame(x=x[[i]], y=y[[i]], col=col[[i]], group=i))
    }
    
    plot1 <- plot1 + geom_point(aes(x, y, colour=col, group=group), shape=pch, size=2*cex, data=df1)
    
    # return plot
    plot1
}
