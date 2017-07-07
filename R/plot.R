
# -----------------------------------
#' plot_nice1
#'
#' Simple plotting function for producing plot similar to base plot(), but in ggplot2. Has the advantage of being able to specify legends outside of plotting area more easily.
#'
#' @export

plot_nice1 <- function(x, y, col=1, xlab=NULL, ylab=NULL, main=NULL, pch=1, cex=1, legend=FALSE, legend_position="right", legend_title="legend_title", legend_label=NULL) {
    
    # if no y variable then first argument is y
    if (missing(y)) {
        y <- x
        x <- 1:length(y)
    }
    
    # set default axis names
    if (is.null(xlab)) {
        xlab <- substitute(x)
    }
    if (is.null(ylab)) {
        ylab <- substitute(y)
    }
    
    # if col is numeric then convert to default colours
    if (is.numeric(col)) {
        col <- (col-1)%%8+1	# recycle colours
        col <- palette()[col]
    }
    
    # create data frame from inputs
    df1 <- data.frame(x=x, y=y, col=col)
    
    # extract levels from colour column (ensures correct ordering)
    col_levels <- levels(df1$col)
    
    # create default legend labels
    legend_labels <- paste0("legend_label", 1:length(col_levels))
    
    # create plot and theme
    plot1 <- ggplot(data=df1)
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
    plot1 <- plot1 + geom_point(aes(x, y, colour=col), shape=pch, size=2*cex)
    
    # add titles
    plot1 <- plot1 + labs(x=xlab, y=ylab)
    plot1 <- plot1 + ggtitle(main)
    
    # add legend
    plot1 <- plot1 + scale_colour_manual(labels=legend_label, values=col_levels, name=legend_title)
    
    # plot
    plot1
}
