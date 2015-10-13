##' append a heatmap of a matrix to right side of phylogenetic tree
##'
##' 
##' @title gheatmap
##' @param p tree view
##' @param data matrix or data.frame
##' @param offset offset of heatmap to tree
##' @param width total width of heatmap, compare to width of tree
##' @param low color of lowest value
##' @param high color of highest value
##' @param color color of heatmap cell border
##' @param colnames logical, add matrix colnames or not
##' @param colnames_position one of 'bottom' or 'top'
##' @param font.size font size of matrix colnames
##' @return tree view
##' @importFrom reshape2 melt
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 geom_text
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 guides
##' @importFrom ggplot2 guide_legend
##' @export
##' @author Guangchuang Yu
gbarmap <- function(p, data, offset=0, width=1,
                     color="black", colnames=TRUE, colnames_position="bottom", font.size=4) {

    colnames_position %<>% match.arg(c("bottom", "top"))
    
    ## if (is.null(width)) {
    ##     width <- (p$data$x %>% range %>% diff)/30
    ## }

    ## convert width to width of each cell
    width <- width * (p$data$x %>% range %>% diff) / ncol(data)
    
    isTip <- x <- y <- variable <- value <- from <- to <- NULL
 
    df=p$data
    df=df[df$isTip,]
    start <- max(df$x) + offset

    dd <- data[df$label[order(df$y)],]
    dd$y <- sort(df$y)

    dd$lab <- rownames(dd)
    dd <- melt(dd, id=c("lab", "y"))
    
    if (any(dd$value == "")) {
        dd$value[dd$value == ""] <- NA
    }

    V2 <- start + as.numeric(dd$variable) * width
    mapping <- data.frame(from=dd$variable, to=V2)
    mapping <- unique(mapping)

    dd$x <- V2

    p2 <- p + geom_bar(data=dd, aes(x, y, fill=value, stat="identity"), color=color, inherit.aes=FALSE)

    #if (is(dd$value,"numeric")) {
    #    p2 <- p2 + scale_fill_gradient(low=low, high=high, na.value="white")
    #} else {
    #    p2 <- p2 + scale_fill_discrete(na.value="white")
    #}
    
    if (colnames) {
        if (colnames_position == "bottom") {
            y <- 0
        } else {
            y <- max(p$data$y) + 1
        }
        p2 <- p2 + geom_text(data=mapping, aes(x=to, label=from), y=y, size=font.size, inherit.aes = FALSE)
    }

    p2 <- p2 + theme(legend.position="right", legend.title=element_blank())
    p2 <- p2 + guides(fill = guide_legend(override.aes = list(colour = NULL)))
    
    attr(p2, "mapping") <- mapping
    return(p2)
}


##' scale x for tree with heatmap
##'
##' 
##' @title scale_x_ggtree
##' @param p tree view
##' @param breaks breaks for tree
##' @param labels lables for corresponding breaks
##' @return tree view
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_x_date
##' @export
##' @author Guangchuang Yu
scale_x_ggtree <- function(p, breaks=NULL, labels=NULL) {
    mrsd <- attr(p, "mrsd")
    if (!is.null(mrsd) && class(p$data$x) == "Date") {
        x <- Date2decimal(p$data$x)
    } else {
        x <- p$data$x
    }

    if (is.null(breaks)) {
        breaks <- hist(x, breaks=5, plot=FALSE)$breaks
    }
    m <- attr(p, "mapping")

    if (!is.null(mrsd) &&class(m$to) == "Date") {
        to <- Date2decimal(m$to)
    } else {
        to <- m$to
    }
    
    idx <- which(sapply(breaks, function(x) any(x > m$to)))
    if (length(idx)) {
        breaks <- breaks[-idx]
    }
    
    if (is.null(labels)) {
        labels <- breaks
    }
    
    breaks <- c(breaks, to)
    labels <- c(labels, gsub("\\.", "", as.character(m$from)))

    if (!is.null(mrsd) && class(p$data$x) == "Date") {
        p <- p + scale_x_date(breaks=decimal2Date(breaks), labels)
    } else {
        p <- p + scale_x_continuous(breaks=breaks, labels=labels)
    }
    return(p)    
}



##' view tree and associated matrix
##'
##' @title gplot
##' @param p tree view
##' @param data matrix
##' @param low low color
##' @param high high color
##' @param widths widths of sub plot
##' @param color color
##' @param font.size font size
##' @return list of figure
##' @importFrom gridExtra grid.arrange
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_y_continuous
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
##' @examples
##' nwk <- system.file("extdata", "sample.nwk", package="ggtree")
##' tree <- read.tree(nwk)
##' p <- ggtree(tree)
##' d <- matrix(abs(rnorm(52)), ncol=4)
##' rownames(d) <- tree$tip.label
##' colnames(d) <- paste0("G", 1:4)
##' gplot(p, d, low="green", high="red")
gplot <- function(p, data, widths=c(0.5, 0.5), color="blacl", font.size=14) {
    ## p <- p + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0.6))
    p1 <- p + scale_y_continuous(expand = c(0, 0.6))
    ## p1 <- p + theme(panel.margin=unit(0, "null"))
    ## p1 <- p1 + theme(plot.margin = unit(c(1, -1, 1.5, 1), "lines"))
    p2 <- gplot.barmap(p, data, color, font.size)
    grid.arrange(p1, p2, ncol=2, widths=widths)
    invisible(list(p1=p1, p2=p2))
}


##' @importFrom grid unit
##' @importFrom ggplot2 scale_fill_gradient
##' @importFrom ggplot2 scale_fill_discrete
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 guides
##' @importFrom ggplot2 guide_legend
##' @importFrom reshape2 melt
gplot.barmap <- function(p, data, color="black", font.size) {
    isTip <- x <- Var1 <- Var2 <- value <- NULL
    dd=melt(as.matrix(data))
    ## p <- ggtree(tree) ## + theme_tree2()
    ## p <- p + geom_text(aes(x = max(x)*1.1, label=label), subset=.(isTip), hjust=0)
    ## p <- p+geom_segment(aes(x=x*1.02, xend=max(x)*1.08, yend=y), subset=.(isTip), linetype="dashed", size=0.4)
    df=p$data
    df=df[df$isTip,]
    
    dd$Var1 <- factor(dd$Var1, levels = df$label[order(df$y)])
    if (any(dd$value == "")) {
        dd$value[dd$value == ""] <- NA
    }
    
    p2 <- ggplot(dd, aes(Var2, Var1, fill=value))+geom_tile(color=color)
    if (is(dd$value,"numeric")) {
        p2 <- p2 + scale_fill_gradient(low=low, high=high, na.value="white")
    } else {
        p2 <- p2 + scale_fill_discrete(na.value="white")
    }
    
    p2 <- p2+xlab("")+ylab("")
    p2 <- p2+theme_tree2() + theme(axis.ticks.x = element_blank(),
                                   axis.line.x=element_blank())
    ## p1 <- p1 + theme(axis.text.x = element_text(size = font.size))
    p2 <- p2 + theme(axis.ticks.margin = unit(0, "lines")) 
    p2 <- p2 + theme(axis.text.x = element_text(size = font.size))
    ## p2 <- p2 + theme(axis.text.y = element_text(size=font.size))
    
    ## plot.margin   margin around entire plot (unit with the sizes of the top, right, bottom, and left margins) 
    ## units can be given in "lines" or  something more specific like "cm"...

    
    p2 <- p2 + theme(panel.margin=unit(0, "null"))
    p2 <- p2 + theme(plot.margin = unit(c(1, 1, .5, -0.5), "lines"))
    p2 <- p2 + theme(legend.position = "right")
    p2 <- p2 + guides(fill = guide_legend(override.aes = list(colour = NULL)))
    ## p2 <- p2 + labs(fill="")
    
    return(p2)
}


coplot <- function(tree1, tree2, hjust=0) {
    x <- y <- label <- isTip <- tree <- NULL
    dx <- fortify(tree1)
    dx$tree <- "A"

    offset <- max(dx$x) * 1.3
    dy <- fortify(tree2)
    dy <- reverse.treeview.data(dy)
    dy$x <- dy$x + offset + hjust
    dy$tree <- "B"

    dd <- rbind(dx, dy)
    p <- ggplot(dd, aes(x, y)) +
        geom_tree(layout="phylogram", subset=.(tree=="A")) +
            geom_tree(layout="phylogram", subset=.(tree=="B")) +
                theme_tree()
 
    p <- p  + geom_text(aes(label=label),
                        subset=.(isTip & tree == "A"),
                        hjust=-offset/40) +
                            geom_text(aes(label=label),
                                      subset=.(isTip & tree == "B"),
                                      hjust = offset/20)
    return(p)
}







