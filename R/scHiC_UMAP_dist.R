
utils::globalVariables(c(
  "umap.defaults"
))

#' UMAP 2D projection of a dissimilarity matrix
#'
#' This function projects a dissimilarity matrix into the two-dimensional (2D) space through UMAP (Uniform Manifold Approximation and Projection).
#'
#' @importFrom umap umap
#' @importFrom stats kmeans
#' @importFrom cluster pam
#' @importFrom factoextra hcut
#' @importFrom scales alpha
#' @importFrom stats dist as.dist
#' @importFrom graphics points legend
#' @param data  An `n_single`x`n_single` symmetric dissimilarity matrix.
#' @param cell_type  A vector of length `n_single` indicating the true cell (sub)type each cell belongs to.
#' @param seed An integer of random seed value. It is recommended to set random seed to ensure reproducibility. Default value is 1234.
#' @param main The title given to the output plot. Default is empty ("").
#' @param clustering A character indicating whether/which type of clustering is performed on the projected space.
#' \describe{
#'   \item{"KMHW"}{Perform K-means (Hartigan-Wong) clustering. This is the default value. }
#'   \item{"KML"}{Perform K-means (Lloyd) clustering.}
#'   \item{"KMMQ"}{Perform K-means (MacQueen) clustering.}
#'   \item{"pam"}{Perform PAM (Partitioning Around Medoids) clustering.}
#'   \item{"hierarchical"}{Perform hierarchical clustering.}
#'   \item{NULL}{No clustering is performed on the projected space.}
#' }
#' @param ncenters An integer specifying the number of clusters. Default is 3. If `clustering = NULL`, then the value of `ncenters` is ignored.
#' @param show.legend A logical value about whether legend is present. Default value is TRUE.
#' @param alpha_value A numeric value ranges from 0 to 1, with lower values corresponding to more transparent colors for the points. Default value is 1.
#' @return A plot that visualizes each cell as a point in the 2D space. If clustering is performed on the projected space, the shape of a point represents the cluster a single cell belongs to, and the color represents the true cell (sub)type (`cell_type`). If clustering is not performed on the projected space (`clustering=NULL`), both point shape and point color represent the true cell (sub)type.
#' @examples
#' n_single <- 50
#' cell.type <- factor(c(rep(1,n_single), rep(2,n_single), rep(3,n_single)), levels=c("1", "2", "3"))
#' scHiC_UMAP_dist(euclid.diss, cell_type=cell.type, seed=1234, main="Euclid", alpha_value=1)
#' @export
scHiC_UMAP_dist <- function(data, cell_type, seed=1234,
                      main="",
                      clustering="KMHW", ncenters=3,
                      show.legend=TRUE, alpha_value=1) {

  custom.config <- umap.defaults
  custom.config$random_state <- 123
  custom.config$n_neighbors <- 5
  custom.config$input <- "dist"
  set.seed(seed)
  x <- umap(data, config=custom.config)
  xymin <- min(x$layout)
  xymax <- max(x$layout)
  xylim <- c(xymin, xymax)
  xylim <- xylim + ((xylim[2]-xylim[1])*0.1)*c(-0.5, 0.5)

  cex <- 1
  cex.lab <- 1
  cex.main <- 2
  cex.legend <- 1.3

  xlab <- "UMAP1"
  ylab <- "UMAP2"
  pch <- c(1, 2, 4, 5)
  colors <- c("#ff7f00", "#17becf", "#e377c2", "#1D7634")

  layout <- x$layout

  plot(xylim, xylim, asp = 1, type="n", axes=T, frame=T, xlab=xlab, ylab=ylab, cex.axis=cex, cex.lab=cex.lab, main=main, cex.main=cex.main)

  my.xy = data.frame(x=layout[, 1], y=layout[, 2])
  euclid.dist4 <- as.matrix(dist(layout))
  set.seed(seed)
  if(is.null(clustering)){
    points(layout[,1], layout[,2], col=alpha(colors[as.integer(cell_type)], alpha_value), cex=cex, pch=pch[as.integer(cell_type)])
  } else if(clustering=="KMHW"){
    res.km <- kmeans(layout, centers = ncenters, nstart = 10, iter.max=100, algorithm = "Hartigan-Wong")
    my.xy$cluster <- factor(res.km$cluster)
    points(layout[,1], layout[,2], col=alpha(colors[as.integer(cell_type)], alpha_value), cex=cex, pch=pch[as.integer(my.xy$cluster)])
  } else if (clustering=="KML"){
    res.km <- kmeans(layout, centers = ncenters, nstart = 10, iter.max=100, algorithm = "Lloy")
    my.xy$cluster <- factor(res.km$cluster)
    points(layout[,1], layout[,2], col=alpha(colors[as.integer(cell_type)], alpha_value), cex=cex, pch=pch[as.integer(my.xy$cluster)])
  } else if (clustering=="KMMQ"){
    res.km <- kmeans(layout, centers = ncenters, nstart = 10, iter.max=100, algorithm = "MacQueen")
    my.xy$cluster <- factor(res.km$cluster)
    points(layout[,1], layout[,2], col=alpha(colors[as.integer(cell_type)], alpha_value), cex=cex, pch=pch[as.integer(my.xy$cluster)])
  } else if (clustering=="pam"){
    res.pam <- pam(layout, nstart=100, k=ncenters, metric="euclidean")
    my.xy$cluster <- factor(res.pam$clustering)
    points(layout[,1], layout[,2], col=alpha(colors[as.integer(cell_type)], alpha_value), cex=cex, pch=pch[as.integer(my.xy$cluster)])
  } else if (clustering=="hierarchical"){
    rownames(euclid.dist4) <- NULL
    colnames(euclid.dist4) <- NULL
    res.hclust <- hcut(as.dist(euclid.dist4), hc_func="hclust", hc_method="complete", k = ncenters, isdiss = TRUE)
    my.xy$cluster <- factor(res.hclust$cluster)
    points(layout[,1], layout[,2], col=alpha(colors[as.integer(cell_type)], alpha_value), cex=cex, pch=pch[as.integer(my.xy$cluster)])
  } else {
    # default is KMHW
    res.km <- kmeans(layout, centers = ncenters, nstart = 10, iter.max=100, algorithm = "Hartigan-Wong")
    my.xy$cluster <- factor(res.km$cluster)
    points(layout[,1], layout[,2], col=alpha(colors[as.integer(cell_type)], alpha_value), cex=cex, pch=pch[as.integer(my.xy$cluster)])
  }

  cell_type.u <- sort(unique(cell_type))
  legend.pos <- "topright"
  legend.text <- as.character(cell_type.u)

  if (show.legend) {
    if(is.null(clustering)){
      legend(legend.pos, legend=legend.text, inset=0.03, title="cell type",
             col=colors[as.integer(cell_type.u)],
             bty="n", pch=pch[as.integer(cell_type.u)], cex=cex.legend)
    } else {

      cluster.u <- sort(unique(my.xy$cluster))
      legend_labels <- c("cell type",
                         paste0("  ", cell_type.u),
                         "cluster",
                         paste0("  ", cluster.u))
      legend_colors <- c(NA, colors[as.integer(cell_type.u)], NA, rep("black", length(cluster.u)))
      legend_pch    <- c(NA, rep(19, length(cell_type.u)), NA, pch[as.integer(cluster.u)])

      # Add combined legend
      legend("topright",
             legend = legend_labels,
             col = legend_colors,
             pch = legend_pch,
             bty = "n",
             title = NULL,
             text.col = "black",
             y.intersp = 1.2)
    }
  }
}

