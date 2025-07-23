
utils::globalVariables(c(
  "x", "y", "type", "cluster"
))

#' t-SNE 2D projection of a dissimilarity matrix
#'
#' This function projects a dissimilarity matrix into the two-dimensional (2D) space through t-SNE (t-Distributed Stochastic Neighbor Embedding).
#'
#' @importFrom Rtsne Rtsne
#' @importFrom stats kmeans
#' @importFrom cluster pam
#' @importFrom factoextra hcut
#' @importFrom scales alpha
#' @importFrom stats dist as.dist
#' @importFrom ggplot2 ggplot geom_point aes xlab ylab ggtitle scale_shape_manual scale_color_manual theme theme_bw coord_fixed xlim ylim guides element_text guide_legend
#' @param data  An `n_single`x`n_single` symmetric dissimilarity matrix.
#' @param cell_type  A vector of length `n_single` indicating the true cell (sub)type each cell belongs to.
#' @param perplexity Perplexity parameter (should not be bigger than `1/3 * (n_single - 1)`). Default value is 10.
#' @param seed An integer of random seed value. It is recommended to set random seed to ensure reproducibility. Default value is 1234.
#' @param title The title given to the output ggplot. Default is NULL.
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
#' @param limit The x- and y-limits of the projected space. Default is c(-55, 55).
#' @param alpha_value A numeric value ranges from 0 to 1, with lower values corresponding to more transparent colors for the points. Default value is 1.
#' @return A ggplot object that visualizes each cell as a point in the 2D space. If clustering is performed on the projected space, the shape of a point represents the cluster a single cell belongs to, and the color represents the true cell (sub)type (`cell_type`). If clustering is not performed on the projected space (`clustering=NULL`), both point shape and point color represent the true cell (sub)type.
#' @examples
#' n_single <- 50
#' cell.type <- factor(c(rep(1,n_single), rep(2,n_single), rep(3,n_single)), levels=c("1", "2", "3"))
#' scHiC_tSNE_dist(euclid.diss, cell_type=cell.type, seed=1234, title="Eucildean distance", alpha_value=1)
#' @export
scHiC_tSNE_dist <- function(data, cell_type, perplexity = 10,
                            seed = 1234, title = NULL, clustering = "KMHW", ncenters = 3,
                            show.legend = TRUE, limit = c(-55, 55), alpha_value = 1){
  mydata <- t(data)
  theta <- 0.0
  dims <- 2  # An integer of the dimensionality of the projected space. Default value is 2.
  check_duplicates <- FALSE  # A logical value checking whether duplicates are present. Default value is FALSE.
  is_distance <- TRUE
  set.seed(seed)
  tsne_dat <- Rtsne(mydata, dims = dims, theta = theta, perplexity = perplexity, check_duplicates = check_duplicates, is_distance = is_distance)
  data = tsne_dat
  my.xy = data.frame(x = data$Y[, 1], y = data$Y[, 2])
  my.xy$type <- c("#ff7f00", "#17becf", "#e377c2", "#1D7634")[as.integer(cell_type)]
  if (!is.null(clustering)) {
    tsne.scale <- my.xy[,1:2]
    euclid.dist4 <- as.matrix(dist(tsne.scale))
    set.seed(seed)

    if(clustering=="KMHW"){
      res.km <- kmeans(tsne.scale, centers = ncenters, nstart = 100, iter.max=1000, algorithm = "Hartigan-Wong")
      # print(table(cell_type, res.km$cluster))
      # cat("Total sum of squares: ", res.km$totss, "\n")
      # cat("Within cluster sum of squares: ", res.km$withinss, "\n", sep=" ")
      # cat("Between cluster sum of squares: ", res.km$betweenss, ", which is ", round(55*res.km$betweenss/res.km$totss, digits=2), "% of total SS", "\n")
      my.xy$cluster <- factor(res.km$cluster)
    } else if (clustering=="KML"){
      res.km <- kmeans(tsne.scale, centers = ncenters, nstart = 100, iter.max=1000, algorithm = "Lloyd")
      my.xy$cluster <- factor(res.km$cluster)
    } else if (clustering=="KMMQ"){
      res.km <- kmeans(tsne.scale, centers = ncenters, nstart = 100, iter.max=1000, algorithm = "MacQueen")
      my.xy$cluster <- factor(res.km$cluster)
    } else if (clustering=="pam"){
      res.pam <- cluster::pam(tsne.scale, nstart=100, k=ncenters, metric="euclidean")
      my.xy$cluster <- factor(res.pam$clustering)
    } else if (clustering=="hierarchical"){
      rownames(euclid.dist4) <- NULL
      colnames(euclid.dist4) <- NULL
      res.hclust <- factoextra::hcut(as.dist(euclid.dist4), hc_func="hclust", hc_method="complete", k = ncenters, isdiss = TRUE)
      my.xy$cluster <- factor(res.hclust$cluster)
    } else {
      # default is KMHW
      res.km <- kmeans(tsne.scale, centers = ncenters, nstart = 100, iter.max=1000, algorithm = "Hartigan-Wong")
      my.xy$cluster <- factor(res.km$cluster)
    }

    return(
      ggplot(my.xy) + geom_point(aes(x = x, y = y, color = type, shape = cluster), alpha=alpha_value, show.legend = show.legend) + ggtitle(title) +
        xlab("tSNE1") + ylab("tSNE2") + coord_fixed() +
        scale_shape_manual(values = c(1, 2, 4, 5, 17, 10, 15, 9, 11)[1:ncenters]) +
        scale_color_manual(name = "cell type", values = c("#ff7f00", "#17becf", "#e377c2", "#1D7634")[1:length(unique(cell_type))], labels = levels(cell_type)) + theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        # theme(legend.position = c(0.9, 0.45),
        #       legend.background = element_rect(fill="white", linetype="soli", colour ="black")) +
        xlim(limit) + ylim(limit) +
        guides(color = guide_legend(override.aes = list(size = 1), order=1),
               shape = guide_legend(override.aes = list(size = 1), order=2))
      )
  } else {
    return(
      ggplot(my.xy) + geom_point(aes(x = x, y = y, color = cell_type, shape = cell_type), alpha=alpha_value, show.legend = show.legend) +
        ggtitle(title) +
        xlab("tSNE1") + ylab("tSNE2") + coord_fixed() +
        scale_shape_manual(name = "cell type", values = c(1, 2, 4, 5, 17, 10, 15, 9, 11)[1:length(unique(cell_type))], labels = levels(cell_type)) +
        scale_color_manual(name = "cell type", values = c("#ff7f00", "#17becf", "#e377c2", "#1D7634")[1:length(unique(cell_type))], labels = levels(cell_type)) + theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        xlim(limit) + ylim(limit)
    )
  }
}

