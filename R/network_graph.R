#' Apply colors to a correlation matrix using ColorBrewer palettes
#'
#' This function is used internally to (a) translate correlation coefficients
#' to colors on a (diverging) ColorBrewer palette and (b) return both the
#' color and the corresponding correlation coefficient.
#'
#' @param R A correlation matrix
#' @param name The name of a ColorBrewer palette. This should be a diverging
#' palette. Default: RdYlBu
#'
#' @return A data frame with columns `color` and `label`. `label` is the
#' correlation coefficient.
apply_colors <- function(R, name = "RdYlBu") {
  ## get database of information about ColorBrewer palettes
  ##
  tmp <- RColorBrewer::brewer.pal.info
  pal <- which(rownames(tmp) == name)
  ## get maximum number of colors in selected palette
  ##
  n_breaks <- tmp$maxcolors[pal]
  my_palette <- RColorBrewer::brewer.pal(n_breaks, name)
  ## This gets even spacing on [0, 1]
  ##
  breaks <- seq(from = 0, to = n_breaks)/n_breaks
  ## Shift it to [-1, 1]
  ##
  breaks <- 2*breaks - 1
  color <- R
  for (i in 1:nrow(R)) {
    for (j in 1:ncol(R)) {
      idx <- which(breaks > R[i, j])
      color[i, j] <- my_palette[idx[1]]
    }
  }
  color_vector <- numeric(0)
  r_vector <- numeric(0)
  n_traits <- nrow(R)
  for (i in 1:(n_traits - 1)) {
    for (j in (i+1):n_traits) {
      color_vector <- c(color_vector, color[i, j])
      r_vector <- c(r_vector, R[i, j])
    }
  }
  dat <- data.frame(color = color_vector,
                    label = r_vector)
  return(dat)
}

#' Visualize a the posterior distribution of a correlation matrix from
#' the output of `brm()`
#'
#' This function plots the posterior correlation matrix from a `brm()` analysis
#' as a series of arcs with the color of each connection representing the
#' posterior mean of the correlation coefficient and the width representing
#' the one-sided probability that the coefficient is different from 0.
#'
#' @param output_model An object of class `brmsfit`
#' @param labels Labels for the nodes. If `NULL` then the names used in the
#' `brms` model are also used here
#' @param style Style of visualization, "arc" (default) or "line"
#'
#' @return A gggraph object that can be further modified
#'
#' @export
plot_posterior_correlation <- function(output_model, labels = NULL,
                                       style = "arc")
{
  Omega <- get_Omega(output_model)
  model_df <- as.data.frame(output_model)
  tmp <- model_df %>%
    dplyr::select(starts_with("b_"))
  traits <- gsub("b_(.*)_Intercept", "\\1", colnames(tmp))

  n_traits <- length(traits)
  R <- matrix(nrow = n_traits, ncol = n_traits)
  if (is.null(labels)) {
    colnames(R) <- traits
    rownames(R) <- traits
  } else {
    colnames(R) <- labels
    rownames(R) <- labels
  }
  weight <- R
  for (i in 1:n_traits) {
    for (j in 1:n_traits) {
      R[i, j] <- mean(Omega[, i, j])
      tmp <- Omega[, i, j]
      if (R[i, j] < 0) {
        count <- length(tmp[tmp < 0])
      } else {
        count <- length(tmp[tmp > 0])
      }
      weight[i, j] <- count/length(tmp)
    }
  }
  weight_vector <- numeric(0)
  for (i in 1:(n_traits - 1)) {
    for (j in (i+1):n_traits) {
      weight_vector <- c(weight_vector, 4*exp(-4*(1 - weight[i, j])))
    }
  }

  graph <- igraph::graph_from_adjacency_matrix(R,
                                               weighted = TRUE,
                                               mode = "undirected",
                                               diag = FALSE)

  igraph::E(graph)$Correlation <- apply_colors(R)$color
  igraph::E(graph)$P_value <- weight_vector

  if (style == "arc") {
    p <- ggraph::ggraph(graph, layout = "linear", circular = TRUE) +
      ggraph::geom_edge_arc(ggplot2::aes(edge_width = P_value,
                                         edge_color = Correlation)) +
      ggraph::geom_node_point(size = 5) +
      ggraph::geom_node_label(ggplot2::aes(label = name),
                              size = 5) +
      ggraph::theme_graph()
  } else if (style == "line") {
    p <- ggraph::ggraph(graph, layout = "linear", circular = TRUE) +
      ggraph::geom_edge_link(ggplot2::aes(edge_width = P_value,
                                          edge_color = Correlation)) +
      ggraph::geom_node_point(size = 5) +
      ggraph::geom_node_label(ggplot2::aes(label = name),
                              size = 5) +
      ggraph::theme_graph()
  } else {
    stop(style, " is not a recognized style!")
  }
  return(p)
}

