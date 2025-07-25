% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_graph.R
\name{plot_posterior_correlation}
\alias{plot_posterior_correlation}
\title{Visualize a the posterior distribution of a correlation matrix from
the output of \code{brm()}}
\usage{
plot_posterior_correlation(
  output_model,
  labels = NULL,
  style = "arc",
  palette = "RdBu"
)
}
\arguments{
\item{output_model}{An object of class \code{brmsfit}}

\item{labels}{Labels for the nodes. If \code{NULL} then the names used in the
\code{brms} model are also used here}

\item{style}{Style of visualization, "arc" (default) or "line"}

\item{palette}{RColorBrewer palette to use. See \code{apply_colors()} for
details. Default: "RdBu"}
}
\value{
A list with two elements

\code{p} A gggraph object that can be further modified

\code{R} The posterior mean correlation matrixtmp
}
\description{
This function plots the posterior correlation matrix from a \code{brm()} analysis
as a series of arcs with the color of each connection representing the
posterior mean of the correlation coefficient and the width representing
the one-sided probability that the coefficient is different from 0.
}
