% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biglm_mapper.R
\name{formula_gen}
\alias{formula_gen}
\title{Helper function to generate formulas}
\usage{
formula_gen(lhs, rhs, icept = T)
}
\arguments{
\item{lhs}{name of left-hand-side variable, as character}

\item{rhs}{name(s) of right-hand-side variables, as character string}

\item{icept}{include an intercept?}
}
\value{
A formula.
}
\description{
Helper function to generate formulas
}
\examples{
formula_gen("y", c("x","z"))
}
