require(tidyverse)
require(biglm)


#' Helper function to generate formulas
#'
#' @param lhs name of left-hand-side variable, as character
#' @param rhs name(s) of right-hand-side variables, as character string
#' @param icept include an intercept?
#' @return A formula.
#' @examples
#' formula_gen("y", c("x","z"))
#' @export

formula_gen <- function(lhs, rhs, icept=T) {
	if(icept) as.formula(paste(lhs, paste(rhs, collapse=" + "),sep=" ~ "))
	else as.formula(paste(paste(lhs, paste(rhs, collapse=" + "),sep=" ~ "), "- 1", sep=""))	
}

### functions to deal with varying factor levels in different data sets

# get list, one element per variable, of union of levels of a set of variables in 
# a set of files
unite_levels <- function(file_list, factor_vars) {
	file_list %>% 
		map(~ . %>% readRDS %>% extract_levels(factor_vars)) %>%
		reduce(list_union)
}

list_union <- function(l1,l2) map2(l1,l2,union)

extract_levels <- function(data, vars) {
	map(vars, ~ data[,.] %>% unique)
}


#' Standardize levels of factor variables.
#'
#' @param data Dataset to standardize
#' @param levels_union Named list. Names should be variable names, values should be character vectors.
#' @return `data` with variables in `levels_union` converted to factors with levels specified in `levels_union` 
#' @examples
#' update_levels(d, list(var1=c("a","b","c"), var2=c("x","y")))
#' @export
update_levels <- function(data, levels_union){
	data_factors <- data %>% 
		ungroup %>% 
		select(one_of(names(levels_union))) %>%
		map2_dfc(levels_union, ~ factor(.x, levels=.y))

	data %>% 
		select(-one_of(names(levels_union))) %>% 
		bind_cols(data_factors)
}

### functions to estimate biglm over series of datasets stored on disk
chunk_biglm <- function(model, file, levels_union = list()) {
	chunk <- readRDS(file)

	if (length(levels_union) > 0) {
		# make levels of factor variables common
		chunk <- chunk %>% update_levels(levels_union)
	}

	update(model, moredata=chunk)
}

#' Map-reduce with `biglm`
#'
#' @param file_list Character vector of data file names. Must be in `.rds` format.
#' @param formula Model formula.
#' @param levels_union Named list. Names should be variable names, values should be character vectors.
#' @param ... Additional arguments to pass to `biglm`
#' @return A fitted `biglm` object.
#' @examples
#' biglm_mapper(data_files, y ~ x + z)
#' @export

biglm_mapper <- function(file_list, formula, levels_union = list(), ...) {
	# init model on first chunk
	chunk0 <- readRDS(file_list[1])
	if (length(levels_union) > 0) {
		# make levels of factor variables common
		chunk0 <- chunk0 %>% update_levels(levels_union)
	}
	m0 <- biglm(formula, chunk0, ...)

	file_list[2:length(file_list)] %>%
		reduce(chunk_biglm, levels_union=levels_union, .init=m0)
}

#### functions to compute cluster se's
make_group_cov <- function(group_data, b) {
	y <- group_data$frac
	x <- model.matrix(b, group_data)
	resid <- y[as.numeric(rownames(x))] - x %*% coef(b)

	list(g_cov=tcrossprod(crossprod(x,resid)), e_hat = mean(resid^2), n_g=length(resid))
}

combine_groups <- function(g1,g2) {
	n_gtot <- g1$n_g + g2$n_g
	list(g_cov = g1$g_cov + g2$g_cov, e_hat = (g1$n_g*g1$e_hat + g2$n_g*g2$e_hat) / n_gtot, n_g = n_gtot)
}

default_g_cov <- function(m) list(g_cov=matrix(0, nrow=length(coef(m)), ncol=length(coef(m))), e_hat=0, n_g=0)

chunk_clusterse <- function(file, fitted_model, clustervar, levels_union=list()) {
	chunk <- readRDS(file)

	if (length(levels_union) > 0) {
		# make levels of factor variables common
		chunk <- chunk %>% update_levels(levels_union)
	}

	chunk %>% 
		group_by_at(.vars=clustervar) %>% 
		nest %>%
		pull(data) %>%
		map(possibly(make_group_cov, 
			         otherwise=default_g_cov(fitted_model)),
		 	b=fitted_model) %>%
		reduce(combine_groups)
}

#' Map-reduce clustered standard errors with `biglm`
#'
#' @param file_list Character vector of data file names. Must be in `.rds` format.
#' @param fitted_model Fitted `biglm` object.
#' @param ... Additional arguments to pass to `chunk_clusterse`
#' * `clustervar`: character name of variable to cluster on
#' * `levels_union`: named list with factor variables to standardize
#' @return A list with elements
#' * `xxinv`: (X'X)^-1
#' * `cluster_vcov`: cluster-robust variance-covariance matrix estimate
#' * `g_cov`: inner matrix in cluster-robust vcov matrix
#' * `e_hat`: mean of squared residuals
#' @examples
#' cluster_se(data_files, m, clustervar="group")
#' @details Computation DOES NOT do small sample DOF adjustment (I assume if you're using `biglm`, the data is large).
#' Also importantly, computation assumes no group overlap across files (all observations from each group are in one file).
#' @export

cluster_se <- function(file_list, fitted_model, ...) {
	
	group_cov_parts <- file_list %>% 
		map(possibly(chunk_clusterse, otherwise=default_g_cov(fitted_model)), fitted_model=fitted_model, ...) %>%
		reduce(combine_groups)

	xxinv <- vcov(m) / group_cov_parts$e_hat
	cluster_vcov <- xxinv %*% group_cov_parts$g_cov %*% xxinv
	
	group_cov_parts$xxinv <- xxinv
	group_cov_parts$cluster_vcov <- cluster_vcov

	group_cov_parts

}