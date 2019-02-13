#################
#' Supervised pseudotime
#'
#' @param x Either SingleCellExperiment class containing all cells and genes required, or matrix of log TPM values.
#' @param y Vector of labels, which should have same length as number of columns in sce / x. Factor levels will be taken as the intended order for training.
#' @param y_labels Alternative ordering of the labels in y. All labels must be present in y. May be a subset of the labels in y.
#' @param sel_genes Subset of genes to be used for training.
#' @return psupertime object
#' @export
psupertime <- function(x, y, y_labels=NULL, 
	sel_genes='hvg', gene_list=NULL, scale=TRUE, smooth=TRUE, min_expression=0.01,
	penalization='1se', method='cumulative', score='x_entropy', 
	n_folds=5, test_propn=0.1, lambdas=NULL, max_iters=1e3, seed=1234) {
	# parse params
	params 			= check_params(x, y, y_labels, 
		sel_genes, gene_list, scale, smooth, 
		min_expression, penalization, method, score, 
		n_folds, test_propn, lambdas, max_iters, seed)

	# select genes, do processing of data
	sel_genes 		= select_genes(x, params)
	x_data 			= make_x_data(x, sel_genes, params)

	# make nice data for ordinal regression
	test_idx 		= get_test_idx(y, params)
	
	# get test data
	y_test 			= y[test_idx]
	x_test 			= x_data[test_idx, ]
	y_train 		= y[!test_idx]
	x_train 		= x_data[!test_idx, ]

	# do tests on different folds
	fold_list 		= get_fold_list(y_train, params)
	scores_dt 		= train_on_folds(x_train, y_train, fold_list, params)

	# find best scoring lambda options, train model on these
	mean_scores 	= calc_mean_scores(scores_dt)
	best_dt 		= calc_best_lambdas(mean_scores)
	best_lambdas 	= get_best_lambdas(best_dt, params)
	glmnet_best 	= get_best_fit(x_train, y_train, params)

	# do projections with this model
	proj_dt 		= calc_proj_dt(glmnet_best, x_data, y, best_lambdas)
	beta_dt 		= make_best_beta(glmnet_best, best_lambdas)

	# make output
	psuper_obj 		= make_psuper_obj(glmnet_best, x_data, y, x_test, y_test, proj_dt, beta_dt, best_lambdas, best_dt, scores_dt, params)

	return(psuper_obj)
}

#' check all parameters
#'
#' @return list of validated parameters
#' @keywords internal
check_params <- function(x, y, y_labels, sel_genes, gene_list, scale, smooth, min_expression, 
	penalization, method, score, n_folds, test_propn, lambdas, max_iters, seed) {
	# check input data looks ok
	if ( class(x)=='SingleCellExperiment') {
		if ( !('logcounts' %in% names(SummarizedExperiment::assays(x))) ) {
			stop('if x is a SingleCellExperiment, it must contain the assay logcounts')
		}
		if (ncol(x)!=length(y)) {
			stop('length of y must be same as number of cells (columns) in SingleCellExperiment x')
		}
	} else if ( is.matrix(x) & is.numeric(x) ) {
		if (nrow(x)!=length(y)) {
			stop('length of y must be same as number of rows in matrix x')
		}
		if ( is.null(colnames(x)) ) {
			stop('column names of x must be given, as gene names')
		}
	} else {
		stop('x must be either a SingleCellExperiment or a matrix of counts')
	}
	if (!is.factor(y)) {
		y 	= factor(y)
	}
	if (!is.null(y_labels)) {
		if ( !all(y_labels %in% levels(y)) ) {
			stop('y_labels must be a subset of the labels for y')
		}
	}

	# check selection of genes is valid
	sel_genes_list 	= c('hvg', 'all', 'TF', 'list')
	span_default 	= 0.1
	if (!is.character(sel_genes)) {
		if ( !is.list(sel_genes) || !all(c('hvg_cutoff', 'bio_cutoff') %in% names(sel_genes)) ) {
			err_message 	= paste0('sel_genes must be one of ', paste(sel_genes_list, collapse=', '), ', or a list containing the following named numeric elements: hvg_cutoff, bio_cutoff')
			stop(err_message)
		}
		hvg_cutoff 	= sel_genes$hvg_cutoff
		bio_cutoff 	= sel_genes$bio_cutoff
		if (is.null(sel_genes$span)) {
			span 		= span_default
		} else {
			span 		= sel_genes$span
		}
		sel_genes 	= 'hvg'
	} else {
		if ( !(sel_genes %in% sel_genes_list) ) {
			err_message 	= paste0('invalid value for sel_genes; please use one of ', paste(sel_genes_list, collapse=', '))
			stop(err_message)
		} else if (sel_genes=='list') {
			if (is.null(gene_list) || !is.character(gene_list)) {
				stop("to use 'list' as sel_genes value, you must also give a character vector as gene_list")
			}
		} else if (sel_genes=='hvg') {
			message('using default parameters to identify highly variable genes')
			hvg_cutoff 		= 0.1
			bio_cutoff		= 0.5
			span 			= span_default
		} else {
			hvg_cutoff 		= NULL
			bio_cutoff		= NULL
			span 			= NULL
		}
	}

	# do smoothing, scaling?
	stopifnot(is.logical(smooth))
	stopifnot(is.logical(scale))
	
	# what proportion of cells must express a gene for it to be included?
	if ( !( is.numeric(min_expression) && ( min_expression>=0 & min_expression<=1) ) ) {
		stop('min_expression must be a number greater than 0 and less than 1')
	}

	# how much regularization to use?
	penalty_list 	= c('1se', 'best')
	penalization 	= match.arg(penalization, penalty_list)

	# which statisical model to use for orginal logistic regression?
	method_list 	= c('cumulative', 'forward', 'backward')
	method 			= match.arg(method, method_list)
	
	# which statistical model to use for orginal logistic regression?
	score_list 		= c('x_entropy', 'accuracy')
	score 			= match.arg(score, score_list)

	# check inputs for training
	if ( !(floor(n_folds)==n_folds) || n_folds<=2 ) {
		stop('n_folds must be an integer greater than 2')
	}
	if ( !( is.numeric(test_propn) && ( test_propn>0 & test_propn<1) ) ) {
		stop('test_propn must be a number greater than 0 and less than 1')
	}
	if (is.null(lambdas)) {
		lambdas 	= 10^seq(from=0, to=-4, by=-0.1)
	} else {
		if ( !( is.numeric(lambdas) && all(lambdas == cummin(lambdas)) ) ) {
			stop('lambdas must be a monotonically decreasing vector')
		}
	}
	if ( !(floor(max_iters)==max_iters) || max_iters<=0 ) {
		stop('max_iters must be a positive integer')
	}
	if ( !(floor(seed)==seed) ) {
		stop('seed must be an integer')
	}

	# put into list
	params 	= list(
		sel_genes 		= sel_genes
		,hvg_cutoff 	= hvg_cutoff
		,bio_cutoff 	= bio_cutoff
		,span 			= span
		,gene_list 		= gene_list
		,smooth 		= smooth
		,scale 			= scale
		,min_expression = min_expression
		,penalization 	= penalization
		,method 		= method
		,score 			= score
		,n_folds 		= n_folds
		,test_propn 	= test_propn
		,lambdas 		= lambdas
		,max_iters 		= max_iters
		,seed 			= seed
		)
	return(params)
}

#' Select genes for use in regression
#'
#' @param x SingleCellExperiment class containing all cells and genes required, or matrix of counts
#' @param params List of all parameters specified.
#' @keywords internal
select_genes <- function(x, params) {
	# unpack
	sel_genes 	= params$sel_genes
	if ( sel_genes=='hvg' ) {
		if ( class(x)=='SingleCellExperiment' ) {
			sce 		= x
		} else if ( class(x)=='matrix' ) {
			sce 		= SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = x))
		} else { stop('class of x must be either matrix or SingleCellExperiment') }

		# calculate selected genes
		sel_genes 	= calc_hvg_genes(sce, params, do_plot=FALSE)

	} else if ( sel_genes=='list' ) {
		sel_genes 	= params$gene_list

	} else {
		if ( sel_genes=='all' ) {
			if (class(x)=='SingleCellExperiment') {
				sel_genes 	= rownames(x)
			} else if (class(x)=='matrix') {
				sel_genes 	= colnames(x)
			} else { stop() }

		} else if ( sel_genes=='TF' ) {
			# sel_genes 	= get_tf_list()
			sel_genes 	= tf_list

		} else {
			stop()
		}
	}

	# restrict to genes which are expressed in at least some proportion of cells
	if ( params$min_expression > 0 ) {
		# calc expressed genes
		expressed_genes 	= calc_expressed_genes(x, params)

		# list missing genes
		missing_g 			= setdiff(sel_genes, expressed_genes)
		if (length(missing_g)>0) {
			message('the following genes have insufficient expression and will not be plotted:')
			message(paste(missing_g, collapse=', '))
		}
		sel_genes 			= intersect(expressed_genes, sel_genes)	
	}
	stopifnot( length(sel_genes)>0 )

	return(sel_genes)	
}

#' Calculates list of highly variable genes (according to approach in scran).
#'
#' @param x SingleCellExperiment class or matrix of log counts
#' @param params List of all parameters specified.
#' @import data.table
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bin2d
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 labs
#' @keywords internal
calc_hvg_genes <- function(sce, params, do_plot=FALSE) {
	message('identifying highly variable genes')
	assay_type 		= 'logcounts'
	# var_fit 		= scran::trendVar(sce, assay.type=assay_type, method="loess", use.spikes=FALSE, span=0.1)
	span 			= ifelse(is.null(params$span), 0.1, params$span)
	var_fit 		= scran::trendVar(SummarizedExperiment::assay(sce, assay_type), method="loess", loess.args=list(span=span))
	var_out 		= scran::decomposeVar(sce, var_fit, assay.type=assay_type)

	# plot trends identified
	var_dt 			= as.data.table(var_out)
	var_dt 			= var_dt[, symbol:=rownames(var_out) ]
	data.table::setorder(var_dt, mean)
	if (do_plot) {
		g 	= ggplot(var_dt[ mean>0.5 ]) +
			aes( x=mean, y=total ) +
			geom_bin2d() +
			geom_point( size=0.1 ) +
			geom_line( aes(y=tech)) +
			scale_fill_distiller( palette='RdBu' ) +
			theme_light() +
			labs(
				x 		= "Mean log-expression",
				y 		= "Variance of log-expression"
				)
		print(g)
	}

	# restrict to highly variable genes
	hvg_dt 			= var_dt[ FDR <= params$hvg_cutoff & bio >= params$bio_cutoff ]
	data.table::setorder(hvg_dt, -bio)
	sel_genes 		= hvg_dt$symbol

	return(sel_genes)
}

#' Restrict to genes with minimum proportion of expression defined in params$min_expression
#'
#' @param x SingleCellExperiment class containing all cells and genes required
#' @param params List of all parameters specified.
#' @keywords internal
calc_expressed_genes <- function(x, params) {
	# check whether necessary
	if (params$min_expression==0) {
		return(rownames(x))
	}

	# otherwise calculate it
	if ( class(x)=='SingleCellExperiment' ) {
		x_mat 		= SummarizedExperiment::assay(x, 'logcounts')
	} else if ( class(x)=='matrix' ) {
		x_mat 		= x
	} else { stop('x must be either SingleCellExperiment or matrix') }
	prop_expressed 	= rowMeans( x_mat>0 )
	expressed_genes = names(prop_expressed[ prop_expressed>params$min_expression ])

	return(expressed_genes)
}

#' Get list of transcription factors
#'
#' @return List of all transcription factors specified.
#' @keywords internal
get_tf_list <- function(dirs) {
	tf_path 		= file.path(dirs$data_root, 'fgcz_annotations', 'tf_list.txt')
	tf_full 		= data.table::fread(tf_path)
	tf_list 		= tf_full$symbol
	return(tf_list)
}

#' @keywords internal
get_go_list <- function(dirs) {
	stop('not implemented yet')
	lookup_path 	= file.path(dirs$data_root, 'fgcz_annotations', 'genes_annotation_byGene.txt')
	ensembl_dt 		= fread(lookup_path)
	old_names 		= names(ensembl_dt)
	new_names 		= stringr::str_replace_all(old_names, ' ', '_')
	data.table::setnames(ensembl_dt, old_names, new_names)
	tf_term 		= 'GO:0003700'
	tf_full 		= ensembl_dt[ stringr::str_detect(GO_MF, tf_term) ]
	tf_list 		= tf_full[, list(symbol=gene_name, description)]

	return(tf_list)
}

#' Process input data
#' @param x SingleCellExperiment or matrix of log counts
#' @param sel_genes Selected genes
#' @param params Full list of parameters
#' @return Matrix of dimension # cells by # selected genes
#' @keywords internal
make_x_data <- function(x, sel_genes, params) {
	message('processing data')
	# get matrix
	if ( class(x)=='SingleCellExperiment' ) {
		x_data 		= t(SummarizedExperiment::assay(x, 'logcounts'))
	} else if ( class(x)=='matrix' ) {
		x_data 		= x
	} else {
		stop('x must be either a SingleCellExperiment or a matrix of counts')
	}

	# check if any genes missing, restrict to selected genes
	all_genes 		= colnames(x_data)
	missing_genes 	= setdiff(sel_genes, all_genes)
	n_missing 		= length(missing_genes)
	if ( n_missing>0 ) {
		message('    ', n_missing, ' genes missing:', sep='')
		message('    ', paste(missing_genes[1:min(n_missing,20)], collapse=', '), sep='')
	}
	x_data 		= x_data[, intersect(sel_genes, all_genes)]

	# exclude any genes with zero SD
	message('    checking for zero SD genes')
	col_sd 		= base::apply(x_data, 2, sd)
	# col_sd 		= matrixStats::colSds(x_data)
	sd_0_idx 	= col_sd==0
	if ( sum(sd_0_idx)>0 ) {
		message('\nthe following genes have zero SD and are removed:')
		message(paste(colnames(x_data[, sd_0_idx]), collapse=', '))
		x_data 		= x_data[, !sd_0_idx]
	}

	# do smoothing
	if (params$smooth) {
		message('    denoising data')
		if ( is.null(params$knn) ) {
			knn 	= 10
		} else {
			knn 	= params$knn
		}

		# calculate correlations between all cells
		x_t 			= t(x_data)
		cor_mat 		= cor(x_t)

		# each column is ranked list of nearest neighbours of the column cell
		nhbr_mat 		= apply(-cor_mat, 1, rank, ties.method='random')
		idx_mat 		= nhbr_mat <= knn
		avg_knn_mat 	= sweep(idx_mat, 2, colSums(idx_mat), '/')
		stopifnot( all(colSums(avg_knn_mat)==1) )
		
		# calculate average over all kNNs
		imputed_mat 	= x_t %*% avg_knn_mat
		x_t 			= imputed_mat
		x_data 			= t(x_t)
	}

	# do scaling
	if (params$scale) {
		message('    scaling data')
		x_data 		= apply(x_data, 2, scale)
	}

	# make all gene names nice
	warning('is replacing hyphens in gene names necessary?')
	old_names 			= colnames(x_data)
	new_names 			= stringr::str_replace_all(old_names, '-', '.')
	colnames(x_data) 	= new_names

	return(x_data)
}

#' Get list of cells to keep aside as test set
#'
#' @param y list of y labels
#' @return Indices for test set
#' @keywords internal
get_test_idx <- function(y, params) {
	set.seed(params$seed)
	n_samples 	= length(y)
	test_idx 	= sample(n_samples, round(n_samples*params$test_propn))
	test_idx 	= 1:n_samples %in% test_idx

	return(test_idx)
}

#' @keywords internal
get_fold_list <- function(y_train, params) {
	set.seed(params$seed + 1)
	n_samples 		= length(y_train)
	fold_labels 	= rep_len(1:params$n_folds, n_samples)
	fold_list 		= fold_labels[ sample(n_samples, n_samples) ]

	return(fold_list)
}

#' @keywords internal
train_on_folds <- function(x_train, y_train, fold_list, params) {
	message(sprintf('cross-validation training, %d folds:', params$n_folds))
	# unpack
	n_folds 		= params$n_folds
	lambdas 		= params$lambdas
	max_iters 		= params$max_iters
	method 			= params$method

	# loop
	scores_dt 		= data.table::data.table()
	for (kk in 1:n_folds) {
		message(sprintf('    fold %d', kk))

		# split the folds
		fold_idx 		= fold_list==kk
		y_valid_k 		= y_train[ fold_idx ]
		x_valid_k 		= x_train[ fold_idx, ]
		y_train_k 		= y_train[ !fold_idx ]
		x_train_k 		= x_train[ !fold_idx, ]

		# train model
		glmnet_fit 		= glmnetcr_cumul(x_train_k, y_train_k, 
			method 	= method
			,lambda = lambdas
			,maxit 	= max_iters
			)

		# validate model
		temp_dt 		= calc_scores_for_one_fit(glmnet_fit, x_valid_k, y_valid_k)
		temp_dt[, fold := kk ]
		scores_dt 		= rbind(scores_dt, temp_dt)
	}

	return(scores_dt)
}

#' @keywords internal
glmnetcr_cumul <- function(x, y, method = "cumulative", weights = NULL, offset = NULL, 
    alpha = 1, nlambda = 100, lambda.min.ratio = NULL, lambda = NULL, 
    standardize = TRUE, thresh = 1e-04, exclude, penalty.factor = NULL, 
    maxit = 100) {
    if (length(unique(y)) == 2) 
        stop("Binary response: Use glmnet with family='binomial' parameter")
    method 	<- c("backward", "forward", "cumulative")[charmatch(method, 
    	c("backward", "forward", "cumulative"))]
    n <- nobs <- dim(x)[1]
    p <- m <- nvars <- dim(x)[2]
    k <- length(unique(y))
    x <- as.matrix(x)
    if (is.null(penalty.factor)) 
        penalty.factor <- rep(1, nvars)
    else penalty.factor <- penalty.factor
    if (is.null(lambda.min.ratio)) 
        lambda.min.ratio <- ifelse(nobs < nvars, 0.01, 1e-04)
    if (is.null(weights)) 
        weights <- rep(1, length(y))
    if (method == "backward") {
        restructure <- cr.backward(x = x, y = y, weights = weights)
    }
    if (method == "forward") {
        restructure <- cr.forward(x = x, y = y, weights = weights)
    }
    if (method == "cumulative") {
        restructure <- restructure_cumul(x = x, y = y, weights = weights)
    }
    glmnet.data <- list(x = restructure[, -c(1, 2)], y = restructure[, 
        "y"], weights = restructure[, "weights"])
    object <- glmnet::glmnet(glmnet.data$x, glmnet.data$y, family = "binomial", 
        weights = glmnet.data$weights, offset = offset, alpha = alpha, 
        nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, 
        lambda = lambda, standardize = standardize, thresh = thresh, 
        exclude = exclude, penalty.factor = c(penalty.factor, rep(0, k - 1)), 
        maxit = maxit, type.gaussian = ifelse(nvars < 500, "covariance", "naive"))
    object$x <- x
    object$y <- y
    object$method <- method
    class(object) <- "glmnetcr"
    object
}

#' @keywords internal
restructure_cumul <- function(x, y, weights) {
	yname 		= as.character(substitute(y))
	if (!is.factor(y)) { y = factor(y, exclude = NA) }
	ylevels 	= levels(y)
	kint 		= length(ylevels) - 1
	y 			= as.numeric(y)
	names 		= dimnames(x)[[2]]
	if (length(names) == 0) { names = paste("V", 1:dim(x)[2], sep = "") }
	expand 		= list()
	for (k in 2:(kint + 1)) {
		expand[[k - 1]] 		= cbind(y, weights, x)
		expand[[k - 1]][, 1] 	= ifelse(expand[[k - 1]][, 1] >= k, 1, 0)
		cp = matrix(
			rep(0, dim(expand[[k - 1]])[1] * kint), 
			ncol = kint
			)
		cp[, k - 1] 	= 1
		dimnames(cp)[[2]] 	= paste("cp", 1:kint, sep = "")
		expand[[k - 1]] 	= cbind(expand[[k - 1]], cp)
		dimnames(expand[[k - 1]])[[2]] = c("y", "weights", names, paste("cp", 1:kint, sep = ""))
	}
	newx = expand[[1]]
	for (k in 2:kint) { newx = rbind(newx, expand[[k]]) }
	newx
}

#' @keywords internal
predict_glmnetcr_cumul <- function(object, newx=NULL, newy=NULL, ...) {
	if (is.null(newx)) {
		newx 		= object$x
		y 			= object$y
	} else {
		if (is.null(newy)) {stop('please include newy')}
		y 			= newy
	}
	method 		= object$method
	method 		= c("backward", "forward", "cumulative")[charmatch(method, c("backward", "forward", "cumulative"))]
	y_levels 	= levels(object$y)
	k 			= length(y_levels)
	if (class(newx) == "numeric") {
		newx 		= matrix(newx, ncol = dim(object$x)[2])
	}
	beta.est 	= object$beta

	# split betas into cutpoints, gene coefficients
	beta_genes 	= rownames(beta.est)
	cut_idx 	= stringr::str_detect(beta_genes, '^cp[0-9]+$')
	coeff_cuts 	= beta.est[cut_idx, ]
	coeff_genes = beta.est[!cut_idx, ]

	# restrict to common genes
	data_genes 	= colnames(newx)
	beta_genes 	= rownames(coeff_genes)
	missing_g 	= setdiff(beta_genes, data_genes)
	if (length(missing_g)>0) {
		# decide which to keep
		message("    these genes are missing from the input data and are not used for projecting:")
		message('        ', paste(missing_g, collapse=', '), sep='')
		message("    this may affect the projection\n")
		both_genes 	= intersect(beta_genes, data_genes)

		# tweak inputs accordingly
		newx 		= newx[, both_genes]
		coeff_genes = coeff_genes[both_genes, ]
		beta.est 	= rbind(coeff_genes, coeff_cuts)
	}
	stopifnot( (ncol(newx) + sum(cut_idx)) == nrow(beta.est))
	
	# carry on
	n 			= dim(newx)[1]
	p 			= dim(newx)[2]
	y.mat 		= matrix(0, nrow = n, ncol = k)
	for (i in 1:n) {
		y.mat[i, which(y_levels==y[i])] 	= 1
	}
	n_lambdas 	= dim(beta.est)[2]
	n_nzero 	= apply(beta.est, 2, function(x) sum(x != 0))
	glmnet.BIC 	= glmnet.AIC = numeric()
	pi 			= array(NA, dim=c(n, k, n_lambdas))
	p.class 	= matrix(NA, nrow=n, ncol=n_lambdas)
	LL_mat 		= matrix(0, nrow=n, ncol=n_lambdas)
	LL 			= rep(NA, length=n_lambdas)
	# cycle through each value of lambda
	for (i in 1:n_lambdas) {
		# get beta estimates
		beta 		= beta.est[, i]
		logit 		= matrix(rep(0, n * (k - 1)), ncol=k - 1)
		# project x for each cutpoint
		b_by_x 		= beta[1:p] %*% t(as.matrix(newx))
		for (j in 1:(k - 1)) {
			cp 			= paste("cp", j, sep = "")
			logit[, j] 	= object$a0[i] + beta[names(beta) == cp] + b_by_x
		}
		# do inverse logit
		delta 		= matrix(rep(0, n * (k - 1)), ncol = k - 1)
		for (j in 1:(k - 1)) {
			delta[, j] 	= exp(logit[, j])/(1 + exp(logit[, j]))
		}
		# delta 		= exp(logit)/(1 + exp(logit))
		minus.delta = 1 - delta
		if (method == "backward") {
			for (j in k:2) {
				if (j == k) {
					pi[, j, i] 	= delta[, k - 1]
				} else if (class(minus.delta[, j:(k - 1)]) == "numeric") {
					pi[, j, i] 	= delta[, j - 1] * minus.delta[, j]
				} else if (dim(minus.delta[, j:(k - 1)])[2] >= 2) {
					pi[, j, i] 	= delta[, j - 1] * apply(minus.delta[, j:(k - 1)], 1, prod)
				}
			}
			if (n == 1) {
				pi[, 1, i] 	= 1 - sum(pi[, 2:k, i])
			} else {
				pi[, 1, i] 	= 1 - apply(pi[, 2:k, i], 1, sum)
			}
		}
		if (method == "forward") {
			for (j in 1:(k - 1)) {
				if (j == 1) {
					pi[, j, i] 	= delta[, j]
				} else if (j == 2) {
					pi[, j, i] 	= delta[, j] * minus.delta[, j - 1]
				} else if (j > 2 && j < k) {
					pi[, j, i] 	= delta[, j] * apply(minus.delta[, 1:(j - 1)], 1, prod)
				}
			}
			if (n == 1) {
				pi[, k, i] 	= 1 - sum(pi[, 1:(k - 1), i])
			} else {
				pi[, k, i] 	= 1 - apply(pi[, 1:(k - 1), i], 1, sum)
			}
		}
		if (method == "cumulative") {
			for (j in 1:k) {
				if (j == 1) {
					pi[, j, i] 	= minus.delta[, j]
				} else if (j > 1 && j < k) {
					pi[, j, i] 	= delta[, j - 1] - delta[, j]
				} else if (j == k) {
					pi[, j, i] 	= delta[, j-1]
				}
			}
			if (n == 1) {
				if ( abs(sum(pi[, , i])-1) > 1e-15 ) {
					warning('some probabilities didn\'t sum to one')
					# pi[, , i] 	= pi[, , i]/sum(pi[, , i])
				}
			} else {
				if ( any(abs( rowSums(pi[, , i]) - 1 ) > 1e-15 ) ) {
					warning('some probabilities didn\'t sum to one')
					# pi[, , i] 	= sweep(pi[, , i], 1, rowSums(pi[, , i]), '/')
				}
			}
		}
		# calculate log likelihoods
		if (method == "backward") {
			for (j in 1:(k - 1)) {
				if (class(y.mat[, 1:j]) == "matrix") {
					ylth 		= apply(y.mat[, 1:j], 1, sum)}
				else {
					ylth 		= y.mat[, 1]
				}
				LL_mat[, i] = LL_mat[, i] + log(delta[, j]) * y.mat[, j + 1] + log(1 - delta[, j]) * ylth
			}
		}
		if (method == "forward") {
			for (j in 1:(k - 1)) {
				if (class(y.mat[, j:k]) == "matrix") {
					ygeh 		= apply(y.mat[, j:k], 1, sum)
				} else {
					ygeh 		= y.mat[, k]
				}
				LL_mat[, i] = LL_mat[, i] + log(delta[, j]) * y.mat[, j] + log(1 - delta[, j]) * ygeh
			}
		}
		if (method == "cumulative") {
			for (j in 1:(k - 1)) {
				if (class(y.mat[, j:k]) == "matrix") {
					ygeh 		= apply(y.mat[, j:k], 1, sum)
				} else {
					ygeh 		= y.mat[, k]
				}
				LL_mat[, i] = LL_mat[, i] + log(delta[, j]) * y.mat[, j] + log(1 - delta[, j]) * ygeh
			}
		}
		LL_temp 		= sum(LL_mat[, i])
		glmnet.BIC[i] 	= -2 * LL_temp + n_nzero[i] * log(n)
		glmnet.AIC[i] 	= -2 * LL_temp + 2 * n_nzero[i]
		if (n == 1) {
			p.class[, i] 	= which.max(pi[, , i])
		} else {
			p.class[, i] 	= apply(pi[, , i], 1, which.max)
		}
	}
	LL 					= colSums(LL_mat)
	class 				= matrix(y_levels[p.class], ncol = ncol(p.class))
	names(glmnet.BIC) 	= names(glmnet.AIC) = dimnames(p.class)[[2]] = dimnames(pi)[[3]] = names(object$a0)
	dimnames(pi)[[2]] 	= y_levels

	# output
	list(BIC=glmnet.BIC, AIC=glmnet.AIC, class=class, probs=pi, LL=LL, LL_mat=LL_mat)
}

#' @keywords internal
calc_scores_for_one_fit <- function(glmnet_fit, x_valid, y_valid) {
	# get predictions
	predictions 	= predict_glmnetcr_cumul(glmnet_fit, x_valid, y_valid)
	pred_classes 	= predictions$class
	probs 			= predictions$probs
	lambdas 		= glmnet_fit$lambda

	class_levels 	= levels(y_valid)
	pred_int 		= apply(pred_classes, c(1,2), function(ij) which(ij==class_levels))
	n_lambdas 		= ncol(pred_classes)

	# calculate various accuracy measures
	scores_mat 		= sapply(
		1:n_lambdas, 
		function(jj) calc_multiple_scores(pred_classes[,jj], probs[,,jj], y_valid, class_levels)
		)
	# store results
	scores_wide 	= data.table::data.table(
		lambda 			= lambdas
		,t(scores_mat)
		)
	scores_dt 		= melt(scores_wide, id='lambda', variable.name='score_var', value.name='score_val')

	return(scores_dt)
}

#' @keywords internal
calc_multiple_scores <- function(pred_classes, probs, y_valid, class_levels) {
	# calculate some intermediate variables
	# y_valid_int 	= as.integer(y_valid)
	# pred_int 		= sapply(pred_classes, function(i) which(i==class_levels))
	# bin_mat 		= t(sapply(pred_classes, function(i) i==class_levels))
	bin_mat 		= t(sapply(y_valid, function(i) i==class_levels))

	# calculate optional scores
	accuracy 		= mean(pred_classes==y_valid)
	# log_cohens_k 	= log10(mean( abs(pred_int - y_valid_int) ))
	# log_cohens_k_2 	= log10(mean( (pred_int - y_valid_int)^2 ))
	x_entropy 		= mean(x_entropy_fn(probs, bin_mat))

	scores_vec 		= c(
		accuracy 		= accuracy, 
		# log_cohens_k 	= log_cohens_k, 
		# log_cohens_k_2 	= log_cohens_k_2, 
		x_entropy 		= x_entropy
		)

	return(scores_vec)
}

#' @keywords internal
x_entropy_fn <- function(p_mat, bin_mat) {
	return( -rowSums(bin_mat * log2(p_mat), na.rm=TRUE) )
}

#' @keywords internal
calc_mean_scores <- function(scores_dt) {
	mean_scores 	= scores_dt[, list(mean=mean(score_val), se=sd(score_val)/sqrt(.N)), by=list(lambda,score_var) ]
	dirn_lookup 	= data.table::data.table(
		score_var 	= c('accuracy', 'log_cohens_k', 'log_cohens_k_2', 'x_entropy')
		,dirn 		= c(1, -1, -1, -1)
		)
	mean_scores 	= dirn_lookup[ mean_scores, on='score_var' ]
	mean_scores[, mean_pos := mean * dirn ]

	return(mean_scores)
}

#' @keywords internal
calc_best_lambdas <- function(mean_scores) {
	# 
	best_dt = mean_scores[, 
		list(
			best_lambda = .SD[ which.max(mean_pos) ]$lambda,
			next_lambda = max(.SD[ mean_pos > max(mean_pos) - se ]$lambda)
			),
		by 	= score_var
		]
	# get indices for lambdas
	lambdas 	= rev(sort(unique(mean_scores$lambda)))
	best_dt[, best_idx := which(lambdas==best_lambda), by = score_var ]
	best_dt[, next_idx := which(lambdas==next_lambda), by = score_var ]

	return(best_dt)
}

#' @keywords internal
get_best_lambdas <- function(best_dt, params) {
	# restrict to selected score, extract values
	sel_score 		= params$score
	best_lambdas 	= list(
			best_idx 		= best_dt[ score_var==sel_score ]$best_idx
			,next_idx 		= best_dt[ score_var==sel_score ]$next_idx
			,best_lambda 	= best_dt[ score_var==sel_score ]$best_lambda
			,next_lambda 	= best_dt[ score_var==sel_score ]$next_lambda
		)
	if (params$penalization=='1se') {
		best_lambdas$which_idx 	= best_lambdas$next_idx
	} else if (params$penalization=='best') {
		best_lambdas$which_idx 	= best_lambdas$best_idx
	} else {
		stop('invalid penalization')
	}
	return(best_lambdas)
}

#' @keywords internal
get_best_fit <- function(x_train, y_train, params) {
	message('fitting best model with all training data')
	glmnet_best 	= glmnetcr_cumul(
		x_train, y_train
		,method = params$method
		,lambda = params$lambdas
		,maxit 	= params$max_iters
		)
	return(glmnet_best)
}

#' @keywords internal
calc_proj_dt <- function(glmnet_best, x_data, y_labels, best_lambdas) {
	# unpack
	which_idx 		= best_lambdas$which_idx

	# get best one
	cut_idx 		= stringr::str_detect(rownames(glmnet_best$beta), '^cp[0-9]+$')
	beta_best 		= glmnet_best$beta[!cut_idx, which_idx]

	# remove any missing genes if necessary
	coeff_genes 	= names(beta_best)
	data_genes 		= colnames(x_data)
	missing_genes 	= setdiff(coeff_genes, data_genes)
	n_missing 		= length(missing_genes)
	if ( n_missing>0 ) {
		message("    these genes are missing from the input data and are not used for projecting:")
		message("        ", paste(missing_genes[1:min(n_missing,20)], collapse=', '), sep='')
		message("    this may affect the projection")
		both_genes 		= intersect(coeff_genes, data_genes)
		x_data 			= x_data[, both_genes]
		beta_best 		= beta_best[both_genes]
	}

	# a0_best 	= glmnet_best$a0[[which_idx]]
	# y_proj 		= a0_best + x_data %*% matrix(beta_best, ncol=1)
	psuper 			= x_data %*% matrix(beta_best, ncol=1)
	predictions 	= predict_glmnetcr_cumul(glmnet_best, x_data, y_labels)
	pred_classes 	= factor(predictions$class[, which_idx], levels=levels(glmnet_best$y))

	# put into data.table
	proj_dt 	= data.table::data.table(
		psuper 			= psuper[, 1]
		,label_input 	= y_labels
		,label_psuper 	= pred_classes
		)
	# proj_dt 	= data.table::data.table(
	# 	y_proj 		= y_proj[, 1]
	# 	,y_label 	= y_labels
	# 	,test 		= ifelse(test_idx, 'Test', 'Train')
	# 	)
	# proj_dt[, test := factor(test, levels=c('Train', 'Test'))]

	return(proj_dt)
}

#' Extracts best coefficients.
#' 
#' @return data.table containing learned coefficients for all genes used as input.
#' @keywords internal
make_best_beta <- function(glmnet_best, best_lambdas) {
	cut_idx 	= stringr::str_detect(rownames(glmnet_best$beta), '^cp[0-9]+$')
	best_beta 	= glmnet_best$beta[!cut_idx, best_lambdas$which_idx]
	beta_dt 	= data.table::data.table( beta=best_beta, symbol=names(best_beta) )
	beta_dt[, abs_beta := abs(beta) ]
	data.table::setorder(beta_dt, -abs_beta)
	beta_dt[, symbol:=factor(symbol, levels=beta_dt$symbol)]
		
	return(beta_dt)
}

#' @keywords internal
make_psuper_obj <- function(glmnet_best, x_data, y, x_test, y_test, proj_dt, beta_dt, best_lambdas, best_dt, scores_dt, params) {
	# make cuts_dt
	which_idx 	= best_lambdas$which_idx
	cut_idx 	= stringr::str_detect(rownames(glmnet_best$beta), '^cp[0-9]+$')
	cuts_dt 	= data.table::data.table(
		psuper 			= c(NA, -(glmnet_best$beta[ cut_idx, which_idx ] + glmnet_best$a0[[ which_idx ]]))
		,label_input 	= factor(levels(proj_dt$label_input), levels=levels(proj_dt$label_input))
		)

	# what do we want here?
	# for both best, and 1se
		# best betas
		# projection of original data
		# probabilities for each label
		# predicted labels
	# which of best / 1se is in use
	psuper_obj 	= list(
		params 			= params
		,glmnet_best 	= glmnet_best
		,x_data 		= x_data
		,y 				= y
		,x_test 		= x_test
		,y_test 		= y_test
		,proj_dt 		= proj_dt
		,cuts_dt 		= cuts_dt
		,beta_dt 		= beta_dt
		,best_lambdas 	= best_lambdas
		,best_dt 		= best_dt
		,scores_dt 		= scores_dt
		)

	return(psuper_obj)
}