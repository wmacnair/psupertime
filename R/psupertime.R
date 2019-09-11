#################
#' Supervised pseudotime
#'
#' @param x Either SingleCellExperiment object containing a matrix of genes * cells required, or a matrix of log TPM values (also genes * cells).
#' @param y Vector of labels, which should have same length as number of columns in sce / x. Factor levels will be taken as the intended order for training.
#' @param y_labels Alternative ordering of the labels in y. All labels must be present in y.
#' @param assay_type If a SingleCellExperiment object is used as input, specifies which assay is to be used.
#' @param sel_genes Method to be used to select interesting genes to be used in psupertime. Must be a string, with permitted values 'hvg', 'all', 'tf_mouse', 'tf_human' and 'list', corresponding to: highly variable genes, all genes, transcription factors in mouse, transcription factors in human, and a user-selected list. If sel_genes='list', then the parameter gene_list must also be specified as input, containing the user-specified list of genes. sel_genes may alternatively be a list, itself, specifying the parameters to be used for selecting highly variable genes via scran, with names 'hvg_cutoff', 'bio_cutoff' (optionally also 'span'). 
#' @param gene_list If sel_genes is specified as 'list', gene_list specifies the list of user-specified genes.
#' @param scale Should the log expression data for each gene be scaled to have mean zero and SD 1? Having the same scale ensures that L1-penalization functions properly; typically you would only set this to FALSE if you have already done your own scaling.
#' @param smooth Should the data be smoothed over neighbours? This is done to denoise the data; if you already done your own denoising, set this to FALSE.
#' @param min_expression Cutoff for excluding genes based on non-zero expression in only a small proportion of cells; default is 1\% of cells. 
#' @param penalization Method of selecting level of L1-penalization. 'best' uses the value of lambda giving the best cross-validation accuracy; '1se' corresponds to largest value of lambda within 1 standard error of the best. This increases sparsity with minimal increased error (and is the default). 
#' @param method Statistical model used for ordinal logistic regression, one of 'proportional', 'forward' and 'backward', corresponding to cumulative proportional odds, forward continuation ratio and backward continuation ratio. 
#' @param score Cross-validated accuracy to be used to select model. May take values 'x_entropy' (default), or 'class_error', corresponding to cross-entropy and classification error respectively. Cross-entropy is a smooth measure, while classification error is based on discrete labels and tends to be a bit 'lumpy'.
#' @param n_folds Number of folds to use for cross-validation; default is 5.
#' @param test_propn Proportion of data to hold out for testing, separate to the cross-validation; default is 0.1 (10\%).
#' @param lambdas User-specified sequence of lambda values. Should be in decreasing order. 
#' @param max_iters Maximum number of iterations to run in glmnet.
#' @param seed Random seed for specifying cross-validation folds and test data
#' @return psupertime object
#' @export
psupertime <- function(x, y, y_labels=NULL, assay_type='logcounts',
	sel_genes='hvg', gene_list=NULL, scale=TRUE, smooth=TRUE, min_expression=0.01,
	penalization='1se', method='proportional', score='xentropy', 
	n_folds=5, test_propn=0.1, lambdas=NULL, max_iters=1e3, seed=1234) {
	# parse params
	params 			= check_params(x, y, y_labels, 
		assay_type, sel_genes, gene_list, scale, smooth, 
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
	scores_train 	= train_on_folds(x_train, y_train, fold_list, params)

	# find best scoring lambda options, train model on these
	mean_train 		= calc_mean_train(scores_train)
	best_dt 		= calc_best_lambdas(mean_train)
	best_lambdas 	= get_best_lambdas(best_dt, params)
	glmnet_best 	= get_best_fit(x_train, y_train, params)
	scores_dt 		= make_scores_dt(glmnet_best, x_test, y_test, scores_train)

	# do projections with this model
	proj_dt 		= calc_proj_dt(glmnet_best, x_data, y, best_lambdas)
	beta_dt 		= make_best_beta(glmnet_best, best_lambdas)

	# make output
	psuper_obj 		= make_psuper_obj(glmnet_best, x_data, y, x_test, y_test, test_idx, fold_list, proj_dt, beta_dt, best_lambdas, best_dt, scores_dt, params)

	return(psuper_obj)
}

#' check all parameters
#'
#' @importFrom SummarizedExperiment assays
#' @return list of validated parameters
#' @keywords internal
check_params <- function(x, y, y_labels, assay_type, sel_genes, gene_list, scale, smooth, min_expression, 
	penalization, method, score, n_folds, test_propn, lambdas, max_iters, seed) {
	# check input data looks ok
	if (!(class(x) %in% c('SingleCellExperiment', 'matrix', 'dgCMatrix', 'dgRMatrix'))) {
		stop('x must be either a SingleCellExperiment or a matrix of log counts')
	}
	if ( class(x)=='SingleCellExperiment') {
		if ( !(assay_type %in% names(assays(x))) ) {
			stop(paste0('SingleCellExperiment x does not contain the specified assay, ', assay_type))
		}
	} else {
		if ( is.null(rownames(x)) ) {
			stop('row names of x must be given, as gene names')
		}
	}
	if (ncol(x)!=length(y)) {
		stop('length of y must be same as number of cells (columns) in SingleCellExperiment x')
	}
	n_genes 		= nrow(x)

	if (!is.factor(y)) {
		y 	= factor(y)
		message('converting y to a factor. label ordering used for training psupertime is:')
		message(paste(levels(y), collapse=', '))
	}
	if (!is.null(y_labels)) {
		stop('y_labels not yet implemented')
		if ( !all(y_labels %in% levels(y)) ) {
			stop('y_labels must be a subset of the labels for y')
		}
	}

	# check selection of genes is valid
	sel_genes_list 	= c('hvg', 'all', 'tf_mouse', 'tf_human', 'list')
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
			hvg_cutoff 		= NULL
			bio_cutoff		= NULL
			span 			= NULL
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

	# give warning about scaling
	if (scale==FALSE) {
		warning("'scale' is set to FALSE. If you are using pre-scaled data, ignore this warning.\nIf not, be aware that for LASSO regression to work properly, the input variables should be on the same scale.")
		if (min_expression > 0) {
			warning("If you are using pre-scaled data, then psupertime cannot tell which are zero values; consider setting 'min_expression' to 0.")
		}
	}

	# what proportion of cells must express a gene for it to be included?
	if ( !( is.numeric(min_expression) && ( min_expression>=0 & min_expression<=1) ) ) {
		stop('min_expression must be a number greater than 0 and less than 1')
	}

	# how much regularization to use?
	penalty_list 	= c('1se', 'best')
	penalization 	= match.arg(penalization, penalty_list)

	# which statisical model to use for orginal logistic regression?
	method_list 	= c('proportional', 'forward', 'backward')
	method 			= match.arg(method, method_list)
	
	# which statistical model to use for orginal logistic regression?
	score_list 		= c('xentropy', 'class_error')
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
		n_genes 		= n_genes
		,assay_type 	= assay_type
		,sel_genes 		= sel_genes
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
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @param x SingleCellExperiment class containing all cells and genes required, or matrix of counts
#' @param params List of all parameters specified.
#' @keywords internal
select_genes <- function(x, params) {
	# unpack
	sel_genes 	= params$sel_genes
	if ( sel_genes=='hvg' ) {
		if ( class(x)=='SingleCellExperiment' ) {
			sce 		= x
		} else if ( class(x) %in% c('matrix', 'dgCMatrix', 'dgRMatrix') ) {
			sce 		= SingleCellExperiment(assays = list(logcounts = x))
		} else { stop('class of x must be either matrix or SingleCellExperiment') }

		# calculate selected genes
		sel_genes 	= calc_hvg_genes(sce, params, do_plot=FALSE)

	} else if ( sel_genes=='list' ) {
		sel_genes 	= params$gene_list

	} else {
		if ( sel_genes=='all' ) {
			sel_genes 	= rownames(x)

		} else if ( sel_genes=='tf_mouse' ) {
			sel_genes 	= tf_mouse

		} else if ( sel_genes=='tf_human' ) {
			sel_genes 	= tf_human

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
		n_missing 			= length(missing_g)
		if (n_missing>0) {
			message(sprintf('%d genes have insufficient expression and will not be used as input to psupertime', n_missing))
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
#' @importFrom data.table setorder
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bin2d
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 labs
#' @importFrom scran trendVar
#' @importFrom scran decomposeVar
#' @importFrom SummarizedExperiment assay
#' @keywords internal
calc_hvg_genes <- function(sce, params, do_plot=FALSE) {
	message('identifying highly variable genes')
	assay_type 		= 'logcounts'
	if (!(assay_type %in% names(assays(sce)))) {
		stop('to calculate highly variable genes (HVGs) with scran, x must contain the assay "logcounts"')
	}
	span 			= ifelse(is.null(params$span), 0.1, params$span)
	var_fit 		= trendVar(assay(sce, assay_type), method="loess", loess.args=list(span=span))
	var_out 		= decomposeVar(sce, var_fit, assay.type=assay_type)

	# plot trends identified
	var_dt 			= as.data.table(var_out)
	var_dt 			= var_dt[, symbol:=rownames(var_out) ]
	setorder(var_dt, mean)
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
	setorder(hvg_dt, -bio)
	sel_genes 		= hvg_dt$symbol

	return(sel_genes)
}

#' Restrict to genes with minimum proportion of expression defined in params$min_expression
#'
#' @param x SingleCellExperiment class containing all cells and genes required
#' @param params List of all parameters specified.
#' @importFrom Matrix rowMeans
#' @importFrom SummarizedExperiment assay
#' @keywords internal
calc_expressed_genes <- function(x, params) {
	# check whether necessary
	if (params$min_expression==0) {
		return(rownames(x))
	}

	# otherwise calculate it
	if ( class(x)=='SingleCellExperiment' ) {
		x_mat 		= assay(x, 'logcounts')
	} else if ( class(x) %in% c('matrix', 'dgCMatrix', 'dgCMatrix') ) {
		x_mat 		= x
	} else { stop('x must be either SingleCellExperiment or (possibly sparse) matrix') }
	prop_expressed 	= rowMeans( x_mat>0 )
	expressed_genes = names(prop_expressed[ prop_expressed>params$min_expression ])

	return(expressed_genes)
}

#' Get list of transcription factors
#'
#' @importFrom data.table fread
#' @return List of all transcription factors specified.
#' @keywords internal
get_tf_list <- function(dirs) {
	tf_path 		= file.path(dirs$data_root, 'fgcz_annotations', 'tf_list.txt')
	tf_full 		= fread(tf_path)
	tf_list 		= tf_full$symbol
	return(tf_list)
}

#' @importFrom data.table setnames
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_detect
#' @keywords internal
get_go_list <- function(dirs) {
	stop('not implemented yet')
	lookup_path 	= file.path(dirs$data_root, 'fgcz_annotations', 'genes_annotation_byGene.txt')
	ensembl_dt 		= fread(lookup_path)
	old_names 		= names(ensembl_dt)
	new_names 		= str_replace_all(old_names, ' ', '_')
	setnames(ensembl_dt, old_names, new_names)
	tf_term 		= 'GO:0003700'
	tf_full 		= ensembl_dt[ str_detect(GO_MF, tf_term) ]
	tf_list 		= tf_full[, list(symbol=gene_name, description)]

	return(tf_list)
}

#' Process input data
#' @param x SingleCellExperiment or matrix of log counts
#' @param sel_genes Selected genes
#' @param params Full list of parameters
#' @importFrom Matrix t
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace_all
#' @importFrom SummarizedExperiment assay
#' @return Matrix of dimension # cells by # selected genes
#' @keywords internal
make_x_data <- function(x, sel_genes, params) {
	message('processing data')
	# get matrix
	if ( class(x)=='SingleCellExperiment' ) {
		x_data 		= assay(x, params$assay_type)
	} else if (class(x) %in% c('matrix', 'dgCMatrix', 'dgRMatrix')) {
		x_data 		= x
	} else {
		stop('x must be either a SingleCellExperiment or a matrix of counts')
	}

	# transpose
	x_data 		= t(x_data)

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
	col_sd 		= apply(x_data, 2, sd)
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
		cor_mat 		= cor(as.matrix(x_t))

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
	old_names 			= colnames(x_data)
	hyphen_idx 			= str_detect(old_names, '-')
	if (any(hyphen_idx)) {
		message('    hyphens detected in the following gene names:')
		message('        ', appendLF=FALSE)
		message(paste(old_names[hyphen_idx], collapse=', '))
		message('    these have been replaced with .s')
		new_names 			= str_replace_all(old_names, '-', '.')
		colnames(x_data) 	= new_names
	}

	message(sprintf('    processed data is %d cells * %d genes', nrow(x_data), ncol(x_data)))

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

#' @importFrom data.table data.table
#' @keywords internal
train_on_folds <- function(x_train, y_train, fold_list, params) {
	message(sprintf('cross-validation training, %d folds:', params$n_folds))
	# unpack
	n_folds 		= params$n_folds
	lambdas 		= params$lambdas
	max_iters 		= params$max_iters
	method 			= params$method

	# loop
	scores_dt 		= data.table()
	for (kk in 1:n_folds) {
		message(sprintf('    fold %d', kk))

		# split the folds
		fold_idx 		= fold_list==kk
		y_valid_k 		= y_train[ fold_idx ]
		x_valid_k 		= x_train[ fold_idx, ]
		y_train_k 		= y_train[ !fold_idx ]
		x_train_k 		= x_train[ !fold_idx, ]

		# train model
		glmnet_fit 		= glmnetcr_propn(x_train_k, y_train_k, 
			method 	= method
			,lambda = lambdas
			,maxit 	= max_iters
			)

		# validate model
		temp_dt 		= calc_scores_for_one_fit(glmnet_fit, x_valid_k, y_valid_k)
		temp_dt[, fold := kk ]
		scores_dt 		= rbind(scores_dt, temp_dt)
	}

	# add label
	scores_dt[, data := 'train' ]

	return(scores_dt)
}

#' @importFrom glmnet glmnet
#' @keywords internal
glmnetcr_propn <- function(x, y, method = "proportional", weights = NULL, offset = NULL, 
    alpha = 1, nlambda = 100, lambda.min.ratio = NULL, lambda = NULL, 
    standardize = TRUE, thresh = 1e-04, exclude, penalty.factor = NULL, 
    maxit = 100) {
    if (length(unique(y)) == 2) 
        stop("Binary response: Use glmnet with family='binomial' parameter")
    method 	<- c("backward", "forward", "proportional")[charmatch(method, 
    	c("backward", "forward", "proportional"))]
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
    if (method == "proportional") {
        restructure <- restructure_propodds(x = x, y = y, weights = weights)
    }
    glmnet.data <- list(x = restructure[, -c(1, 2)], y = restructure[, 
        "y"], weights = restructure[, "weights"])
    object <- glmnet(glmnet.data$x, glmnet.data$y, family = "binomial", 
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
restructure_propodds <- function(x, y, weights) {
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

#' @importFrom stringr str_detect
#' @keywords internal
predict_glmnetcr_propodds <- function(object, newx=NULL, newy=NULL, ...) {
	if (is.null(newx)) {
		newx 		= object$x
		y 			= object$y
	} else {
		if (is.null(newy)) {stop('please include newy')}
		y 			= newy
	}
	method 		= object$method
	method 		= c("backward", "forward", "proportional")[charmatch(method, c("backward", "forward", "proportional"))]
	y_levels 	= levels(object$y)
	k 			= length(y_levels)
	if (class(newx) == "numeric") {
		newx 		= matrix(newx, ncol = dim(object$x)[2])
	}
	beta.est 	= object$beta

	# split betas into cutpoints, gene coefficients
	beta_genes 	= rownames(beta.est)
	cut_idx 	= str_detect(beta_genes, '^cp[0-9]+$')
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
			exp_val 	= exp(logit[, j])
			delta[, j] 	= exp_val/(1 + exp_val)
			# check no infinite values
			if ( any(exp_val==Inf) ) {
				delta[exp_val==Inf, j] 	= 1 - 1e-16
			}
		}
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
		if (method == "proportional") {
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
		if (method == "proportional") {
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

#' @importFrom data.table data.table
#' @keywords internal
calc_scores_for_one_fit <- function(glmnet_fit, x_valid, y_valid) {
	# get predictions
	predictions 	= predict_glmnetcr_propodds(glmnet_fit, x_valid, y_valid)
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
	scores_wide 	= data.table(
		lambda 			= lambdas
		,t(scores_mat)
		)
	scores_dt 		= melt(scores_wide, id='lambda', variable.name='score_var', value.name='score_val')

	# put scores in nice order
	scores_dt[, score_var := factor(score_var, levels=c('xentropy', 'class_error')) ]

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
	class_error 	= mean(pred_classes!=y_valid)
	xentropy 		= mean(xentropy_fn(probs, bin_mat))

	scores_vec 		= c(
		class_error 	= class_error, 
		xentropy 		= xentropy
		)

	return(scores_vec)
}

#' @keywords internal
xentropy_fn <- function(p_mat, bin_mat) {
	# calculate standard xentropy values
	xentropy 	= -rowSums(bin_mat * log2(p_mat), na.rm=TRUE)

	# deal with any -Inf values
	inf_idx 	= xentropy==Inf
	if ( sum(inf_idx) ) {
		xentropy[inf_idx] 	= -log2(1e-16)
	}

	return( xentropy )
}

#' @keywords internal
calc_mean_train <- function(scores_train) {
	# calculate mean scores
	mean_train 			= scores_train[, 
		list(
			mean 	= mean(score_val), 
			se 		= sd(score_val)/sqrt(.N)
			), 
		by = list(lambda, score_var)
		]

	# set up levels
	mean_train$score_var 	= factor(mean_train$score_var, levels=levels(scores_train$score_var))

	return(mean_train)
}

#' @keywords internal
calc_best_lambdas <- function(mean_train) {
	# 
	best_dt = mean_train[, 
		list(
			best_lambda = .SD[ which.min(mean) ]$lambda,
			next_lambda = max(.SD[ mean < min(mean) + se ]$lambda)
			),
		by 	= score_var
		]
	# get indices for lambdas
	lambdas 	= rev(sort(unique(mean_train$lambda)))
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
		best_lambdas$which_idx 		= best_lambdas$next_idx
		best_lambdas$which_lambda 	= best_lambdas$next_lambda
	} else if (params$penalization=='best') {
		best_lambdas$which_idx 		= best_lambdas$best_idx
		best_lambdas$which_lambda 	= best_lambdas$best_lambda
	} else {
		stop('invalid penalization')
	}
	return(best_lambdas)
}

#' @keywords internal
get_best_fit <- function(x_train, y_train, params) {
	message('fitting best model with all training data')
	glmnet_best 	= glmnetcr_propn(
		x_train, y_train
		,method = params$method
		,lambda = params$lambdas
		,maxit 	= params$max_iters
		)
	return(glmnet_best)
}

make_scores_dt <- function(glmnet_best, x_test, y_test, scores_train) {
	scores_test 	= calc_scores_for_one_fit(glmnet_best, x_test, y_test)
	scores_test[, fold := NA]
	scores_test[, data := 'test' ]
	scores_dt 		= rbind(scores_train, scores_test)
	scores_dt[, data := factor(data, levels=c('train', 'test'))]

	return(scores_dt)
}

#' @importFrom data.table data.table
#' @importFrom stringr str_detect
#' @keywords internal
calc_proj_dt <- function(glmnet_best, x_data, y_labels, best_lambdas) {
	# unpack
	which_idx 		= best_lambdas$which_idx

	# get best one
	cut_idx 		= str_detect(rownames(glmnet_best$beta), '^cp[0-9]+$')
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
	predictions 	= predict_glmnetcr_propodds(glmnet_best, x_data, y_labels)
	pred_classes 	= factor(predictions$class[, which_idx], levels=levels(glmnet_best$y))

	# put into data.table
	proj_dt 	= data.table(
		psuper 			= psuper[, 1]
		,label_input 	= y_labels
		,label_psuper 	= pred_classes
		)

	return(proj_dt)
}

#' Extracts best coefficients.
#' 
#' @importFrom data.table data.table
#' @importFrom data.table setorder
#' @importFrom stringr str_detect
#' @return data.table containing learned coefficients for all genes used as input.
#' @keywords internal
make_best_beta <- function(glmnet_best, best_lambdas) {
	cut_idx 	= str_detect(rownames(glmnet_best$beta), '^cp[0-9]+$')
	best_beta 	= glmnet_best$beta[!cut_idx, best_lambdas$which_idx]
	beta_dt 	= data.table( beta=best_beta, symbol=names(best_beta) )
	beta_dt[, abs_beta := abs(beta) ]
	setorder(beta_dt, -abs_beta)
	beta_dt[, symbol:=factor(symbol, levels=beta_dt$symbol)]
		
	return(beta_dt)
}

#' @importFrom data.table data.table
#' @importFrom stringr str_detect
#' @keywords internal
make_psuper_obj <- function(glmnet_best, x_data, y, x_test, y_test, test_idx, fold_list, proj_dt, beta_dt, best_lambdas, best_dt, scores_dt, params) {
	# make cuts_dt
	which_idx 	= best_lambdas$which_idx
	cut_idx 	= str_detect(rownames(glmnet_best$beta), '^cp[0-9]+$')
	cuts_dt 	= data.table(
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
		,test_idx 		= test_idx
		,fold_list 		= fold_list
		,proj_dt 		= proj_dt
		,cuts_dt 		= cuts_dt
		,beta_dt 		= beta_dt
		,best_lambdas 	= best_lambdas
		,best_dt 		= best_dt
		,scores_dt 		= scores_dt
		)

	# make psupertime object
	class(psuper_obj) 	= c('psupertime', class(psuper_obj))

	return(psuper_obj)
}

#' Text to summarize psupertime object
#' @keywords internal
psummarize <- function(psuper_obj) {
	# what trained on
	n_cells 	= dim(psuper_obj$x_data)[1]
	n_genes 	= psuper_obj$params$n_genes
	n_sel 		= dim(psuper_obj$x_data)[2]
	sel_genes 	= psuper_obj$params$sel_genes

	# labels used
	label_order = paste(levels(psuper_obj$y), collapse=', ')

	# accuracy + sparsity
	sel_lambda 	= psuper_obj$best_lambdas$which_lambda
	mean_acc_dt = psuper_obj$scores_dt[ score_var=='class_error', list(mean_acc=mean(score_val)), by=list(lambda, data) ]
	acc_train 	= 1 - mean_acc_dt[ lambda==sel_lambda & data=='train' ]$mean_acc
	acc_test 	= 1 - mean_acc_dt[ lambda==sel_lambda & data=='test' ]$mean_acc
	n_nzero 	= sum(psuper_obj$beta_dt$abs_beta>0)
	sparse_prop = n_nzero / n_sel

	# define outputs
	line_1 		= sprintf('psupertime object using %d cells * %d genes as input\n', n_cells, n_genes)
	line_2 		= sprintf('    label ordering used for training: %s\n', label_order)
	line_3 		= sprintf('    genes selected for input: %s\n', sel_genes)
	line_4 		= sprintf('    # genes taken forward for training: %d\n', n_sel)
	line_5 		= sprintf('    # genes identified as relevant: %d (= %.0f%% of training genes)\n', n_nzero, 100*sparse_prop)
	line_6 		= sprintf('    mean training accuracy: %.0f%%\n', 100*acc_train)
	line_7 		= sprintf('    mean test accuracy: %.0f%%\n', 100*acc_test)

	# join lines together
	psummary 	= paste(line_1, line_2, line_3, line_4, line_5, line_6, line_7, sep = "")
	return(psummary)
}

#' @export
print.psupertime <- function(psuper_obj) {
	psummary 	= psummarize(psuper_obj)
	cat(psummary)
}

#' @importFrom knitr asis_output
#' @keywords internal
knit_print.psupertime = function(psuper_obj, ...) {
	psummary 	= psummarize(psuper_obj)
	asis_output(psummary)
}
