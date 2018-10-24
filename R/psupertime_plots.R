# psupertime_plots.R

#' Convenience function to do multiple plots
#'
#' @param sce SingleCellExperiment class containing all cells and genes required
#' @param x Matrix of log TPM values. Either this or sce should be provided, not both.
#' @param y Vector of labels, which should have same length as number of columns in sce / x. Factor levels will be taken as the intended order for training.
#' @export
psupertime_plot_all <- function(psuper_obj, output_dir='.', tag='', ext='png') {
	# validate model
	cat('plotting results\n')
	g 			= plot_train_results(psuper_obj)
	plot_file 	= file.path(output_dir, sprintf('%s training results.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=6)

	g 			= plot_labels_over_psupertime(psuper_obj)
	plot_file 	= file.path(output_dir, sprintf('%s labels over psupertime.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=12)

	g 			= plot_identified_gene_coefficients(psuper_obj)
	plot_file 	= file.path(output_dir, sprintf('%s identified genes.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=8)

	g 			= plot_identified_genes_over_psupertime(psuper_obj)
	plot_file 	= file.path(output_dir, sprintf('%s identified genes over psupertime.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=8, width=12)

	g 				= plot_predictions_against_classes(psuper_obj)
	plot_file 		= file.path(output_dir, sprintf('%s predictions over psupertime, original data.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=10)
}

#' Plot results of training
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @return ggplot2 object showing test and training performance of classifier.
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_linerange
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 theme_bw
plot_train_results <- function(psuper_obj) {
	# unpack
	params 			= psuper_obj$params
	glmnet_best 	= psuper_obj$glmnet_best
	x_test 			= psuper_obj$x_test
	y_test 			= psuper_obj$y_test
	scores_dt 		= psuper_obj$scores_dt
	mean_scores 	= calc_mean_scores(scores_dt)

	# where should vertical lines go?
	lines_best 		= copy(psuper_obj$best_dt)
	lines_best[, selected := score_var==params$score ]

	# calc test results
	test_scores_dt 	= calc_scores_for_one_fit(glmnet_best, x_test, y_test)

	# put together for plot
	plot_dt 		= rbind(
		test_scores_dt[ , list(lambda, score_var, mean_score=score_val, Data='test') ],
		mean_scores[ 	, list(lambda, score_var, mean_score=mean, Data='train') ]
		)
	plot_dt[, Data := factor(Data, levels=c('train', 'test'))]

	# set up
	g = ggplot() +
		aes( x=log10(lambda) )
	
	# plot each fold
	g = g + geom_point(data=scores_dt, aes(fill=factor(fold), y=score_val), colour='transparent', shape=21 ) +
		scale_fill_brewer( palette='Set1' )

	# plot test and training data
	g = g + geom_linerange(data=mean_scores, aes(ymin=mean-se, ymax=mean+se), colour='grey' ) +
		geom_point(data=plot_dt, aes(y=mean_score, colour=Data) ) +
		geom_line(data=plot_dt, aes(y=mean_score, colour=Data) ) +
		scale_colour_manual( values=c('grey', 'black') )

	# annotate with best lambdas, tidy up
	g = g + geom_vline(data=lines_best, aes(xintercept=log10(best_lambda), size=selected), colour='grey', linetype='solid' ) +
		geom_vline(data=lines_best, aes(xintercept=log10(next_lambda), size=selected), colour='grey', linetype='dashed' ) +
		scale_size_manual( values = c(0.5, 1) ) +
		guides( size=FALSE )

	# label nicely
	g = g + 
		facet_grid( score_var ~ ., scales='free_y' ) +
		theme_bw() +
		labs(
			x 		= 'log10( lambda )'
			,y 		= 'Accuracy measure'
			,colour = 'Data'
			,fill 	= 'Fold'
			,title 	= 'Grey is mean over training data with line indicating SE; black is test data'
			)

	return(g)
}

#' Define RColorBrewer palette to use; default is RdBu.
#'
#' @param y_labels List of labels used for training
#' @return Colour values
#' @keywords internal
make_col_vals <- function(y_labels, palette='RdBu') {
	n_labels 	= length(levels(y_labels))
	max_col 	= 11
	if (n_labels==1) {
		col_vals 	= RColorBrewer::brewer.pal(3, palette)
		col_vals 	= col_vals[1]
	} else if (n_labels==2) {
		col_vals 	= RColorBrewer::brewer.pal(3, palette)
		col_vals 	= col_vals[-2]
	} else if (n_labels<=max_col) {
		col_vals 	= RColorBrewer::brewer.pal(n_labels, palette)
	} else {
		col_pal 	= RColorBrewer::brewer.pal(max_col, palette)
		col_vals 	= rev(grDevices::colorRampPalette(col_pal)(n_labels))
	}

	return(col_vals)
}

#' Plots labels over their projected values on psupertime.
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @return ggplot2 object
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_x_continuous
plot_labels_over_psupertime <- function(psuper_obj, palette='RdBu') {
	# unpack
	proj_dt 		= psuper_obj$proj_dt
	glmnet_best 	= psuper_obj$glmnet_best
	best_lambdas 	= psuper_obj$best_lambdas

	# make nice colours
	col_vals 		= make_col_vals(proj_dt$y_label, palette)

	# find cutpoints
	which_idx 	= best_lambdas$which_idx
	cut_idx 	= stringr::str_detect(rownames(glmnet_best$beta), '^cp[0-9]+$')
	cuts_dt 	= data.table::data.table(
		y_proj 		= c(NA, -(glmnet_best$beta[ cut_idx, which_idx ] + glmnet_best$a0[[ which_idx ]]))
		,y_label 	= factor(levels(proj_dt$y_label), levels=levels(proj_dt$y_label))
		)

	g = ggplot(proj_dt) +
		aes( x=y_proj, fill=y_label) +
		geom_density( alpha=0.5, colour=NA ) +
		scale_fill_manual( values=col_vals ) +
		geom_vline( data=cuts_dt, aes(xintercept=y_proj, colour=y_label) ) +
		scale_colour_manual( values=col_vals ) +
		guides(
			fill 	= guide_legend(override.aes = list(alpha=1))
			,colour = FALSE
			) +
		scale_x_continuous( breaks=scales::pretty_breaks() ) +
		labs(
			x 		= 'Pseudotime'
			,y 		= 'Density'
			,fill 	= 'Ordered labels'
			) +
		theme_bw()

	return(g)
}

#' Plots top coefficients
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @return ggplot2 object
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
plot_identified_gene_coefficients <- function(psuper_obj, n=20, abs_cutoff=0.05) {
	# prepare plot
	plot_dt 	= psuper_obj$beta_dt[ abs_beta > abs_cutoff ]
	plot_dt 	= plot_dt[ 1:min(n, nrow(plot_dt)) ]
	max_val 	= ceiling(max(plot_dt$abs_beta)*10)/10

	# plot
	g 	= ggplot(plot_dt) +
		aes( x=symbol, xend=symbol, y=beta, yend=0 ) +
		geom_segment( colour='black' ) +
		geom_point( colour='blue', size=5 ) +
		# geom_hline( yintercept=0, colour='grey' ) +
		scale_y_continuous( breaks=scales::pretty_breaks(), limits=c(-max_val, max_val) ) +
		theme_bw() +
		theme(
			axis.text.x 	= element_text( angle=-45, hjust=0 )
			) +
		labs(
			x 	= 'Gene',
			y 	= 'Coefficient value'
			)
	return(g)
}

#' Plots profiles of identified genes against psupertime.
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @return ggplot2 object
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
plot_identified_genes_over_psupertime <- function(psuper_obj, n_to_plot=25, palette='RdBu') {
	# unpack
	proj_dt 	= psuper_obj$proj_dt
	beta_dt 	= psuper_obj$beta_dt
	x_data 		= psuper_obj$x_data
	params 		= psuper_obj$params

	# aset
	beta_nzero 	= beta_dt[ abs_beta > 0 ]
	n_nzero 	= nrow(beta_nzero)
	top_genes 	= as.character(beta_nzero[1:min(n_to_plot, nrow(beta_nzero))]$symbol)

	# set up data for plotting
	plot_wide 	= cbind(proj_dt, data.table::data.table(x_data[, top_genes]))
	plot_dt 	= data.table::melt.data.table(plot_wide, id=c('y_proj', 'y_label'), variable.name='symbol')
	plot_dt[, symbol := factor(symbol, levels=top_genes)]

	# get colours
	col_vals 	= make_col_vals(plot_dt$y_label, palette)

	# plot
	g =	ggplot(plot_dt) +
		aes( x=y_proj, y=value) +
		geom_point( size=1, aes(colour=y_label) ) +
		geom_smooth(se=FALSE, colour='black') +
		scale_colour_manual( values=col_vals ) +
		scale_shape_manual( values=c(1, 16) ) +
		scale_x_continuous( breaks=scales::pretty_breaks() ) +
		scale_y_continuous( breaks=scales::pretty_breaks() ) +
		facet_wrap( ~ symbol, scales='free_y' ) +
		theme_bw() +
		theme(
			axis.text.x = element_blank()
			) +
		labs(
			x 		= 'Pseudotime'
			,y 		= 'z-scored log2 expression'
			,colour = 'Condition'
			)
	return(g)
}

#' Plots profiles of hand-selected genes against psupertime.
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @return ggplot2 object
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
plot_specified_genes_over_psupertime <- function(psuper_obj, extra_genes, palette='RdBu') {
	stop('not implemented yet')
	# get smoothed data across all genes
	x_all 		= make_x_data(sce, rownames(sce), y_labels, params)

	# restrict to just this set
	extra_genes = intersect(extra_genes, colnames(x_all))
	plot_wide 	= cbind(proj_dt, data.table(x_all[, extra_genes]))
	plot_dt 	= data.table::melt.data.table(plot_wide, id=c('y_proj', 'y_label', 'test'), variable.name='symbol')
	corrs_dt 	= plot_dt[, list(abs_cor = abs(cor(y_proj, value))), by=symbol]
	data.table::setorder(corrs_dt, -abs_cor)
	plot_dt[, symbol := factor(symbol, levels=corrs_dt$symbol)]

	# set up plot
	col_vals 	= make_col_vals(plot_dt$y_label, palette)
	n_genes 	= length(extra_genes)
	plot_ratio 	= 5/4
	ncol 		= ceiling(sqrt(n_genes*plot_ratio))
	nrow 		= ceiling(n_genes/ncol)
	plot_unit 	= 2.5

	# plot
	g =	ggplot(plot_dt) +
		aes( x=y_proj, y=value ) +
		geom_point( size=1, aes(colour=y_label) ) +
		geom_smooth(se=FALSE, colour='black') +
		scale_colour_manual( values=col_vals ) +
		scale_shape_manual( values=c(1, 16) ) +
		scale_x_continuous( breaks=scales::pretty_breaks() ) +
		scale_y_continuous( breaks=scales::pretty_breaks() ) +
		facet_wrap( ~ symbol, scales='free_y', nrow=nrow, ncol=ncol ) +
		theme_bw() +
		theme(
			axis.text.x = element_blank()
			) +
		labs(
			x 		= 'Pseudotime'
			,y 		= 'z-scored log2 expression'
			,colour = 'Condition'
			# ,shape 	= 'Test / train data'
			)
	return(g)
}

#' @keywords internal
process_new_data <- function(psuper_obj, new_x) {
	# process new_x
	params_copy 				= psuper_obj$params
	params_copy$sel_genes 		= 'list'
	params_copy$gene_list 		= colnames(psuper_obj$x_data)
	params_copy$min_expression 	= 0
	sel_genes 					= select_genes(new_x, params_copy)
	new_data 					= make_x_data(new_x, sel_genes, params_copy)
	return(new_data)
}

#' Gives projection of data onto psupertime (either using original data, or new data)
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @param new_x, new_y Optional pair of new data and labels
#' @return data.table with projection and labels
#' @export
project_onto_psupertime <- function(psuper_obj, new_x=NULL, new_y=NULL) {
	# unpack
	glmnet_best 	= psuper_obj$glmnet_best
	best_lambdas 	= psuper_obj$best_lambdas

	# project new data
	if ( is.null(new_x) & is.null(new_y) ) {
		x_in 			= psuper_obj$x_data
		y_in 			= psuper_obj$y
	} else if ( !is.null(new_x) & !is.null(new_y) ) {
		x_in 			= process_new_data(psuper_obj, new_x)
		if (!is.factor(new_y)) {
			new_y 			= factor(new_y)
			message('converting new_y into factor, with the following ordered values:')
			message(paste(levels(new_y), ', '))
			message('(define new_y as a factor if you prefer a different ordering)')
		}
		y_in 			= factor(new_y)
	} else {
		stop('either both of new_x and new_y must be given, or neither')
	}

	proj_dt 		= calc_proj_dt(glmnet_best, x_in, y_in, best_lambdas)

	return(proj_dt)
}

#' Plots profiles of hand-selected genes against psupertime.
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @param new_x,new_y Optional data to predict with psuper_obj
#' @return ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_raster
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom ggplot2 expand_limits
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
plot_new_data_over_psupertime <- function(psuper_obj, new_x, new_y, palette='BrBG') {
	# project new data
	proj_new 		= project_onto_psupertime(psuper_obj, new_x, new_y)

	# make nice colours
	col_vals 		= make_col_vals(proj_new$y_label, palette)

	# find cutpoints
	glmnet_best 	= psuper_obj$glmnet_best
	which_idx 		= psuper_obj$best_lambdas$which_idx
	cut_idx 		= stringr::str_detect(rownames(glmnet_best$beta), '^cp[0-9]+$')
	cuts_dt 		= data.table::data.table(
		y_proj 			= -(glmnet_best$beta[ cut_idx, which_idx ] + glmnet_best$a0[[ which_idx ]])
		)

	# do plot
	g = ggplot(proj_new) +
		aes( x=y_proj, fill=y_label) +
		geom_density( alpha=0.5, colour=NA ) +
		scale_fill_manual( values=col_vals ) +
		geom_vline( data=cuts_dt, aes(xintercept=y_proj), colour='grey' ) +
		scale_colour_manual( values=col_vals ) +
		guides(
			fill 	= guide_legend(override.aes = list(alpha=1))
			,colour = FALSE
			) +
		scale_x_continuous( breaks=scales::pretty_breaks() ) +
		labs(
			x 		= 'Pseudotime'
			,y 		= 'Density'
			,fill 	= 'Ordered labels'
			) +
		theme_bw()

	return(g)
}


#' Plots confusion matrix of true labels against predicted labels.
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @param new_x,new_y Optional data to predict with psuper_obj
#' @param plot_var Variable to plot: prop_true is proportion of true labels, prop_predict is proportion of predicted labels, N is # of cells
#' @return ggplot2 object
#' @export
#' @keywords internal
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 expand_limits
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_raster
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme_bw
plot_predictions_against_classes <- function(psuper_obj, new_x=NULL, new_y=NULL, plot_var='prop_true') {
	# decide what to plot
	plot_var_list 	= c('prop_true', 'N', 'prop_predict')
	plot_var 		= match.arg(plot_var, plot_var_list)
	labels_list 	= c(prop_true='Proportion\nof labelled\nclass\n', N='# of cells', prop_predict='Proportion\nof predicted\nclass\n')
	plot_label 		= labels_list[[plot_var]]

	# unpack
	which_idx 		= psuper_obj$best_lambdas$which_idx
	glmnet_best 	= psuper_obj$glmnet_best

	# define fn to handle y
	get_y_in <- function(new_y) {
		if (is.null(new_x)) {
			if ( length(new_y) != length(psuper_obj$y) ) {
				stop('when no new_x given, new_y must be same length as original y')
			}
		}
		if (!is.factor(new_y)) {
			new_y 			= factor(new_y)
			message('converting new_y into factor, with the following ordered values:')
			message(paste(levels(new_y), ', '))
			message('(define new_y as a factor if you prefer a different ordering)')
		}
		y_in 			= factor(new_y)
		return(y_in)
	}

	# what inputs to use?
	if ( is.null(new_x) & is.null(new_y) ) {
		x_in 			= psuper_obj$x_data
		y_in 			= psuper_obj$y

	} else if ( is.null(new_x) & !is.null(new_y) ) {
		x_in 			= psuper_obj$x_data
		y_in 			= get_y_in(new_y)

	} else if ( !is.null(new_x) & !is.null(new_y) ) {
		x_in 			= process_new_data(psuper_obj, new_x)
		y_in 			= get_y_in(new_y)

	} else if ( !is.null(new_x) & is.null(new_y) ) {
		stop('to use new_x, new_y must also be given')

	} else {
		stop('aargh some unexpected error')
	}

	# get predicted classes for each thing
	predictions 	= predict_glmnetcr_cumul(glmnet_best, x_in, y_in)
	pred_classes 	= factor(predictions$class[, which_idx], levels=levels(psuper_obj$y))
	predict_dt 		= data.table::data.table( predicted=pred_classes, true=y_in )

	# count and average
	counts_dt 		= predict_dt[, .N, by=list(predicted, true)]
	counts_dt[, prop_true 		:= N / sum(N), by=true ]
	counts_dt[, prop_predict 	:= N / sum(N), by=predicted ]

	# define where borders should be
	borders_dt 		= counts_dt[as.character(true)==as.character(predicted), list(true, predicted) ]

	# plot grid
	g = ggplot(counts_dt) +
		aes( y=true, x=predicted ) +
		geom_tile( aes_string(fill=plot_var) ) +
		geom_tile(data=borders_dt, aes(y=true, x=predicted), fill=NA, colour='black', size=0.5) +
		geom_text( aes(label=N) ) +
		scale_x_discrete( drop=FALSE ) +
		scale_fill_distiller( palette='BuPu', direction=1, breaks=scales::pretty_breaks() )
	if (plot_var=='N') {
		g = g + expand_limits( fill=0 )
	} else {
		g = g + expand_limits( fill=c(0,1) )
	}
	g = g + labs(
			x 		= 'Predicted class'
			,y 		= 'Labelled class'
			,fill 	= plot_label
			) +
		theme_bw()
	return(g)
}
