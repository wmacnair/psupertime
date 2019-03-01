# psupertime_plots.R

#' Convenience function to do multiple plots
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @param output_dir Directory to save to
#' @param tag Label for all files
#' @param label_name Description for the ordered labels in the legend (e.g. 'Developmental stage (days)')
#' @param ext Image format for outputs, compatible with ggsave (eps, ps, tex, pdf, jpeg, tiff, png, bmp, svg, wmf)
#' @export
psupertime_plot_all <- function(psuper_obj, output_dir='.', tag='', label_name='Ordered labels', ext='png') {
	# validate model
	cat('plotting results\n')
	g 			= plot_train_results(psuper_obj)
	plot_file 	= file.path(output_dir, sprintf('%s training results.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=6)

	g 			= plot_labels_over_psupertime(psuper_obj, label_name)
	plot_file 	= file.path(output_dir, sprintf('%s labels over psupertime.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=12)

	g 			= plot_identified_gene_coefficients(psuper_obj)
	plot_file 	= file.path(output_dir, sprintf('%s identified genes.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=8)

	g 			= plot_identified_genes_over_psupertime(psuper_obj, label_name)
	plot_file 	= file.path(output_dir, sprintf('%s identified genes over psupertime.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=8, width=12)

	g 			= plot_predictions_against_classes(psuper_obj)
	plot_file 	= file.path(output_dir, sprintf('%s predictions over psupertime, original data.%s', tag, ext))
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

	# add nice labels for accuracy measures
	measure_list 	= c(xentropy='Cross entropy', class_error='Classification error')
	mean_scores[, 	nice_score_var := factor(measure_list[score_var], levels=measure_list) ]
	plot_dt[, 		nice_score_var := factor(measure_list[score_var], levels=measure_list) ]
	lines_best[, 	nice_score_var := factor(measure_list[score_var], levels=measure_list) ]

	# set up
	g = ggplot() +
		aes( x=log10(lambda) )
	
	# # plot each fold
	# g = g + geom_point(data=scores_dt, aes(fill=factor(fold), y=score_val), colour='transparent', shape=21 ) +
	# 	scale_fill_brewer( palette='Set1' )

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
		facet_grid( nice_score_var ~ ., scales='free_y' ) +
		theme_bw() +
		labs(
			x 		= 'log10( lambda )'
			,y 		= 'Accuracy measure'
			,colour = 'Data'
			# ,fill 	= 'Fold'
			,title 	= sprintf('%s used for model selection', measure_list[params$score])
			) +
		theme(
			plot.title 	= element_text( size=10, hjust=1 )
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
		col_vals 	= grDevices::colorRampPalette(col_pal)(n_labels)
	}
	col_vals 	= rev(col_vals)

	return(col_vals)
}

#' Plots labels over their projected values on psupertime.
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @param label_name Description for the ordered labels in the legend (e.g. 'Developmental stage (days)')
#' @param palette RColorBrewer palette to use
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
plot_labels_over_psupertime <- function(psuper_obj, label_name='Ordered labels', palette='RdBu') {
	# unpack
	proj_dt 		= psuper_obj$proj_dt
	cuts_dt 		= psuper_obj$cuts_dt

	# make nice colours
	col_vals 		= make_col_vals(proj_dt$label_input, palette)

	# plot
	g = ggplot(proj_dt) +
		aes( x=psuper, fill=label_input ) +
		geom_density( alpha=0.5, colour='black' ) +
		scale_fill_manual( values=col_vals ) +
		geom_vline( data=cuts_dt, aes(xintercept=psuper, colour=label_input) ) +
		scale_colour_manual( values=col_vals ) +
		guides(
			fill 	= guide_legend(override.aes = list(alpha=1))
			,colour = FALSE
			) +
		scale_x_continuous( breaks=scales::pretty_breaks() ) +
		labs(
			x 		= 'psupertime'
			,y 		= 'Density'
			,fill 	= label_name
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
#' @param label_name Description for the ordered labels in the legend (e.g. 'Developmental stage (days)')
#' @param n_to_plot Maximum number of genes to plot (default 20)
#' @param palette RColorBrewer palette to use
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
plot_identified_genes_over_psupertime <- function(psuper_obj, label_name='Ordered labels', n_to_plot=20, palette='RdBu') {
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
	plot_dt 	= data.table::melt.data.table(plot_wide, id=c('psuper', 'label_input', 'label_psuper'), variable.name='symbol')
	plot_dt[, symbol := factor(symbol, levels=top_genes)]

	# get colours
	col_vals 	= make_col_vals(plot_dt$label_input, palette)

	# plot
	g =	ggplot(plot_dt) +
		aes( x=psuper, y=value) +
		geom_point( size=1, aes(colour=label_input) ) +
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
			x 		= 'psupertime'
			,y 		= 'z-scored log2 expression'
			,colour = label_name
			)
	return(g)
}

#' Plots profiles of hand-selected genes against psupertime.
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @param extra_genes 
#' @param label_name Description for the ordered labels in the legend (e.g. 'Developmental stage (days)')
#' @param palette RColorBrewer palette to use
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
plot_specified_genes_over_psupertime <- function(psuper_obj, extra_genes, label_name='Ordered labels', palette='RdBu') {
	stop('not implemented yet')
	# get smoothed data across all genes
	x_all 		= make_x_data(sce, rownames(sce), y_labels, params)

	# restrict to just this set
	extra_genes = intersect(extra_genes, colnames(x_all))
	plot_wide 	= cbind(proj_dt, data.table(x_all[, extra_genes]))
	plot_dt 	= data.table::melt.data.table(plot_wide, id=c('psuper', 'label_input', 'label_psuper'), variable.name='symbol')
	corrs_dt 	= plot_dt[, list(abs_cor = abs(cor(psuper, value))), by=symbol]
	data.table::setorder(corrs_dt, -abs_cor)
	plot_dt[, symbol := factor(symbol, levels=corrs_dt$symbol)]

	# set up plot
	col_vals 	= make_col_vals(plot_dt$label_input, palette)
	n_genes 	= length(extra_genes)
	plot_ratio 	= 5/4
	ncol 		= ceiling(sqrt(n_genes*plot_ratio))
	nrow 		= ceiling(n_genes/ncol)
	plot_unit 	= 2.5

	# plot
	g =	ggplot(plot_dt) +
		aes( x=psuper, y=value ) +
		geom_point( size=1, aes(colour=label_input) ) +
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
			x 		= 'psupertime'
			,y 		= 'z-scored log2 expression'
			,colour = label_name
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
#' @param palette RColorBrewer palette to use
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
	col_vals 		= make_col_vals(proj_new$label_input, palette)

	# get cutpoints
	cuts_dt 		= psuper_obj$cuts_dt

	# do plot
	g = ggplot(proj_new) +
		aes( x=psuper, fill=label_input) +
		geom_density( alpha=0.5, colour=NA ) +
		scale_fill_manual( values=col_vals ) +
		geom_vline( data=cuts_dt, aes(xintercept=psuper, colour=label_input) ) +
		scale_colour_manual( values=col_vals ) +
		guides(
			fill 	= guide_legend(override.aes = list(alpha=1))
			,colour = FALSE
			) +
		scale_x_continuous( breaks=scales::pretty_breaks() ) +
		labs(
			x 		= 'psupertime'
			,y 		= 'Density'
			,fill 	= 'Ordered labels'
			) +
		theme_bw()

	return(g)
}


#' Check variables for confusion matrices
#'
#' @param plot_var Variable to plot: prop_true is proportion of true labels, prop_predict is proportion of predicted labels, N is # of cells
#' @return list with checked plot_var, and nice label
#' @internal
check_conf_params <- function(plot_var) {
	plot_var_list 	= c('prop_true', 'N', 'prop_predict')
	plot_var 		= match.arg(plot_var, plot_var_list)
	labels_list 	= c(prop_true='Proportion\nof labelled\nclass\n', N='# of cells', prop_predict='Proportion\nof predicted\nclass\n')
	plot_label 		= labels_list[[plot_var]]

	return( list(plot_var=plot_var, plot_label=plot_label) )
}

#' Plots confusion matrix of true labels against predicted labels.
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @param new_x,new_y Optional data to predict with psuper_obj
#' @param plot_var Variable to plot: prop_true is proportion of true labels, prop_predict is proportion of predicted labels, N is # of cells
#' @param palette RColorBrewer palette to use
#' @return ggplot2 object
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 expand_limits
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_raster
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme_bw
plot_predictions_against_classes <- function(psuper_obj, new_x=NULL, new_y=NULL, plot_var='prop_true', palette='BuPu') {
	# decide what to plot
	conf_params 	= check_conf_params(plot_var)
	plot_var 		= conf_params$plot_var
	plot_label 		= conf_params$plot_label

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
		scale_fill_distiller( palette=palette, direction=1, breaks=scales::pretty_breaks() )
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

#' Projects two different psupertimes onto each other
#'
#' @param psuper_1, psuper_2 Two previously calculated psupertime objects
#' @param labels Character vector of length two, labelling the psupertime inputs
#' @return data.table containing projections in both directions
#' @export
double_psupertime <- function(psuper_1, psuper_2, run_names=NULL) {
	# check run_names
	if ( is.null(run_names) ) {
		run_names 	= c('1','2')
		message('using default values for run_names:', paste(run_names, sep=', '))
	} else {
		if ( !is.character(run_names) | length(unique(run_names))!=2 ) {
			stop('run_names must be character vector of length two with no repeated values')
		}
	}

	# repack
	psuper_list 	= list(psuper_1, psuper_2)
	n_psupers 		= length(psuper_list)
	# names(psuper_list) 	= run_names

	# loop through projections on both
	doubles_dt 		= data.table()
	for (ii in 1:n_psupers) {
		# unload
		psuper_ii 		= psuper_list[[ii]]
		label_ii 		= run_names[[ii]]

		for (jj in 1:n_psupers) {
			# unload
			psuper_jj 		= psuper_list[[jj]]
			label_jj 		= run_names[[jj]]

			# get appropriate projection
			if (ii == jj) {
				proj_ii_on_jj 	= psuper_ii$proj_dt
			} else {
				proj_ii_on_jj 	= project_onto_psupertime(psuper_jj, psuper_ii$x_data, psuper_ii$y)
			}

			# label
			proj_ii_on_jj[, input 		:= label_ii ]
			proj_ii_on_jj[, projection 	:= label_jj ]
			n_digits = ceiling(log10(nrow(psuper_ii$x_data)))
			proj_ii_on_jj[, cell_id 	:= sprintf(sprintf('%%s_%%0%dd', n_digits), label_ii, 1:nrow(psuper_ii$x_data)) ]

			# store
			doubles_dt 		= rbind(doubles_dt, proj_ii_on_jj)
		}
	}

	# sort out levels
	lvls_all 		= c()
	for (ii in 1:n_psupers) {
		lvls_temp 	= setdiff(levels(psuper_list[[ii]]$y), lvls_all)
		lvls_all 	= c(lvls_all, lvls_temp)
	}
	doubles_dt[, label_input := factor(label_input, levels=lvls_all) ]

	# make wide, sort out levels?
	doubles_wide 	= dcast(doubles_dt, input + cell_id + label_input ~ projection, value.var=c('psuper', 'label_psuper'))
	for (ii in 1:n_psupers) {
		label 	= run_names[[ii]]
		doubles_wide[[ paste0('label_psuper_', label) ]] 	= fct_drop(doubles_wide[[ paste0('label_psuper_', label) ]])
		# levels(doubles_wide[[ paste0('label_psuper_', label) ]]) 	= levels(psuper_list[[ii]]$y)
	}

	# make lists of levels
	levels_list 		= lapply(psuper_list, function(p) levels(p$y) )
	names(levels_list) 	= run_names

	# put into list
	double_obj 	= list(
		run_names 		= run_names
		,levels_list 	= levels_list
		,doubles_dt 	= doubles_dt
		,doubles_wide 	= doubles_wide
		)
	return(double_obj)
}

#' Projects two different psupertimes onto each other, using points, side by side
#'
#' @param double_obj Result of applying double_psupertime to two previously calculated psupertime objects
#' @param psuper_1, psuper_2 Two previously calculated psupertime objects
#' @param run_names Character vector of length two, labelling the psupertime inputs
#' @return ggplot object plotting the two against each other
#' @export
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 theme_bw
plot_double_psupertime <- function(double_obj, psuper_1=NULL, psuper_2=NULL, run_names=NULL) {
	# check inputs
	if (is.null(double_obj)) {
		if ( is.null(psuper_1) | is.null(psuper_2) ) {
			stop('either a double_obj must be given, or psuper_1 and psuper_2 must both be given')
			double_obj 		= double_psupertime(psuper_1, psuper_2, run_names)
		}
	}

	# unpack
	run_names 		= double_obj$run_names
	label_x 		= run_names[[1]]
	label_y 		= run_names[[2]]
	doubles_wide 	= double_obj$doubles_wide

	# make colours
	col_vals 		= make_col_vals(doubles_wide$label_input)

	# add facet labels
	plot_dt 		= copy(doubles_wide)
	plot_dt[, input_label := paste0('Input data: ', input) ]

	# do some plotting
	g = ggplot(plot_dt) +
		aes_string(
			x 			= paste0('psuper_', label_x)
			,y 			= paste0('psuper_', label_y)
			,colour 	= paste0('label_input')
			) +
		geom_point() +
		scale_colour_manual( values=col_vals ) +
		facet_grid( . ~ input_label) +
		theme_bw() +
		labs(
			x 			= paste0('Psupertime trained on ', label_x)
			,y 			= paste0('Psupertime trained on ', label_y)
			,colour 	= 'Known\nlabels'
			)

	return(g)
}

#' Projects two different psupertimes on top of each other
#'
#' @param double_obj Result of applying double_psupertime to two previously calculated psupertime objects
#' @param psuper_1, psuper_2 Two previously calculated psupertime objects
#' @param run_names Character vector of length two, labelling the psupertime inputs
#' @return ggplot object plotting the two against each other
#' @export
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_density2d
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_colour_brewer
#' @importFrom ggplot2 theme_bw
plot_double_psupertime_contour <- function(double_obj, psuper_1=NULL, psuper_2=NULL, run_names=NULL) {
	# check run_names
	if ( is.null(run_names) ) {
		run_names 	= c('1','2')
		message('using default values for run_names:', paste(run_names, sep=', '))
	} else {
		if ( !is.character(run_names) | length(unique(run_names))!=2 ) {
			stop('run_names must be character vector of length two with no repeated values')
		}
	}
	# check inputs
	if (is.null(double_obj)) {
		if ( is.null(psuper_1) | is.null(psuper_2) ) {
			stop('either a double_obj must be given, or psuper_1 and psuper_2 must both be given')
			double_obj 		= double_psupertime(psuper_1, psuper_2, run_names)
		}
	}

	# unpack
	run_names 		= double_obj$run_names
	label_x 		= run_names[[1]]
	label_y 		= run_names[[2]]
	doubles_wide 	= double_obj$doubles_wide

	# do some plotting
	g = ggplot(doubles_wide) +
		aes_string(
			x 			= paste0('psuper_', label_x)
			,y 			= paste0('psuper_', label_y)
			,colour 	= 'input'
			) +
		geom_density2d() +
		scale_colour_brewer( palette='Set1' ) +
		theme_bw() +
		labs(
			x 			= paste0('Psupertime run on ', label_x)
			,y 			= paste0('Psupertime run on ', label_y)
			,colour 	= 'Input\ndata'
			)

	return(g)
}

#' Compares coefficients for genes learned from different psupertimes
#'
#' @param psuper_1, psuper_2 Two previously calculated psupertime objects
#' @param run_names Character vector of length two, labelling the psupertime inputs
#' @return ggplot object plotting the two sets of coefficients
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
plot_double_psupertime_genes <- function(psuper_1, psuper_2, run_names=NULL) {
	# check run_names
	if ( is.null(run_names) ) {
		run_names 	= c('1','2')
		message('using default values for run_names:', paste(run_names, sep=', '))
	} else {
		if ( !is.character(run_names) | length(unique(run_names))!=2 ) {
			stop('run_names must be character vector of length two with no repeated values')
		}
	}

	# get genes from both
	old_names 	= c('beta', 'abs_beta')
	genes_1_dt 	= psuper_1$beta_dt[ abs_beta > 0 ]
	data.table::setnames(genes_1_dt, old_names, paste0(old_names, '_1'))
	genes_2_dt 	= psuper_2$beta_dt[ abs_beta > 0 ]
	data.table::setnames(genes_2_dt, old_names, paste0(old_names, '_2'))

	# join together, tidy up
	genes_dt 	= merge(genes_1_dt, genes_2_dt, by='symbol', all=TRUE, )
	genes_dt[ is.na(beta_1), beta_1 := 0 ]
	genes_dt[ is.na(abs_beta_1), abs_beta_1 := 0 ]
	genes_dt[ is.na(beta_2), beta_2 := 0 ]
	genes_dt[ is.na(abs_beta_2), abs_beta_2 := 0 ]

	# plot
	g = ggplot(genes_dt) +
		aes( x=beta_1, y=beta_2 ) +
		geom_point( alpha=0.5 ) +
		theme_bw() +
		labs(
			x 		= paste0('Coefficient for ', run_names[[1]])
			,y 		= paste0('Coefficient for ', run_names[[2]])
			)

	return(g)
}

#' Plots the confusion matrices of two psupertime objects against each other
#'
#' @param double_obj Result of applying double_psupertime to two previously calculated psupertime objects
#' @param psuper_1, psuper_2 Two previously calculated psupertime objects
#' @param run_names Character vector of length two, labelling the psupertime inputs
#' @param palette RColorBrewer palette to use
#' @return cowplot plot_grid object, showing known and predicted labels for each dataset, and each set of predictions
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 expand_limits
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme_bw
plot_double_psupertime_confusion <- function(double_obj, psuper_1=NULL, psuper_2=NULL, run_names=NULL, plot_var='prop_true', palette='BuPu') {
	if ( !requireNamespace("cowplot", quietly=TRUE) ) {
		message('cowplot not installed; not plotting confusion matrix')
		return()
	}

	# decide what to plot
	conf_params 	= check_conf_params(plot_var)
	plot_var 		= conf_params$plot_var
	plot_label 		= conf_params$plot_label

	# check inputs
	if (is.null(double_obj)) {
		if ( is.null(psuper_1) | is.null(psuper_2) ) {
			stop('either a double_obj must be given, or psuper_1 and psuper_2 must both be given')
			double_obj 		= double_psupertime(psuper_1, psuper_2, run_names)
		}
	}

	# unpack
	run_names 		= double_obj$run_names
	label_x 		= run_names[[1]]
	label_y 		= run_names[[2]]
	doubles_dt 		= double_obj$doubles_dt

	# set up
	input_list 		= unique(doubles_dt$input)
	n_inputs 		= length(input_list)
	proj_list 		= unique(doubles_dt$projection)
	n_projs 		= length(proj_list)
	g_list 			= list()
		
	# get factor lists
	levels_list 	= double_obj$levels_list

	# do multiple plots
	for (ii in 1:n_inputs) {
		for (jj in 1:n_projs) {
			# restrict to this combo of inputs/predictions
			input_ii 		= input_list[[ii]]
			psuper_jj 		= proj_list[[jj]]
			counts_dt 		= doubles_dt[ input==input_ii & projection==psuper_jj, .N, by=list(label_input, label_psuper) ]

			# calculate proportions
			counts_dt[, prop_true 		:= N / sum(N), 		by=label_input ]
			counts_dt[, prop_predict 	:= N / sum(N), 	by=label_psuper ]

			# tidy up labels
			counts_dt[, label_input 	:= forcats::fct_drop(label_input) ]
			counts_dt[, label_input 	:= factor(label_input, levels=levels_list[[input_ii]])]
			counts_dt[, label_psuper 	:= forcats::fct_drop(label_psuper) ]
			counts_dt[, label_psuper 	:= factor(label_psuper, levels=levels_list[[psuper_jj]])]

			# define where borders should be
			borders_dt 		= counts_dt[as.character(label_input)==as.character(label_psuper), list(label_input, label_psuper) ]

			# plot grid
			g = ggplot(counts_dt) +
				aes( y=label_input, x=label_psuper ) +
				geom_tile( aes_string(fill=plot_var) ) +
				geom_tile(data=borders_dt, aes(y=label_input, x=label_psuper), fill=NA, colour='black', size=0.5) +
				geom_text( aes(label=N) ) +
				scale_x_discrete( drop=FALSE ) +
				scale_fill_distiller( palette=palette, direction=1, breaks=scales::pretty_breaks(), guide=FALSE ) +
				theme_bw()

			# colouring for tiles
			if (plot_var=='N') {
				g = g + expand_limits( fill=0 )
			} else {
				g = g + expand_limits( fill=c(0,1) )
			}

			# x, y labels
			if ( ii==n_inputs ) {
				g = g + labs( x=paste0('Predicted: ', run_names[[jj]]) )
			} else {
				g = g + labs( x=NULL )
			}
			if ( jj==1 ) {
				g = g + labs( y=paste0('Known: ', run_names[[ii]]) )
			} else {
				g = g + labs( y=NULL )
			}

			g_list[[ (ii - 1)*n_inputs + jj ]] 	= g
		}
	}

	g_grid 			= cowplot::plot_grid(plotlist=g_list, labels=NULL, nrow=n_inputs, ncol=n_projs, align='h', axis='b')

	return(g_grid)
}

#' GO enrichment analysis for genes learned from different psupertimes
#'
#' @param psuper_obj A previously calculated psupertime object
#' @param org_mapping Organism to use for annotations (e.g. 'org.Mm.eg.db', 'org.Hs.eg.db')
#' @return data.table containing results of GO enrichment analysis
#' @internal
psupertime_go_analysis_old <- function(psuper_obj, org_mapping) {
	# can we do this?
	if ( !requireNamespace("topGO", quietly=TRUE) ) {
		message('topGO not installed; not doing GO analysis')
		return()
	}
	library('topGO')

	# unpack
	psuper 			= scale(psuper_obj$proj_dt$psuper)
	n_obs 			= length(psuper)
	x_data 			= psuper_obj$x_data

	# calculate correlations
	corrs 			= as.vector(matrix(psuper, nrow=1) %*% x_data) / n_obs
	names(corrs) 	= colnames(x_data)

	# calculate p values for these
	t_stat 			= (corrs*sqrt(n_obs-2))/sqrt(1-corrs^2)
	p_vals 			= 2*(1 - pt(abs(t_stat),(n_obs-2)))

	# do GO in various ways
	go_dt 			= data.table::data.table()
	for (up_or_down in c('both', 'up', 'down')) {
		# do ranking
		if (up_or_down=='both') {
			scores 		= abs(corrs)

		} else if (up_or_down=='up') {
			scores 		= corrs

		} else if (up_or_down=='down') {
			scores 		= -corrs

		}
		scores[ scores < 0 ] 	= 0
		scores 		= sort(scores, decreasing=TRUE)
		if ( sum(scores > 0)==0 ) {
			next
		}

		# make topGO object
		topGO_data = new("topGOdata", 
			description 		= up_or_down, 
			allGenes 			= scores, 
			geneSel 			= function(x) {x>0.1},
			annot 				= topGO::annFUN.org, 
			mapping 			= org_mapping, 
			ontology 			= 'BP',
			ID 					= 'symbol'
			)

		# run enrichment tests on these, extract results
		go_weight 	= topGO::runTest(topGO_data, algorithm = "weight01", statistic = "fisher")
		go_temp 	= data.table::data.table(topGO::GenTable(topGO_data, 
			p_go 	= go_weight, 
			orderBy 		= 'p_go', 
			ranksOf 		= 'p_go', 
			topNodes 		= 1000
			))
		data.table::setnames(go_temp, 'p_go', 'temp')
		go_temp[, p_go := as.numeric(temp) ]
		go_temp[ temp == '< 1e-30', p_go := 9e-31 ]
		go_temp[, temp := NULL ]
		go_temp[, direction := up_or_down]
		go_temp[, rank := 1:nrow(go_temp)]

		# store
		go_dt 		= rbind(go_dt, go_temp)
	}

	# change column order
	setcolorder(go_dt, c('direction', 'rank'))
	# print top terms
	p_cutoff 		= 5e-2
	n_terms_cutoff 	= 5
	print_dt 		= go_dt[ p_go < p_cutoff & Significant>n_terms_cutoff ]
	if (nrow(print_dt)==0) {
		message(sprintf('no GO terms met the cutoffs (p-value < %.1e and at least %d genes significant)', p_cutoff, n_terms_cutoff))
	} else {
		message('Significant GO terms:')
		print(print_dt)
	}

	return(go_dt)
}

#' GO enrichment analysis for genes learned from different psupertimes
#'
#' @param psuper_obj A previously calculated psupertime object
#' @param org_mapping Organism to use for annotations (e.g. 'org.Mm.eg.db', 'org.Hs.eg.db')
#' @return data.table containing results of GO enrichment analysis
#' @export
psupertime_go_analysis <- function(psuper_obj, org_mapping, k=5, sig_cutoff=5) {
	if ( !requireNamespace("topGO", quietly=TRUE) ) {
		message('topGO not installed; not doing GO analysis')
		return()
	}
	if ( !requireNamespace("fastcluster", quietly=TRUE) ) {
		message('fastcluster not installed; not doing GO analysis')
		return()
	}
	library('topGO')
	library('fastcluster')

	# unpack
	glmnet_best 	= psuper_obj$glmnet_best
	best_lambdas 	= psuper_obj$best_lambdas
	proj_dt 		= copy(psuper_obj$proj_dt)
	x_data 			= copy(psuper_obj$x_data)
	beta_dt 		= psuper_obj$beta_dt
	cuts_dt 		= psuper_obj$cuts_dt

	# put cells in nice order, label projections
	rownames(x_data) 	= sprintf('cell_%04d', 1:nrow(x_data))
	proj_dt[, cell_id := rownames(x_data) ]
	setorder(proj_dt, psuper)
	proj_dt[, cell_id := factor(cell_id) ]

	# do clustering on symbols
	message('clustering genes')
	hclust_obj 		= fastcluster::hclust(dist(t(x_data)), method='complete')

	# extract clusters from them
	clusters_dt 	= calc_clusters_dt(hclust_obj, x_data, proj_dt, k)
	go_results 		= do_topgo_for_cluster(clusters_dt, sig_cutoff, org_mapping)

	# make plot_dt
	plot_dt 		= make_plot_dt(x_data, hclust_obj, proj_dt, clusters_dt)

	# assemble outputs
	go_list = list(
		clusters_dt 	= clusters_dt,
		go_results 		= go_results,
		plot_dt 		= plot_dt,
		cuts_dt 		= copy(psuper_obj$cuts_dt)
		)

	return(go_list)
}

#' make nice data.table of hierarchical clusters
#'
#' @param hclust_obj Result of hclust
#' @param x_data Data used to calculate psuper_obj
#' @param proj_dt Projection of cells onto psupertime
#' @return data.table containing clusters of genes, ordered according to correlation with psupertime
#' @internal
calc_clusters_dt <- function(hclust_obj, x_data, proj_dt, k=5) {
	# make thing
	clusters_dt 	= data.table( h_clust=cutree(hclust_obj, k=k), symbol=colnames(x_data))
	# add clustering
	clusters_dt[, N:=.N, by=h_clust ]

	# order by correlation with psupertime
	temp_dt 			= data.table(melt(x_data, varnames=c('cell_id', 'symbol')))
	temp_dt 			= clusters_dt[ temp_dt, on='symbol' ]
	means_dt 			= temp_dt[, list(mean=mean(value)), by=list(cell_id, h_clust) ]
	means_dt 			= proj_dt[ means_dt, on='cell_id' ]
	corrs_dt 			= means_dt[, list( cor=cor(mean, psuper) ), by=h_clust]
	setorder(corrs_dt, cor)
	corrs_dt[, clust := 1:.N ]
	corrs_dt[, clust := factor(clust)]
	
	# add clusters ordered by size back in
	clusters_dt			= corrs_dt[ clusters_dt, on='h_clust' ]
	clusters_dt[, clust_label := factor(sprintf('%02d (%d genes)', clust, N)) ]
	clusters_dt[, h_clust := NULL ]
	setorder(clusters_dt, clust, symbol)

	return(clusters_dt)
}

#' Calculate GO enrichment for each cluster vs all other genes
#'
#' @param clusters_dt 
#' @param sig_cutoff How many genes should be in the cluster for us to consider a GO term?
#' @return data.table with GO term results
#' @internal
do_topgo_for_cluster <- function(clusters_dt, sig_cutoff, org_mapping) {
	# set up
	all_clusters 	= unique(clusters_dt[N>=sig_cutoff]$clust)
	go_results 		= data.table()

	# loop through clusters
	message(sprintf('calculating GO enrichments for %d clusters:', length(all_clusters)))
	for (c in all_clusters) {
		message('.', appendLF=FALSE)
		gene_list 			= factor( as.integer(clusters_dt$clust == c) )
		names(gene_list) 	= clusters_dt$symbol

		# make topGO object
		suppressMessages({
			topGO_data 	= new("topGOdata", 
			description 		= c, 
			allGenes 			= gene_list, 
			# geneSelectionFun 	= function(x) {x==TRUE},
			annot 				= annFUN.org, 
			mapping 			= org_mapping, 
			ontology 			= 'BP',
			ID 					= 'symbol'
			)
		})
		# run enrichment tests on these, extract results
		suppressMessages({go_weight 	= runTest(topGO_data, algorithm = "weight", statistic = "fisher")})
		n_terms 		= length(go_weight@score)
		temp_results 	= data.table(GenTable(topGO_data, 
			p_go 	= go_weight, 
			orderBy 		= 'p_go', 
			ranksOf 		= 'p_go', 
			topNodes 		= n_terms
			))
		temp_results[ , cluster := c ]

		# store
		go_results 	= rbind(go_results, temp_results)
	}
	message('')

	# tidy up
	setnames(go_results, 'p_go', 'tmp')
	go_results[, p_go := as.numeric(tmp) ]
	go_results[ tmp == '< 1e-30', p_go := 9e-31 ]
	go_results[ , tmp := NULL ]
	go_results[ , cluster := factor(cluster, levels=all_clusters) ]

	return(go_results)
}

#' Internal function
#'
#' @param x_data 
#' @param hclust_obj 
#' @param proj_dt 
#' @param clusters_dt 
#' @return data.table for plotting
#' @internal
make_plot_dt <- function(x_data, hclust_obj, proj_dt, clusters_dt) {
	# plot
	plot_dt 			= data.table(melt(x_data, varnames=c('cell_id', 'symbol')))

	# nice ordering
	symbol_order 		= colnames(x_data)[hclust_obj$order]
	plot_dt[, symbol 	:= factor(symbol, levels=symbol_order)]
	plot_dt[, cell_id 	:= factor(cell_id, levels=proj_dt$cell_id)]

	# put this into plotting 
	plot_dt 			= clusters_dt[ plot_dt, on='symbol' ]
	plot_dt 			= proj_dt[ plot_dt, on='cell_id' ]

	return(plot_dt)
}

#' Plots the significant GO terms for each cluster
#'
#' @param go_results Output from GO analysis
#' @param sig_cutoff What is the minimum number of annotated genes to display a GO term?
#' @param p_cutoff What is the maximum p-value to display a GO term?
#' @return bar plot
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
plot_go_results <- function(go_list, sig_cutoff=5, p_cutoff=0.1) {
	# unpack
	go_results 	= go_list$go_results

	# set up
	plot_dt 	= go_results[ Significant>=sig_cutoff & p_go<p_cutoff ]
	data.table::setorder(plot_dt, cluster, -p_go)
	plot_dt[, N := .N, by=Term ]
	plot_dt[	  , term_n := Term ]
	plot_dt[ N > 1, term_n := paste0(Term, '_', 1:.N), by=Term ]
	plot_dt[, term_n := factor(term_n, levels=plot_dt$term_n) ]

	# plot
	g = ggplot(plot_dt) +
		aes( x=term_n, y=-log10(p_go) ) +
		geom_col() +
		scale_y_continuous( breaks=scales::pretty_breaks() ) +
		facet_grid( cluster ~ ., scales='free_y', space='free_y') +
		coord_flip() +
		labs(
			x 	= NULL
			,y 	= '-log10( p-value )'
			) +
		theme_bw()
	return(g)
}

#' Plot heatmap of gene clusters
#'
#' @param go_list Output from GO analysis
#' @return ggplot object
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
plot_heatmap_of_gene_clusters <- function(go_list) {
	# unpack
	plot_dt 	= go_list$plot_dt

	# plot
	g = ggplot(plot_dt) +
		aes( x=cell_id, y=symbol, fill=value ) +
		geom_tile() +
		scale_fill_distiller( palette='RdBu', limits=c(-3, 3) ) +
		facet_grid( clust_label ~ ., scale='free_y', space='free_y' ) +
		theme_bw() +
		theme(
			axis.text 	 = element_blank()
			,axis.ticks  = element_blank()
			) +
		labs(
			x 		= 'Cell'
			,y 		= 'Symbol'
			,fill 	= 'z-scored gene\nexpression'
			)
	return(g)
}

#' Plot heatmap of gene clusters
#'
#' @param go_list Output from GO analysis
#' @param label_name Description for the ordered labels in the legend (e.g. 'Developmental stage (days)')
#' @param palette RColorBrewer palette to use
#' @return ggplot object
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
plot_profiles_of_gene_clusters <- function(go_list, label_name='Ordered labels', palette='RdBu') {
	# unpack
	plot_dt 	= go_list$plot_dt
	cuts_dt 	= go_list$cuts_dt

	# set up what to plot
	means_dt 	= plot_dt[, list(value=mean(value)), by=list(psuper, clust_label)]

	# make nice colours
	col_vals 	= make_col_vals(cuts_dt$label_input, palette)

	# plot
	g = ggplot(means_dt) +
		aes( x=psuper, y=value ) +
		geom_vline(data=cuts_dt, aes(xintercept=psuper, colour=label_input)) +
		scale_colour_manual( values=col_vals ) +
		geom_smooth( colour='black', span=0.2, method='loess' ) +
		facet_grid( clust_label ~ ., scales='free_y' ) +
		theme_bw() +
		theme(
			axis.text 	= element_blank()
			) +
		labs(
			x 		= 'psupertime'
			,y 		= 'z-scored gene expression'
			,colour = label_name
			)
	return(g)
}
