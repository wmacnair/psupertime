# psupertime_plots.R

#' Convenience function to do multiple plots
#'
#' @param psuper_obj Psupertime object, output from psupertime
#' @param output_dir Directory to save to
#' @param tag Label for all files
#' @param ext Image format for outputs, compatible with ggsave (eps, ps, tex, pdf, jpeg, tiff, png, bmp, svg, wmf)
#' @param org_mapping Organism to use for annotation in GO enrichment analysis
#' @export
psupertime_plot_all <- function(psuper_obj, output_dir='.', tag='', ext='png', org_mapping='org.Mm.eg.db') {
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

	g 			= plot_predictions_against_classes(psuper_obj)
	plot_file 	= file.path(output_dir, sprintf('%s predictions over psupertime, original data.%s', tag, ext))
	ggplot2::ggsave(plot_file, g, height=6, width=10)

	# do GO term analysis if we can
	go_dt 		= psupertime_go_analysis(psuper_obj, org_mapping=org_mapping)
	go_file 	= file.path(output_dir, sprintf('%s go analysis.txt', tag))
	data.table::fwrite(go_dt, file=go_file)
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
			,title 	= sprintf('%s used for model selection', params$score)
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
		levels(doubles_dt[[ paste0('label_psuper_', label) ]]) 	= levels(psuper_list$y)
	}

	# put into list
	double_obj 	= list(
		run_names 		= run_names
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

#' GO enrichment analysis for genes learned from different psupertimes
#'
#' @param psuper_obj A previously calculated psupertime object
#' @param org_mapping Organism to use for annotations
#' @return data.table containing results of GO enrichment analysis
#' @export
psupertime_go_analysis <- function(psuper_obj, org_mapping='org.Mm.eg.db') {
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
		go_temp 	= data.table(topGO::GenTable(topGO_data, 
			weightFisher 	= go_weight, 
			orderBy 		= 'weightFisher', 
			ranksOf 		= 'weightFisher', 
			topNodes 		= 1000
			))
		data.table::setnames(go_temp, 'weightFisher', 'temp')
		go_temp[, weightFisher := as.numeric(temp) ]
		go_temp[ temp == '< 1e-30', weightFisher := 9e-31 ]
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
	print_dt 		= go_dt[ weightFisher < p_cutoff & Significant>n_terms_cutoff ]
	if (nrow(print_dt)==0) {
		message(sprintf('no GO terms met the cutoffs (p-value < %.1e and at least %d genes significant)', p_cutoff, n_terms_cutoff))
	} else {
		message('Significant GO terms:')
		print(print_dt)
	}

	return(go_dt)
}
