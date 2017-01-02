# first part: this is how this is done for the RNAcontext paper:
# A) using predictions on training data:
# * select best parameters for each combination type,protein,set
# * compute mean apr for each combination type,protein
# * select full/weak type according to mean apr
# B) predictions on test data:
# * use selected models from selected type to report prediction score on test data
#
# second part: don't do selection of input data on training data, just select a best model from within full and from within weak

library(plyr)

# best for my own computations:
d_apr <- read.table("results_aucpr.csv", col.names=c('pred','protein','type','set','length','apr'));
d_correlation <- read.table("results_correlation.csv", col.names=c('pred','protein','type','set','length','correlation','pval'));

evaluate_apr <- function(d) {
	# select best models using training data
	selection_best_model_train <- subset(ddply(subset(d, pred=='train'), .(protein,type,set), transform, max_apr=max(apr)), apr==max_apr)[,c(2,3,4,5,6)]

	# select best input type using training data
	selection_best_type_train <- subset(ddply(ddply(selection_best_model_train, .(protein,type), transform, mean_apr=mean(apr)), .(protein), transform, max_mean_apr=max(mean_apr)), mean_apr==max_mean_apr)[,c(1,2,3,4)]

	# use selected models of selected type to extract results from test sets
	selection_models_test <- join(selection_best_type_train, subset(d, pred=='test'))

	# compute final results
	selection_models_test_result <- ddply(selection_models_test, .(protein), summarize, mean_apr=mean(apr))

	# save
	write.table(selection_models_test, "selection_models_test.csv", quote=F, sep="\t", row.names=F)
	write.table(selection_models_test_result, "selection_models_test_result.csv", quote=F, sep="\t", row.names=F)


	# second part: use all sets

	# get best models using training data
	all_selection_best_models_train <- subset(ddply(subset(d, pred=='train'), .(protein,type,set), transform, max_apr=max(apr)), apr==max_apr)[,c(2,3,4,5)]
	# use selected models to extract results from test sets
	all_selection_models_test <- join(all_selection_best_models_train, subset(d, pred=='test'))
	# compute results
	all_selection_models_test_result <- ddply(all_selection_models_test, .(protein,type), summarize, mean_apr=mean(apr))
	# save
	write.table(all_selection_models_test, 'all_selection_models_test.csv', quote=F, sep="\t", row.names=F)
	write.table(all_selection_models_test_result, 'all_selection_models_test_result.csv', quote=F, sep="\t", row.names=F)
}

evaluate_correlation <- function(d) {
 	# select best models using training data
	selection_best_model_train <- subset(ddply(subset(d, pred=='train'), .(protein,type,set), transform, max_correlation=max(correlation)), correlation==max_correlation)[,c(2,3,4,5,6)]

	# select best input type using training data
	selection_best_type_train <- subset(ddply(ddply(selection_best_model_train, .(protein,type), transform, mean_correlation=mean(correlation)), .(protein), transform, max_mean_correlation=max(mean_correlation)), mean_correlation==max_mean_correlation)[,c(1,2,3,4)]

	# use selected models of selected type to extract results from test sets
	selection_models_test <- join(selection_best_type_train, subset(d, pred=='test'))

	# compute final results
	selection_models_test_result <- ddply(selection_models_test, .(protein), summarize, mean_correlation=mean(correlation))

	# save
	write.table(selection_models_test, "selection_models_test_correlation.csv", quote=F, sep="\t", row.names=F)
	write.table(selection_models_test_result, "selection_models_test_result_correlation.csv", quote=F, sep="\t", row.names=F)


	# second part: use all sets

	# get best models using training data
	all_selection_best_models_train <- subset(ddply(subset(d, pred=='train'), .(protein,type,set), transform, max_correlation=max(correlation)), correlation==max_correlation)[,c(2,3,4,5)]
	# use selected models to extract results from test sets
	all_selection_models_test <- join(all_selection_best_models_train, subset(d, pred=='test'))
	# compute results
	all_selection_models_test_result <- ddply(all_selection_models_test, .(protein,type), summarize, mean_correlation=mean(correlation))
	# save
	write.table(all_selection_models_test, 'all_selection_models_test_correlation.csv', quote=F, sep="\t", row.names=F)
	write.table(all_selection_models_test_result, 'all_selection_models_test_result_correlation.csv', quote=F, sep="\t", row.names=F)
}

evaluate_apr(d_apr)
evaluate_correlation(d_correlation)
