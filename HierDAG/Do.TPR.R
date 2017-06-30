# High level function to compute hierarchical correction according to TPR-threshold algorithm. 

##*********************************************************************##
# Copyright (C) 2017 Marco Notaro

# This file is part of HierDAG library. 

# HierDAG is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
##*********************************************************************##

# March 2016	  
# April 2016: all high level function modified: some input parameters of find.best.f function added
# April 2016: all high level function modified: normaliation method added
# July 2016: added PRC computed by precrec package (only in holdout)
# July 2016: added Do.tpr.'tpr_variants'.holdout
# November 2016: added b.per.example as input parameter 
# December 2016: added PRC computed by precrec pkg to each Do.tpr.'tpr_variants'.cv
# December 2016: added AUPRC.single.over.classes and AUROC.single.over.classes for computing AUC/PRC by precrec also with class with zero annotations,
# These function are useful to compute AUC/PRC in a contest of k-fold cross-validation, where a fold can have a class with zero annotations..
# March 2017: figure RAM memory prloblem out!! see AUPROC.R source file for explanation..

##***************************************************************************************##
# library to call 
library("PerfMeas"); 				## compute PXR
library("precrec");					## compute AUPRC and AUROC
library("preprocessCore"); 			## Qnorm
source("tpr.R");					## TPR-DAG algorithm
source("flat.score.norm.R");		## Maxnorm
source("graph.utils.R");			## graph utility functions 
source("F.hier.R");					## compute Kiritchenko-like multi-label F-scores
source("Do.flat.normalization.R");	## high level function to compute Maxnorm e Qnorm
source("AUPROC.R");					## functions to compute AUPRC and AUROC 
##***************************************************************************************##

# It perform a k fold cross-validation to find the best threshold maximizing on F.max measure.
# INPUT: 
# threshold: range of threshold values to be tested in order to find the best threshold. The denser the range is, the higher the probability to find the best
# 			 theshold, but the execution time will be increasing (def: from:0.1, to:0.9, by:0.1 step).
# kk:  number of folds of the cross validation (def: 5).
# seed:  intialization seed for the random generator to create folds (def:0). If NULL (default) no initialization is performed.
# norm: boolean value: 1.TRUE means that the flat scores matrix has been already normalized in according to a normalization method (def);
#       			   2.FALSE means that the flat scores matrix has NOT been normalized yet.
# norm.type: this variable can assume three values: MaxNorm, Qnorm, NONE. We have two case respect to norm:
#			 1.if norm==FALSE, two kind of normalizations are possible:	1. MaxNorm: each score is divided w.r.t. the max of each class 
#																		2. Qnorm: quantile normalization is applied. Library preprocessCore is used.
#			 2.if norm==TRUE, set norm.type=="NONE" (def);
# flat.file: name of flat scores matrix already normalized or to be normalized in according to norm.type (without rda extension); 
# ann.file: name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples (def: ann.file)
# dag.file: name of the graph that represents the hierarchy of the classes. 
# flat.dir: relative path to folder where flat normalized scores matrix is stored
# ann.dir: relative path to folder where annotation matrix is stored
# dag.dir: relative path to folder where graph is stored
# flat.norm.dir: 1.if norm=FALSE, relative path where flat normalized scores matrix is strored;
#				 2.if norm=TRUE, the flat scores matrix is already normalized, than it is set to NULL (def)
# b.per.example : if TRUE (def.) precision, recall, F-measure, specificity and accuracy are returned for each example, 
#                 otherwise only the averages are returned.
# n.round: 	number of rounding digits to be applied to the hierarchical scores matrix (def. 3). 
#			It's used for choosing the best threshold on the basis of the best F.measure (see f.criterion parameter).
# f.criterion: character. Type of F-measure to be used to select the best F.measure. There are 2 possibilities: 
#              1. "F" (default) corresponds to the harmonic mean between the average precision and recall; 
#              2. "avF" corresponds to the per-example F-score averaged across all the examples.
# hierScore.dir: relative path to folder where the matrix with the scores of the classes corrected according to TPR-threshold algorithm is stored 
# perf.dir: relative path to folder where the macro-centric measures (i.e. AUC, PRC and PxR across classes) and 
# 			where example-centric measures (i.e. Precision, Recall, Specificity, F-measure, Accuracy across example) are stored
# OUTPUT:
# 5 rda files stored in the corresponding output directory:
# - Matrix with examples on rows and classes on colums representing the hierarchical scores of the classes computed with TPR-threshold algorithm.
#	Stored in hierScore.dir folder
# - Example-centric measures computed through find.best.f from F-hier.R file. Stored in perf.dir folder
# - AUC (average and per classes) computed through AUC.single.over.classes from package PerfMeas. Stored in perf.dir
# - Precision at fixed recall levels (average and per classes) computed through precision.at.multiple.recall.level.over.classes from package PerfMeas.
# 	Stored in perf.dir folder 
# - PRC (average and per.class) computed by precrec package. Stored in perf.dir
Do.tpr.threshold.cv <- function(	threshold=seq(from=0.1, to=0.9, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type= "NONE",
									flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, 
									ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion ="F", 
									b.per.example=TRUE, hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"
								){
	## Loading Data ############
	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); 
	
	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); #saving PRC result in the same format of package PerfMeas

	## Hierarchical Correction #########################
	## TPR-Threshold correction using k-fold cross-validation to find best threshold 

	## splitting data in k unstratified fold
	folds <- do.unstratified.cv.data(S, kk=kk, seed=seed); 
	
	## storing the best F.max and the best threshold for each of k training set
	training.top.Fmax <- vector(mode="list", length=kk);
	names(training.top.Fmax) <- paste0(rep("fold",kk), 1:kk);

	## storing all the best average protein-centric measures (average and per.example) for each of k test set 
	best.avg.meas.test <- vector(mode="list", length=kk);
	names(best.avg.meas.test) <- paste0(rep("fold",kk), 1:kk);
	FMM.per.example <- c();

	## storing average macro AUC for each of k test set 
	AUC.average.test <- vector(mode="list", length=kk);
	names(AUC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	AUC.class.test <- vector(mode="list", length=kk);
	names(AUC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average PxR for each of k test set 
	PXR.average.test <- vector(mode="list", length=kk);
	names(PXR.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class PxR for each of k test set 
	PXR.class.test <- vector(mode="list", length=kk);
	names(PXR.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average macro AUC for each of k test set 
	PRC.average.test <- vector(mode="list", length=kk);
	names(PRC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	PRC.class.test <- vector(mode="list", length=kk);
	names(PRC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## variable for hosting the k sub-matrix assembled 
	S.tpr <- c();

	## Let's start k-fold crossing validation for choosing best threshold maximizing on F.max...
	for(k in 1:kk){
		test <- S[folds[[k]],];		# test set: 1 out of k folds (e.g.: if k=5, testing set is 1/5 of data)
		training <- S[!rownames(S) %in% rownames(test),];		# training set: (k-1) out of k folds (e.g.: if k=5, training set is 4/5 of data)
		
		top.Fmax <- 0;
		best.Fmaxt <- 0;
		for(t in threshold){
			pred.training <- tpr.threshold(training, g, root=root, t=t);	## training set hierarchical correction...
			target.training <- ann[rownames(pred.training),colnames(pred.training)];
			training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=FALSE); 
			training.Fmax <- training.avg.meas[4];	## F.max maximization...
			if(training.Fmax > top.Fmax){
				top.Fmax <- training.Fmax;
				best.Fmaxt <- t;
				training.top.Fmax[[k]] <- c(best.Fmax=top.Fmax, best.thres=best.Fmaxt);
				cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, sep="\t", "\n");
			}
		}
		pred.test <- tpr.threshold(test, g, root=root, t=training.top.Fmax[[k]][2]);
		target.test <- ann[rownames(pred.test),colnames(pred.test)];
		test.avg.meas <- find.best.f(target.test, pred.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);	
		best.avg.meas.test[[k]] <- test.avg.meas$average;
		FMM.per.example <- rbind(FMM.per.example,test.avg.meas$per.example);

		## AUC (average and per.class) computed with PerfMeas package
		AUC.average.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$average;	## storing average AUC of each k testing set
		AUC.class.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$per.class;	## storing AUC per.class of each k testing set

		## PxR at fixed recall levels (average and per.class) computed with PerfMeas package
		PXR.average.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$avgPXR;	## storing average PxR of each k testing set
		PXR.class.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$PXR;	## storing PxR per.class of each k testing set

		## PRC (average and per.class) computed with precrec package
		PRC.average.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$average;	## storing average PRC of each k testing set
		PRC.class.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$per.class;	## storing PRC per.class of each k testing set

		## assembling of all the hierarchical scores of each k sub-matrix to build the full matrix 
		S.tpr <- rbind(S.tpr, pred.test);
	}
	## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the full beginning matrix
	S.tpr <- S.tpr[rownames(S),];

	## remove no longer useful variables..
	rm(S, folds, pred.test); gc();

	## Averaging Performances Measures across k testing set ###################
	## averaging protein-centric measures across k testing sets
	F.cv <- Reduce("+", best.avg.meas.test)/kk;
	F.meas <-  apply(FMM.per.example,2,mean);
	names(F.meas)[4] <- "avF";
	F.max <- 2*(F.meas[["P"]] * F.meas[["R"]])/((F.meas[["P"]] + F.meas[["R"]]));
	names(F.max) <- "F";
	FMM.avg <- append(F.meas, F.max, after=3);
	FMM.avg <- append(FMM.avg,F.cv["T"]);
	FMM.tpr <- list(average=FMM.avg, per.example=FMM.per.example);

	## averaging AUC (average and per.class) across k testing sets 
	AUC.average.tpr.over.test <- Reduce("+", AUC.average.test)/kk;
	AUC.class.tpr.over.test <- Reduce("+", AUC.class.test)/kk;
	AUC.tpr <- list(average=AUC.average.tpr.over.test, per.class=AUC.class.tpr.over.test);

	## averaging PxR (average and per.class) across k testing sets 
	PXR.average.tpr.over.test <- Reduce("+", PXR.average.test)/kk;
	PXR.class.tpr.over.test <- Reduce("+", PXR.class.test)/kk;
	PXR.tpr <- list(average=PXR.average.tpr.over.test, per.class=PXR.class.tpr.over.test);

	## averaging PRC (average and per.class) across k testing sets 
	PRC.average.tpr.over.test <- Reduce("+", PRC.average.test)/kk;
	PRC.class.tpr.over.test <- Reduce("+", PRC.class.test)/kk;
	PRC.tpr <- list(average=PRC.average.tpr.over.test, per.class=PRC.class.tpr.over.test);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprT.rda"), compress=TRUE);
	}
}

# TPR-threshold Hold-Out version
# High level function to correct the scores with a hierarchy according to TPR-threshold algorithm performing a classical holdout procedure
# All input parameters are the same of the above corresponding function that implement a kfcv, except the following:
# ind.test.set: vector of integer. Indices refer to the examples of the adjancency matrix to be used in the test set. 
# ind.dir: relative path to folder where ind.test.set is stored
Do.tpr.threshold.holdout <- function(	threshold=seq(from=0.1, to=0.9, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type= "NONE",
										flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set,
										ind.dir=ind.dir,flat.dir=flat.dir, ann.dir=ann.dir,dag.dir=dag.dir, flat.norm.dir=NULL, 
										n.round=3, f.criterion ="F", b.per.example=TRUE, hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"
									){
	## Loading Data ############
	## loading examples indices of the test set
	ind.set <- paste0(ind.dir, ind.test.set, ".rda");
	ind.test <- get(load(ind.set));

	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## scores flat matrix shrinked to test and training test respectively
	S.test <- S[ind.test,];
	S.training <- S[-ind.test,];
	rm(S);

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table and shrinking annotation table to tet set
	ann.test <- ann[ind.test,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann.test, S.test); 

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann.test, S.test); 

	## Hierarchical Correction #########################
	## TPR-Threshold correction using k-fold cross-validation to find best threshold on training set
	folds <- do.unstratified.cv.data(S.training, kk=kk, seed=seed);
			
	## Let's start k-fold crossing validation for choosing best threshold maximizing on the basis of F.max value...
	for(k in 1:kk){
		training <- S.training[folds[[k]],];
		top.Fmax <- 0;
		best.Fmaxt <- 0;
		for(t in threshold){
			pred.training <- tpr.threshold(training, g, root=root, t=t);	## training set hierarchical correction...
			target.training <- ann[rownames(pred.training),colnames(pred.training)];
			training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE);
			training.Fmax <- training.avg.meas[4];	## F.max maximization...
			if(training.Fmax > top.Fmax){
				top.Fmax <- training.Fmax;
				best.Fmaxt <- t;
				cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, sep="\t", "\n");
			}
		}
	}
	S.test <- tpr.threshold(S.test, g, root=root, t=best.Fmaxt);	## testing set hierarchical correction..
		
	## AUC (average and per.class) computed with PerfMeas package on the test set
	AUC.tpr <- AUROC.single.over.classes(ann.test, S.test);

	## PxR at fixed recall levels (average and per.class) computed with PerfMeas package on the test set
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);	

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann.test, S.test); 

	## storing the hierarchical matrix
	S.tpr <- S.test;
	rm(S.test, S.training); gc();

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.holdout.tprT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.holdout.tprT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.holdout.tprT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.holdout.tprT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.holdout.tprT.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.holdout.tprT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.holdout.tprT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.holdout.tprT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.holdout.tprT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.holdout.tprT.rda"), compress=TRUE);
	}
}


# High level function to compute hierarchical correction according to TPR-weighted-threshold-free algorithm. 
# It perform a k fold cross-validation to find the best weight maximizing on F.max measure.
# INPUT: 
# weight: range of weight values to be tested in order to find the best weight. The denser the range is, the higher the probability to find the best
# 		  weigth, but the execution time will be increasing (def: from:0.1, to:1, by:0.1 step).
# kk:  number of folds of the cross validation (def: 5).
# seed:  intialization seed for the random generator to create folds (def:0). If NULL (default) no initialization is performed.
# norm: boolean value: 1.TRUE means that the flat scores matrix has been already normalized in according to a normalization method (def);
#       			   2.FALSE means that the flat scores matrix has NOT been normalized yet.
# norm.type: this variable can assume three values: MaxNorm, Qnorm, NONE. We have two case respect to norm:
#			 1.if norm==FALSE, two kind of normalizations are possible:	1. MaxNorm: each score is divided w.r.t. the max of each class 
#																		2. Qnorm: quantile normalization is applied. Library preprocessCore is used.
#			 2.if norm==TRUE, set norm.type=="NONE" (def);
# flat.file: name of flat scores matrix already normalized or to be normalized in according to norm.type (without rda extension); 
# ann.file: name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples (def: ann.file)
# dag.file: name of the graph that represents the hierarchy of the classes. 
# flat.dir: relative path to folder where flat normalized scores matrix is stored
# ann.dir: relative path to folder where annotation matrix is stored
# dag.dir: relative path to folder where graph is stored
# flat.norm.dir: 1.if norm=FALSE, relative path where flat normalized scores matrix is strored;
#				 2.if norm=TRUE, the flat scores matrix is already normalized, than it is set to NULL (def)
# b.per.example : if TRUE (def.) precision, recall, F-measure, specificity and accuracy are returned for each example, 
#                 otherwise only the averages are returned.
# n.round: 	number of rounding digits to be applied to the hierarchical scores matrix (def. 3). 
#			It's used for choosing the best threshold on the basis of the best F.measure (see f.criterion parameter).
# f.criterion: character. Type of F-measure to be used to select the best F.measure. There are 2 possibilities: 
#              1. "F" (default) corresponds to the harmonic mean between the average precision and recall; 
#              2. "avF" corresponds to the per-example F-score averaged across all the examples.
# hierScore.dir: relative path to folder where the matrix with the scores of the classes corrected in according to TPR-weighted-threshold-free algorithm is stored 
# perf.dir: relative path to folder where the macro-centric measures (i.e. AUC, PRC and PxR across classes) and 
# 			where example-centric measures (i.e. Precision, Recall, Specificity, F-measure, Accuracy across example) are stored
# OUTPUT:
# 5 rda files stored in the rispective output directory:
# - Matrix with examples on rows and classes on colums representing the hierarchical scores of the classes computed with TPR-weighted-threshold-free algorithm.
#	Stored in hierScore.dir folder
# - Example-centric measures computed through find.best.f from F-hier.R file. Stored in perf.dir folder
# - AUC (average and per classes) computed through AUC.single.over.classes from package PerfMeas. Stored in perf.dir
# - Precision at fixed recall levels (average and per classes) computed through precision.at.multiple.recall.level.over.classes from package PerfMeas.
# 	Stored in perf.dir folder 
# - PRC (average and per.class) computed by precrec package. Stored in perf.dir
Do.tpr.weighted.threshold.free.cv <- function(	weight=seq(from=0.1, to=1, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type= "NONE",
												flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, 
												ann.dir=ann.dir,dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion ="F", 
												b.per.example=TRUE, hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"
											){
	## Loading Data ############
	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); gc();

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S); gc();

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example); gc();

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); gc();

	## Hierarchical Correction #########################
	## TPR-Weighted-Threshold-Free correction using k-fold cross-validation to find best weight 

	## splitting data in k stratified fold
	folds <- do.unstratified.cv.data(S, kk=kk, seed=seed); 

	## storing the best F.max and the best threshold for each of k training set
	training.top.Fmax <- vector(mode="list", length=kk);
	names(training.top.Fmax) <- paste0(rep("fold",kk), 1:kk);

	## storing all the best protein-centric measures for each of k test set 
	best.avg.meas.test <- vector(mode="list", length=kk);
	names(best.avg.meas.test) <- paste0(rep("fold",kk), 1:kk);
	FMM.per.example <- c();

	## storing average macro AUC for each of k test set 
	AUC.average.test <- vector(mode="list", length=kk);
	names(AUC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	AUC.class.test <- vector(mode="list", length=kk);
	names(AUC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average PxR for each of k test set 
	PXR.average.test <- vector(mode="list", length=kk);
	names(PXR.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class PxR for each of k test set 
	PXR.class.test <- vector(mode="list", length=kk);
	names(PXR.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average macro AUC for each of k test set 
	PRC.average.test <- vector(mode="list", length=kk);
	names(PRC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	PRC.class.test <- vector(mode="list", length=kk);
	names(PRC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## variable for hosting the k sub-matrix assembled 
	S.tpr <- c();

	## Let's start k-fold crossing validation for choosing best threshold maximizing on F.max...
	for(k in 1:kk){
		test <- S[folds[[k]],];		# test set: 1 out of k folds (e.g.: if k=5, testing set is 1/5 of data)
		training <- S[!rownames(S) %in% rownames(test),];		# training set: (k-1) out of k folds (e.g.: if k=5, training set is 4/5 of data)
		gc();

		top.Fmax <- 0;
		best.Fmaxw <- 0;
		for(w in weight){
			pred.training <- tpr.weighted.threshold.free(training, g, root=root, w=w);	## training set hierarchical correction...
			target.training <- ann[rownames(pred.training),colnames(pred.training)];
			training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE);	 
			training.Fmax <- training.avg.meas[4];	## F.max maximization...
			if(training.Fmax > top.Fmax){
				top.Fmax <- training.Fmax;
				best.Fmaxw <- w;
				training.top.Fmax[[k]] <- c(best.Fmax=top.Fmax, best.weigh=best.Fmaxw);
				cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best weight:",best.Fmaxw, sep="\t", "\n");
			}
		}
		pred.test <- tpr.weighted.threshold.free(test, g, root=root, w=training.top.Fmax[[k]][2]);	## testing set hierarchical correction...
		target.test <- ann[rownames(pred.test),colnames(pred.test)];
		test.avg.meas <- find.best.f(target.test, pred.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);	
		best.avg.meas.test[[k]] <- test.avg.meas$average;
		FMM.per.example <- rbind(FMM.per.example,test.avg.meas$per.example);

		## Hierarchical AUC (average and per.class) computed with precrec package
		AUC.average.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$average;	## storing average AUC of each k testing set
		AUC.class.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$per.class;	## storing AUC per.class of each k testing set

		## Hierarchical PxR at fixed recall levels (average and per.class) computed with PerfMeas package
		PXR.average.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$avgPXR;	## storing average PxR of each k testing set
		PXR.class.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$PXR;	## storing PxR per.class of each k testing set

		## PRC (average and per.class) computed with precrec package
		PRC.average.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$average;	## storing average PRC of each k testing set
		PRC.class.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$per.class;	## storing PRC per.class of each k testing set

		## assembling of all the hierarchical scores of each k sub-matrix to build the full matrix 
		S.tpr <- rbind(S.tpr, pred.test);
	}
	## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the full beginning matrix
	S.tpr <- S.tpr[rownames(S),];

	## remove no longer useful variables..
	rm(S, folds, pred.test); gc();

	## Averaging Performances Measures across k testing set ###################
	## averaging protein-centric measures across k testing sets
	F.cv <- Reduce("+", best.avg.meas.test)/kk;
	F.meas <-  apply(FMM.per.example,2,mean);
	names(F.meas)[4] <- "avF";
	F.max <- 2*(F.meas[["P"]] * F.meas[["R"]])/((F.meas[["P"]] + F.meas[["R"]]));
	names(F.max) <- "F";
	FMM.avg <- append(F.meas, F.max, after=3);
	FMM.avg <- append(FMM.avg,F.cv["T"]);
	FMM.tpr <- list(average=FMM.avg, per.example=FMM.per.example);

	## averaging macro-AUC (average and per.class) across k testing sets 
	AUC.average.tpr.over.test <- Reduce("+", AUC.average.test)/kk;
	AUC.class.tpr.over.test <- Reduce("+", AUC.class.test)/kk;
	AUC.tpr <- list(average=AUC.average.tpr.over.test, per.class=AUC.class.tpr.over.test);

	## averaging PxR (average and per.class) across k testing sets 
	PXR.average.tpr.over.test <- Reduce("+", PXR.average.test)/kk;
	PXR.class.tpr.over.test <- Reduce("+", PXR.class.test)/kk;
	PXR.tpr <- list(average=PXR.average.tpr.over.test, per.class=PXR.class.tpr.over.test);

	## averaging PRC (average and per.class) across k testing sets 
	PRC.average.tpr.over.test <- Reduce("+", PRC.average.test)/kk;
	PRC.class.tpr.over.test <- Reduce("+", PRC.class.test)/kk;
	PRC.tpr <- list(average=PRC.average.tpr.over.test, per.class=PRC.class.tpr.over.test);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
	}
}

# High level function to correct the scores with a hierarchy according to TPR-weighted-threshold-free algorithm performing a classical holdout procedure
# All input parameters are the same of the above corresponding function that implement a kfcv, except the following:
# ind.test.set: vector of integer. Indices refer to the examples of the adjancency matrix to be used in the test set. 
# ind.dir: relative path to folder where ind.test.set is stored
Do.tpr.weighted.threshold.free.holdout <- function(	weight=seq(from=0.1, to=1, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type= "NONE",
													flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set,
													ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir,dag.dir=dag.dir, flat.norm.dir=NULL, 	 
													n.round=3, f.criterion ="F", b.per.example=TRUE, hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"
												){
	## Loading Data ############
	## loading examples indices of the test set
	ind.set <- paste0(ind.dir, ind.test.set, ".rda");
	ind.test <- get(load(ind.set));

	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## scores flat matrix shrinked to test and training test respectively
	S.test <- S[ind.test,];
	S.training <- S[-ind.test,];
	rm(S);

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table and shrinking annotation table to tet set
	ann.test <- ann[ind.test,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann.test, S.test); 
	
	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann.test, S.test); 
	
	## Hierarchical Correction #########################
	## TPR-Threshold correction using k-fold cross-validation to find best threshold on training set
	folds <- do.unstratified.cv.data(S.training, kk=kk, seed=seed);
			
	## Let's start k-fold crossing validation for choosing best threshold maximizing on the basis of F.max value...
	for(k in 1:kk){
		training <- S.training[folds[[k]],];	
		top.Fmax <- 0;
		best.Fmaxw <- 0;
		for(w in weight){
			pred.training <- tpr.weighted.threshold.free(training, g, root=root, w=w);	## training set hierarchical correction...
			target.training <- ann[rownames(pred.training),colnames(pred.training)];
			training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE); 
			training.Fmax <- training.avg.meas[4];	## F.max maximization...
			if(training.Fmax > top.Fmax){
				top.Fmax <- training.Fmax;
				best.Fmaxw <- w;
				cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxw, sep="\t", "\n");
			}
		}
	}
	S.test <- tpr.weighted.threshold.free(S.test, g, root=root, w=best.Fmaxw);	## testing set hierarchical correction..
		
	## AUC (average and per.class) computed with PerfMeas package on the test set
	AUC.tpr <- AUROC.single.over.classes(ann.test, S.test);	

	## PxR at fixed recall levels (average and per.class) computed with PerfMeas package on the test set
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann.test, S.test); 

	## storing the hierarchical matrix
	S.tpr <- S.test;
	rm(S.test, S.training);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprW.rda"), compress=TRUE);
	}
}


# High level function to compute hierarchical correction according to TPR-weighted-threshold algorithm. 
# It perform a k fold cross-validation to find the best threshold and weight maximizing on F.max measure.
# INPUT: 
# threshold: range of threshold values to be tested in order to find the best threshold. The denser the range is, the higher the probability to find the best
# 			 theshold, but the execution time will be increasing (def: from:0.1, to:0.9, by:0.1 step).
# weight: range of weight values to be tested in order to find the best weight. The denser the range is, the higher the probability to find the best
# 		  weigth, but the execution time will be increasing (def: from:0.1, to:1, by:0.1 step).
# kk:  number of folds of the cross validation (def: 5).
# seed:  intialization seed for the random generator to create folds (def:0). If NULL (default) no initialization is performed.
# norm: boolean value: 1.TRUE means that the flat scores matrix has been already normalized in according to a normalization method (def);
#       			   2.FALSE means that the flat scores matrix has NOT been normalized yet.
# norm.type: this variable can assume three values: MaxNorm, Qnorm, NONE. We have two case respect to norm:
#			 1.if norm==FALSE, two kind of normalizations are possible:	1. MaxNorm: each score is divided w.r.t. the max of each class 
#																		2. Qnorm: quantile normalization is applied. Library preprocessCore is used.
#			 2.if norm==TRUE, set norm.type=="NONE" (def);
# flat.file: name of flat scores matrix already normalized or to be normalized in according to norm.type (without rda extension); 
# ann.file: name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples (def: ann.file)
# dag.file: name of the graph that represents the hierarchy of the classes. 
# flat.dir: relative path to folder where flat normalized scores matrix is stored
# ann.dir: relative path to folder where annotation matrix is stored
# dag.dir: relative path to folder where graph is stored
# flat.norm.dir: 1.if norm=FALSE, relative path where flat normalized scores matrix is strored;
#				 2.if norm=TRUE, the flat scores matrix is already normalized, than it is set to NULL (def)
# b.per.example : if TRUE (def.) precision, recall, F-measure, specificity and accuracy are returned for each example, 
#                 otherwise only the averages are returned.
# n.round: 	number of rounding digits to be applied to the hierarchical scores matrix (def. 3). 
#			It's used for choosing the best threshold on the basis of the best F.measure (see f.criterion parameter).
# f.criterion: character. Type of F-measure to be used to select the best F.measure. There are 2 possibilities: 
#              1. "F" (default) corresponds to the harmonic mean between the average precision and recall; 
#              2. "avF" corresponds to the per-example F-score averaged across all the examples.
# hierScore.dir: relative path to folder where the matrix with the scores of the classes corrected in according to TPR-weighted-threshold algorithm is stored 
# perf.dir: relative path to folder where the macro-centric measures (i.e. AUC, PRC and PxR across classes) and 
# 			where example-centric measures (i.e. Precision, Recall, Specificity, F-measure, Accuracy across example) are stored
# OUTPUT:
# 5 rda files stored in the rispective output directory:
# - Matrix with examples on rows and classes on colums representing the hierarchical scores of the classes computed with TPR-weighted-threshold algorithm.
#	Stored in hierScore.dir folder
# - Example-centric measures computed through find.best.f from F-hier.R file. Stored in perf.dir folder
# - AUC (average and per classes) computed through AUC.single.over.classes from package PerfMeas. Stored in perf.dir
# - Precision at fixed recall levels (average and per classes) computed through precision.at.multiple.recall.level.over.classes from package PerfMeas.
# 	Stored in perf.dir folder 
# - PRC (average and per.class) computed by precrec package. Stored in perf.dir
Do.tpr.weighted.threshold.cv <- function(	threshold=seq(from=0.1, to=0.9, by=0.1), weight=seq(from=0.1, to=1, by=0.1), 
											kk=5, seed=NULL,  norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, 
											dag.file=dag.file, flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL,
											n.round=3,  b.per.example=TRUE, f.criterion ="F", hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"
										){
	## Loading Data ############
	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); 

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); 

	## Hierarchical Correction #########################
	## TPR-Weighted-Threshold correction using k-fold cross-validation to find best threshold and weight 

	## splitting data in k unstratified folds 
	folds <- do.unstratified.cv.data(S, kk=kk, seed=seed); 
	
	## storing the best F.max and the best threshold and weight for each of k training set
	training.top.Fmax <- vector(mode="list", length=kk);
	names(training.top.Fmax) <- paste0(rep("fold",kk), 1:kk);

	## storing all the best protein-centric measures for each of k test set 
	best.avg.meas.test <- vector(mode="list", length=kk);
	names(best.avg.meas.test) <- paste0(rep("fold",kk), 1:kk);
	FMM.per.example <- c();

	## storing average macro AUC for each of k test set 
	AUC.average.test <- vector(mode="list", length=kk);
	names(AUC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	AUC.class.test <- vector(mode="list", length=kk);
	names(AUC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average PxR for each of k test set 
	PXR.average.test <- vector(mode="list", length=kk);
	names(PXR.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class PxR for each of k test set 
	PXR.class.test <- vector(mode="list", length=kk);
	names(PXR.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average macro PRCC for each of k test set 
	PRC.average.test <- vector(mode="list", length=kk);
	names(PRC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro PRC for each of k test set 
	PRC.class.test <- vector(mode="list", length=kk);
	names(PRC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## variable for hosting the k sub-matrix assembled 
	S.tpr <- c();

	## Let's start k-fold crossing validation for choosing best threshold and weight maximizing on F.max...
	for(k in 1:kk){
		test <- S[folds[[k]],];	# test set: 1 out of k folds (e.g.: if k=5, testing set is 1/5 of data)
		training <- S[!rownames(S) %in% rownames(test),];		# training set: (k-1) out of k folds (e.g.: if k=5, training set is 4/5 of data)
		gc();

		top.Fmax <- 0;
		best.Fmaxt <- 0;
		best.Fmaxw <- 0;
		for(t in threshold){
			for(w in weight){
				pred.training <- tpr.weighted.threshold(training, g, root=root, w=w, t=t);	## training set hierarchical correction...
				target.training <- ann[rownames(pred.training),colnames(pred.training)];
				training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE);	 
				training.Fmax <- training.avg.meas[4];	## F.max maximization...
				if(training.Fmax > top.Fmax){
					top.Fmax <- training.Fmax;
					best.Fmaxt <- t;
					best.Fmaxw <- w;
					training.top.Fmax[[k]] <- c(best.Fmax=top.Fmax, best.thres=best.Fmaxt, best.weight=best.Fmaxw);
					cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, "best weight:", best.Fmaxw,sep="\t", "\n");
				}
			}
		}
		pred.test <- tpr.weighted.threshold(test, g, root=root, t=training.top.Fmax[[k]][2], w=training.top.Fmax[[k]][3]);
		target.test <- ann[rownames(pred.test),colnames(pred.test)];
		test.avg.meas <- find.best.f(target.test, pred.test, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=b.per.example);	
		best.avg.meas.test[[k]] <- test.avg.meas$average;
		FMM.per.example <- rbind(FMM.per.example,test.avg.meas$per.example);

		## Hierarchical AUC (average and per.class) computed with PerfMeas package
		AUC.average.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$average;	## storing average AUC of each k testing set
		AUC.class.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$per.class;	## storing AUC per.class of each k testing set

		## Hierarchical PxR at fixed recall levels (average and per.class) computed with PerfMeas package
		PXR.average.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$avgPXR;	## storing average PxR of each k testing set
		PXR.class.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$PXR;	## storing PxR per.class of each k testing set

		## PRC (average and per.class) computed with precrec package
		PRC.average.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$average;	## storing average AUC of each k testing set
		PRC.class.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$per.class;	## storing AUC per.class of each k testing set

		## assembling of all the hierarchical scores of each k sub-matrix to build the full matrix 
		S.tpr <- rbind(S.tpr, pred.test);
	}
	## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the full beginning matrix
	S.tpr <- S.tpr[rownames(S),];

	## remove no longer useful variables..
	rm(S, folds, pred.test);

	## Averaging Performances Measures across k testing set ###################
	## averaging protein-centric measures across k testing sets
	F.cv <- Reduce("+", best.avg.meas.test)/kk;
	F.meas <-  apply(FMM.per.example,2,mean);
	names(F.meas)[4] <- "avF";
	F.max <- 2*(F.meas[["P"]] * F.meas[["R"]])/((F.meas[["P"]] + F.meas[["R"]]));
	names(F.max) <- "F";
	FMM.avg <- append(F.meas, F.max, after=3);
	FMM.avg <- append(FMM.avg,F.cv["T"]);
	FMM.tpr <- list(average=FMM.avg, per.example=FMM.per.example);

	## averaging macro-AUC (average and per.class) across k testing sets 
	AUC.average.tpr.over.test <- Reduce("+", AUC.average.test)/kk;
	AUC.class.tpr.over.test <- Reduce("+", AUC.class.test)/kk;
	AUC.tpr <- list(average=AUC.average.tpr.over.test, per.class=AUC.class.tpr.over.test);

	## averaging PxR (average and per.class) across k testing sets 
	PXR.average.tpr.over.test <- Reduce("+", PXR.average.test)/kk;
	PXR.class.tpr.over.test <- Reduce("+", PXR.class.test)/kk;
	PXR.tpr <- list(average=PXR.average.tpr.over.test, per.class=PXR.class.tpr.over.test);

	## averaging PRC (average and per.class) across k testing sets 
	PRC.average.tpr.over.test <- Reduce("+", PRC.average.test)/kk;
	PRC.class.tpr.over.test <- Reduce("+", PRC.class.test)/kk;
	PRC.tpr <- list(average=PRC.average.tpr.over.test, per.class=PRC.class.tpr.over.test);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
	}
}

# Hold-Out version
# High level function to correct the scores with a hierarchy according to TPR-weighted-threshold algorithm performing a classical holdout procedure
# All input parameters are the same of the corresponding function above, except the following:
# ind.test.set: vector of integer. Indices refer to the examples of the adjancency matrix to be used in the test set. 
# ind.dir: relative path to folder where ind.test.set is stored
Do.tpr.weighted.threshold.holdout <- function(	threshold=seq(from=0.1, to=0.9, by=0.1), weight=seq(from=0.1, to=1, by=0.1),
												kk=5, seed=NULL, norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, 
												dag.file=dag.file, ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, 
												ann.dir=ann.dir,dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion ="F", 
												b.per.example=TRUE, hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"										
											){
	## Loading Data ############
	## loading examples indices of the test set
	ind.set <- paste0(ind.dir, ind.test.set, ".rda");
	ind.test <- get(load(ind.set));

	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## scores flat matrix shrinked to test and training test respectively
	S.test <- S[ind.test,];
	S.training <- S[-ind.test,];
	rm(S);

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table and shrinking annotation table to tet set
	ann.test <- ann[ind.test,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann.test, S.test); 

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann.test, S.test); 

	## Hierarchical Correction #########################
	## TPR-Threshold correction using k-fold cross-validation to find best threshold on training set
	folds <- do.unstratified.cv.data(S.training, kk=kk, seed=seed);
			
	## Let's start k-fold crossing validation for choosing best threshold maximizing on the basis of F.max value...
	for(k in 1:kk){
		training <- S.training[folds[[k]],];	
		
		top.Fmax <- 0;
		best.Fmaxt <- 0;
		best.Fmaxw <- 0;
		for(t in threshold){
			for(w in weight){
				pred.training <- tpr.weighted.threshold(training, g, root=root, w=w, t=t); ## training set hierarchical correction...
				target.training <- ann[rownames(pred.training),colnames(pred.training)];
				training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE); 
				training.Fmax <- training.avg.meas[4];	## F.max maximization...
				if(training.Fmax > top.Fmax){
					top.Fmax <- training.Fmax;
					best.Fmaxt <- t;
					best.Fmaxw <- w;
					cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, "best weight:",best.Fmaxw, sep="\t","\n");
				}
			}
		}
	}
	S.test <- tpr.weighted.threshold(S.test, g, root=root, t=best.Fmaxt, w=best.Fmaxw);	## testing set hierarchical correction..
		
	## AUC (average and per.class) computed with PerfMeas package on the test set
	AUC.tpr <- AUROC.single.over.classes(ann.test, S.test);

	## PxR at fixed recall levels (average and per.class) computed with PerfMeas package on the test set
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);	

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann.test, S.test); 

	## storing the hierarchical matrix
	S.tpr <- S.test;
	rm(S.test, S.training);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprWT.rda"), compress=TRUE);
	}
}

# High level function to compute hierarchical correction according to TPR-threshold-free algorithm. 
# INPUT: 
# norm: boolean value: 1.TRUE means that the flat scores matrix has been already normalized in according to a normalization method (def);
#       			   2.FALSE means that the flat scores matrix has NOT been normalized yet.
# norm.type: this variable can assume three values: MaxNorm, Qnorm, NONE. We have two case respect to norm:
#			 1.if norm==FALSE, two kind of normalizations are possible:	1. MaxNorm: each score is divided w.r.t. the max of each class 
#																		2. Qnorm: quantile normalization is applied. Library preprocessCore is used.
#			 2.if norm==TRUE, set norm.type=="NONE" (def);
# flat.file: name of flat scores matrix already normalized or to be normalized in according to norm.type (without rda extension);
# ann.file: name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples (def: ann.file)
# dag.file: name of the graph that represents the hierarchy of the classes (def dag.file) 
# flat.dir: relative path to folder where flat scores matrix (already normalized or to normalize) is stored (def flat.dir)
# ann.dir: relative path to folder where annotation matrix is stored
# dag.dir: relative path to folder where graph is stored
# flat.norm.dir: 1.if norm=FALSE, relative path where flat normalized scores matrix is strored;
#				 2.if norm=TRUE, the flat scores matrix is already normalized, than it is set to NULL (def)
# b.per.example : if TRUE (def.) precision, recall, F-measure, specificity and accuracy are returned for each example, 
#                 otherwise only the averages are returned.
# n.round: 	number of rounding digits to be applied to the hierarchical scores matrix (def. 3). 
#			It's used for choosing the best threshold on the basis of the best F.measure (see f.criterion parameter).
# f.criterion: character. Type of F-measure to be used to select the best F.measure. There are 2 possibilities: 
#              1. "F" (default) corresponds to the harmonic mean between the average precision and recall; 
#              2. "avF" corresponds to the per-example F-score averaged across all the examples.
# hierScore.dir: relative path to folder where the matrix with the scores of the classes corrected in according to TPR-threshold-free algorithm is stored 
# perf.dir: relative path to folder where the macro-centric measures (i.e. AUC, PRC and PxR across classes) and 
# 			where example-centric measures (i.e. Precision, Recall, Specificity, F-measure, Accuracy across example) are stored
# OUTPUT:
# 5 rda files stored in the rispective output directory:
# - Matrix with examples on rows and classes on colums representing the hierarchical scores of the classes computed with TPR-threshold-free algorithm.
#	Stored in hierScore.dir folder
# - Example-centric measures computed through find.best.f from F-hier.R file. Stored in perf.dir folder
# - AUC (average and per classes) computed through AUC.single.over.classes from package PerfMeas. Stored in perf.dir
# - Precision at fixed recall levels (average and per classes) computed through precision.at.multiple.recall.level.over.classes from package PerfMeas.
# 	Stored in perf.dir folder 
# - PRC (average and per.class) computed by precrec package. Stored in perf.dir
Do.tpr.threshold.free <- function(	norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file,
									flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, 
									b.per.example=TRUE, f.criterion ="F", hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"
								){
	## Loading Data ############
	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## Computing Flat Perfromances
	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); 
	
	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); 

	## Hierarchical Correction #########################
	## TPR-Threshold-free correction 
	S <- tpr.threshold.free(S, g, root=root);

	## Computing Hierarchical Performances
	## Hierarchical AUC (average and per.class) computed by PerfMeas package
	AUC.tpr <- AUROC.single.over.classes(ann, S);
	
	## Hierarchical PxR at fixed recall levels (average and per.class) computed by PerfMeas package
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann, S);

	## Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann, S); 

	## storing the hierarchical matrix
	S.tpr <- S;
	rm(S);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
	}
}

# Hold-Out version
# High level function to correct the scores with a hierarchy according to TPR-threshold-free algorithm performing a classical holdout procedure
# All input parameters are the same of the corresponding function above, except the following:
# ind.test.set: vector of integer. Indices refer to the examples of the adjancency matrix to be used in the test set. 
# ind.dir: relative path to folder where ind.test.set is stored
Do.tpr.threshold.free.holdout <- function(	norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file,
											flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL,
											ind.test.set=ind.test.set, ind.dir=ind.dir, n.round=3, f.criterion ="F",											 
											b.per.example=TRUE, hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"
										){
	## Loading Data ############
	# loading examples indices of the test set
	ind.set <- paste0(ind.dir, ind.test.set, ".rda");
	ind.test <- get(load(ind.set));

	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm==TRUE){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## shrinking the size of S to the examples of test set
	S <- S[ind.test,];

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[ind.test,-which(colnames(ann)==root)];

	## Computing Flat Perfromances
	## FLAT AUC computed by PerfMeas package
	AUC.flat <- AUROC.single.over.classes(ann, S); 

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); 

	## Hierarchical Correction #########################
	## TPR-Threshold-free correction 
	S <- tpr.threshold.free(S, g, root=root);

	## Computing Hierarchical Performances
	## Hierarchical AUC (average and per.class) computed by PerfMeas package
	AUC.tpr <- AUROC.single.over.classes(ann, S);
	
	## Hierarchical PxR at fixed recall levels (average and per.class) computed by PerfMeas package
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann, S);

	## Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann, S, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann, S); 

	## storing the hierarchical matrix
	S.tpr <- S;
	rm(S);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprTF.rda"), compress=TRUE);
	}
}


##**************************************************##
## TPR version with DESCENDANTS instead of children ##
##**************************************************##
# High level function to compute hierarchical correction according to TPR-threshold-descendants algorithm. 
# It perform a k fold cross-validation to find the best threshold maximizing on F.max measure.
# INPUT: 
# threshold: range of threshold values to be tested in order to find the best threshold. The denser the range is, the higher the probability to find the best
# 			 theshold, but the execution time will be increasing (def: from:0.1, to:0.9, by:0.1 step).
# kk:  number of folds of the cross validation (def: 5).
# seed:  intialization seed for the random generator to create folds (def:0). If NULL (default) no initialization is performed.
# flat.norm.file: name of the flat scores matrix already normalized in according to a normalization method (without rda extension). 
# ann.file: name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples (def: ann.file)
# dag.file: name of the graph that represents the hierarchy of the classes. 
# b.per.example : if TRUE (def.) precision, recall, F-measure, specificity and accuracy are returned for each example, 
#                 otherwise only the averages are returned.
# n.round: 	number of rounding digits to be applied to the hierarchical scores matrix (def. 3). 
#			It's used for choosing the best threshold on the basis of the best F.measure (see f.criterion parameter).
# f.criterion: character. Type of F-measure to be used to select the best F.measure. There are 2 possibilities: 
#              1. "F" (default) corresponds to the harmonic mean between the average precision and recall; 
#              2. "avF" corresponds to the per-example F-score averaged across all the examples.
# flat.dir: relative path to folder where flat normalized scores matrix is stored
# ann.dir: relative path to folder where annotation matrix is stored
# dag.dir: relative path to folder where graph is stored
# hierScore.dir: relative path to folder where the matrix with the scores of the classes corrected in according to TPR-threshold-descendants algorithm is stored 
# perf.dir: relative path to folder where the macro-centric measures (i.e. AUC, PRC and PxR across classes) and 
# 			where example-centric measures (i.e. Precision, Recall, Specificity, F-measure, Accuracy across example) are stored
# OUTPUT:
# 5 rda files stored in the rispective output directory:
# - Matrix with examples on rows and classes on colums representing the hierarchical scores of the classes computed with TPR-threshold-descendants algorithm.
#	Stored in hierScore.dir folder
# - Example-centric measures computed through find.best.f from F-hier.R file. Stored in perf.dir folder
# - AUC (average and per classes) computed through AUC.single.over.classes from package PerfMeas. Stored in perf.dir
# - Precision at fixed recall levels (average and per classes) computed through precision.at.multiple.recall.level.over.classes from package PerfMeas.
# 	Stored in perf.dir folder 
# - PRC (average and per.class) computed by precrec package. Stored in perf.dir
Do.tpr.threshold.descendants.cv <- function(	threshold=seq(from=0.1, to=0.9, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type= "NONE",
												flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, flat.dir=flat.dir, 
												ann.dir=ann.dir,dag.dir=dag.dir, flat.norm.dir=NULL, b.per.example=TRUE,										
												n.round=3, f.criterion ="F", hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"
											){
	## Loading Data ############
	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## FLAT AUC computed by PerfMeas package
	AUC.flat <- AUROC.single.over.classes(ann, S); 
	
	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); 

	## Hierarchical Correction #########################
	## TPR-threshold-descendants correction using k-fold cross-validation to find best threshold 

	## splitting data in k stratified fold
	folds <- do.unstratified.cv.data(S, kk=kk, seed=seed); 

	## storing the best F.max and the best threshold for each of k training set
	training.top.Fmax <- vector(mode="list", length=kk);
	names(training.top.Fmax) <- paste0(rep("fold",kk), 1:kk);

	## storing all the best protein-centric measures for each of k test set 
	best.avg.meas.test <- vector(mode="list", length=kk);
	names(best.avg.meas.test) <- paste0(rep("fold",kk), 1:kk);
	FMM.per.example <- c();

	## storing average macro AUC for each of k test set 
	AUC.average.test <- vector(mode="list", length=kk);
	names(AUC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro AUC for each of k test set 
	AUC.class.test <- vector(mode="list", length=kk);
	names(AUC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average PxR for each of k test set 
	PXR.average.test <- vector(mode="list", length=kk);
	names(PXR.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class PxR for each of k test set 
	PXR.class.test <- vector(mode="list", length=kk);
	names(PXR.class.test) <- paste0(rep("fold",kk), 1:kk);

	## storing average macro PRC for each of k test set 
	PRC.average.test <- vector(mode="list", length=kk);
	names(PRC.average.test) <- paste0(rep("fold",kk), 1:kk);

	## storing per class macro PRC for each of k test set 
	PRC.class.test <- vector(mode="list", length=kk);
	names(PRC.class.test) <- paste0(rep("fold",kk), 1:kk);

	## variable for hosting the k sub-matrix assembled 
	S.tpr <- c();

	## Let's start k-fold crossing validation for choosing best threshold maximizing on F.max...
	for(k in 1:kk){
		test <- S[folds[[k]],];		# test set: 1 out of k folds (e.g.: if k=5, testing set is 1/5 of data)
		training <- S[!rownames(S) %in% rownames(test),];		# training set: (k-1) out of k folds (e.g.: if k=5, training set is 4/5 of data)
		gc();

		top.Fmax <- 0;
		best.Fmaxt <- 0;
		for(t in threshold){
			pred.training <- tpr.threshold.descendants(training, g, root=root, t=t);	## training set hierarchical correction...
			target.training <- ann[rownames(pred.training),colnames(pred.training)];
			training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE);	 
			training.Fmax <- training.avg.meas[4];	## F.max maximization...
			if(training.Fmax > top.Fmax){
				top.Fmax <- training.Fmax;
				best.Fmaxt <- t;
				training.top.Fmax[[k]] <- c(best.Fmax=top.Fmax, best.thres=best.Fmaxt);
				cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, sep="\t", "\n");
			}
		} 
		pred.test <- tpr.threshold.descendants(test, g, root=root, t=training.top.Fmax[[k]][2]);	## testing set hierarchical correction...
		target.test <- ann[rownames(pred.test),colnames(pred.test)];
		test.avg.meas <- find.best.f(target.test, pred.test, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=b.per.example);
		best.avg.meas.test[[k]] <- test.avg.meas$average;
		FMM.per.example <- rbind(FMM.per.example,test.avg.meas$per.example);

		## Hierarchical AUC (average and per.class) computed with PerfMeas package
		AUC.average.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$average;	## storing average AUC of each k testing set
		AUC.class.test[[k]] <- AUROC.single.over.classes(target.test,pred.test)$per.class;	## storing AUC per.class of each k testing set

		## Hierarchical PxR at fixed recall levels (average and per.class) computed with PerfMeas package
		PXR.average.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$avgPXR;	## storing average PxR of each k testing set
		PXR.class.test[[k]] <- precision.at.multiple.recall.level.over.classes(target.test,pred.test)$PXR;	## storing PxR per.class of each k testing set

		## PRC (average and per.class) computed with precrec package
		PRC.average.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$average;	## storing average AUC of each k testing set
		PRC.class.test[[k]] <- AUPRC.single.over.classes(target.test,pred.test)$per.class;	## storing AUC per.class of each k testing set

		## assembling of all the hierarchical scores of each k sub-matrix to build the full matrix 
		S.tpr <- rbind(S.tpr, pred.test);
	}
	## put the rows (i.e. genes) of assembled k sub-matrix in the same order of the full beginning matrix
	S.tpr <- S.tpr[rownames(S),];

	## remove no longer useful variables..
	rm(S, folds, pred.test);

	## Averaging Performances Measures across k testing set ###################
	## averaging protein-centric measures across k testing sets
	F.cv <- Reduce("+", best.avg.meas.test)/kk;
	F.meas <-  apply(FMM.per.example,2,mean);
	names(F.meas)[4] <- "avF";
	F.max <- 2*(F.meas[["P"]] * F.meas[["R"]])/((F.meas[["P"]] + F.meas[["R"]]));
	names(F.max) <- "F";
	FMM.avg <- append(F.meas, F.max, after=3);
	FMM.avg <- append(FMM.avg,F.cv["T"]);
	FMM.tpr <- list(average=FMM.avg, per.example=FMM.per.example);

	## averaging macro-AUC (average and per.class) across k testing sets 
	AUC.average.tpr.over.test <- Reduce("+", AUC.average.test)/kk;
	AUC.class.tpr.over.test <- Reduce("+", AUC.class.test)/kk;
	AUC.tpr <- list(average=AUC.average.tpr.over.test, per.class=AUC.class.tpr.over.test);

	## averaging PxR (average and per.class) across k testing sets 
	PXR.average.tpr.over.test <- Reduce("+", PXR.average.test)/kk;
	PXR.class.tpr.over.test <- Reduce("+", PXR.class.test)/kk;
	PXR.tpr <- list(average=PXR.average.tpr.over.test, per.class=PXR.class.tpr.over.test);

	## averaging PRC (average and per.class) across k testing sets 
	PRC.average.tpr.over.test <- Reduce("+", PRC.average.test)/kk;
	PRC.class.tpr.over.test <- Reduce("+", PRC.class.test)/kk;
	PRC.tpr <- list(average=PRC.average.tpr.over.test, per.class=PRC.class.tpr.over.test);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
	}
}

# High level function to correct the scores with a hierarchy according to TPR-threshold-descendants algorithm performing a classical holdout procedure
# All input parameters are the same of the above corresponding function that implement a kfcv, except the following:
# ind.test.set: vector of integer. Indices refer to the examples of the adjancency matrix to be used in the test set. 
# ind.dir: relative path to folder where ind.test.set is stored
Do.tpr.threshold.descendants.holdout <- function(	threshold=seq(from=0.1, to=0.9, by=0.1), kk=5, seed=NULL, norm=TRUE, norm.type= "NONE",
													flat.file=flat.file, ann.file=ann.file, dag.file=dag.file, ind.test.set=ind.test.set,
													ind.dir=ind.dir,flat.dir=flat.dir, ann.dir=ann.dir,dag.dir=dag.dir, flat.norm.dir=NULL, 
													n.round=3, b.per.example=TRUE, f.criterion ="F", hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"		
												){
	## Loading Data ############
	## loading examples indices of the test set
	ind.set <- paste0(ind.dir, ind.test.set, ".rda");
	ind.test <- get(load(ind.set));

	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.file,".rda");
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir);
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## scores flat matrix shrinked to test and training test respectively
	S.test <- S[ind.test,];
	S.training <- S[-ind.test,];
	rm(S);

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table and shrinking annotation table to tet set
	ann.test <- ann[ind.test,-which(colnames(ann)==root)];

	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann.test, S.test); 

	# FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann.test, S.test); 

	## Hierarchical Correction #########################
	## splitting data in k stratified fold
	folds <- do.unstratified.cv.data(S.training, kk=kk, seed=seed); 
			
	## Let's start k-fold crossing validation for choosing best threshold maximizing on the basis of F.max value...
	for(k in 1:kk){
		training <- S.training[folds[[k]],];	
		top.Fmax <- 0;
		best.Fmaxt <- 0;
		for(t in threshold){
			pred.training <- tpr.threshold.descendants(training, g, root=root, t=t);	## training set hierarchical correction...
			target.training <- ann[rownames(pred.training),colnames(pred.training)];
			training.avg.meas <- find.best.f(target.training, pred.training, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=FALSE); 
			training.Fmax <- training.avg.meas[4];	## F.max maximization...
			if(training.Fmax > top.Fmax){
				top.Fmax <- training.Fmax;
				best.Fmaxt <- t;
				cat("training fold:",k, "better Fmax.avg found:",top.Fmax, "best threshold:",best.Fmaxt, sep="\t", "\n");
			}
		}
	}
	S.test <- tpr.threshold.descendants(S.test, g, root=root, t=best.Fmaxt);	## testing set hierarchical correction..
		
	## AUC (average and per.class) computed with PerfMeas package on the test set
	AUC.tpr <- AUROC.single.over.classes(ann.test, S.test);
	
	## PxR at fixed recall levels (average and per.class) computed with PerfMeas package on the test set
	PXR.tpr <- precision.at.multiple.recall.level.over.classes(ann.test, S.test);

	## F.measure: Computing Hierarchical Examples-Measures 
	FMM.tpr <- find.best.f(ann.test, S.test, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example);	

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.tpr <- AUPRC.single.over.classes(ann.test, S.test); 

	## storing the hierarchical matrix
	S.tpr <- S.test;
	rm(S.test, S.training);

	## Storing Results #########
	if(norm){
		save(S.tpr, file=paste0(hierScore.dir, flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
	}else{
		save(S.tpr, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(AUC.flat, AUC.tpr, file=paste0(perf.dir, "AUC.", norm.type, ".", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(PXR.flat, PXR.tpr, file=paste0(perf.dir, "PXR.", norm.type, ".", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(FMM.flat, FMM.tpr, file=paste0(perf.dir, "PCM.", norm.type, ".", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
		save(PRC.flat, PRC.tpr, file=paste0(perf.dir, "PRC.", norm.type, ".", flat.file, ".hierScores.tprTD.rda"), compress=TRUE);
	}
}


