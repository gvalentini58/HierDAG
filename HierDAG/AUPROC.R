## High-level function to compute AUROC and AUPRC through precrec package, more precise than PerfMeas (especially for AUPRC)
## Note: PerfMeas computes the PRC in a *wrong* way whereas the AUC not, but precrec should be slightly more accurate.

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

library(precrec);
## NOTE: classes that have zero annotations are impossible to predict by definition. If we have a class with zero annotation the package 
## precrec stop and return the following error message: 'Error: length(unique(labels)) not equal to 2L'. We have to avoid this problem, especially 
## if we need to compute AUC and PRC (by precrec pkg) in a contest of k-fold cross-validation, where is highly likely that a fold can have
## classes with zero annotations..

# Function to compute PRC for a single class
# Input
# target: vector of the true labels (0 negative, 1 positive examples)
# pred: numeric vector of  the values of the predicted labels
# Output:
# a numeric value corresponding to the AUC
AUPRC.single.class <- function(target, pred){
	if(sum(target)==0){
		PRC <- 0;
		return(PRC)
	}
	else{
		res <- evalmod(scores=pred, labels=target);
		aucs <- auc(res);
		prc <- subset(aucs, curvetypes == "PRC");
		PRC <- prc$aucs;
		return(PRC);
	}
}


# Function to compute AUC for a single class
# Input
# target: vector of the true labels (0 negative, 1 positive examples)
# pred: numeric vector of  the values of the predicted labels
# Output:
# a numeric value corresponding to the AUC
AUROC.single.class <- function(target, pred){
	if(sum(target)==0){
		AUC <- 0.5;
		return(AUC)
	}
	else{
		res <- evalmod(scores=pred, labels=target);
		aucs <- auc(res);
		roc <- subset(aucs, curvetypes == "ROC");
		AUC <- roc$aucs;
		return(AUC);
	}
}


# Function that computes Precision and Recall Curves (PRC) for each class. Both the target e predicted matrices have a 
# number of rows equal to the number of examples and a number of columns equal to the number of the classes.
# Input 
# target: matrix with the target multilabels
# pred: a numeric matrix with predicted values
# Output 
# a list with two elements:
# average: the average AUC across classes.              
# per.class: a named vector with AUC for each class. Names correspond to classes
AUPRC.single.over.classes <- function(target, pred){
	## if there are classes with zero annotations, we remove them..
	target.names <- colnames(target);
	class.ann <- apply(target,2,sum);
	class.noann <- which(class.ann==0);
	check <- length(class.noann)!=0;
	if(check){
		target <- target[,-class.noann];
		pred <- pred[,-class.noann];
	}
	
	## compute PRC considering only those class with non-zero annotations
	PRC.class <- rep(0,ncol(pred));
	names(PRC.class) <- colnames(pred);
	for(i in 1:ncol(pred)){
		PRC.class[i] <- AUPRC.single.class(target[,i],pred[,i]);
	}
	
	## if there are classes with zero annotations, set the prc of those classes to zero and restore the start classes order 
	if(check){
		PRC.class <- PRC.class[target.names];
		PRC.class[is.na(PRC.class)] <- 0;
		names(PRC.class) <- target.names; 
	}

	#saving PRC result in the same format of package PerfMeas
	PRC.mean <- mean(PRC.class);
	PRC.res <- list(average=PRC.mean, per.class=PRC.class); 
	return(PRC.res);
}

# Function that computes Area Under Curves (AUC) for each class. Both the target e predicted matrices have a 
# number of rows equal to the number of examples and a number of columns equal to the number of the classes.
# Input
# target: matrix with the target multilabels
# pred: a numeric matrix with predicted values
# Output 
# a list with two elements:
# average: the average AUC across classes.              
# per.class: a named vector with AUC for each class. Names correspond to classes
AUROC.single.over.classes <- function(target, pred){
	## if there are classes with zero annotations, we remove them..
	target.names <- colnames(target);
	class.ann <- apply(target,2,sum);
	class.noann <- which(class.ann==0);
	check <- length(class.noann)!=0;
	if(check){
		target <- target[,-class.noann];
		pred <- pred[,-class.noann];
	}
	
	## compute AUC considering only those class with non-zero annotations
	AUC.class <- rep(0,ncol(pred));
	names(AUC.class) <- colnames(pred);
	for(i in 1:ncol(pred)){
		AUC.class[i] <- AUROC.single.class(target[,i],pred[,i]); 
	}

	## if there are classes with zero annotations, set the AUC of those classes to zero and restore the start classes order 
	if(check){
		AUC.class <- AUC.class[target.names];
		AUC.class[is.na(AUC.class)] <- 0.5;
		names(AUC.class) <- target.names; 
	}
		
	#saving AUC result in the same format of package PerfMeas
	AUC.mean <- mean(AUC.class);
	AUC.res <- list(average=AUC.mean, per.class=AUC.class); 
	return(AUC.res);
}


######################################################
#### OLD FUNCTIONs: too RAM memory with BigMatrix ####
######################################################
## NOTE: the old function takes use of "join_scores" (i.e. list of double scores) that fill in flash the RAM if we are dealing with big
# size matrix... 
# Function that computes Precision and Recall Curves (PRC) for each class. Both the target e predicted matrices have a
# number of rows equal to the number of examples and a number of columns equal to the number of the classes.
# Input
# target: matrix with the target multilabels
# pred: a numeric matrix with predicted values
# Output 
# a list with two elements:
# average: the average AUC across classes.              
# per.class: a named vector with AUC for each class. Names correspond to classes
OLD_AUPRC.single.over.classes <- function(target, pred){
	## if there are classes with zero annotations, we remove them..
	target.names <- colnames(target);
	class.ann <- apply(target,2,sum);
	class.noann <- which(class.ann==0);
	check <- length(class.noann)!=0;
	if(check){
		target <- target[,-class.noann];
		pred <- pred[,-class.noann];
	}
	
	## compute PRC considering only those class with non-zero annotations
	## NOTE: the functions join_label/scores fill in a flash the ram memory with matrix of big size!!
	## NOTE: of course we get the same result if we use lapply instead...
	## SOLUTION: use for loop: process one vector of scores and labels at time, even if slowly it's the best solutions
	labels <- join_labels(target);	## lapply(seq_len(ncol(target)), function(i) target[,i])	
	scores <- join_scores(pred);   	## lapply(seq_len(ncol(pred)), function(i) pred[,i])
	res <- evalmod(scores=scores, labels=labels, dsids=1:ncol(pred), modnames=colnames(pred));
	aucs <- auc(res);
	prc <- subset(aucs, curvetypes == "PRC");
	PRC.class <- prc$aucs;
	names(PRC.class) <- prc$modnames;
		
	## if there are classes with zero annotations, set the prc of those classes to zero and restore the start classes order 
	if(check){
		PRC.class <- PRC.class[target.names];
		PRC.class[is.na(PRC.class)] <- 0;
		names(PRC.class) <- target.names; 
	}
		
	#saving PRC result in the same format of package PerfMeas
	PRC.mean <- mean(PRC.class);
	PRC.res <- list(average=PRC.mean, per.class=PRC.class); 
	return(PRC.res);
}
# Function that computes Area Under Curves (AUC) for each class. Both the target e predicted matrices have a 
# number of rows equal to the number of examples and a number of columns equal to the number of the classes.
# Input
# target: matrix with the target multilabels
# pred: a numeric matrix with predicted values
# Output 
# a list with two elements:
# average: the average AUC across classes.              
# per.class: a named vector with AUC for each class. Names correspond to classes
OLD_AUROC.single.over.classes <- function(target, pred){
	## if there are classes with zero annotations, we remove them..
	target.names <- colnames(target);
	class.ann <- apply(target,2,sum);
	class.noann <- which(class.ann==0);
	check <- length(class.noann)!=0;
	if(check){
		target <- target[,-class.noann];
		pred <- pred[,-class.noann];
	}
	
	## compute AUC considering only those class with non-zero annotations
	labels <- join_labels(target);	
	scores <- join_scores(pred);
	res <- evalmod(scores=scores, labels=labels, dsids=1:ncol(pred), modnames=colnames(pred));
	aucs <- auc(res);
	roc <- subset(aucs, curvetypes == "ROC");
	AUC.class <- roc$aucs;
	names(AUC.class) <- roc$modnames;
		
	## if there are classes with zero annotations, set the AUC of those classes to zero and restore the start classes order 
	if(check){
		AUC.class <- AUC.class[target.names];
		AUC.class[is.na(AUC.class)] <- 0.5;
		names(AUC.class) <- target.names; 
	}
		
	#saving AUC result in the same format of package PerfMeas
	AUC.mean <- mean(AUC.class);
	AUC.res <- list(average=AUC.mean, per.class=AUC.class); 
	return(AUC.res);
}



