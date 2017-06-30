##########################################################################
# Functions to compute Kiritchenko-like multi-label F-scores
# March 2016
##########################################################################

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

# Generic function for computing, precision, recall, specificity and F-measure for multiclass multilabel classification
setGeneric("F.measure.multilabel", 
                 function(target, predicted, b.per.example=FALSE) standardGeneric("F.measure.multilabel"));

# Method that computes precision, recall, specificity and F-measure for multiclass multilabel classification
# Both the target e predicted matrices have a number of rows equal to the number of examples
# and a number of columns equal to the number of the classes.
# Input:
# target : matrix with the target multilabels. 1 entries correspond to postives, 0 to negatives
# predicted : matrix with the predicted multilabels. 1 entries correspond to postives, 0 to negatives
# b.per.example : if TRUE (def: FALSE) precision recall and F-measure are returned for each example, otherwise only the averages are returned.
# Output :
# if  b.per.example == FALSE the function returns a list with a single element average: a named vector with average  precision (P), recall (R), specificity (S)  F-measure (F), average F-measure (avF) and accuracy (A) across examples;
# otherwise it returns a list with  two elements:
# average :  a named vector with the average precision, recall, specificity F-measure and av.F-measure across examples.
#            The elements correspond respectively to the average precision, recall, specificity; F-measure is the 
#            F-measure computed as the harmonic mean between the average precision and recall;
#            av.f.measure is the f.measure computed as the average across examples.
# per.example : a named matrix with the precision, recall, specificity and F-measure for each example.
#                       Named rows correspond to examples,
#                       named columns correspond respectively to precision, recall (sensitivity), specificity and F-measure.
setMethod("F.measure.multilabel", signature(target="matrix", predicted="matrix"),
  function(target, predicted, b.per.example=FALSE) { 
       n.examples <- nrow(target);
       n.classes <- ncol(target);
       if ((n.examples!=nrow(predicted)) || (n.classes!=ncol(predicted)))
          stop ("F.measure.multilabel: number of rows or columns do not match between target and predicted classes");
	  
       z <- target + predicted;
       TP <- apply(z, 1, function(x)  {
                            return(sum(x==2));
                          });
       TN <- apply(z, 1, function(x)  {
                            return(sum(x==0));
                          });
       z <- predicted - target;
       FP <- apply(z, 1, function(x)  {
                            return(sum(x==1));
                          });
       FN <- apply(z, 1, function(x)  {
                            return(sum(x== -1));
                          });
       rm(z);
       n <- sum(TP)+sum(TN)+sum(FN)+sum(FP);
       if ( n != (n.examples*n.classes)) { 
	       cat("n = ", n, "\n n.examples = ", n.examples, "\n n.classes = ", n.classes, "\n");
		   cat (" sum(TP) = ", sum(TP), "\n sum(TN) = ", sum(TN), "\n sum(FN) = ", sum(FN), "\n sum(FP) = ", sum(FP), "\n");
           warning("F.measure.multilabel: Something went wrong in F-measure)");
	   }
	   
       P <- TP+FP;
       P[which(P==0)] <- 1;  # to avoid division by 0 in precision
       
       sum.TP.FN <- TP+FN;
       sum.TN.FP <- TN+FP;
       
        sum.TP.FN[which(sum.TP.FN==0)] <- 1;  # to avoid division by 0 in recall
        sum.TN.FP[which(sum.TN.FP==0)] <- 1;  # to avoid division by 0 in specificity
          
       precision <- TP/P;
       recall <- TP/sum.TP.FN;
       specificity <- TN/sum.TN.FP;
       
       prec.rec <- precision+recall;
       prec.rec[which(prec.rec==0)] <- 1;  # to avoid division by 0 for f.measure
       f.measure <- (2*precision*recall)/prec.rec;
       accuracy <- (TP+TN)/n.classes;
       
       av.precision <- sum(precision)/n.examples; 
       av.recall <- sum(recall)/n.examples; 
       av.specificity <- sum(specificity)/n.examples; 
	     av.prec.rec <- av.precision+av.recall;
	   if (av.prec.rec == 0)  av.prec.rec <- 1;
	     overall.av.f.measure <- (2*av.precision*av.recall)/av.prec.rec;
       av.f.measure <- sum(f.measure)/n.examples; 
       av.accuracy  <- sum(accuracy)/n.examples; 
       
       average <- c(av.precision, av.recall, av.specificity, overall.av.f.measure, av.f.measure,av.accuracy);
       names(average) <- c("P", "R", "S", "F", "avF", "A");
       
       if (b.per.example)  {
          per.example <- cbind(precision, recall, specificity, f.measure, accuracy);
    colnames(per.example) <- c("P", "R", "S", "F","A");
          return (list(average=average, per.example=per.example))
       } else
          return (list(average=average));
   } 
)


# Function to select the best hierarchical F-score by choosing an appropriate threshold in the scores
# N.B. All the examples having no positive annotations are discarded
#  Arguments:
# target : matrix with the target multilabels. 1 stands for positive, 0 for negative
# pred : matrix with the predicted scores. Values are assumed to be postive 
# n.round : number of rounding digits to be applied to pred (default=3)
# f.criterion : character. Type of F-measure to be used to select the best F.  There are 2 possibilities: 
#               1. "F" (default) corresponds to the harmonic mean between the average precision and recall; 
#               2. "avF" corresponds to the per-example F-score averaged across all the examples.
# verbose : boolean. If TRUE (def) the number of iterations are printed on stdout
# The pred matrix is rounded according to n.round and all the values of pred are divided by max(pred). 
# Then all the thresholds corresponding to all the different values included in pred are attempted, and the threshold 
# leading to the maximum f.measure is selected.
# b.per.example : if TRUE (def: FALSE) precision, recall, F-measure, specificity and accuracy are returned for each example, 
#                 otherwise only the averages are returned.
# Output:
# if b.per.example == FALSE (def.) the function returns a vector with 7 elements relative to the best result in terms of the f.measure:  
# precision, recall, specificity, f.measure, av.f.measure, accuracy, thresh, where thresh thresh is the selected best threshold,  
# av.f.measure is the f.measure averaged across examples and  f.measure is the f-score computed as the harmonic mean between the 
# average precision and the average recall   
# otherwise (if b.per.example == TRUE) the function returns a list with  two elements:
# 0) average: the same vector with 7 elements aforementioned
# 1) per.example: a named matrix with the precision, recall, specificity and F-measure for each example.
#                 Named rows correspond to examples, named columns correspond respectively to precision, recall (sensitivity), specificity and F-measure.
find.best.f <- function(target, pred, n.round=3, f.criterion ="F", verbose=TRUE, b.per.example=FALSE)  {
  
  x<- apply(target,1,sum);
  selected <- which(x>0);
  target <- target[selected,];
  pred <- pred[selected,];
  pred <- pred/max(pred);
  pred <- round(pred,n.round);
  n.examples <- nrow(pred);
  n.classes <- ncol(pred);
  
  thresh <- unique(as.numeric(pred));
  thresh <- sort(thresh);
  best.res <- best <- best.thresh <- 0;
  i=0;
  for (t in thresh) {
    pred.labels <- matrix(numeric(n.examples*n.classes), nrow=n.examples);
    pred.labels[pred>=t] <-1;
    res <- F.measure.multilabel(target, pred.labels, b.per.example);
    if (res$average[f.criterion] > best) {
       best <- res$average[f.criterion];
       best.res <- res;  
       best.thresh <- t;
    }
    i <- i+1;
    if (i%%100 == 0  && verbose) 
      cat("iteration ", i,  "\n");
  }
  if(b.per.example){
    best.res$average <- c(best.res$average, best.thresh);
    names(best.res$average)[7] <- "T"; 
    return(best.res);
  }else{
    best.res <- c(best.res$average, best.thresh);
    names(best.res)[7] <- "T";
    return(best.res);
  }
}


# Function to select the best hierarchical F-score by choosing an appropriate threshold in the scores
# N.B. All the examples having no positive annotations are discarded
#  Arguments:
# target : matrix with the target multilabels. 1 stands for positive, 0 for negative
# pred : matrix with the predicted scores. Values are assumed to be postive 
# n.round : number of rounding digits to be applied to pred (default=3)
# f.criterion : character. Type of F-measure to be used to select the best F.  There are 2 possibilities: 
#               1. "F" (default) corresponds to the harmonic mean between the average precision and recall; 
#               2. "avF" corresponds to the per-example F-score averaged across all the examples.
# verbose : boolean. If TRUE (def) the number of iterations are printed on stdout
# The pred matrix is rounded according to n.round and all the values of pred are divided by max(pred). 
# Then all the thresholds corresponding to all the different values included in pred are attempted, and the threshold 
# leading to the maximum f.measure is selected.
# Output:
# a vector with 7 elements relative to the best result in terms of the f.measure:   precision, recall, specificity, f.measure, av.f.measure, accuracy, thresh
# Note: thresh is the selected best threshold   
# av.f.measure if the f.measure averaged across examples
# f.measure is the f-score computed as the harmonic mean between the average precision and the average recall   
find.best.f.OLD <- function(target, pred, n.round=3, f.criterion ="F", verbose=TRUE)  {
  
  x<- apply(target,1,sum);
  selected <- which(x>0);
  target <- target[selected,];
  pred <- pred[selected,];
  pred <- pred/max(pred);
  pred <- round(pred,n.round);
  n.examples <- nrow(pred);
  n.classes <- ncol(pred);
  
  thresh <- unique(as.numeric(pred));
  thresh <- sort(thresh);
  best.res <- best <- best.thresh <- 0;
  i=0;
  for (t in thresh) {
    pred.labels <- matrix(numeric(n.examples*n.classes), nrow=n.examples);
    pred.labels[pred>=t] <-1;
    res <- F.measure.multilabel(target, pred.labels, b.per.example=FALSE);
    if (res$average[f.criterion] > best) {
       best <- res$average[f.criterion];
       best.res <- res$average;  
       best.thresh <- t;
    }
    i <- i+1;
    if (i%%100 == 0  && verbose) 
      cat("iteration ", i,  "\n");
  }
  best.res <- c(best.res, best.thresh);
  names(best.res)[7] <- "T"; 
  return(best.res);
}

# Function to select the best hierarchical F-score by choosing an appropriate threshold in the scores. This is the old version.
#  Arguments:
# target : matrix with the target multilabels
# pred : matrix with the predicted scores
# attempts : number of different thresholds to be attempted to find the best F-score. If attempts=0 all the scores in pred are used as possible threshold
# Output:
# a vector with 5 elements relative to the best results:   precision, recall, specificity, f.measure, accuracy      
find.best.f.old <- function(target, pred, attempts=500.0)  {

  thresh <- unique(as.numeric(pred));
  thresh <- sort(thresh);
  len.values <- length(thresh);
  best.res <- best <- 0;
  i=0;
  if (attempts==0) {
     for (t in thresh) {
        pred.labels <- matrix(numeric(nrow(pred)*ncol(pred)), nrow=nrow(pred));
        pred.labels[pred>t] <-1;
        res <- F.measure.multilabel(target, pred.labels, b.per.example=FALSE);
        if (res$average[4] > best) {
            best <- res$average[4];
            best.res <- res$average;  
        }
        i <- i+1;
        if (i%%10 == 0) 
           cat("attempt ", i,  "\n");
     }   
  } else {
     if (len.values>attempts) 
   	step <- len.values/attempts  else   step <- 1; 
     n.step <- 0;
     while(n.step<attempts) {
       n.step <- n.step + 1;
       thresh.value <- thresh[floor(n.step*step)];
       pred.labels <- matrix(numeric(nrow(pred)*ncol(pred)), nrow=nrow(pred));
       pred.labels[pred>thresh.value] <-1;
       pred.labels[pred<=thresh.value]<-0;
       res <- F.measure.multilabel(target, pred.labels, b.per.example=FALSE);
       if (res$average[4] > best) {
   	  best <- res$average[4];
   	  best.res <- res$average;  
       }
       i <- i+1;
       if (i%%10 == 0) 
   	 cat("attempt ", i,  "\n");
     }
  }
  return(best.res);
}


# Usage example:
# Here target is a matrix with 0/1 entries of the target multilabels, and pred a numeric matrix with the predicted scores. Attempts is the number of different thresholds tested: larger is the number most accurate is the result.
#res <- find.best.f(target, pred, n.round=3);
