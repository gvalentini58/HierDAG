#############################################################
## Hierarchical True Path Rule for Directed Acyclic Graphs ##
#############################################################

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

# July 2014: tpr.threshold, tpr.threshold.free, 
# Semptemper 2014: added tpr.weighted.threshold and tpr.weighted.threshold.free
# April 2016: added tpr.threshold.descendants, tpr.threshold.free.descendants

## library to call
library("graph");
library("RBGL");
library("Rgraphviz");
source("graph.utils.R");


##----------------------------------------------------------------------------------------------------------------------------##
# Different variants of TPR algorithm are implemented. They are a two-step hierchical correction methods:
# 1. in the first step they compute a complete bottom-up step from the leaves to the root to propagate positive 
#	 predictions across the hierarchy;
# 2. in the second step they compute a complete top-down step from the root to the leaves according to the HTD algorithm 
#	 in order to assure the hierarchical consistency of the predictions.
# All the TPR variants can be applied to any DAG structured hierarchy of classes
##-----------------------------------------------------------------------------------------------------------------------------##


# True Path Rule with Threshold (TPR-T)
# Input:
# S : the flat scores matrix relative to a set of examples. 
#	  The colnames must correspond to the names of the classes of the hierarchy (except the root node)
# g : a graph of class graphNEL. It represents the hierarchy of the classes. The names of the classes must correspond to the name of the matrix of scores.
# root: name of the class that it is the top-level (root) of the hierarchy (def:00)
# t: threshold for the choice of positive children (def. = 0.5)
# Output:
# a matrix with the scores of the classes corrected according to the TPR-T algorithm.
tpr.threshold <- function(S,g, root="00",t=0.5){
	levels <- levels.graph(g,root);
	# si aggiunge una "dummy root column" se non presente nella matrice degli scores flat S
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	# visita bottom-up del grafo
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
		if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
				child.set <- children[j,] > t;    # scelta dei positive children in base alla soglia t
				child.pos <- children[j,][child.set];
				parent[j] <- (parent[j] + sum(child.pos))/(1+length(child.pos));   # correzione predizione flat
			}
			S[,names(chd.bup[i])] <- parent;
		}
	}
	# visita top-down del grafo
	par.tod <- get.parents.top.down(g,levels,root)
	for(i in 1:length(par.tod)){
		child <- S[,names(par.tod[i])];
		parents <- as.matrix(S[,par.tod[[i]]]);
		# colnames(parents) <- par.tod[[i]]
		# Note: the version with an apply and an ifelse statement is slower ...
		for(j in 1:length(child)){
			x <- min(parents[j,]);
			if(x < child[j]){
				child[j] <- x;    # hierarchy correction
			}
		}
		S[,names(par.tod[i])] <- child;
	}
	# the dummy root node is removed 
	S <- S[,-which(colnames(S)==root)];
	return(S);
}



# True Path Rule Threshold-Free (TPR-TF)
# Input:
# S : the flat scores matrix relative to a set of examples. 
#	  The colnames must correspond to the names of the classes of the hierarchy (except the root node)
# g : a graph of class graphNEL. It represents the hierarchy of the classes. The names of the classes must correspond to the name of the matrix of scores.
# root: name of the class that it is the top-level (root) of the hierarchy (def:00)
# Output:
# a matrix with the scores of the classes corrected according to the TPR-TF algorithm.
tpr.threshold.free <- function(S,g, root="00"){
	levels <- levels.graph(g,root)
	# si aggiunge una "dummy root column" se non presente nella matrice degli scores flat S
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	# visita bottom-up del grafo
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
		if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
				child.set <- children[j,] > parent[j];  # scelta dei positive children
				child.pos <- children[j,][child.set];
				parent[j] <- (parent[j] + sum(child.pos))/(1+length(child.pos));  # correzione predizione flat
			}
		S[,names(chd.bup[i])] <- parent;
	  } 
	}
	# visita top-down del grafo
	par.tod <- get.parents.top.down(g,levels,root);
	for(i in 1:length(par.tod)){
		child <- S[,names(par.tod[i])];
		parents <- as.matrix(S[,par.tod[[i]]]);
		# colnames(parents) <- par.tod[[i]]
		# Note: the version with an apply and an ifelse statement is slower ...
		for(j in 1:length(child)){
			x <- min(parents[j,]);
			if(x < child[j]){
				child[j] <- x;   # hierarchy correction
			}
		}
		S[,names(par.tod[i])] <- child;
	}
	# the dummy root node is removed 
	S <- S[,-which(colnames(S)==root)];
	return(S);
}



# Weighted True Path Rule with Threshold (TPR-WT)
# Input:
# S : the flat scores matrix relative to a set of examples. 
#	  The colnames must correspond to the names of the classes of the hierarchy (except the root node)
# g : a graph of class graphNEL. It represents the hierarchy of the classes. The names of the classes must correspond to the name of the matrix of scores.
# root: name of the class that it is the top-level (root) of the hierarchy (def:00)
# t: threshold for the choice of positive children (def. = 0.5)
# w: weight to balance the contribution between the node i and its positive children (def. = 0.5)
# Output:
# a matrix with the scores of the classes corrected according to the TPR-WT algorithm.
tpr.weighted.threshold <- function(S,g, root="00", t=0.5,w=0.5){
	levels <- levels.graph(g,root)
	# si aggiunge una "dummy root column" se non presente nella matrice degli scores flat S
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	# visita bottom-up del grafo
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
	  	if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
			  	child.set <- children[j,] > t;    # scelta dei positive children in base alla soglia t
		  		child.pos <- children[j,][child.set];
		  		if(length(child.pos)!=0){
		 			parent[j] <- w*parent[j] + (1-w)*sum(child.pos)/length(child.pos);  # correzione predizione flat
		  		}
			}
			S[,names(chd.bup[i])] <- parent;
	  	}
	}
	# visita top-down del grafo
	par.tod <- get.parents.top.down(g,levels,root)
	for(i in 1:length(par.tod)){
		child <- S[,names(par.tod[i])];
	  	parents <- as.matrix(S[,par.tod[[i]]]);
		# colnames(parents) <- par.tod[[i]]
	  	# Note: the version with an apply and an ifelse statement is slower ...
	  	for(j in 1:length(child)){
			x <- min(parents[j,]);
			if(x < child[j]){
		  		child[j] <- x;    # hierarchy correction
			}
		}
		S[,names(par.tod[i])] <- child;
	}
	# the dummy root node is removed 
	S <- S[,-which(colnames(S)==root)];
	return(S);
}



# Weighted True Path Rule with Threshold-Fee (TPR-W)
# Input:
# S : the flat scores matrix relative to a set of examples. 
#	  The colnames must correspond to the names of the classes of the hierarchy (except the root node)
# g : a graph of class graphNEL. It represents the hierarchy of the classes. The names of the classes must correspond to the name of the matrix of scores.
# root: name of the class that it is the top-level (root) of the hierarchy (def:00)
# w: weight to balance the contribution between the node i and its positive children (def. = 0.5)
# Output:
# a matrix with the scores of the classes corrected according to the TPR-W algorithm.
tpr.weighted.threshold.free <- function(S,g, root="00", w=0.5){
	levels <- levels.graph(g,root);
	# si aggiunge una "dummy root column" se non presente nella matrice degli scores flat S
	if(!(root %in% colnames(S))){
	   max.score <- max(S);
	   z <- rep(max.score,nrow(S));
	   S <- cbind(z,S);
	   colnames(S)[1] <- root;
	}
	# visita bottom-up del grafo
	chd.bup <- get.children.bottom.up(g,levels);
	for(i in 1:length(chd.bup)){
	  	if(length(chd.bup[[i]])!=0){
			parent <- S[,names(chd.bup[i])];
			children <- as.matrix(S[,chd.bup[[i]]]);
			# colnames(children) <- chd.bup[[i]]
			for(j in 1:length(parent)){
			 	child.set <- children[j,] > parent[j];    # scelta dei positive children
			  	child.pos <- children[j,][child.set];
		  		if(length(child.pos)!=0){
			 		parent[j] <- w*parent[j] + (1-w)*sum(child.pos)/length(child.pos);  # correzione predizione flat
		  		}
			}
			S[,names(chd.bup[i])] <- parent;
	  	}
	}
	# visita top-down del grafo
	par.tod <- get.parents.top.down(g,levels,root)
	for(i in 1:length(par.tod)){
		child <- S[,names(par.tod[i])];
	  	parents <- as.matrix(S[,par.tod[[i]]]);
		# colnames(parents) <- par.tod[[i]]
	  	# Note: the version with an apply and an ifelse statement is slower ...
	  	for(j in 1:length(child)){
			x <- min(parents[j,]);
			if(x < child[j]){
		   		child[j] <- x;    # hierarchy correction
			}
	 	}
		S[,names(par.tod[i])] <- child;
	}
	# the dummy root node is removed 
	S <- S[,-which(colnames(S)==root)];
	return(S);
}


##***********************************************##
## TPR variant with DESCENDANTS instead children ##
##***********************************************##

# True Path Rule with Threshold -- VARIANT WITH DESCENDANTS (TPR-T-desc) 
# Input:
# S : the flat scores matrix relative to a set of examples. 
#	  The colnames must correspond to the names of the classes of the hierarchy (except the root node)
# g : a graph of class graphNEL. It represents the hierarchy of the classes. The names of the classes must correspond to the name of the matrix of scores.
# root: name of the class that it is the top-level (root) of the hierarchy (def:00)
# t: threshold for the choice of positive children (def. = 0.5)
# Output:
# a matrix with the scores of the classes corrected according to the TPR-T-desc algorithm.
tpr.threshold.descendants <- function(S,g, root="00", t=0.5){
	levels <- levels.graph(g,root)
	# si aggiunge una "dummy root column" se non presente nella matrice degli scores flat S
	if(!(root %in% colnames(S))){
	   max.score <- max(S);
	   z <- rep(max.score,nrow(S));
	   S <- cbind(z,S);
	   colnames(S)[1] <- root;
	}
	# visita bottom-up del grafo
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			node.curr <- nodes[i];
		  	if(length(desc.bup[[i]])!=0){
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				# if(i<=10) cat(node.curr, length(tmp), length(desc.bup[[i]]), "****************", sep="\n");
				desc <- as.matrix(S[,tmp]);
				# colnames(desc) <- desc.bup[[i]]
				for(j in 1:length(parent)){
			  		desc.set <- desc[j,] > t;    # scelta dei "positive" descendants in base alla soglia t
					desc.pos <- desc[j,][desc.set];
					parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # correzione predizione flat
				}
				S[,names(desc.bup[i])] <- parent;
		  	}
		}
		# visita top-down del grafo
		par.tod <- get.parents.top.down(g,levels,root)
		for(i in 1:length(par.tod)){
		  	child <- S[,names(par.tod[i])];
		  	parents <- as.matrix(S[,par.tod[[i]]]);
		  	# colnames(parents) <- par.tod[[i]]
		  	# Note: the version with an apply and an ifelse statement is slower ...
		  	for(j in 1:length(child)){
				x <- min(parents[j,]);
				if(x < child[j]){
			   		child[j] <- x;    # hierarchy correction
				}
		 	}
		  	S[,names(par.tod[i])] <- child;
		}
	# the dummy root node is removed if dummynode=TRUE
	S <- S[,-which(colnames(S)==root)];
	return(S);
}


# True Path Rule Threshold-Free -- VARIANT WITH DESCENDANTS (TPR-T-desc) 
# Input:
# S : the flat scores matrix relative to a set of examples. 
#	  The colnames must correspond to the names of the classes of the hierarchy (except the root node)
# g : a graph of class graphNEL. It represents the hierarchy of the classes. The names of the classes must correspond to the name of the matrix of scores.
# root: name of the class that it is the top-level (root) of the hierarchy (def:00)
# Output:
# a matrix with the scores of the classes corrected according to the TPR-TF algorithm.
tpr.threshold.free.descendants <- function(S,g, root="00"){
	levels <- levels.graph(g,root)
	# si aggiunge una "dummy root column" se non presente nella matrice degli scores flat S
	if(!(root %in% colnames(S))){
	   max.score <- max(S);
	   z <- rep(max.score,nrow(S));
	   S <- cbind(z,S);
	   colnames(S)[1] <- root;
	}
	# visita bottom-up del grafo
	desc.bup <- build.descendants.bottom.up(g,levels);
	nodes <- names(desc.bup);
		for(i in 1:length(desc.bup)){
			node.curr <- nodes[i];
		  	if(length(desc.bup[[i]])!=0){
				parent <- S[,names(desc.bup[i])];
				tmp <- setdiff(desc.bup[[i]],node.curr);
				# if(i<=10) cat(node.curr, length(tmp), length(desc.bup[[i]]), "****************", sep="\n");
				desc <- as.matrix(S[,tmp]);
				# colnames(desc) <- desc.bup[[i]]
				for(j in 1:length(parent)){
			  		desc.set <- desc[j,] > parent[j];
					desc.pos <- desc[j,][desc.set];
					parent[j] <- (parent[j] + sum(desc.pos))/(1+length(desc.pos));   # correzione predizione flat
				}
				S[,names(desc.bup[i])] <- parent;
		  	}
		}
		# visita top-down del grafo
		par.tod <- get.parents.top.down(g,levels,root)
		for(i in 1:length(par.tod)){
		  	child <- S[,names(par.tod[i])];
		  	parents <- as.matrix(S[,par.tod[[i]]]);
		  	# colnames(parents) <- par.tod[[i]]
		  	# Note: the version with an apply and an ifelse statement is slower ...
		  	for(j in 1:length(child)){
				x <- min(parents[j,]);
				if(x < child[j]){
			   		child[j] <- x;    # hierarchy correction
				}
		 	}
		  	S[,names(par.tod[i])] <- child;
		}
	# the dummy root node is removed if dummynode=TRUE
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

