########################################################
## Hierarchical Top-Down for Directed Acyclic Graphs) ## 
########################################################

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

# April 2014
# July 2014

## library to call
library("graph");
library("RBGL");
#library(Rgraphviz);
source("graph.utils.R"); 



# It implements a top-down procedure to correct the computed scores in the hierarchy according to the constraints that the score 
# of a  node cannot be greater than a score of its parents. 
# The algorithm can be applied to any DAG structured hierarchy of classes
# Input:
# S : the flat scores matrix relative to a set of examples. 
#	  The colnames must correspond to the names of the classes of the hierarchy (except the root node)
# g : a graph of class graphNEL. It represents the hierarchy of the classes. The names of the classes must correspond to the name of the matrix of scores.
# root: name of the class that it is the top-level (root) of the hierarchy (def:00)
# Output:
# a matrix with the scores of the classes corrected according to the HTD algorithm.
htd <- function(S,g, root="00"){
	levels <- levels.graph(g,root);
	# si aggiunge una "dummy root column" se non presente nella matrice degli scores flat S
	if(!(root %in% colnames(S))){
		max.score <- max(S);
		z <- rep(max.score,nrow(S));
		S <- cbind(z,S);
		colnames(S)[1] <- root;
	}
	# nodes are scanned from top to bottom: a list par.tod with the parents for each node (ordered from top to bottom) is obtained	
	par.tod <- get.parents.top.down(g,levels,root)
	for(i in 1:length(par.tod)){  								    
		child <- S[,names(par.tod[i])]; 	# score del nodo figlio					    
		parents <- as.matrix(S[,par.tod[[i]]]);	# score dei nodi genitori  			    
		# colnames(parents) <- par.tod[[i]];
		# Note: the version with an apply and an ifelse statement is slower ...						    
		for(j in 1:length(child)){ 							    
			x <- min(parents[j,]);								    
			if(x < child[j])									    
				child[j] <- x;									    
		}  													    
		S[,names(par.tod[i])] <- child;						    
	}													    
	# the dummy root column is removed
	S <- S[,-which(colnames(S)==root)];
	return(S);
}

