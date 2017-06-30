# Function to normalize the scores of a flat scores matrix per class

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

# The scores of each class are normalized by dividing the score values for the maximum score of that class.
# If the max score of a class is zero, no normalization is needed, otherwise NaN value will be printed as results of 0 out of 0 division.
# Input:
# S : matrix with the raw non normalized scores. Rows are examples (proteins), columns are the classes (GO terms)
# Output:
# A score matrix with the same dimensions of S, but with scores max/normalized separately for each class.
normalize.max <- function(S){
	classes <- colnames(S);
	maximum <- apply(S,2,max);
	for(class in classes){
		if(maximum[class] != 0){
	 		S[,class] <- S[,class]/maximum[class];
	 	}
	}
	return(S);
}

###########
## BUG: the function below does not take in account if all the values of a specific class of the flat scores matrix are equal to ZERO (then the max of 
## the that class will be ZERO). Then NaN value will be returned as result of 0 out of 0 division. 
## This is likely to happen with FILTERED flat scores matrix without removing singleton nodes (e.g. STRING10_FILTERED) 
normalize.max.BUG <- function(S){
	max.class <- apply(S,2,max);
	for(i in 1:ncol(S)){
	 	S[,i] <- S[,i]/max.class[i]
	}
	return(S);
}

