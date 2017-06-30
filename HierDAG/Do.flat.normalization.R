# High level function to normalize flat scores matrix w.r.t. MaxNorm or Quantile Normalization (Qnorm)

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

#April 2016

## library to call
source("graph.utils.R")
source("flat.score.norm.R")

# input: 
# norm.type: it must be a normalization method:
#			 - MaxNorm: each score is divided w.r.t. the max of each class
# 			 - Qnorm: a quantile normalization is applied. Library preprocessCore is used.
# flat.file: name of the flat scores matrix (without rda extension). 
# output:
# the matrix of the scores flat normalized w.r.t. MaxNorm or Qnorm
Do.FLAT.scores.normalization <- function(norm.type= "MaxNorm", flat.file=flat.file, flat.norm.dir=flat.norm.dir){
	## normalization
	if(norm.type=="MaxNorm"){
		## Max Normalization
		S <- normalize.max(S);		
	}else{
		## Quantile Normalization
		S <- normalize.quantiles(S);
		dimnames(S) <- list(rownames(S),colnames(S));	
	}
	## Storing results
	save(S, file=paste0(flat.norm.dir, norm.type, ".", flat.file, ".rda"));
}

