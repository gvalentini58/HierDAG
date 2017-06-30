# High-level function to get the *FULL* HPO annotation matrix (i.e. labels matrix) w.r.t a given weighted adiacency matrix

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

## library to call:
source("graph.utils.R")

# Input:
# anc.file.name: name of the list of ancestors to load (without rda extension). 
#				 It must be an rda file containing for each node the list of all its ancestor (def: anc.name) 
# anc.dir: relative path to directory where the ancestor file is stored (def: anc.dir)
# net : name of the dataset to load (without rda extension). It must be an .rda file containing the weighted adjiacency matrix of the graph. (def: net)
# net.dir: relative path to directory where the weighted adjiacency matrix is stored (def: net.dir)
# ann.file.name: name of the most specific annotations matrix (without rda extension)
# ann.dir: relative path to directory where the most specific annotation matrix is stored (def: ann.dir)
# output.name : name of the output file without rda extension (def: output.name)
# output.dir: relative path to directory where the final HPO annotation matrix is stored  (def: output)
# Output:
# a matrix with annotations gene-phenotype stored in an .rda file in "/output" directory. 
# Rows correspond to genes and columns to phenotpyes. Let's denote "ann" the annotation matrix, the entries ann[i,j] 
# correspond to the i-th gene and j-th phenotype. ann[i,j]=1 means that gene i is annotated with phenotype j, O otherwise.
Do.full.annotation.matrix <- function(	anc.file.name=anc.file.name, anc.dir="/anc.dir" ,net=net, net.dir="/net.dir", 
										ann.file.name=ann.file.name, ann.dir="/ann.dir", output.name=output.name, output.dir="/output"){
	## loading list of ancestors
	anc.path <- paste0(anc.dir, anc.file.name, ".rda");
	anc.name <- load(anc.path);
	anc <- eval(parse(text=anc.name));  

	## loading wadj matrix
	net.path <- paste0(net.dir, net, ".rda");
	net.name <- load(net.path);
	W <- eval(parse(text=net.name));  

	## loading the specific annotation matrix
	ann.path <- paste0(ann.dir, ann.file.name, ".rda");
	ann.name <- load(ann.path);
	ann.spec <- eval(parse(text=ann.name));

	## costruction of full HPO annotation matrix
	ann <- full.annotation.matrix(W=W, anc=anc, ann.spec=ann.spec);
	
	## saving labels matrix
	ann.file <- paste0(output.dir, output.name, ".rda");
	save(ann, file=ann.file, compress=TRUE);
}

