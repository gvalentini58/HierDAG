# High level function to select the best F-score by choosing an appropriate threshold on scores matrix. 

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

# April 2016

## library to call 
source("graph.utils.R");
source("F-hier.R");

# input: 
# score.file: name of the scores matrix (without rda extension). 
# ann.file: name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples (def: ann.file)
# dag.file: name of the graph that represents the hierarchy of the classes. 
# n.round: 	number of rounding digits to be applied to the hierarchical scores matrix (def. 3). 
#			It's used for choosing the best threshold on the basis of the best F.measure (see f.criterion parameter).
# f.criterion: character. Type of F-measure to be used to select the best F.measure. There are 2 possibilities: 
#              1. "F" (default) corresponds to the harmonic mean between the average precision and recall; 
#              2. "avF" corresponds to the per-example F-score averaged across all the examples.
# verbose: boolean. If TRUE the number of iterations are printed on stdout, FALSE (def) otherwise.
# score.dir: relative path to folder where normalized scores matrix is stored
# ann.dir: relative path to folder where annotation matrix is stored
# dag.dir: relative path to folder where graph is stored
# Fmeas.dir: relative path to folder where example-centric measures (i.e. Precision, Recall, Specificity, F-measure, Accuracy across example) are stored
# output:
# an rda files stored in macro.dir folder contains Example-centric measures computed through find.best.f from F-hier.R file
Do.best.F.score <- function(score.file=score.file, ann.file=ann.file, dag.file=dag.file,
							n.round=3, f.criterion ="F", verbose=FALSE, b.per.example=TRUE,
							score.dir=score.dir, ann.dir=ann.dir, dag.dir=dag.dir, Fmeas.dir="Fmeas.dir/"
							){
	
	## Loading Data ############
	## loading dag
	dag.path <- paste0(dag.dir, dag.file,".rda");
	g <- get(load(dag.path));
	
	##root node
	root <- root.node(g);

	## loading flat scores matrix relative to a specific subontology
	flat.path <- paste0(flat.dir, flat.filescore.file,".rda");
	S <- get(load(flat.path));
	gc();	##in order to save ram memory..
	
	## removing root node from flat matrix if it exists
	if(root %in% colnames(S)){
		S <- S[,-which(colnames(S)==root)];
	}

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## Computing Hierarchical Examples-Measures 
	F.meas <- find.best.f(ann, S, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=b.per.example);	

	## Storing Results #########
	save(F.meas, file=paste0(Fmeas.dir,"PCM.",score.file,".hierScoreS.rda"), compress=TRUE);
}

