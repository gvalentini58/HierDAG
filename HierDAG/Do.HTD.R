# High level function to compute hierarchical correction according to HTD algorithm. 

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
# April 2016: modified Do.HTD: normalization methods added
# July 2016: added PRC computed by precrec package; 
# July 2016: added Do.HTD.holdout
# November 2016: added AUC computed by precrec pkg; 
# November 2016: added b.per.example as input parameter
# January 2017: RAM memory problem: see lines 109-111, after 'Hierarchical Top Down Correction': 
# March 2017: figure RAM memory prloblem out!! see AUPROC.R source file for explanation..

##***************************************************************************************##
## library to call 
library("PerfMeas"); 				## compute PXR
library("precrec");					## compute AUPRC and AUROC
library("preprocessCore"); 			## Qnorm
source("htd.R");					## HTD-DAG algorithm
source("flat.score.norm.R");		## Maxnorm
source("graph.utils.R");			## graph utility functions 
source("F-hier.R");					## compute Kiritchenko-like multi-label F-scores
source("Do.flat.normalization.R");	## high level function to compute Maxnorm e Qnorm
source("AUPROC.R");					## functions to compute AUPRC and AUROC 
##***************************************************************************************##

# INPUT: 
# norm: boolean value: 1.TRUE means that the flat scores matrix has been already normalized in according to a normalization method (def);
#       			   2.FALSE means that the flat scores matrix has NOT been normalized yet.
# norm.type: this variable can assume three values: MaxNorm, Qnorm, NONE. We have two case respect to norm:
#			 1.if norm==FALSE, two kind of normalizations are possible:	1. MaxNorm: each score is divided w.r.t. the max of each class 
#																		2. Qnorm: quantile normalization is applied. PreprocessCore library is used.
#			 2.if norm==TRUE, set norm.type=="NONE" (def);
# flat.file: name of flat scores matrix already normalized or to be normalized in according to norm.type (without rda extension);
# ann.file: name of the target labels (without rda extension). It must be  an .rda file containing the label matrix of the examples (def: ann.file)
# dag.file: name of the graph that represents the hierarchy of the classes (def dag.file) 
# flat.dir: relative path to folder where flat scores matrix (already normalized or to normalize) is stored (def flat.dir)
# ann.dir: relative path to folder where annotation matrix is stored
# dag.dir: relative path to folder where graph is stored
# flat.norm.dir: 1.if norm=FALSE, relative path where flat normalized scores matrix will be strored;
#				 2.if norm=TRUE, the flat scores matrix is already normalized, than it is set to NULL (def)
# b.per.example : if TRUE (def.) precision, recall, F-measure, specificity and accuracy are returned for each example, 
#                 otherwise only the averages are returned.
# n.round: 	number of rounding digits to be applied to the hierarchical scores matrix (def. 3). 
#			It's used for choosing the best threshold on the basis of the best F.measure (see f.criterion parameter).
# f.criterion: character. Type of F-measure to be used to select the best F.measure. There are 2 possibilities: 
#              1. "F" (default) corresponds to the harmonic mean between the average precision and recall; 
#              2. "avF" corresponds to the per-example F-score averaged across all the examples.
# hierScore.dir: relative path to folder where the matrix with the scores of the classes corrected in according to HTD algorithm is stored 
# perf.dir: relative path to folder where the macro-centric measures (i.e. AUC, PRC and PxR across classes) and 
# 			where example-centric measures (i.e. Precision, Recall, Specificity, F-measure, Accuracy across example) are stored
# OUTPUT:
# 5 rda files stored in the rispective output directory:
# - Matrix with examples on rows and classes on colums representing the hierarchical scores of the classes computed with HTD algorithm.
#	Stored in hierScore.dir folder
# - Example-centric measures computed through find.best.f from F-hier.R file. Stored in perf.dir folder
# - AUC (average and per classes) computed through AUC.single.over.classes from package PerfMeas. Stored in perf.dir
# - Precision at fixed recall levels (average and per classes) computed through precision.at.multiple.recall.level.over.classes from package PerfMeas.
# 	Stored in perf.dir folder
# - PRC (average and per.class) computed by precrec package. Stored in perf.dir
Do.HTD <- function	(	norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file,
						flat.dir=flat.dir, ann.dir=ann.dir, dag.dir=dag.dir, flat.norm.dir=NULL, b.per.example=TRUE,
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
		# gc();	##in order to save ram memory..
		
		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir)
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table 
	ann <- ann[,-which(colnames(ann)==root)];

	## Computing FLAT Performances
	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); gc();
	
	## FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S); gc();

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example); gc();

	## FLAT PRC computed by precrec package (more precise and accurate than PerfMeas)
	PRC.flat <- AUPRC.single.over.classes(ann, S); gc();

	## Hierarchical Top Down Correction ####################
	## in this way we fill memory because we store two double-float matrix. Solution overwrite!! we have already calculated the flat performances..
	# S.htd <- htd(S,g,root);
	S <- htd(S, g, root);
	
	## Computing Hier Performances
	## Hierarchical AUC (average and per.class) computed by precrec package
	AUC.htd <- AUROC.single.over.classes(ann, S); gc();

	## Hierarchical PxR at fixed recall levels 
	PXR.htd <- precision.at.multiple.recall.level.over.classes(ann, S); gc();

	## Computing Hierarchical Examples-Measures 
	FMM.htd <- find.best.f(ann, S, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=b.per.example); gc();

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.htd <- AUPRC.single.over.classes(ann, S); gc();

	## storing the hierarchical matrix
	S.htd <- S;
	rm(S);

	## Storing Results #########
	if(norm){
		save(S.htd, file=paste0(hierScore.dir, flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(AUC.flat, AUC.htd, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PXR.flat, PXR.htd, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(FMM.flat, FMM.htd, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.htd, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}else{
		save(S.htd, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(AUC.flat, AUC.htd, file=paste0(perf.dir, "AUC.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(PXR.flat, PXR.htd, file=paste0(perf.dir, "PXR.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(FMM.flat, FMM.htd, file=paste0(perf.dir, "PCM.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.htd, file=paste0(perf.dir, "PRC.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}
}

# High level function to correct the scores with a hierarchy according to HTD algorithm performing a classical holdout procedure
# All input paramenters are the same of the function above, except the following:
# ind.test.set: vector of integer. Indices refer to the examples of the adjancency matrix to be used in the test set. 
# ind.dir: relative path to folder where ind.test.set is stored
Do.HTD.holdout <- function	(	norm=TRUE, norm.type= "NONE", flat.file=flat.file, ann.file=ann.file, dag.file=dag.file,
								ind.test.set=ind.test.set, ind.dir=ind.dir, flat.dir=flat.dir, ann.dir=ann.dir,  
								dag.dir=dag.dir, flat.norm.dir=NULL, n.round=3, f.criterion ="F", b.per.example=TRUE,
								hierScore.dir="hierScore.dir/", perf.dir="perf.dir/"
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
	if(norm){
		S <- get(load(flat.path));
		gc();	##in order to save ram memory..

		## removing root node from flat norm matrix if it exists
		if(root %in% colnames(S)){
			S <- S[,-which(colnames(S)==root)];
		}
	}else{
		Do.FLAT.scores.normalization(norm.type= norm.type, flat.file=flat.file, flat.norm.dir=flat.norm.dir)
		flat.path <- paste0(flat.norm.dir, norm.type,".",flat.file,".rda");
		S <- get(load(flat.path));
	}

	## shrinking the size of S to the examples of test set
	S <- S[ind.test,];

	## loading annotation matrix
	ann.path <- paste0(ann.dir, ann.file,".rda");
	ann <- get(load(ann.path));
	gc();

	## removing root node from annotation table and shrinking the size of annotation table to the examples of test set
	ann <- ann[ind.test,-which(colnames(ann)==root)];

	## Computing FLAT Performances
	## FLAT AUC computed by precrec package
	AUC.flat <- AUROC.single.over.classes(ann, S); 	gc();
			
	## FLAT PxRs computed by PerfMeas pacakge
	PXR.flat <- precision.at.multiple.recall.level.over.classes(ann, S); gc();

	## F.measure: Computing Flat Examples-Measures 
	FMM.flat <- find.best.f(ann, S, n.round=n.round, f.criterion=f.criterion, verbose=FALSE, b.per.example=b.per.example); gc();

	## FLAT PRC computed by precrec package 
	PRC.flat <- AUPRC.single.over.classes(ann, S);  gc();

	## Hierarchical Top Down Correction ####################
	## in this way we fill memory because we store two double-float matrix. Solution overwrite!! we have already calculated the flat performances..
	# S.htd <- htd(S,g,root);
	S <- htd(S,g,root);
	
	## Computing Hier Performances
	## Hierarchical AUC (average and per.class) computed by precrec package
	AUC.htd <- AUROC.single.over.classes(ann, S); gc();
		
	## Hierarchical PxR at fixed recall levels 
	PXR.htd <- precision.at.multiple.recall.level.over.classes(ann, S); gc();

	## Computing Hierarchical Examples-Measures 
	FMM.htd <- find.best.f(ann, S, n.round=n.round, f.criterion =f.criterion, verbose=FALSE, b.per.example=b.per.example);	gc();

	## Hierarchical PRC (average and per.class) computed by precrec package
	PRC.htd <- AUPRC.single.over.classes(ann, S);  gc();

	## storing the hierarchical matrix
	S.htd <- S;
	rm(S);
	
	## Storing Results #########
	if(norm){
		save(S.htd, file=paste0(hierScore.dir, flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(AUC.flat, AUC.htd, file=paste0(perf.dir, "AUC.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PXR.flat, PXR.htd, file=paste0(perf.dir, "PXR.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(FMM.flat, FMM.htd, file=paste0(perf.dir, "PCM.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.htd, file=paste0(perf.dir, "PRC.", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}else{
		save(S.htd, file=paste0(hierScore.dir, norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(AUC.flat, AUC.htd, file=paste0(perf.dir, "AUC.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(PXR.flat, PXR.htd, file=paste0(perf.dir, "PXR.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);	
		save(FMM.flat, FMM.htd, file=paste0(perf.dir, "PCM.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
		save(PRC.flat, PRC.htd, file=paste0(perf.dir, "PRC.", norm.type,".", flat.file, ".hierScores.htd.rda"), compress=TRUE);
	}
}
