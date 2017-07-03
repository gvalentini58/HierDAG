
## What is HierDAG library?

**HierDAG** is an R software library implementing two hierarchical ensemble methods for Directed Acyclic Graphs (DAGs):<br>
	1. Hierarchical Top-Down for DAG (**HTD-DAG**);<br>
	2. True-Path-Rule for DAG (**TPR-DAG**);<br>

## Files contained in the **HierDAG** library: brief description

1. Functions implementing Hierarchical Ensemble Methods for DAGs:
	* _htd.R_: 	implementation of the HTD algorithm;
	* _tpr.R_: 	implementation of the TPR algorithm and its variants;
	* _Do.HTD.R_: 	high level function to compute hierarchical correction according to HTD algorithm. 
	* _Do.TPR.R_:	high level functions to compute hierarchical correction according to TPR-DAG algorithm and its variants;

2. Utility Functions:
	* _graph.utils.R_:	utility functions to process and analyze a graph;
	* _IOgraph.R_:		IO functions to store and build graph both in plain text and in rda compressed format;
	* _flat.score.norm.R_:		function to normalize the flat scores according to the maximum score of each class;
	* _Do.full.annotations.table.R_:	high level functions to compute the full annotation table;
	* _Do.best.F.score_:	high level functions to select the best F-score by choosing an appropriate threshold in a scores matrix;
	* _Do.FLAT.scores.normalization_:	high level functions to normalize flat scores matrix w.r.t. MaxNorm or Qnorm;
	* _F.hier.R_:	function to compute precision, recall, F-measure, specificity and accuracy for multi-class multi-label classification task
	(Kiritchenko-like multi-label F-scores)
	* _AUPROC.R_:	function to compute AUROC and AUPRC through the R package precrec

## Loading 
Even if HTD-DAG and TPR-DAG algorithm depend on several other source codes, loading them in the R environment is straightforward. For instance to load the high-level function implementing the HTD-DAG algorithm, just open the R environment in the same folder where you download the **HierDAG** library and type:
```
source("Do.HTD.R");
```

**NB**: to rightly load the **HierDAG** library, the following R libraries are required:
```
library("graph");
library("RBGL");
library("Rgraphviz");
library("PerfMeas"); 		
library("precrec");			
library("preprocessCore"); 	
```



