# Utility functions to process and analyze graphs

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

# July 2014: added levels.graph, compute.flipped.graph, get.parents.top.down, get.children.bottom.up, get.parents, build.children
# August 2014: added build.descendants, build.descendants.per.level, build.ancestors, build.ancestors.per.level
# September 2014: added get.parents.bottom.up, get.children.top.down, build.ancestors.bottom.up and build.descendants.bottom.up
# 				  added check.hierarchy.single.sample, constraints.matrix
# February 2016: added weighted.adjacency.matrix, specific.annotation.list, specific.annotation.matrix, full.annotation.matrix, 
# 				 check.annotation.matrix.integrity, check.DAG.integrity, do.submatrix, do.subgraph 
# April 2016: correct bug in weighted.adjacency.matrix, add do.unstratified.cv.data
# April 2016: added transitive.closure.annotations and modified full.annotation.matrix: might happen that there are same HPO IDs that 
#			  are classified as "obsolete" in obo file, but that still exist in the annotation file (e.g. hpo-DAG Build #1700 and hpo-ann Build #110)
# May 2016: added do.my.subgraph. N.B.: do.my.subgraph is 0.6 sec slower than do.subgraph ( Tested on 2016-04-01 HPO obo and annotation release; 
#			see the note below the function). 
# May 2016: root.node and find.leaves up to date and new function root.or.leaf added. Briefly, leaves are nodes with outdegree equal to zero (i.e. nodes 
#			having no out-coming edges); root is the node with indegree equal to zero (i.e. node having no in-coming edge).
# December: do.unstratified.cv.data, do.stratified.cv.data.single.class, do.stratified.cv.fold

## library to call
library("graph");
library("RBGL");
#library("Rgraphviz");


# Function to group a set of nodes in according to their maximum depth in the graph. It simply inverts the weights of the graph 
# and then applies the Bellman Ford algorithm to find the shortest path, achieving in this way the longest path 
# Input: 
# g: a graph of class graphNEL 
# root : name of
# the root node (def. 00) 
# Output: 
# a list of the nodes grouped w.r.t. the distance from the root: the first element of the list corresponds to the root node (level 0),
# the second to nodes at maximum distance 1 (level 1), the third to the node at maximum distance 3 (level 2) and so on.
levels.graph <- function(g,root="00") {
	ed <- edges(g);
	ew <- edgeWeights(g);
	for(i in 1:length(ed)){
	  l <- length(ew[[i]]);
	  if(l!=0)
		ew[[i]][1:l] <- -1;
	}
	edL <- vector(mode="list", length=length(ed));
	names(edL) <- names(ed);
	for(i in 1:length(ed)){
		edL[[i]] <- list(edges=ed[[i]], weights=ew[[i]]);
	}
	G <- graphNEL(nodes=nodes(g), edgeL=edL, edgemode="directed");  # grafo con i pesi degli archi negativi
	depth.G <- bellman.ford.sp(G,root)$distance;
	depth.G <- -depth.G  # si settano le distanze positive
	levels <- vector(mode="list", length=max(depth.G)+1);
	names(levels) <- paste(rep("level", max(depth.G)+1), 0:max(depth.G), sep="_");
	for(i in 1:(max(depth.G)+1)){
		levels[[i]] <- names(which(depth.G==i-1));  # il livello 0 corrisponde al nodo root
	}
	return(levels);
}



# It computes a directed graph with edges in the opposite direction.
# Input:
# g : a graphNEL directed graph
# Output:
# a graph with edges in the opposite direction w.r.t. to g.
compute.flipped.graph <- function(g){
	ed <- edges(g);
	ndL <- vector(mode="list", length=length(ed));
	names(ndL) <- names(ed);
	for(i in 1:length(ed)){
		children <- ed[[i]];
		parent   <- names(ed[i]);
		if(length(children)!=0){
			for(j in 1:length(children)){
				ndL[[children[j]]] <- c(ndL[[children[j]]],parent); # memorizzo per ogni nodo i rispettivi genitori
			}
		}
	}
	for (i in 1:length(ndL)){
		ndL[[i]] <- list(edges=ndL[[i]]);
	}
	og <- graphNEL(nodes=nodes(g), edgeL=ndL, edgemode="directed");
	return(og);
}


# N.B.: computa la stessa cosa di compute.flipped.graph: da valutare qual e' la versione piu' veloce ...
# It computes a directed graph with edges in the opposite direction.
# Input:
# g : a graphNEL directed graph
# Output:
# a graph with edges in the opposite direction w.r.t. to g.
compute.opposite.graph  <- function(g){
   # if ( !identical(g@edgemode,"directed"))
   #     stop("compute.opposite.graph:  the graph must be directed");
	v <- nodes(g);
	edge.list <- vector("list", length = length(v));
	names(edge.list) <- v;
	#browser();    
	l <- adj(g, v);
	for (i in 1:length(l)) {
	   x <- l[[i]];
	   if (length(x)>0){
		   for (j in 1:length(x)){
				edge.list[[x[j]]] <- c(edge.list[[x[j]]], names(l[i]));
			}
		}
	}
	for (i in 1:length(edge.list)){
		edge.list[[i]] <- list(edges=edge.list[[i]]);
	}
	g.opposed <- graphNEL(nodes = v, edgeL = edge.list, edgemode = "directed");
	return(g.opposed);
}




# Function to compute the parents of each node of a graph
# Input:
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# root : name of the root node (def. 00)
# Output:
# a named list of character vectors. Each component corresponds to a node x of the graph (i.e. child node) and its vector is 
# the set of its parents (the root node is not included)
get.parents <- function(g, root="00") {
	nd <- nodes(g)
	ndL <- vector(mode="list", length=length(nd));
	names(ndL) <- nd;
	ed <- edges(g);
	for(i in 1:length(ed)){
		children <- ed[[i]];
		parent   <- names(ed[i]);
		if(length(children)!=0){
			for(j in 1:length(children)){
				ndL[[children[j]]] <- c(ndL[[children[j]]],parent); # memorizzo per ogni nodo i rispettivi genitori
			}
		}
	}
	ndL <- ndL[-which(names(ndL)==root)]; # si elimina il nodo root
	return(ndL);
}



# Function to compute the parents of each node of a graph by a per-level visiting from top to bottom
# Input:
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# levels: a list of character vectors. Each component represents a graph level and the elements of any component correspond to nodes.
#	      The level 0 coincides with the root node. It can be obtained by calling levels.graph function.
# root : name of the root node (def. 00)
# Output:
# a named list of character vectors. Each component corresponds to a node x of the graph (i.e. child node) and its vector is 
# the set of its parents. The nodes order follows the levels of the graph from root (excluded) to leaves. 
get.parents.top.down <- function(g,levels, root="00") {
	ord.nd <- unlist(levels); # ordered character vector con i nomi dei nodi ordinati dalla root (inclusa) ai nodi foglia
	ndL <- vector(mode="list", length=length(ord.nd));
	names(ndL) <- ord.nd;
	ed <- edges(g);
	for(i in 1:length(ed)){
		children <- ed[[i]];
		parent   <- names(ed[i]);
		if(length(children)!=0){
			for(j in 1:length(children)){
				ndL[[children[j]]] <- c(ndL[[children[j]]],parent); # memorizzo per ogni nodo i rispettivi genitori
			}
		}
	}
	ndL <- ndL[-which(names(ndL)==root)]; # si elimina il nodo root
	return(ndL);
}



# Function to compute the parents of each node of a graph by a per-level visiting from bottom to top
# Input:
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# levels: a list of character vectors. Each component represents a graph level and the elements of any component correspond to nodes.
#	      The level 0 coincides with the root node. It can be obtained by calling levels.graph function.
# root : name of the root node (def. 00)
# Output:
# a named list of character vectors. Each component corresponds to a node x of the graph (i.e. child node) and its vector is 
# the set of its parents. The nodes are ordered from leaves to root (excluded). 
get.parents.bottom.up <- function(g,levels, root="00") {
	flip.ord.nd <- rev(unlist(levels)); # ordered character vector con i nomi dei nodi ordinati dai nodi foglia alla root (inclusa)
	ndL <- vector(mode="list", length=length(flip.ord.nd));
	names(ndL) <- flip.ord.nd;
	ed <- edges(g);
	for(i in 1:length(ed)){
		children <- ed[[i]];
		parent   <- names(ed[i]);
		if(length(children)!=0){
			for(j in 1:length(children)){
				ndL[[children[j]]] <- c(ndL[[children[j]]],parent); # memorizzo per ogni nodo i rispettivi genitori
			}
		}
	}
	ndL <- ndL[-which(names(ndL)==root)]; # si elimina il nodo root
	return(ndL);
}



# Function to compute the descendants of each node of a graph
# Input :
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# Output:
# a named list of vectors. Each component corresponds to a node x of the graph, and its vector is the set of its descendants including also x.
build.descendants <- function(g) {
	name.nodes <- nodes(g);
	g2 <- transitive.closure(g);
	desc <- edges(g2);
	for(x in name.nodes){
		desc[[x]] <- c(desc[[x]],x);
	}
	return(desc);
}



# Function to compute the descendants of each node of a graph by a per-level visiting, from top to bottom
# Input:
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# levels: a list of character vectors. Each component represents a graph level and the elements of any component correspond to nodes.
#	      The level 0 coincides with the root node. It can be obtained by calling levels.graph function.
# Output:
# a named list of vectors. Each component corresponds to a node x of the graph and its vector is the set of its descendants including also x.
# The nodes are ordered from root (included) to leaves.
build.descendants.per.level <- function(g,levels) {
	ord.nd <- unlist(levels);
	g2 <- transitive.closure(g);
	desc <- edges(g2)[ord.nd];
	for(x in ord.nd){
		desc[[x]] <- c(desc[[x]],x);
	}
	return(desc);
}



# Function to compute the descendants of each node of a graph by a per-level visiting, from bottom to top
# Input:
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# levels: a list of character vectors. Each component represents a graph level and the elements of any component correspond to nodes.
#	  The level 0 coincides with the root node. It can be obtained by calling levels.graph function.
# Output:
# a named list of vectors. Each component corresponds to a node x of the graph and its vector is the set of its descendants including also x.
# The nodes are ordered from root (included) to leaves.
build.descendants.bottom.up <- function(g,levels) {
	flip.ord.nd <- rev(unlist(levels));
	g2 <- transitive.closure(g);
	desc <- edges(g2)[flip.ord.nd];
	for(x in flip.ord.nd)
	  desc[[x]] <- c(desc[[x]],x);
	return(desc);
}



# Function to compute the children of each node of a graph
# Input :
# g : a graph of class graphNEL. It represents the hierarchy of the classes.
# Output:
# a named list of vectors. Each component corresponds to a node x of the graph, and its vector is the set of its children
build.children <- function(g) {
	return(edges(g));
}



# Function to compute the children of each node of a graph by a per-level visiting, from top to bottom
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# levels: a list of character vectors. Each component represents a graph level and the elements of any component correspond to nodes.
#	  	  The level 0 coincides with the root node. It can be obtained by calling levels.graph function.
# Output:
# a named list of character vectors. Each component corresponds to a node x of the graph (i.e. parent node) and its vector is 
# the set of its children. The nodes are ordered from root (included) to leaves.
get.children.top.down <- function(g,levels){
	child <- build.children(g)
	nd <- c();
	for(i in 1:length(levels)){
		level.nodes <- levels[[i]];
		nd <- append(nd,child[level.nodes]);
	}
	return(nd);
}



# Function to compute the children of each node of a graph by a per-level visiting, from bottom to top
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# levels: a list of character vectors. Each component represents a graph level and the elements of any component correspond to nodes.
#	  	  The level 0 coincides with the root node. It can be obtained by calling levels.graph function.
# Output:
# a named list of character vectors. Each component corresponds to a node x of the graph (i.e. parent node) and its vector is 
# the set of its children. The nodes are ordered from leaves (included) to root.
get.children.bottom.up <- function(g,levels){
	ed <- edges(g);  
	nd <- c();
	n.levels <- length(levels);
	for(i in n.levels:1){
	  level.nodes <- levels[[i]];
	  nd <- append(nd,ed[level.nodes]);
	}
	return(nd);
}



# Function to compute the ancestor of each node of a graph
# Input:
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# Output:
# a named list of vectors. Each component corresponds to a node x of the graph,
# and its vector is the set of its ancestors including also x.
build.ancestors <- function(g){
	og <- compute.flipped.graph(g);
	names.nodes <- nodes(og);
	og2 <- transitive.closure(og);
	anc <- edges(og2);
	for(x in names.nodes){
		anc[[x]] <- c(anc[[x]],x);
	}
	return(anc);
}

# MN: the same of above, but the input parameter 'itself' is added. 
# It's a boolean value: if TRUE in the output list the node x is included in the ancestors of node x, FALSE otherwise 
build.ancestors2 <- function(g,itself=TRUE){
	og <- compute.flipped.graph(g);
	names.nodes <- nodes(og);
	og2 <- transitive.closure(og);
	anc <- edges(og2);
	if(itself){
		for(x in names.nodes){
			anc[[x]] <- c(anc[[x]],x);
		}
		return(anc);
	}else{
		return(anc);
	}	
}

# Function to compute the ancestors of each node of a graph per level
# Input :
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# levels: a list of character vectors. Each component represents a graph level and the elements of any component correspond to nodes.
#         The level 0 coincides with the root node. It can be obtained by calling levels.graph function.
# root: name of the class that is on the top-level of the hierarchy (def:"00")
# Output:
# a named list of vectors. Each component corresponds to a node x of the graph and its vector is the set of its ancestors including also x.
# The nodes are ordered from root (included) to leaves.
build.ancestors.per.level <- function(g,levels, root="00"){
	og <- compute.flipped.graph(g);
	ord.nd <- unlist(levels);
	og2 <- transitive.closure(og);
	anc <- edges(og2)[ord.nd];
	for(x in ord.nd)
		anc[[x]] <- c(anc[[x]],x);
	return(anc);
}



# Function to compute the ancestors of each node of a graph from bottom to top
# Input :
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# levels: a list of character vectors. Each component represents a graph level and the elements of any component correspond to nodes.
#	  	  The level 0 coincides with the root node. It can be obtained by calling levels.graph function.
# root: name of the class that is on the top-level of the hierarchy (def:"00")
# Output:
# a named list of vectors. Each component corresponds to a node x of the graph and its vector is the set of its ancestors including also x.
# The nodes are ordered from leaves to root (included).
build.ancestors.bottom.up <- function(g,levels, root="00"){
	og <- compute.flipped.graph(g);
	flip.ord.nd <- rev(unlist(levels));
	og2 <- transitive.closure(og);
	anc <- edges(og2)[flip.ord.nd];
	for(x in flip.ord.nd)
	  anc[[x]] <- c(anc[[x]],x);
	return(anc);
}



# Function to find the node that is on the top-level of hierarchy (i.e. root node)
# Input:
# g : a graph of class graphNEL. It represents the hierarchy of the classes.
# Output:
# name of the root node 
root.node <- function(g){
	d <- degree(g);
	root <- names(which(d$inDegree==0));
	return(root);
}

root.node.old <- function(g){
	og <- compute.flipped.graph(g);
	root <- find.leaves(og);
	return(root);
}

# Function to compute the leaves of a graph
# Input :
# g : a graph of class graphNEL. It represents the hierarchy of the classes.
# Output:
# a vector with the names of the leaves of g
find.leaves <- function(g){
	d <- degree(g);
	leaves <- names(which(d$outDegree==0));
	return(leaves);
}

find.leaves.old <- function(g){
	child <- edges(g);
	leaves <- c();
	for(i in 1:length(child)){
	  if(length(child[[i]])==0)
		leaves <- append(leaves,names(child[i]));
	}
	return(leaves);
}


# Function to find root and leaves using in/out degree. The graph must be direct.
# Leaves are nodes with outdegree equal to zero; root is the node with indegree equal to zero.
# Input :
# g : a graph of class graphNEL. It represents the hierarchy of the classes. It must be a direct graph.
# deg: can assume two values: 	1. "indeg": nodes have no in coming edges are returned: root node is found.
#								2. "outdeg": nodes have no out coming edges are returned: leaves nodes are found.
# Output:
# a vector with the names of the leaves or the root name of g
root.or.leaf <- function(g,deg){
	d <- degree(g);
	if(deg=="indeg"){
		root <- names(which(d$inDegree==0));
		return(root);
	}
	if(deg=="outdeg"){
		leaves <- names(which(d$outDegree==0));
		return(leaves);
	}
}

# Function to compute the distances of each node of the graph from the leaves
# It returns the minimum distance of each node from one of the leaves of the graph
# Input :
# g : a graph of class graphNEL. It represents the hierarchy of the classes.
# Output:
# a named vector. The names are the names of the nodes of the graph g, and their values represent the distance from the leaves.
# A value equal to 0 is assigned to the leaves, 1 to nodes with distance 1 from a leaf and so on.
distances.from.leaves <- function(g){
	leaves <- find.leaves(g);
	n.leaves <- length(leaves);
	og <- compute.flipped.graph(g);
	og <- addNode("root", og)
	#   for (x in leaves)
	#    og = addEdge("root", x, og, 1);
	og <- addEdge(rep("root",n.leaves), leaves, og, rep(1,n.leaves));
	dist <- acc(og,"root")[[1]]-1;
	return(dist);
}


# It computes the pairwise constraints matrix which needs to have 2 columns and as many rows as constraints. Referring to DAG this matrix must
# define a partial order. The specification is always as follow: element column 2 \leq element column 1, e.g the first row states x_1 \leq x_2.
# It follows that the first column includes the set of the children nodes whereas the second column includes the set of parents nodes.
# In the constraints matrix the children and the parents nodes are specified by a number and not by their names.
# Input:
# g: a graph of class graphNELL. Represent the hierarchy of the class
# Output:
# a constraints matrix w.r.t the graph g
constraints.matrix <- function(g){
	eM <- t(edgeMatrix(g));
	eM <- cbind(eM[,2],eM[,1]);
	nd <- nodes(g);
	dimnames(eM) <- list(nd[eM[,2]], c("child","parent"))
	return(eM);
}

# It is equivalent to the previous one but no transposition is needed ... This could be faster on big graphs.
constraints.matrix2 <- function(g){
	eM <- edgeMatrix(g);
	eM <- cbind(eM[2,],eM[1,]);
	nd <- nodes(g);
	dimnames(eM) <- list(nd[eM[,2]], c("child","parent"))
	return(eM);
}

constraints.matrix.old <- function(g){
	edL <- edgeL(g);
	parent <- c();
	for(i in 1:length(edL)) {
	  x <- rep(i,length(edL[[i]][[1]]));
	  parent <- append(parent,x);
	}
	child <-c();
	for(i in 1:length(edL)){
	  child <- append(child,edL[[i]][[1]]);
	}
	eM <- cbind(child,parent);
	nd <- nodes(g)
	rownames(eM) <- nd[parent];
	return(eM);
}


# Conctruction of Weighted Adjacency Matrix
# Input:
# file : name of the text file to be read (def. edges). The format of the file is a sequence of rows. 
# Each row corresponds to an edge represented through a pair of vertices separated by blanks and the weight of the edge e.g.:
#  nodeX  nodeY   score
#  nodeZ  nodeK   score   
#  and so on... 
# compressed: boolean value: TRUE (def.) the file is in a .gz compressed format, otherwise is in a plain text format
# nodename: boolean value: TRUE (def.) if the names of nodes are gene symbol (i.e. characters)
#                          FALSE if the names of the nodes are entrez gene ID (i.e. integer numbers)
# 			Note: this boolean value is needed to get names of rows and columns of wadj matrix sorted in ascending order, both if they are 
#			characters or numbers. Indeed if names of rows (or columns) are number but they are stored as characters, the default "sort" function in R
#			sorts the characters-like-numbers by considering only the first charachter (a kind of alphabetic order). For instances:
#			x <- as.character(c(10,1,9)) --> sort(x) = "1" "10" "9". 
#			This can be easily solved: 	1.transforming numbers-like-characters as integers 
#										2.sorting integers 
#										3.transforming sorted integers numbers again as characters 
# Output:
# a named symmetric adjacency weighted matrix of the graph
weighted.adjacency.matrix <- function(file="edges", compressed=TRUE, nodename=TRUE) {
	if(compressed){
		m <- read.table(gzfile(file), colClasses="character", stringsAsFactors=FALSE);
	}else{
		m <- as.matrix(read.table(file, colClasses="character", stringsAsFactors=FALSE));
	}
	if(nodename){
		nodes <- sort(unique(as.vector(as.matrix(m[,1:2])))); ##NB:df must be converted as matrix to make as.vector workig..
	}else{
		nodes <- as.character(sort(as.numeric(unique(as.vector(m[,1:2]))))); 
	}
	n.nodes <- length(nodes);
	# building the adjacency matrix
	W <- matrix(0, nrow=n.nodes, ncol=n.nodes);
	dimnames(W) <- list(nodes,nodes);
	W[cbind(m[,1], m[,2])] <- as.numeric(m[,3]);
	W[cbind(m[,2], m[,1])] <- as.numeric(m[,3]);
	return(W);
}


# NOTA MN: *BUG* della funzione: definita una una coppia di vertici come intA e intB (i.e. interattore A e B) la funzione sottostante considera 
# come nodi univoci solo quelli di intA. Questa funzione ritorna un output corretto solo se il numero univoco di intA eguaglia quello di intB.
# Ad esempio nel database String (9.1) il numero univoco di intA = intB, mentre nel db Biogrid (3.2.106) il numero univoco di intA ≠ intB. 
# Di conseguenza si perderà informazione sul numero totale univoco dei nodi. 
# Conctruction of Weighted Adjacency Matrix
# Input:
# file : name of the text file to be read. The format of the file is a sequence of rows. 
# Each row corresponds to an edge represented through a pair of vertices separated by blanks and the weight of the edge e.g.:
#  nodeX  nodeY   score
#  nodeZ  nodeK   score   
#  and so on... 
# nodename: boolean value: if TRUE the names of nodes are characters (i.e. gene symbol), 
#                          if FALSE the names of the nodes are integer number (i.e. entrez gene ID)   
# Output:
# a named symmetric adjacency weighted matrix of the graph
weighted.adjacency.matrix.BUG <- function(file="edges.txt", nodename=TRUE) {
  m <- read.table(file, colClasses=c("character","character","numeric"));
  if(nodename){
	nodes <-sort(unique(m[,1])); # unique Nodes ##MN: BUG is here!! we consider only the unique node of intA!!
  }else{
	nodes <-sort(unique(as.numeric(m[,1]))); # unique nodes
	nodes <- as.character(nodes);
  }
  n.nodes <- length(nodes);
  # building the adjacency matrix
  W <- matrix(0, nrow=n.nodes, ncol=n.nodes);
  dimnames(W) <- list(nodes,nodes);
  W[cbind(m[,1], m[,2])] <- m[,3];
  W[cbind(m[,2], m[,1])] <- m[,3];
  return(W);
}


# Function to construct the matrix of the most specific annotations
# Input: 
# file: text file representing the most specific associations gene-HPO term (def: "gene2pheno.txt"). 
# 		file format: The file must be written as sequence of rows. Each row represent a gene and *all* its
#			   		 associations with abnormal phenotype separated by blanks, e.g.:
# 			   		 gene_1<tab>phen1<tab>...phen_N
#				  			.					.
#				 			.					.	
# 			   		 gene_M<tab>phen1<tab>...phen_N
# 		Note: this input file can be get parsing the file "ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt" with the perl subroutine "Do.HPO.edges.pl"
# geneame: boolean value:  TRUE (def.) if the names of genes are gene_symbol (i.e. characters)
#                          FALSE if the names of gene are entrez gene ID (i.e. integer numbers);
# 		   NOTE: see "weighted.adjacency.matrix" for a detailed explanation (see input parameter 'nodename' explanation)
# Output:
# the annotation matrix of the most specific annotations (0/1): rows are genes and columns are HPO terms. 
specific.annotation.matrix <- function(file="gene2pheno.txt", genename="TRUE"){
	line <- readLines(file);
	tmp <- strsplit(line, split="[ +\t]");

	gene.names <- c();
	for(i in 1:length(tmp)) gene.names <- c(gene.names,tmp[[i]][1]);

	ann.list <- list();
	for(i in 1:length(tmp)) ann.list[[i]] <- tmp[[i]][-1];
	names(ann.list) <- gene.names;

	hpoID <- unique(unlist(ann.list));
 
	n.genes <- length(gene.names);
	n.hpoID <- length(hpoID);
	m <- matrix(integer(n.genes * n.hpoID), nrow=n.genes);
	rownames(m) <- gene.names;
	colnames(m) <- hpoID;

	for (i in gene.names){
		spec.ann <- ann.list[[i]]; 
	    m[i, spec.ann] <- 1;  
	}

	if(genename){
		m <- m[sort(rownames(m)),sort(colnames(m))];
	}else{
		rname <- as.character(sort(as.numeric(gene.names)));
		m <- m[rname, sort(colnames(m))];
	}

	return(m);
}

## older and slight slower version..
specific.annotation.matrix.OLD <- function(file="gene2pheno.txt"){
	line <- readLines(file);
	tmp <- strsplit(line, split="[ +\t]");

	gene <- c();
	for(i in 1:length(tmp)) gene <- c(gene,tmp[[i]][1]);

	ann.list <- list();
	for(i in 1:length(tmp)) ann.list[[i]] <- tmp[[i]][-1];
	names(ann.list) <- gene;

	hpo <- c();
	for(i in 1:length(ann.list)) hpo <- c(hpo,ann.list[[i]]);
	hpo <- unique(hpo);

	ann <- matrix(0,length(gene),length(hpo));
	dimnames(ann) <- list(gene,hpo);

	for(gn in rownames(ann)) ann[gn,ann.list[[gn]]] <- 1;
	ann <- ann[,sort(colnames(ann))];

	return(ann);
}


# Function to construct a list of the most specific annotations form the table of the most specific annotations
# Input:
# ann : annotation matrix. Rows are examples and columns are most specific terms. It must be a named matrix
# Output:
# a list : names correspond to examples (genes). Elements are the most specific classes associated to that genes.
specific.annotation.list <- function(ann){
 ann.list <- apply(ann, 1, function(gene){
		terms <- which(gene==1);
		return(names(gene[terms])) 
	});
 	return(ann.list);
}

# Function to perform the transitive closure of the annotations using the most specific annotations and its ancestors. 
# The annotations are propagated from bottom to top, enriching the most specific annotations table. 
# The rows of the matrix correspond to the genes of the most specific annotation table and the columns to terms/classes
# Input:
# ann.spec: the annotation matrix of the most specific annotations (0/1): rows are genes and columns are HPO terms.
#			It can be obtained by calling specific.annotation.matrix function
# anc: list of the ancestors of the ontology. It can be obtained by calling build.ancestors function
# Output:
# an annotation table T: rows correspond to genes and columns to HPO terms. T[i,j]=1 means that gene i is annotated for the term j,
#             		  	 T[i,j]=0 means that gene i is not annotated for the term j.
transitive.closure.annotations <- function(ann.spec, anc){
	## costructiion of annotation list
	ann.list <- specific.annotation.list(ann.spec);

	## cotruction the full empty annotation matrix
	entrezIDs <- rownames(ann.spec);
	n.genes <- length(entrezIDs);
	HPOIDs <- names(anc);
	n.HPOID <- length(anc);
	hpo.ann <- matrix(numeric(n.HPOID * n.genes), nrow=n.genes, ncol=n.HPOID);	#empty label matrix
	dimnames(hpo.ann) <- list(entrezIDs,HPOIDs);
	
	## fill the full empty annotation matrix with the most specific annotation 	
	hpo.spec.term <- colnames(ann.spec); # the most specific hpo terms
	# might happen that there are same HPO IDs that are classified as "obsolete" in obo file, but that still exist in the annotation file 
	hpo.spec.term.sel <- HPOIDs[HPOIDs  %in% hpo.spec.term]; # removing obsolete HPO terms...
	hpo.ann[entrezIDs,hpo.spec.term.sel] <- ann.spec[,hpo.spec.term.sel];	

	## transitive closure: annotation propagation from the most specific nodes to all its ancestors
	for (i in entrezIDs){
		spec.ann <- ann.list[[i]];
		all.anc <- lapply(spec.ann, function(x) return(anc[[x]]));
		all.anc <- unique(unlist(all.anc));
		hpo.ann[i, all.anc] <- 1;  # setting the annotations derived by transitive closure
	}

	## remove HPO empty terms 
	hpo.ann <- hpo.ann[,colSums(hpo.ann)!=0];
	#x <- apply(hpo.ann,2,sum);
	#ind.hpo.id <- which(x==0);
	#hpo.ann <- hpo.ann[,-ind.hpo.id];

	return(hpo.ann);
}


# Function to construct the full annotations table using the most specific annotations and its ancestors w.r.t. a given weighted adjacency matrix. 
# The rows of the full annotations matrix correspond to all the examples of the given weighted adjacency matrix and the columns to the class/terms.
# The transitive closure of the annotations is performed.
# VERY IMPORTANT NOTE: examples that are present in the annotation matrix (ann.spec) but not in the symmetric adjacency weighted matrix (W) are purged
# Input:
# W: symmetric adjacency weighted matrix of the graph 
# anc: list of the ancestors of the ontology. It can be obtained by calling build.ancestors function
# ann.spec: the annotation matrix of the most specific annotations (0/1): rows are genes and columns are HPO terms.
#			It can be obtained by calling specific.annotation.matrix function
# Output:
# a full annotation table T: rows correspond to genes and columns to HPO terms. T[i,j]=1 means that gene i is annotated for the term j,
#             				 T[i,j]=0 means that gene i is not annotated for the term j.
full.annotation.matrix <- function(W, anc, ann.spec){
	## construction of annotation list
	ann.list <- specific.annotation.list(ann.spec);

	## construction the full empty annotation matrix
	entrezIDs <- rownames(W);
	n.genes <- length(entrezIDs);
	HPOIDs <- names(anc);
	n.HPOID <- length(anc);
	hpo.ann <- matrix(numeric(n.HPOID * n.genes), nrow=n.genes, ncol=n.HPOID);	#empty label matrix
	dimnames(hpo.ann) <- list(entrezIDs,HPOIDs);
	
	## fill the full empty annotation matrix with the most specific annotation 
	entrezIDs2hpo <- rownames(ann.spec);								# all genes that are associated with hpo terms
	entrezIDs.sel <- entrezIDs[entrezIDs %in% entrezIDs2hpo];			# genes 2 hpo terms 2 entrez id of wadj
	hpo.spec.term <- colnames(ann.spec);								# the most specific hpo terms
	#might happen that there are same HPO IDs that are classified as "obsolete" in obo file, but that still exist in the annotation file (e.g. build 1233)
	hpo.spec.term.sel <- HPOIDs[HPOIDs  %in% hpo.spec.term]; # removing obsolete HPO terms...
	hpo.ann[entrezIDs.sel,hpo.spec.term.sel] <- ann.spec[entrezIDs.sel,hpo.spec.term.sel];	# setting the most specific annotations

	## transitive closure: annotation propagation from the most specific nodes to all its ancestors
	for (i in entrezIDs){
		spec.ann <- ann.list[[i]];
		all.anc <- lapply(spec.ann, function(x) return(anc[[x]]));
		all.anc <- unique(unlist(all.anc));
		hpo.ann[i, all.anc] <- 1;  # setting the annotations derived by transitive closure
	}

	## remove HPO empty terms 
	hpo.ann <- hpo.ann[,colSums(hpo.ann)!=0];	# faster way than below
	# x <- apply(hpo.ann,2,sum);
	# ind.hpo.id <- which(x==0);
	# hpo.ann <- hpo.ann[,-ind.hpo.id];

	return(hpo.ann);
}


# Function to build an annotation matrix with only those terms having more than n annotations (n included).
# In other words terms having less than n annotations are pruned.
# Input:
# hpo.ann: the full annotations matrix (0/1). It can be obtained by calling 'full.annotation.matrix' function 
# n: number of annotations to be pruned (n is not included)
# Output:
# Matrix of annotations having only those terms with more than n annotations (n included)
do.submatrix <- function(hpo.ann,n){
	hpo.ann.sel <- hpo.ann[,colSums(hpo.ann)>n];
	return(hpo.ann.sel);
}


# Function to build a subgraph using only those classes present in an annotation matrix
# Input:
# nodes: a named vector with only the nodes for which we want to build the subgraph
# hpo.univ: a universal graph of class graphNELL. Represent the hierarchy of the whole class 
# Output:
# a subgraph with only the HPO terms included in the annotation matrix closed by transitive closure
do.subgraph <- function(nodes, hpo.univ){
	hpo <- subGraph(nodes,hpo.univ);
	return(hpo);
}

##MN: computa la stessa cosa di do.subgraph, ma in modo leggermente più lento...
## Test: dal grafo universo si costruisce il sotto-grafo avente come nodi soli quelli che hanno almeno un associazione ad 
##  	un fentipo patologico (termini della full annotation matrix).
##		versione file obo: data-version: releases/2016-04-01
##		versione file annotazioni: Build #110 (Jan 25, 2016 7:00:22 PM) 
## tempi di calcolo: 
## system.time(do.subgraph(hpo.ann,hpo.univ)): 		0.574 sec
## system.time(do.my.subgraph(hpo.ann,hpo.univ)):	1.154 sec
do.my.subgraph <- function(nd, g){
	ed <- edges(g);
	ed.sel <- ed[nd];

	ndL <- vector(mode="list", length=length(ed.sel));
	names(ndL) <- names(ed.sel);

	for(i in 1:length(ed.sel)){ 
		parent   <- names(ed.sel[i]);
		children <- ed.sel[[i]];
		if(length(children!=0)){ 
			children.map <- children[children %in% nd]
			ndL[[i]] <- append(ndL[[i]],children.map);
		} 
	}

	for (i in 1:length(ndL)){
		ndL[[i]] <- list(edges=ndL[[i]]);
	}

	G <- graphNEL(nodes=nd, edgeL=ndL, edgemode="directed");
	return(G);
}


# Function to assess the integrity of a full annotation table
# Input:
# anc: list of the ancestors of the ontology. It can be obtained by calling build.ancestors function
# ann.spec: the annotation matrix of the most specific annotations (0/1): rows are genes and columns are HPO terms. 
#			It can be obtained by calling 'specific.annotation.matrix' function
# hpo.ann: the full annotations matrix (0/1). It can be obtained by calling 'full.annotation.matrix' function
# Output:
# If the transitive closure of the annotations is performed well (i.e. if the full annotation matrix is right) "OK" is printed,
# otherwise a message error is printed
check.annotation.matrix.integrity <- function(anc, ann.spec, hpo.ann){
	## construction of annotation list
	ann.list <- specific.annotation.list(ann.spec);

	genes <- rownames(hpo.ann);

	check <- c();
	for (i in genes){
		spec.ann <- which(hpo.ann[i,]==1);
		len.ann <- length(spec.ann);
		all.anc <- lapply(ann.list[[i]], function(x) return(anc[[x]]));
		all.anc <- unique(unlist(all.anc));
		len.anc <- length(all.anc);
		cmp <- len.anc == len.ann;
		if(cmp==TRUE){
			check <- c(check,"OK")
		} else {
			check <- c(check,"NOTOK")
		}
	}
	names(check) <- genes;

	violated <- any(check!="OK");
	if(violated){
		n <- names(check)[check=="NOTOK"];
		cat("check.annotation.matrix: NOT_OK. Transitive closure NOT RESPECTED", "\n");
	}else{
		cat("check.annotation.matrix: OK", "\n");	
	}
}


# Function to assess the integrity of a DAG
# If there are nodes not accessible from the root node a message error is displayed
# Input:
# g : a graph of class graphNEL. It represents the hierarchy of the classes.
# root: name of the class that is on the top-level of the hierarchy (def:"00")
# Output:
# If there are no nodes not accessible from the root, then "OK" is printed,
# otherwise a message error and the list of the not accessible nodes is printed.
check.DAG.integrity <- function(g, root="00") {
  all.nodes <- nodes(g);
  acc.nodes <- names(acc(g,root)[[1]]);
  if((length(all.nodes) - length(acc.nodes)) > 1) {
		n <- setdiff(all.nodes,c(acc.nodes,root));
		cat("check.GO.integrity: not all nodes accessible from root", "\n");
		cat("Nodes not accessible from root: \n");
		cat(n,"\n");
	}else{ 
		cat("OK \n")
	};
}


# Function to check the hierarchical constraints of the prediction. It checks whether the predictions verify the hierarchy constraints, that is
# the score of a child node must not be greater than the score of its parents. 
# The algorithm can be applied to any DAG-structured hierarchy of classes
# Input:
# S.hier: the matrix with the scores of the classes corrected in according to hierarchy.
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# root: name of the class that is on the top-level of the hierarchy (def:"00")
# Output:
# a list of 3 elements:
# - Status: a single value: OK if none hierarchical constraints have bee broken, NOTOK if there is at least one hierarchical constraints broken;
# - Hierarchy_Constraints_Broken: boolean values vector: for each terms/classes are possible two value: TRUE if there is some example did not  
#  	respect the hierarchical constraints; FALSE if none sample broke the hierarchical constraints.
# - Hierarchy_costraints_satisfied: how many terms satisfied the hierarchical constraint.
check.hierarchy <- function(S.hier,g, root="00"){
	# si aggiunge una "dummy root column" se non presente nella matrice S.hier
	if(!(root %in% colnames(S.hier))){
	  	max.score <- max(S.hier);
	  	z <- rep(max.score,nrow(S.hier));
	  	S.hier <- cbind(z,S.hier);
	  	colnames(S.hier)[1] <- root;
	}
	par <- get.parents(g,root);
	v <- c()
	for(i in 1:length(par)){
		child <- S.hier[,names(par[i])];
		parents <- S.hier[,par[[i]]]
		x <- parents >= child   # in alternativa: child <= parents
		y <- any(x==0)   # check hierarchy constraints
		v <- append(v,y)
	}
	names(v) <- names(par)
	violated <- any(v==TRUE);
	if(violated)
	  	Status = "NOTOK"
	else
	  	Status = "OK";
	h <- as.factor(v);
	k <- summary(h);
	l <- list(Status=Status, hierarchy.constraints.broken=v, hierarchy.costraints.satisfied=k);
	return(l);
}

#MN: The same function of above, but instead to get in input the whole hierarchical score matrix, receives in input 
# the vector of scores relative to a single example
# Input:
# y.hier: vector of scores relative to a single example. This must be a named numeric vector. 
#		  The names must correspond to the names of the classes of the hierarchy
# g: a graph of class graphNEL. It represents the hierarchy of the classes.
# root: name of the class that is on the top-level of the hierarchy (def:"00")
# Output:
# a list of 3 elements:
# - Status: a single value: OK if none hierarchical constraints have bee broken, NOTOK if there is at least one hierarchical constraints broken;
# - Hierarchy_Constraints_Broken: boolean values vector: for each terms/classes are possible two value: TRUE if there is some example did not  
#  	respect the hierarchical constraints; FALSE if none sample broke the hierarchical constraints.
# - Hierarchy_costraints_satisfied: how many terms satisfied the hierarchical constraint.
check.hierarchy.single.sample <- function(y.hier,g, root="00"){
	# si aggiunge il root node se non presente nel vettore y.hier
	if(!(root %in% names(y.hier))){
		max.score <- max(y.hier);
		y.hier <- c(max.score,y.hier);
		names(y.hier)[1] <- root;
	}
	par <- get.parents(g,root);
	v <- c()
	for(i in 1:length(par)){
	  	child <- y.hier[names(par[i])];
	  	parents <- y.hier[par[[i]]]
	  	x <- parents >= child   # in alternativa: child <= parents
	  	y <- any(x==0)   # check hierarchy constraints
	  	v <- append(v,y)
	}
	names(v) <- names(par)
	violated <- any(v==TRUE);
	if(violated)
	  	Status = "NOTOK"
	else
	  	Status = "OK";
	h <- as.factor(v);
	k <- summary(h);
	l <- list(Status=Status, hierarchy.constraints.broken=v, hierarchy.costraints.satisfied=k);
	return(l);
}


# Function to split a dataset randomly in k fold through an unstratified way. 
# This function is useful to perform experiment with k-fold cross-validation in the contest of hierarchical correction, where split data 
# in stratified fold is not needed.
# This function is called inside the high-level function Do.'tpr-variant'.cv.
# Input:
# S: matrix of the flat scores. It must be a named matrix, where rows are example (e.g. genes) and columns are classes/terms.
# kk: number of folds in which to split the dataset (def: 5)
# seed : seed for the random generator. If NULL (default) no initialization is performed.
# Output:
# a list with k=kk components (folds). Each component is a vector with name of the examples. Following this strategy is highly likely that a fold 
# won't contain equal amount of positive and negative examples (especially for unbalanced dataset).
do.unstratified.cv.data <- function(S, kk=5, seed=NULL){
	set.seed(seed);
	examples <- rownames(S);
	n <- nrow(S);
	size <- c();
	folds <- vector(mode="list", length=kk)
	names(folds) <- paste0(rep("fold",kk), 1:kk)
	for (k in 1:kk) {
		first <- ((k - 1) * n) %/% kk
		last <- (k * n) %/% kk
		size <- last-first;
		x	<- sample(examples,size);
		folds[[k]] <- x;
		examples <- setdiff(examples,x);
	}
	return(folds);
}

# Function to generate data for the stratified cross-validation w.r.t a single class/term.
# Input:
# examples : indices or names of the examples. Can be either a vector of integers or a vector of names. 
# positives: vector of integers or vector of names. The indices (or names) refer to the indices (or names) of 'positive' examples	
# kk : number of folds (def=5)
# seed : seed of the random generator (def=NULL). If is set to NULL no initialization is performed
# Output:
# a list with 2 two components
#   - fold.non.positives : a list with k components. Each component is a vector with the indices (or names) of the non positive elements of the fold
#   - fold.positives : a list with k components. Each component is a vector with the indices (or names) of the positive elements of the fold
# N.B.: in both elements indices (or names) refer to row numbers (or names) of the data matrix	 
do.stratified.cv.data.single.class <- function(examples, positives, kk=5, seed=NULL){
	set.seed(seed);

	negatives <- setdiff(examples,positives); 	

	n <- length(positives);		
	m <- length(negatives);		

	set.pos <- list();
	set.neg <- list();

	for (k in 1:kk) {
		#fold positives 
		last.pos <- (k * n) %/% kk;
		first.pos  <- ((k - 1) * n) %/% kk;
		size.pos <-  last.pos - first.pos;				
		subset.pos <- 	sample(positives, size.pos);			
		set.pos[[k]] <- subset.pos;						
		positives <- setdiff(positives, subset.pos);	
		
		#fold non positives
		last.neg <- (k * m) %/% kk;
		first.neg  <- ((k - 1) * m) %/% kk;
		size.neg <-  last.neg - first.neg;				
		subset.neg <- sample(negatives, size.neg);					
		set.neg[[k]] <- subset.neg;						
		negatives <- setdiff(negatives, subset.neg);	
	}
	return(list(fold.positives=set.pos, fold.negatives=set.neg));
}

# Function to generate stratified fold w.r.t a single class/term.
# Input:
# examples : indices or names of the examples. Can be either a vector of integers or a vector of names. 
# positives: vector of integers or vector of names. The indices (or names) refer to the indices (or names) of 'positive' examples	
# kk : number of folds (def=5)
# seed : seed of the random generator (def=NULL). If is set to NULL no initialization is performed
# Output:
# a list with k components: each element is a fold that contains equal amount of positives and negatives examples 
do.stratified.cv.fold <- function(examples, positives, kk=5, seed=NULL){
	set.seed(seed);
	folds <- vector(mode="list", length=kk);
	names(folds) <- paste0(rep("fold",kk), 1:kk);
	strfold <- do.stratified.cv.data.single.class(examples,positives, kk=kk, seed=seed)
	for(k in 1:kk){
		folds[[k]] <- append(strfold$fold.positives[[k]], strfold$fold.negatives[[k]]);
	}
	return(folds);
}

# Function to split a dataset randomly in k fold through an unstratified way. It is useful to perform experiment with k-fold cross-validation 
# in the contest of hierarchical correction. This function is called inside the high-level function Do.tpr.threshold.cv.
# Input:
# S.flat: matrix of the flat scores. It must be a named matrix, where rows are example (e.g. genes) and columns are classes/terms.
# kk: number of folds in which to split the dataset (def: 5)
# seed : seed for the random generator. If NULL (default) no initialization is performed.
# Output:
# a list with k components where each element is a submatrix: rows are a random subset of genes (i.e. row names of the S.flat matrix) 
# and columns all terms.
do.unstratified.cv.data.OLD <- function(S.flat, kk=5, seed=NULL){
	if (!is.null(seed)) set.seed(seed);
	samples <- rownames(S.flat);
	n.samples <- nrow(S.flat);
	size <- c();
	folds <- vector(mode="list", length=kk)
	names(folds) <- paste0(rep("fold",kk), 1:kk)
	for (k in 1:kk) {
		first <- ((k - 1) * n.samples) %/% kk
		last <- (k * n.samples) %/% kk
		size <- c(size, last-first)
		x	<- sample(samples,size[k]);
		folds[[k]] <- S.flat[x,];
		samples <- setdiff(samples,x);
	}
	return(folds);
}

