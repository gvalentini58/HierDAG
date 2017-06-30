# Utility functions to read and write graphs to files

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

# September 2013
# December 2013: added a function to parse the HPO obo file

## library to call
library("graph");
library("RBGL");
library("Rgraphviz");

###########################################
# Utility functions
###########################################

# Function to extract the edges of the DAG from a HPO obo file
# Input:
# file: a HPO OBO file
# file : name of the file of the edges to be written
# Output:
# a text file representing the edges in the format:
# source  destination
# (i.e. one row for each edge)

do.edges.from.HPO.obo <- function(file="hp.obo", output.file="edge.file") {  
  line <- readLines(file);
  n.lines <- length(line);
  m <- matrix(character(1000000*2), ncol=2);
  colnames(m)=c("source", "destination");
  i=1;
  j=0; # number of edges;
  #browser();
  while(i<=n.lines) {
    while((i<=n.lines) && (line[i]!="[Term]")) {
	  i <- i + 1;
	}
	if (i>=n.lines)
	  break();
	i <- i + 1; # id
	destination <- strsplit(line[i], split="[ +\t]")[[1]][2];
	while( (line[i]!="") && (strsplit(line[i], split="[ +\t]")[[1]][1]!="is_a:") )  # checking first is_a entry
	  i <- i + 1;
	if (line[i] == "") next();  # we are at the end of the record and is_a has been found
	source <- strsplit(line[i], split="[ +\t]")[[1]][2];
	j <- j + 1;
	i <- i + 1;
	m[j,]<-c(source,destination);
	# while(line[i]=="") i <- i+1;
	while( (line[i]!="") && (strsplit(line[i], split="[ +\t]")[[1]][1]=="is_a:") )  {# checking successive is_a entry
      source <- strsplit(line[i], split="[ +\t]")[[1]][2];
      i <- i + 1;
	  j <- j + 1;
	  m[j,]<-c(source,destination);
	} 	  
  }
  m <- m[1:j,];
  write.table(m, file=output.file, quote=FALSE, row.names=FALSE, col.names=FALSE);
}



# Writing a graph to file.
# The graph is written to a file as a sequence of rows. Each row corresponds to an edge
# represented through a pair of vertices separated by blanks, e.g.:
#  GO:0006732  GO:0006733 
#  GO:0006769  GO:0009108
# Input:
# g: a graph of class graphNEL
# file : name of the file to be written
# Output:
# a file representing the graph
write.graph <- function(g, file="graph.txt") {  
  num.edges <- length(unlist(edges(g)));
  num.v <- numNodes(g);
  m <- matrix(character(num.edges*2), ncol=2);
  res <- edges(g);
  count=0;
  node1 <- names(res);
  for (i in 1:num.v) {
    x <- res[[i]];
	len.x <- length(x);
	if (len.x!=0)
	  for (j in 1:len.x) {
	    count <- count + 1;
		m[count,] <- c(node1[i],x[j]);
	  }
  }
  write.table(m, file=file, quote=FALSE, row.names=FALSE, col.names=FALSE);
  return();	    
}



# Reading a graph from file.
# The graph is read from a file and a graphNEL object is built
# The format of the file is a sequence of rows. Each row corresponds to an edge
# represented through a pair of vertices separated by blanks, e.g.:
#  GO:0006732  GO:0006733 
#  GO:0006769  GO:0009108
# Input:
# file : name of the file to be read
# Output:
# a graph of class graphNEL
read.graph <- function(file="graph.txt"){  
  m <- as.matrix(read.table(file, colClasses="character"));
  thenodes<-sort(unique(as.vector(m))); # nodes
  n.nodes <- length(thenodes);
  n.edges <- nrow(m);
  # building the graph
  edL <- vector("list", length=n.nodes);
  names(edL) <- thenodes;
  for(i in 1:n.nodes)
    edL[[i]]<-list(edges=NULL);
  g <- graphNEL(nodes=thenodes, edgeL=edL, edgemode="directed");
  g <- addEdge(m[1:n.edges,1], m[1:n.edges,2], g, rep(1,n.edges));
  return(g);
}


read.graph.old <- function(file="graph.txt") {  
  m <- as.matrix(read.table(file, colClasses="character"));
  thenodes<-sort(unique(as.vector(m))); # nodes
  n.nodes <- length(thenodes);
  n.edges <- nrow(m);
  # building the graph
  edL <- vector("list", length=n.nodes);
  names(edL) <- thenodes;
  for(i in 1:n.nodes)
    edL[[i]]<-list(edges=NULL);
  g <- graphNEL(nodes=thenodes, edgeL=edL, edgemode="directed");

  for (i in 1:n.edges) {
    g <- addEdge(m[i,1], m[i,2], g, 1);
  }
  return(g);
}


# Reading an undirected graph with weights from a file.
# The graph is read from a file and a graphNEL object is built
# The format of the file is a sequence of rows. Each row corresponds to an edge
# represented through a pair of vertices separated by blanks, and the weight of the edge e.g.:
#  GO:0006732  GO:0006733   0.234
#  GO:0006769  GO:0009108   0.106
# Input:
# file : name of the file to be read
# Output:
# a graph of class graphNEL
read.undirected.graph <- function(file="graph.txt") {  
  m <- as.matrix(read.table(file, colClasses="character"));
  thenodes<-sort(unique(as.vector(m[,1:2]))); # nodes
  n.nodes <- length(thenodes);
  n.edges <- nrow(m);
  # building the graph
  edL <- vector("list", length=n.nodes);
  names(edL) <- thenodes;
  for(i in 1:n.nodes)
    edL[[i]]<-list(edges=NULL);
  g <- graphNEL(nodes=thenodes, edgeL=edL, edgemode="undirected");
  g <- addEdge(m[1:n.edges,1], m[1:n.edges,2], g, as.numeric(m[1:n.edges,3]));
  return(g);
}

