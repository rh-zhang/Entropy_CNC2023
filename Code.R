################################################################
########### Simulating Weak Attacks in a New Duplication #######
#################### Divergence Model with Node Loss ###########
###################### Data pre-processing ###################
################################################################


# Created: 04/02/24
# Last edited: 04/02/24
# Edit history:
### 04/02/24: Copied from original scripts. 

# Description: Script for data pre-processing (for the two STRING networks)
# Where files are written to or accessed, the exact path name has been removed.

#install libraries
library(netcom)
library(netdiffuseR)
library(igraph)
library(foreach)
library(data.table)

yeast <- read.table("4932.protein.physical.links.v12.0.txt.gz", header=TRUE)
colnames(yeast)[3] <- "weight"
ppi_yeast <- graph_from_data_frame(yeast, directed = TRUE, vertices = NULL)
ppi_yeast <-simplify(ppi_yeast, remove.multiple=TRUE, remove.loops=TRUE)
ppi_yeast <- delete.edges(ppi_yeast, E(ppi_yeast)[E(ppi_yeast)$weight < 400])
sum(degree(ppi_yeast)==0) ##number of isolated nodes: 1100
ppi_yeast #5925nodes, 140402 edges
#assign equal initial weight
edge.attributes(ppi_yeast)$weight <- rep(1,gsize(ppi_yeast))
#plot
plot(ppi_yeast,layout=layout_nicely)

ecoli <- read.table("511145.protein.physical.links.v12.0.txt.gz", header=TRUE)
colnames(ecoli)[3] <- "weight"
ppi_ecoli <- graph_from_data_frame(ecoli, directed = TRUE, vertices = NULL)
ppi_ecoli <-simplify(ppi_ecoli, remove.multiple=TRUE, remove.loops=TRUE)
ppi_ecoli <- delete.edges(ppi_ecoli, E(ppi_ecoli)[E(ppi_ecoli)$weight < 400])
sum(degree(ppi_ecoli)==0) ##number of isolated nodes: 833
ppi_ecoli #3043nodes, 59182 edges
#assign equal initial weight
edge.attributes(ppi_ecoli)$weight <- rep(1,gsize(ppi_ecoli))
#plot
plot(ppi_ecoli,layout=layout_nicely)


###################### Simulation of weak attacks ##############
################################################################

# Created: 04/02/24
# Last edited: 04/02/24
# Edit history:
### 04/02/24: Copied from original scripts. 

# Description: Script for weak attacks simulations on two STRING networks,
# DD model and the new DD model with node loss

## network efficiency
node_eff <- function(graph,nodeIndex){
  return(sum(1/(distances(graph, V(graph)[nodeIndex])[-nodeIndex])))
}

NE <- function(graph){
  eff <- foreach(node = seq_along(V(graph)), .combine = c, 
                 .packages = "igraph") %do% node_eff(graph, node)
  return(list(sum(eff),eff))
}

## Attack strategies
# Complete attack
attackA <- function(graph, node) {
  graph <- delete_edges(graph,c(E(graph)[from(V(graph)[node])],E(graph)[to(V(graph)[node])]))
  graph
}

# Partial attack
attackB1 <- function(graph, node){
  edges <- E(graph)[c(E(graph)[from(V(graph)[node])],E(graph)[to(V(graph)[node])])]
  deleteEdges <- sample(edges,size = floor(length(edges)/2))
  graph <- delete_edges(graph,deleteEdges)
  graph
}

# Partial attenuation
attackB2 <- function(graph, node){
  edges <- E(graph)[c(E(graph)[from(V(graph)[node])],E(graph)[to(V(graph)[node])])]
  graph <- set.edge.attribute(graph, "weight", index=edges, edges$weight*2)  
  graph
}

# Distributed attack
attackC1 <- function(graph){
  deleteEdges <- rbinom(0.2*length(E(graph)),length(E(graph)),0.2)
  graph <- delete_edges(graph,E(graph)[deleteEdges])
  graph
}

# Distributed attenuation
attackC2 <- function(graph){
  attenuatedEdges <- rbinom(0.2*length(E(graph)),length(E(graph)),0.2)
  graph <- set.edge.attribute(graph, "weight", index=E(graph)[attenuatedEdges], E(graph)[attenuatedEdges]$weight*2)  
  graph 
}

# Node removal by successive maximal damage strategy
greedy_node <- function(graph, attackStrategy, Niter, m){
  nodeList <- V(graph)
  efficiency <- vector()
  efficiency[1] <- 1
  edgeNumber <- vector()
  nodeNumber <- vector()
  nisolated <- vector()
  local_centrality <- list()
  
  n <- 1
  NetEff <- NE(graph)[[1]]
  edgeNumber[1] <- gsize(graph)
  nodeNumber[1] <- gorder(graph)
  nisolated[1] <- sum(degree(graph)==0)
  newgraph <- graph
  
  while (n < Niter){
    i <- 1
    tmp_graph <- newgraph
    removed <- vector()
    
    while (i <= m ){
      new_NE <- foreach(node = seq_along(V(tmp_graph)), .combine = c, 
                        .packages = "igraph") %do% NE(delete_vertices(tmp_graph, node))[[1]]
      damage <- NetEff-new_NE
      removed[i] <-V(tmp_graph)[which.max(damage)][1]
      tmp_graph <- delete_vertices(tmp_graph,removed[i])
      i <- i+1
    }
    
    newgraph <- attackStrategy(newgraph, removed)
    neweff <- NE(newgraph)
    efficiency[n+1] <-  neweff[[1]]/NetEff
    edgeNumber[n+1] <- gsize(newgraph)
    nodeNumber[n+1] <- gorder(newgraph)
    local_centrality[[n]] <- neweff[[2]]
    nisolated[n+1] <- sum(degree(newgraph)==0)
    n <- n+1 
    print(n)
  }
  return(list(efficiency,edgeNumber,nodeNumber, nisolated, local_centrality))
}

# Edge removal by successive maximal damage strategy
greedy_edge <- function(graph, attackStrategy, Niter){
  nodeList <- V(graph)
  efficiency <- vector()
  efficiency[1] <- 1
  edgeNumber <- vector()
  nodeNumber <- vector()
  nisolated <- vector()
  local_centrality <- list()
  
  n <- 1
  NetEff <- NE(graph)[[1]]
  edgeNumber[1] <- gsize(graph)
  nodeNumber[1] <- gorder(graph)
  nisolated[1] <- sum(degree(graph)==0)
  
  while (n < Niter){
    newgraph <- attackStrategy(graph)
    neweff <- NE(newgraph)
    efficiency[n+1] <-  neweff[[1]]/NetEff
    edgeNumber[n+1] <- gsize(newgraph)
    nodeNumber[n+1] <- gorder(newgraph)
    local_centrality[[n]] <- neweff[[2]]
    nisolated[n+1] <- sum(degree(newgraph)==0)
    n <- n+1 
    print(n)
    graph <- newgraph
  }
  return(list(efficiency,edgeNumber,nodeNumber, nisolated, local_centrality))
}

#Simulation of weak attacks on real yeast networks, damage measured by network efficiency
simulationA_yeast <- list()
simulationB1_1node <- list()
simulationB1_2nodes <- list()
simulationB1_5nodes <- list()
simulationC1 <- list()
simulationC1 <- list()

for (i in 1:10){
  simulationA_yeast[[i]] <-  greedy_node(ppi_yeast,attackA,35,1)
  simulationB1_1node[[i]] <- greedy_node(ppi_yeast,attackB1,35, 1)
  simulationB1_2nodes[[i]] <- greedy_node(ppi_yeast,attackB1,35, 2)
  simulationB1_5nodes[[i]] <- greedy_node(ppi_yeast,attackB1,35, 5)
  simulationC1[[i]] <-  greedy_edge(ppi_yeast,attackC1,35)
  simulationC2[[i]] <-  greedy_edge(ppi_yeast,attackC2,35)
}

# Note: the simulations could take long, for a quicker option, can try the following 
# with a smaller number of nodes.
#for (i in 1:10){
#  simulationA_yeast[[i]] <-  greedy_node(ppi_yeast,attackA,10,1)
#  simulationB1_1node[[i]] <- greedy_node(ppi_yeast,attackB1,10, 1)
#  simulationB1_2nodes[[i]] <- greedy_node(ppi_yeast,attackB1,10, 2)
#  simulationB1_5nodes[[i]] <- greedy_node(ppi_yeast,attackB1,10, 5)
#  simulationC1[[i]] <-  greedy_edge(ppi_yeast,attackC1,10)
#  simulationC2[[i]] <-  greedy_edge(ppi_yeast,attackC2,10)
#}

plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), simulationA_yeast)))[1]),
     type="b",lty=5,col="blue",pch=22,
     xlab="Number of attacks", ylab="Network efficiency")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), simulationB1_1node)))[1]),
      type="b",lty=5,col="red",pch=17)
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), simulationB1_2nodes)))[1]),
      type="b",lty=5,col="green",pch=17)
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), simulationB1_5nodes)))[1]),
      type="b",col="purple",pch=12)
legend("bottomleft",legend=c("Knockout","1 node halved","2 nodes halved","5 nodes halved"),
       col=c("blue","red","green","purple"),
       pch=c(22,17,12,5),lty=5,cex=0.8)


par(xpd=FALSE)
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), simulationC1)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), simulationA_yeast)))[2]))
abline(v=10)

par(xpd=FALSE)
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), simulationC2)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), simulationA_yeast)))[2]))
abline(v=10)


## Node Loss Model
NodeLossModel1 <- function(size, divergence, directed = FALSE, q1) {
  ## Start with pair of connected nodes
  matrix <- matrix(0, nrow = 3, ncol = 3)
  matrix[1,2] = 1
  matrix[2,1] = 1
  i <- 1
  node <- 3
  ## Start with node three because first two nodes are part of the initial network
  while (i <= size) {
    duplication <- sample(1:(node-1), 1)
    matrix[node, ] = matrix[duplication, ]
    matrix[, node] = matrix[node, ]
    
    ## Change each edge with probability divergence
    edges <- which(matrix[node, ] != 0)
    for (e in edges) {
      if (stats::runif(1) <= divergence) {
        matrix[node, e] = 0
        
        if (directed == FALSE) {
          matrix[e, node] = 0
        }
      }
    }
    
    ## Delete isolated nodes with probability q1
    g1 <- graph_from_adjacency_matrix(matrix)
    isolated <- which(degree(g1)==0)
    if (length(isolated) != 0){
      g1 <- delete.vertices(g1, sample(isolated,round(length(isolated)*q1)))
    }
    matrix <- cbind(as_adjacency_matrix(g1),0)
    matrix <- rbind(matrix,0)
    node <- nrow(matrix)
    i <- i+1
  }
  
  return(g1)}

#q1=0
result0 <- NodeLossModel1(100, 0.2, directed = FALSE, 0)
E(result0)$weight <- 1
model0afterA <- list()
model0afterB12 <- list()
model0afterB15 <- list()
model0afterB110 <- list()
model0afterB22 <- list()
model0afterB25 <- list()
model0afterB210 <- list()
model0afterC1 <- list()
model0afterC2 <- list()

for (i in 1:10){
  model0afterA[[i]] <- greedy_node(result0, attackA,25,1) 
  model0afterB12[[i]] <- greedy_node(result0,attackB1,25,2)
  model0afterB15[[i]] <- greedy_node(result0,attackB1,25,5)
  model0afterB110[[i]] <- greedy_node(result0,attackB1,25,10)
  model0afterB22[[i]] <- greedy_node(result0,attackB2,25,2)
  model0afterB25[[i]] <- greedy_node(result0,attackB2,25,5)
  model0afterB210[[i]] <- greedy_node(result0,attackB2,25,10)
  model0afterC1[[i]] <- greedy_edge(result0,attackC1,25,1)
  model0afterC2[[i]] <- greedy_edge(result0,attackC2,25,1)
}

plot.new()
par(mfrow=c(2,2)) 
#Knockout
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterA)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=2,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Duplication divergence model with p=0.2, r=0, Knockout")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterB12)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=2,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterB15)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=2,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterB110)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("Knockout","2 nodes halved","5 nodes halved","10 nodes halved"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17,15),lty=2,cex=0.5)

#Attenuation
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterA)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Duplication divergence model with p=0.2, r=0, Attenuation")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterB22)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=5,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterB25)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=5,col="green")
lines(s.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterB210)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("Knockout","2 nodes halved","5 nodes halved","10 nodes halved"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17),lty=2,cex=0.5)

#Distributed knockout
par(xpd=FALSE)
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterC1)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Duplication divergence model with p=0.2, r=0, Distributed Knockout")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterA)))[2]))
abline(v=10)

#Distributed attenuation
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterC2)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Duplication divergence model with p=0.2, r=0, Distributed Attenuation")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model0afterA)))[2]))
abline(v=10)

## r=0
#q1=0.2
result1 <- NodeLossModel1(100, 0.2, directed = FALSE, 0.2)
E(result1)$weight <- 1
model1afterA <- list()
model1afterB12 <- list()
model1afterB15 <- list()
model1afterB110 <- list()
model1afterB22 <- list()
model1afterB25 <- list()
model1afterB210 <- list()
model1afterC1 <- list()
model1afterC2 <- list()

for (i in 1:10){
  model1afterA[[i]] <- greedy_node(result1, attackA,25,1) 
  model1afterB12[[i]] <- greedy_node(result1,attackB1,25,2)
  model1afterB15[[i]] <- greedy_node(result1,attackB1,25,5)
  model1afterB110[[i]] <- greedy_node(result1,attackB1,25,10)
  model1afterB22[[i]] <- greedy_node(result1,attackB2,25,2)
  model1afterB25[[i]] <- greedy_node(result1,attackB2,25,5)
  model1afterB210[[i]] <- greedy_node(result1,attackB2,25,10)
  model1afterC1[[i]] <- greedy_edge(result1,attackC1,25,1)
  model1afterC2[[i]] <- greedy_edge(result1,attackC2,25,1)
}


plot.new()
par(mfrow=c(2,2)) 
#Knockout
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterA)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=2,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.2, p=0.2, r=0, Knockout")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterB12)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=2,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterB15)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=2,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterB110)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("Knockout","2 nodes halved","5 nodes halved","10 nodes halved"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17,15),lty=2,cex=0.5)

#Attenuation
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterA)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.2, p=0.2, r=0, Attenuation")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterB22)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=5,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterB25)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=5,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterB210)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("Knockout","2 nodes halved","5 nodes halved","10 nodes halved"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17),lty=2,cex=0.5)

#Distributed knockout
par(xpd=FALSE)
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterC1)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.2, p=0.2, r=0, Distributed knockout")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterA)))[2]))
abline(v=10)

#Distributed attenuation
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterC2)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.2, p=0.2, r=0, Distributed attenuation")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterA)))[2]))
abline(v=10)


#q=0.4
result2 <- NodeLossModel1(100, 0.2, directed = FALSE, 0.4)
E(result2)$weight <- 1
model12afterA <- list()
model12afterB12 <- list()
model12afterB15 <- list()
model12afterB110 <- list()
model12afterB22 <- list()
model12afterB25 <- list()
model12afterB210 <- list()
model12afterC1 <- list()
model12afterC2 <- list()

for (i in 1:10){
  model12afterA[[i]] <- greedy_node(result2, attackA,25,1) 
  model12afterB12[[i]] <- greedy_node(result2,attackB1,25,2)
  model12afterB15[[i]] <- greedy_node(result2,attackB1,25,5)
  model12afterB110[[i]] <- greedy_node(result2,attackB1,25,10)
  model12afterB22[[i]] <- greedy_node(result2,attackB2,25,2)
  model12afterB25[[i]] <- greedy_node(result2,attackB2,25,5)
  model12afterB210[[i]] <- greedy_node(result2,attackB2,25,10)
  model12afterC1[[i]] <- greedy_edge(result2,attackC1,25,1)
  model12afterC2[[i]] <- greedy_edge(result2,attackC2,25,1)
}

plot.new()
par(mfrow=c(2,2)) 
#Knockout
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterA)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=2,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.4, p=0.2, r=0, Knockout")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterB12)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=2,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterB15)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=2,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterB110)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("Knockout","2 nodes halved","5 nodes halved","10 nodes halved"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17,15),lty=2,cex=0.5)

#Attenuation
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterA)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.4, p=0.2, r=0, Attenuation")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterB22)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=5,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterB25)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=5,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterB210)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("Knockout","2 nodes halved","5 nodes halved","10 nodes halved"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17),lty=2,cex=0.5)

#Distributed knockout
par(xpd=FALSE)
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterC1)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.4, p=0.2, r=0, Distributed knockout")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterA)))[2]))
abline(v=11)

#Distributed attenuation
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterC2)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.4, p=0.2, r=0, Distributed attenuation")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterA)))[2]))
abline(v=11)


#q=0.6
result3 <- NodeLossModel1(100, 0.2, directed = FALSE, 0.6)
E(result3)$weight <- 1
model13afterA <- list()
model13afterB12 <- list()
model13afterB15 <- list()
model13afterB110 <- list()
model13afterB22 <- list()
model13afterB25 <- list()
model13afterB210 <- list()
model13afterC1 <- list()
model13afterC2 <- list()

for (i in 1:10){
  model13afterA[[i]] <- greedy_node(result3, attackA,25,1) 
  model13afterB12[[i]] <- greedy_node(result3,attackB1,25,2)
  model13afterB15[[i]] <- greedy_node(result3,attackB1,25,5)
  model13afterB110[[i]] <- greedy_node(result3,attackB1,25,10)
  model13afterB120[[i]] <- greedy_node(result3,attackB1,25,20)
  model13afterB22[[i]] <- greedy_node(result3,attackB2,25,2)
  model13afterB25[[i]] <- greedy_node(result3,attackB2,25,5)
  model13afterB210[[i]] <- greedy_node(result3,attackB2,25,10)
  model13afterC1[[i]] <- greedy_edge(result3,attackC1,25,1)
  model13afterC2[[i]] <- greedy_edge(result3,attackC2,25,1)
}


plot.new()
par(mfrow=c(2,2)) 
#Knockout
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterA)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=2,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.6, p=0.2, r=0, Knockout")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterB12)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=2,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterB15)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=2,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterB110)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("Knockout","2 nodes halved","5 nodes halved","10 nodes halved"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17,15),lty=2,cex=0.5)

#Attenuation
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterA)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.6, p=0.2, r=0, Attenuation")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterB22)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=5,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterB25)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=5,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterB210)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("Knockout","2 nodes halved","5 nodes halved","10 nodes halved"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17),lty=2,cex=0.5)

#Distributed knockout
par(xpd=FALSE)
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterC1)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.6, p=0.2, r=0, Distributed knockout")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterA)))[2]))
abline(v=12)

#Distributed attenuation
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterC2)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.6, p=0.2, r=0, Distributed attenuation")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterA)))[2]))
abline(v=12)

#q=0.8
result4 <- NodeLossModel1(100, 0.2, directed = FALSE, 0.8)
E(result4)$weight <- 1
model14afterA <- list()
model14afterB12 <- list()
model14afterB15 <- list()
model14afterB110 <- list()
model14afterB22 <- list()
model14afterB25 <- list()
model14afterB210 <- list()
model14afterC1 <- list()
model14afterC2 <- list()

for (i in 1:10){
  model14afterA[[i]] <- greedy_node(result4, attackA,25,1) 
  model14afterB12[[i]] <- greedy_node(result4,attackB1,25,2)
  model14afterB15[[i]] <- greedy_node(result4,attackB1,25,5)
  model14afterB110[[i]] <- greedy_node(result4,attackB1,25,10)
  model14afterB22[[i]] <- greedy_node(result4,attackB2,25,2)
  model14afterB25[[i]] <- greedy_node(result4,attackB2,25,5)
  model14afterB210[[i]] <- greedy_node(result4,attackB2,25,10)
  model14afterC1[[i]] <- greedy_edge(result4,attackC1,25,1)
  model14afterC2[[i]] <- greedy_edge(result4,attackC2,25,1)
}


#Knockout
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterA)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=2,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.8, p=0.2, r=0, Knockout")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterB12)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=2,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterB15)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=2,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterB110)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("Knockout","2 nodes halved","5 nodes halved","10 nodes halved"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17,15),lty=2,cex=0.5)

#Attenuation
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterA)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.8, p=0.2, r=0, Attenuation")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterB22)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=5,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterB25)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=5,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterB210)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("Knockout","2 nodes halved","5 nodes halved","10 nodes halved"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17),lty=2,cex=0.5)

#Distributed knockout
par(xpd=FALSE)
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterC1)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.8, p=0.2, r=0, Distributed knockout")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterA)))[2]))
abline(v=12)

#Distributed attenuation
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterC2)))[1]),
     type="b",lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model with q1=0.8, p=0.2, r=0, Distributed attenuation")
abline(h = as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterA)))[2]))
abline(v=12)


#Effect of q
plot.new()
plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterB12)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model, Partial knockout two halved")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterB12)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=5,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterB12)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=5,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterB12)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("q=0.2","q=0.4","q=0.6","q=0.8"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17,15),lty=2,cex=0.5)

plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterB15)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model, Partial knockout five halved")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterB15)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=5,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterB15)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=5,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterB15)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("q=0.2","q=0.4","q=0.6","q=0.8"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17,15),lty=2,cex=0.5)

plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterB22)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model, Partial attenuation")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterB22)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=5,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterB22)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=5,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterB22)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("q=0.2","q=0.4","q=0.6","q=0.8"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17,15),lty=2,cex=0.5)

plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterC1)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model, Distributed knockout")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterC1)))[1]),
     ylim=c(0,1),type="b",pch=22,lty=5,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterC1)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=5,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterC1)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("q=0.2","q=0.4","q=0.6","q=0.8"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17,15),lty=2,cex=1)

plot(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model1afterC2)))[1]),
     ylim=c(0,1),type="b",pch=5,lty=5,col="blue",
     xlab="Number of attacks", ylab="Network efficiency",
     main="Node Loss Model, Distributed attenuation")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model12afterC2)))[1]),
      ylim=c(0,1),type="b",pch=22,lty=5,col="red")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model13afterC2)))[1]),
      ylim=c(0,1),type="b",pch=17,lty=5,col="green")
lines(as.numeric(colMeans(rbindlist(Map(function(x) data.frame(t(x)), model14afterC2)))[1]),
      ylim=c(0,1),type="b",pch=15,lty=2,col="orange")
legend("bottomleft",legend=c("q=0.2","q=0.4","q=0.6","q=0.8"),
       col=c("blue","red","green","orange"),
       pch=c(5,22,17,15),lty=2,cex=1)
