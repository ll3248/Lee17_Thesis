# README file for StatComps2017-02
# This folder immediately contains all of the components of knitting the file. 
# Simply hit the Knit button on `Lee_Thesis.Rmd` in order to generate a compiled `.pdf`, titles `Lee_Thesis.pdf`

# Abstract, Chapter 0 (Introduction), Chapter 1, Chapter 2, Chapter 3, and Chapter 4 (Conclusion), and Bibliography are self-explanatory
# They correspond to `abstract.Rmd`, `chap0.Rmd`, `chap1.Rmd`, `chap2.Rmd`, `chap3.Rmd`, `conclusion.Rmd`, `bibliography.Rmd`, respectively

# Saved datasets are in the `data` folder; reproducible code to get the `g.er.fbel.Rda` and `g.ws.fbel.Rda` files are shown below

# The `facebook_comined.txt.gx` from comes from the following source: https://snap.stanford.edu/data/

# Name	        Type	      Nodes	  Edges	      Description
# ego-Facebook	Undirected	4,039	  88,234	    Social circles from Facebook (anonymized)

# WARNING: Expected run time for simulation study is about 2-4 hours


### Required Packages
setwd("~/STAT495-Lee/Lee Thesis Data")

library(sand)
library(igraph)
library(network)
library(sna)
library(statnet)
library(ergm)

### FACEBOOK

# this is an edge list
setwd("~/STAT495-Lee/StatComps2017-02 /data")
facebookcombined <- read.table(gzfile("facebook_combined.txt.gz"), header=F)

# characteristics of the Facebook component
nrow(facebookcombined) #88234 edges
length(unique(c(facebookcombined$V2, facebookcombined$V1))) #4039 nodes
no.clusters(fbel) #1 cluster
sum(count_triangles(fbel))/3 #1612010 unique triangles (up to ordering)

# creating `igraph` object
fbel <- graph.data.frame(facebookcombined)

# creating `network` object
A_fb <- get.adjacency(fbel, sparse = FALSE)
fbg <- network::as.network.matrix(A_fb)

# some network statisitcs
transitivity(fbel, type="localaverage")  #0.6170038
diameter(fbel) #17
average.path.length(fbel) #4.337744

# centrality measures 
mean(igraph::degree(fbel)) #43.69101
mean(igraph::betweenness(fbel)) #2072.642
mean(igraph::closeness(fbel)) #8.881448e-08
mean(igraph::eigen_centrality(fbel)$vector) #0.04047316

### Modeling the Facebook Network 

#### Erdős-Rényi(-Gilbert)

# set the seed for reproducible analysis 
set.seed(499)

numsim <- 1000

g.er.fbel.ecount <- rep(NA, numsim)
g.er.fbel.coef <- rep(NA, numsim)
g.er.fbel.apl <- rep(NA, numsim)
g.er.fbel.dia <- rep(NA, numsim)
g.er.fbel.avgdeg <- rep(NA, numsim)
g.er.fbel.avgbtwcen <- rep(NA, numsim)
g.er.fbel.avgclocen <- rep(NA, numsim)
g.er.fbel.avgeigveccen <- rep(NA, numsim)

# value needed for simulation
p = ecount(fbel)/choose(vcount(fbel), 2)

# WARNING: Expected run time for simulation is about 2-4 hours
for (i in 1:numsim) {
  # n = number of vertices, p = probability of a link
  # p is caluclated as number of edges over number of possible edges
  # p is then (4039 choose 2) since there is the possibility of a link between
  # igraph object
  g.er.fbel <- erdos.renyi.game(n = vcount(fbel), p)
  
  # observe the network statistics for this simulated graph
  g.er.fbel.ecount[i] <- ecount(g.er.fbel)
  g.er.fbel.coef[i] <- transitivity(g.er.fbel, type="localaverage")
  g.er.fbel.apl[i] <- average.path.length(g.er.fbel)
  g.er.fbel.dia[i] <- diameter(g.er.fbel)
  g.er.fbel.avgdeg[i] <- mean(igraph::degree(g.er.fbel))
  g.er.fbel.avgbtwcen[i] <- mean(igraph::betweenness(g.er.fbel))
  g.er.fbel.avgclocen[i] <- mean(igraph::closeness(g.er.fbel))
  g.er.fbel.avgeigveccen[i] <- mean(igraph::eigen_centrality(g.er.fbel)$vector)
}

g.er.fbel <- as.data.frame(cbind(g.er.fbel.ecount, g.er.fbel.coef, 
                                 g.er.fbel.apl, g.er.fbel.dia, 
                                 g.er.fbel.avgdeg, g.er.fbel.avgbtwcen, 
                                 g.er.fbel.avgclocen, g.er.fbel.avgeigveccen))

# store all the data in the folder below
# setwd("~/STAT495-Lee/StatComps2017-02 /data")

# save the data frame for furture use
# save(g.er.fbel,file="g.er.fbel.Rda")

# load the data frame for future use
# load("g.er.fbel.Rda")

# attach and detach previously saved data frame
# attach(g.er.fbel)

# do not forget to detach data frame below

# average network statistics using Erdős-Rényi model

avg.g.er.fbel.coef <- mean(g.er.fbel.coef); avg.g.er.fbel.coef #0.01082603
sd(g.er.fbel.coef) #0.000101139
confint(g.er.fbel.coef)

avg.g.er.fbel.apl <- mean(g.er.fbel.apl); avg.g.er.fbel.apl #2.605739
sd(g.er.fbel.apl) #0.001976112

avg.g.er.fbel.dia <- mean(g.er.fbel.dia); avg.g.er.fbel.dia #3.956
sd(g.er.fbel.dia) #0.2051977

avg.g.er.fbel.avgdeg <- mean(g.er.fbel.avgdeg); avg.g.er.fbel.avgdeg #43.69406
sd(g.er.fbel.avgdeg) #0.143238

avg.g.er.fbel.avgbtwcen <- mean(g.er.fbel.avgbtwcen); avg.g.er.fbel.avgbtwcen #3241.987
sd(g.er.fbel.avgbtwcen) #3.989771

avg.g.er.fbel.avgclocen <- mean(g.er.fbel.avgclocen); avg.g.er.fbel.avgclocen #9.507038e-05
sd(g.er.fbel.avgclocen) #7.229988e-08

avg.g.er.fbel.avgeigveccen <- mean(g.er.fbel.avgeigveccen); avg.g.er.fbel.avgeigveccen #0.6202335
sd(g.er.fbel.avgeigveccen) #0.02240826

# detach previously saved data frame
# detach(g.er.fbel)

#### Watts-Strogatz

# set the seed again for reproducible analysis
set.seed(499)

numsim <- 1000

g.ws.fbel.ecount <- rep(NA,numsim)
g.ws.fbel.coef <- rep(NA, numsim)
g.ws.fbel.apl <- rep(NA, numsim)
g.ws.fbel.dia <- rep(NA, numsim)
g.ws.fbel.avgdeg <- rep(NA, numsim)
g.ws.fbel.avgbtwcen <- rep(NA, numsim)
g.ws.fbel.avgclocen <- rep(NA, numsim)
g.ws.fbel.avgeigveccen <- rep(NA, numsim)

# WARNING: Expected run time for simulation is about 2-4 hours
for (i in 1:numsim) {
  #create a lattice with the same number of vertices as 'fbel'
  #let the starting number of neighbors be equal to the minimum vertex degree of 'fbel'
  g.ws.fbel <- watts.strogatz.game(dim = 1, size = vcount(fbel), nei = min(degree(fbg)), p = 0)
  
  #generate list of random vertex values to add edges to lattice
  #each pair in the list represents one random edge
  randomedgepairs <- sample(1:vcount(fbel), 2*(ecount(fbel)-vcount(fbel)), replace=TRUE)
  g.ws.fbel.pre <- add_edges(g.ws.fbel, randomedgepairs)
  
  
  #igraph object
  g.ws.fbel <- simplify(g.ws.fbel.pre)
  
  #observe the network statistics for this simulated graph
  g.ws.fbel.ecount[i] <- ecount(g.ws.fbel)
  g.ws.fbel.coef[i] <- transitivity(g.ws.fbel, type="localaverage")
  g.ws.fbel.apl[i] <- average.path.length(g.ws.fbel)
  g.ws.fbel.dia[i] <- diameter(g.ws.fbel)
  g.ws.fbel.avgdeg[i] <- mean(igraph::degree(g.ws.fbel))
  g.ws.fbel.avgbtwcen[i] <- mean(igraph::betweenness(g.ws.fbel))
  g.ws.fbel.avgclocen[i] <- mean(igraph::closeness(g.ws.fbel))
  g.ws.fbel.avgeigveccen[i] <- mean(igraph::eigen_centrality(g.ws.fbel)$vector)
}

g.ws.fbel <- as.data.frame(cbind(g.ws.fbel.ecount, g.ws.fbel.coef, 
                                 g.ws.fbel.apl, g.ws.fbel.dia, 
                                 g.ws.fbel.avgdeg, g.ws.fbel.avgbtwcen, 
                                 g.ws.fbel.avgclocen, g.ws.fbel.avgeigveccen))

# store all the data in the folder below
# setwd("~/STAT495-Lee/StatComps2017-02 /data")

# save the data frame for furture use
# save(g.ws.fbel,file="g.ws.fbel.Rda")

# load the data frame for future use
# load("g.ws.fbel.Rda")

# attach previously saved data frame
# attach(g.ws.fbel)

# do not forget to detach data frame below

# average network statistics using Watts-Strogatz model

avg.g.ws.fbel.coef <- mean(g.ws.fbel.coef); avg.g.ws.fbel.coef #0.01073166
sd(g.ws.fbel.coef) #9.135127e-05

avg.g.ws.fbel.apl <- mean(g.ws.fbel.apl); avg.g.ws.fbel.apl #2.609323
sd(g.ws.fbel.apl) #0.000193598

avg.g.ws.fbel.dia <- mean(g.ws.fbel.dia); avg.g.ws.fbel.dia #3.95
sd(g.ws.fbel.dia) #0.218054

avg.g.ws.fbel.avgdeg <- mean(g.ws.fbel.avgdeg); avg.g.ws.fbel.avgdeg #43.44555
sd(g.ws.fbel.avgdeg) #0.01105388

avg.g.ws.fbel.avgbtwcen <- mean(g.ws.fbel.avgbtwcen); avg.g.ws.fbel.avgbtwcen #3249.223
sd(g.ws.fbel.avgbtwcen) # 0.3908744

avg.g.ws.fbel.avgclocen <- mean(g.ws.fbel.avgclocen); avg.g.ws.fbel.avgclocen #9.493807e-05
sd(g.ws.fbel.avgclocen) #7.319204e-09

avg.g.ws.fbel.avgeigveccen <- mean(g.ws.fbel.avgeigveccen); avg.g.ws.fbel.avgeigveccen #0.6234877
sd(g.ws.fbel.avgeigveccen) #0.02269783

# detach previously saved data frame
# detach(g.ws.fbel)

#### ERGMs

# needs an object of class network
# uses `fbg` object

# this should be the same as Erdos-Renyi
# one paramter: edges
g.ergm1.fbg <- ergm(fbg ~ edges)

# set the seed for reproducible analysis 
set.seed(499)
g.ergm1.fbg.sim <- simulate(g.ergm1.fbg, numsim=1000)

# two parameters: edges and triangles 
g.ergm2.fbg <- ergm(fbg ~ edges+triangle)

# set the seed for reproducible analysis 
set.seed(499)
g.ergm2.fbg.sim <- simulate(g.ergm2.fbg, numsim=1000)
