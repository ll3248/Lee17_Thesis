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
#setwd("~/STAT495-Lee/Lee Thesis Data")
setwd("~/")

library(sand)
library(igraph)
library(network)
library(sna)
library(statnet)
library(ergm)
library(intergraph)

### FACEBOOK

#setwd("~/Lee17Thesis/Thesis/data")
setwd("~/")
# this is an edge list
facebookcombined <- read.table(gzfile("facebook_combined.txt.gz"), header=F)


# creating `igraph` object
fbel <- graph.data.frame(facebookcombined)

# creating `network` object
A_fb <- get.adjacency(fbel, sparse = FALSE)
fbg <- network::as.network.matrix(A_fb, directed = FALSE)


# characteristics of the Facebook component
nrow(facebookcombined) #88234 edges
length(unique(c(facebookcombined$V2, facebookcombined$V1))) #4039 nodes
no.clusters(fbel) #1 cluster
sum(count_triangles(fbel))/3 #1612010 unique triangles (up to ordering)

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


#### ERGMs

# needs an object of class network
# uses `fbg`, a `network` object

## one paramter: edges
# this should be the same as Erdos-Renyi
# set the seed for reproducible analysis 
set.seed(499)

# model 1a
g.ergm1a.fbg <- ergm(fbg ~ edges)

numsim = 1000
set.seed(499)

g.ergm1a.fbg.sim.ecount <- rep(NA, numsim)
g.ergm1a.fbg.sim.coef <- rep(NA, numsim)
g.ergm1a.fbg.sim.apl <- rep(NA, numsim)
g.ergm1a.fbg.sim.dia <- rep(NA, numsim)
g.ergm1a.fbg.sim.avgdeg <- rep(NA, numsim)
g.ergm1a.fbg.sim.avgbtwcen <- rep(NA, numsim)
g.ergm1a.fbg.sim.avgclocen <- rep(NA, numsim)
g.ergm1a.fbg.sim.avgeigveccen <- rep(NA, numsim)

for (i in 1:numsim) {
  #simulate graphs from ergms model one at a time
  #convert to igraph object
  g.ergm.fbg.sim.convert <- asIgraph(simulate(g.ergm1a.fbg, nsim = 1))
  
  #observe the network statistics for this simulated graph
  g.ergm1a.fbg.sim.ecount[i] <- ecount(g.ergm.fbg.sim.convert)
  g.ergm1a.fbg.sim.coef[i] <- transitivity(g.ergm.fbg.sim.convert, type="localaverage")
  g.ergm1a.fbg.sim.apl[i] <- average.path.length(g.ergm.fbg.sim.convert)
  g.ergm1a.fbg.sim.dia[i] <- diameter(g.ergm.fbg.sim.convert)
  g.ergm1a.fbg.sim.avgdeg[i] <- mean(igraph::degree(g.ergm.fbg.sim.convert))
  g.ergm1a.fbg.sim.avgbtwcen[i] <- mean(igraph::betweenness(g.ergm.fbg.sim.convert))
  g.ergm1a.fbg.sim.avgclocen[i] <- mean(igraph::closeness(g.ergm.fbg.sim.convert))
  g.ergm1a.fbg.sim.avgeigveccen[i] <- mean(igraph::eigen_centrality(g.ergm.fbg.sim.convert)$vector)
}

g.ergm1a.fbg.sim.df <- as.data.frame(cbind(g.ergm1a.fbg.sim.ecount, g.ergm1a.fbg.sim.coef, 
                                           g.ergm1a.fbg.sim.apl, g.ergm1a.fbg.sim.dia, 
                                           g.ergm1a.fbg.sim.avgdeg, g.ergm1a.fbg.sim.avgbtwcen, 
                                           g.ergm1a.fbg.sim.avgclocen, g.ergm1a.fbg.sim.avgeigveccen))

# store all the data in the folder below
# setwd("~/Lee17Thesis/Thesis/data")

# save the data frame for future use
save(g.ergm1a.fbg.sim.df, file = "g.ergm1a.fbg.sim.Rda")

# load the data frame for future use
# load("g.ergm1a.fbg.sim.Rda")

# attach previously saved data frame
attach(g.ergm1a.fbg.sim.df)

# do not forget to detach data frame below

# average network statistics using ERGM 1a 

mean(g.ergm1a.fbg.sim.coef) #0.3696413
sd(g.ergm1a.fbg.sim.coef) #0.001246206

mean(g.ergm1a.fbg.sim.apl) #2.885094
sd(g.ergm1a.fbg.sim.apl) #0.004316369

mean(g.ergm1a.fbg.sim.dia) #5.216
sd(g.ergm1a.fbg.sim.dia) #0.4117202

mean(g.ergm1a.fbg.sim.avgdeg) #43.66728
sd(g.ergm1a.fbg.sim.avgdeg) #0.05240351

mean(g.ergm1a.fbg.sim.avgbtwcen) #3805.712
sd(g.ergm1a.fbg.sim.avgbtwcen) #8.757621

mean(g.ergm1a.fbg.sim.avgclocen) #8.285989e-05
sd(g.ergm1a.fbg.sim.avgclocen) #8.353233e-06

mean(g.ergm1a.fbg.sim.avgeigveccen) #0.04179492
sd(g.ergm1a.fbg.sim.avgeigveccen) #0.0008307558

# detach previously saved data frame
detach(g.ergm1a.fbg.sim.df)



## two parameters: edges and triangles 
# set the seed for reproducible analysis 
set.seed(499)

# model 2a
# WARNINGS: takes time (20 iterations)
g.ergm2a.fbg <- ergm(fbg ~ edges + triangle)

numsim = 1000
set.seed(499)

g.ergm2a.fbg.sim.ecount <- rep(NA, numsim)
g.ergm2a.fbg.sim.coef <- rep(NA, numsim)
g.ergm2a.fbg.sim.apl <- rep(NA, numsim)
g.ergm2a.fbg.sim.dia <- rep(NA, numsim)
g.ergm2a.fbg.sim.avgdeg <- rep(NA, numsim)
g.ergm2a.fbg.sim.avgbtwcen <- rep(NA, numsim)
g.ergm2a.fbg.sim.avgclocen <- rep(NA, numsim)
g.ergm2a.fbg.sim.avgeigveccen <- rep(NA, numsim)

for (i in 1:numsim) {
  #simulate graphs from ergms model one at a time
  #convert to igraph object
  g.ergm.fbg.sim.convert <- asIgraph(simulate(g.ergm2a.fbg, nsim = 1))
  
  #observe the network statistics for this simulated graph
  g.ergm2a.fbg.sim.ecount[i] <- ecount(g.ergm.fbg.sim.convert)
  g.ergm2a.fbg.sim.coef[i] <- transitivity(g.ergm.fbg.sim.convert, type="localaverage")
  g.ergm2a.fbg.sim.apl[i] <- average.path.length(g.ergm.fbg.sim.convert)
  g.ergm2a.fbg.sim.dia[i] <- diameter(g.ergm.fbg.sim.convert)
  g.ergm2a.fbg.sim.avgdeg[i] <- mean(igraph::degree(g.ergm.fbg.sim.convert))
  g.ergm2a.fbg.sim.avgbtwcen[i] <- mean(igraph::betweenness(g.ergm.fbg.sim.convert))
  g.ergm2a.fbg.sim.avgclocen[i] <- mean(igraph::closeness(g.ergm.fbg.sim.convert))
  g.ergm2a.fbg.sim.avgeigveccen[i] <- mean(igraph::eigen_centrality(g.ergm.fbg.sim.convert)$vector)
}

g.ergm2a.fbg.sim.df <- as.data.frame(cbind(g.ergm2a.fbg.sim.ecount, g.ergm2a.fbg.sim.coef, 
                                           g.ergm2a.fbg.sim.apl, g.ergm2a.fbg.sim.dia, 
                                           g.ergm2a.fbg.sim.avgdeg, g.ergm2a.fbg.sim.avgbtwcen, 
                                           g.ergm2a.fbg.sim.avgclocen, g.ergm2a.fbg.sim.avgeigveccen))

# store all the data in the folder below
# setwd("~/Lee17Thesis/Thesis/data")

# save the data frame for future use
# save(g.ergm2a.fbg.sim.df, file="g.ergm2a.fbg.sim.Rda")

# load the data frame for future use
# load("g.ergm2a.fbg.sim.Rda")

# attach previously saved data frame
# attach(g.ergm2a.fbg.sim.df)

# do not forget to detach data frame below

# average network statistics using ERGM 2a 

mean(g.ergm2a.fbg.sim.coef) #0.4823232
sd(g.ergm2a.fbg.sim.coef) #0.001965885

mean(g.ergm2a.fbg.sim.apl) #3.052008
sd(g.ergm2a.fbg.sim.apl) #0.008639291

mean(g.ergm2a.fbg.sim.dia) #6.098
sd(g.ergm2a.fbg.sim.dia) #0.3008097

mean(g.ergm2a.fbg.sim.avgdeg) #44.54161
sd(g.ergm2a.fbg.sim.avgdeg) #0.03491666

mean(g.ergm2a.fbg.sim.avgbtwcen) #4140.19
sd(g.ergm2a.fbg.sim.avgbtwcen) #17.65641

mean(g.ergm2a.fbg.sim.avgclocen) #6.007678e-05
sd(g.ergm2a.fbg.sim.avgclocen) #1.513646e-05

mean(g.ergm2a.fbg.sim.avgeigveccen) #0.04103567
sd(g.ergm2a.fbg.sim.avgeigveccen) #3.577786e-05

# detach previously saved data frame
# detach(g.ergm2a.fbg.sim.df)






## two parameters: edges and k-stars
# WARNINGS: takes time (18 iterations)
# set the seed for reproducible analysis 
set.seed(499)

# model 2b
g.ergm2b.fbg <- ergm(fbg ~ edges + kstar(3))

numsim = 1000
set.seed(499)

g.ergm2b.fbg.sim.ecount <- rep(NA, numsim)
g.ergm2b.fbg.sim.coef <- rep(NA, numsim)
g.ergm2b.fbg.sim.apl <- rep(NA, numsim)
g.ergm2b.fbg.sim.dia <- rep(NA, numsim)
g.ergm2b.fbg.sim.avgdeg <- rep(NA, numsim)
g.ergm2b.fbg.sim.avgbtwcen <- rep(NA, numsim)
g.ergm2b.fbg.sim.avgclocen <- rep(NA, numsim)
g.ergm2b.fbg.sim.avgeigveccen <- rep(NA, numsim)

for (i in 1:numsim) {
  #simulate graphs from ergms model one at a time
  #convert to igraph object
  g.ergm.fbg.sim.convert <- asIgraph(simulate(g.ergm2b.fbg, nsim = 1))
  
  #observe the network statistics for this simulated graph
  g.ergm2b.fbg.sim.ecount[i] <- ecount(g.ergm.fbg.sim.convert)
  g.ergm2b.fbg.sim.coef[i] <- transitivity(g.ergm.fbg.sim.convert, type="localaverage")
  g.ergm2b.fbg.sim.apl[i] <- average.path.length(g.ergm.fbg.sim.convert)
  g.ergm2b.fbg.sim.dia[i] <- diameter(g.ergm.fbg.sim.convert)
  g.ergm2b.fbg.sim.avgdeg[i] <- mean(igraph::degree(g.ergm.fbg.sim.convert))
  g.ergm2b.fbg.sim.avgbtwcen[i] <- mean(igraph::betweenness(g.ergm.fbg.sim.convert))
  g.ergm2b.fbg.sim.avgclocen[i] <- mean(igraph::closeness(g.ergm.fbg.sim.convert))
  g.ergm2b.fbg.sim.avgeigveccen[i] <- mean(igraph::eigen_centrality(g.ergm.fbg.sim.convert)$vector)
}

g.ergm2b.fbg.sim.df <- as.data.frame(cbind(g.ergm2b.fbg.sim.ecount, g.ergm2b.fbg.sim.coef, 
                                           g.ergm2b.fbg.sim.apl, g.ergm2b.fbg.sim.dia, 
                                           g.ergm2b.fbg.sim.avgdeg, g.ergm2b.fbg.sim.avgbtwcen, 
                                           g.ergm2b.fbg.sim.avgclocen, g.ergm2b.fbg.sim.avgeigveccen))

# store all the data in the folder below
# setwd("~/Lee17Thesis/Thesis/data")

# save the data frame for future use
# save(g.ergm2b.fbg.sim.df, file="g.ergm2b.fbg.sim.Rda")

# load the data frame for future use
# load("g.ergm2b.fbg.sim.Rda")

# attach previously saved data frame
# attach(g.ergm2b.fbg.sim.df)

# do not forget to detach data frame below

# average network statistics using ERGM 2b

mean(g.ergm2b.fbg.sim.coef) #0.3786565
sd(g.ergm2b.fbg.sim.coef) #0.001279854

mean(g.ergm2b.fbg.sim.apl) #2.850693
sd(g.ergm2b.fbg.sim.apl) #0.003085883

mean(g.ergm2b.fbg.sim.dia) #5.095
sd(g.ergm2b.fbg.sim.dia) #0.2933617

mean(g.ergm2b.fbg.sim.avgdeg) #44.0645
sd(g.ergm2b.fbg.sim.avgdeg) #0.05715313

mean(g.ergm2b.fbg.sim.avgbtwcen) #3736.415
sd(g.ergm2b.fbg.sim.avgbtwcen) #6.259627

mean(g.ergm2b.fbg.sim.avgclocen) #8.558362e-05
sd(g.ergm2b.fbg.sim.avgclocen) #6.024968e-06

mean(g.ergm2b.fbg.sim.avgeigveccen) #0.03923439
sd(g.ergm2b.fbg.sim.avgeigveccen) #0.0002047637

# detach previously saved data frame
# detach(g.ergm2b.fbg.sim.df)







## three parameters: edges and triangles and k-stars
# set the seed for reproducible analysis 
set.seed(499)

# model 3a
g.ergm3a.fbg <- ergm(fbg ~ edges + triangle + kstar(3))

numsim = 1000
set.seed(499)

g.ergm3a.fbg.sim.ecount <- rep(NA, numsim)
g.ergm3a.fbg.sim.coef <- rep(NA, numsim)
g.ergm3a.fbg.sim.apl <- rep(NA, numsim)
g.ergm3a.fbg.sim.dia <- rep(NA, numsim)
g.ergm3a.fbg.sim.avgdeg <- rep(NA, numsim)
g.ergm3a.fbg.sim.avgbtwcen <- rep(NA, numsim)
g.ergm3a.fbg.sim.avgclocen <- rep(NA, numsim)
g.ergm3a.fbg.sim.avgeigveccen <- rep(NA, numsim)

for (i in 1:numsim) {
  #simulate graphs from ergms model one at a time
  #convert to igraph object
  g.ergm.fbg.sim.convert <- asIgraph(simulate(g.ergm2a.fbg, nsim = 1))
  
  #observe the network statistics for this simulated graph
  g.ergm3a.fbg.sim.ecount[i] <- ecount(g.ergm.fbg.sim.convert)
  g.ergm3a.fbg.sim.coef[i] <- transitivity(g.ergm.fbg.sim.convert, type="localaverage")
  g.ergm3a.fbg.sim.apl[i] <- average.path.length(g.ergm.fbg.sim.convert)
  g.ergm3a.fbg.sim.dia[i] <- diameter(g.ergm.fbg.sim.convert)
  g.ergm3a.fbg.sim.avgdeg[i] <- mean(igraph::degree(g.ergm.fbg.sim.convert))
  g.ergm3a.fbg.sim.avgbtwcen[i] <- mean(igraph::betweenness(g.ergm.fbg.sim.convert))
  g.ergm3a.fbg.sim.avgclocen[i] <- mean(igraph::closeness(g.ergm.fbg.sim.convert))
  g.ergm3a.fbg.sim.avgeigveccen[i] <- mean(igraph::eigen_centrality(g.ergm.fbg.sim.convert)$vector)
}

g.ergm3a.fbg.sim.df <- as.data.frame(cbind(g.ergm3a.fbg.sim.ecount, g.ergm3a.fbg.sim.coef, 
                                           g.ergm3a.fbg.sim.apl, g.ergm3a.fbg.sim.dia, 
                                           g.ergm3a.fbg.sim.avgdeg, g.ergm3a.fbg.sim.avgbtwcen, 
                                           g.ergm3a.fbg.sim.avgclocen, g.ergm3a.fbg.sim.avgeigveccen))

# store all the data in the folder below
# setwd("~/Lee17Thesis/Thesis/data")

# save the data frame for future use
# save(g.ergm3a.fbg.sim.df, file="g.ergm3a.fbg.sim.Rda")

# load the data frame for future use
# load("g.ergm3a.fbg.sim.Rda")

# attach previously saved data frame
# attach(g.ergm3a.fbg.sim.df)

# do not forget to detach data frame below

# average network statistics using ERGM 3a

mean(g.ergm3a.fbg.sim.coef) #0.4823232
sd(g.ergm3a.fbg.sim.coef) #0.001965885

mean(g.ergm3a.fbg.sim.apl) #3.052008
sd(g.ergm3a.fbg.sim.apl) #0.008639291

mean(g.ergm3a.fbg.sim.dia) #6.098
sd(g.ergm3a.fbg.sim.dia) #0.3008097

mean(g.ergm3a.fbg.sim.avgdeg) #44.54161
sd(g.ergm3a.fbg.sim.avgdeg) #0.03491666

mean(g.ergm3a.fbg.sim.avgbtwcen) #4140.19
sd(g.ergm3a.fbg.sim.avgbtwcen) #17.65641

mean(g.ergm3a.fbg.sim.avgclocen) #6.007678e-05
sd(g.ergm3a.fbg.sim.avgclocen) #1.513646e-05

mean(g.ergm3a.fbg.sim.avgeigveccen) #0.04103567
sd(g.ergm3a.fbg.sim.avgeigveccen) #3.577786e-05

# detach previously saved data frame
# detach(g.ergm3a.fbg.sim.df)






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
# setwd("~/Lee17Thesis/Thesis/data")

# save the data frame for future use
# save(g.er.fbel, file="g.er.fbel.Rda")

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
# setwd("~/Lee17Thesis/Thesis/data")

# save the data frame for future use
# save(g.ws.fbel, file="g.ws.fbel.Rda")

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

