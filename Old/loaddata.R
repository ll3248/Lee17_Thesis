# loads Twitter, Facebook, and Google+ data from .gz and .txt.gx files 
# Source: https://snap.stanford.edu/data/

#Name	        Type	      Nodes	  Edges	      Description
#ego-Facebook	Undirected	4,039	  88,234	    Social circles from Facebook (anonymized)
#ego-Gplus	  Directed	  107,614	13,673,453	Social circles from Google+
#ego-Twitter	Directed	  81,306	1,768,149	  Social circles from Twitter

#setwd("~/")
setwd("~/Thesis Data")

library(igraph)
library(sand)
library(network)
library(sna)
library(ggplot2)
library("gridExtra")


## FACEBOOK
#facebook <- untar("facebook.tar.gz",list=TRUE) 
##facebook <- untar("facebook.tar.gz") ##outputs 0L
#head(facebook, 20)

facebookcombined <- read.table(gzfile("facebook_combined.txt.gz"), header=F)
head(facebookcombined, 500)

nrow(facebookcombined)
length(unique(c(facebookcombined$V2, facebookcombined$V1)))

fbel <- graph.data.frame(facebookcombined)
no.clusters(fbel)
summary(fbel)

transitivity(fbel, type="global")
transitivity(fbel, type="localaverage") #<- USE THIS ONE FOR ANALYSIS 
transitivity(fbel, type="weighted")

sum(count_triangles(fbel))/3

diameter(fbel)

A_fb <- get.adjacency(fbel, sparse = FALSE)
fbg <- network::as.network.matrix(A_fb)

hist(degree(fbg))

ggplot(degree(fbg), aes(x=degree(fbg))) + 
  geom_histogram(aes(y = ..density..), colour = "darkgreen", fill = "white", binwidth = 50) + 
  geom_density() + 
  ggtitle("Degree Distribution for Facebook Network")


### Modeling the Facebook Network 

#### Erdos-Renyi-Gilbert 

#```{r, echo=FALSE, fig.align='center', fig.height=2, fig.width=5, warning=FALSE}

set.seed(499)


numsim <- 5

g.er.fbel.coef <- rep(NA, numsim)
g.er.fbel.apl <- rep(NA, numsim)
g.er.fbel.dia <- rep(NA, numsim)

p = ecount(fbel)/choose(vcount(fbel), 2)

for (i in 1:numsim) {
  # n = number of vertices, p = probability of a link
  #p is caluclated as number of edges over number of possible edges
  #p is then (34 choose 2) since there is the possibility of a link between
  g.er.fbel <- erdos.renyi.game(n = vcount(fbel), p)
  g.er.fbel.coef[i] <- transitivity(g.er.fbel, type="global")
  g.er.fbel.apl[i] <- average.path.length(g.er.fbel)
  g.er.fbel.dia[i] <- diameter(g.er.fbel)
}


gplot1 <- ggplot2::ggplot(data = as.data.frame(g.er.fbel.coef), 
                 aes(x = g.er.fbel.coef)) +
  geom_histogram(binwidth=0.0065, 
                 colour="blue", 
                 fill="lavender") + 
  geom_vline(xintercept = transitivity(fbel, type="global"), colour="red") + 
  labs(title="Simulated 
       Clustering Coefficients") + 
  xlab("Clustering Coefficient") + 
  ylab("Frequency") + 
  theme(legend.position="none", 
        plot.title = element_text(size=6),
        axis.text = element_text(size=4),
        axis.title = element_text(size=4)) 


gplot2 <- ggplot(data = as.data.frame(g.er.fbel.apl), aes(x = g.er.fbel.apl)) +
  geom_histogram(binwidth=0.031,
                 colour="blue",
                 fill = "lavender") + 
  geom_vline(xintercept = average.path.length(fbel), colour="red") + 
  labs(title="Simulated 
       Average Path Lengths") + 
  xlab("Average Path Length") +
  ylab("") + 
  theme(legend.position="none", 
        plot.title = element_text(size=6),
        axis.text = element_text(size=4),
        axis.title = element_text(size=4)) 

gplot3 <- ggplot(data = as.data.frame(g.er.fbel.dia), aes(x = g.er.fbel.dia)) +
  geom_histogram(binwidth=0.12,
                 colour="blue",
                 fill = "lavender") + 
  geom_vline(xintercept = diameter(fbel), colour="red") + 
  labs(title="Simulated Diameters") + 
  xlab("Diameter") + 
  ylab("") + 
  theme(legend.position="none", 
        plot.title = element_text(size=6),
        axis.text = element_text(size=4),
        axis.title = element_text(size=4))

arrangeGrob(gplot1, gplot2, gplot3, ncol=3)

#```


#### Watts-Strogatz

numsim <- 100

g.er.fbel.coef <- rep(NA, numsim)
g.er.fbel.apl <- rep(NA, numsim)
g.er.fbel.dia <- rep(NA, numsim)














#sna::gplot.target(fbg, degree(fbg), main = "Degree", 
#                  circ.lab = FALSE, cir.col = "skyblue", 
#                  usearrows = FALSE, 
#                  edge.col = "darkgray")


## TWITTER 
#twitter <- untar("twitter.tar.gz",list=TRUE) 
##twitter <- untar("twitter.tar.gz") #outputs 0L
#head(twitter, 20)

#twittercombined <- read.table(gzfile("twitter_combined.txt.gz"), header=F)

#head(twittercombined, 20)
# outputs an edge list?

#nrow(twittercombined)
#length(unique(c(twittercombined$V2, twittercombined$V1)))

#twel <- graph.data.frame(twittercombined)
#no.clusters(twel)





## GOOGLE PLUS
#gplus <- untar("gplus.tar.gz",list=TRUE) 
#gplus <- untar("gplus.tar.gz") #outputs 0L
#head(gplus, 20)

#gpluscombined <- read.table(gzfile("gplus_combined.txt.gz"), header=F)
#head(gpluscombined, 20)

#gpel <- graph.data.frame(gpluscombined)
#no.clusters(gpel)

#length(unique(c(gpluscombined$V2, gpluscombined$V1)))
