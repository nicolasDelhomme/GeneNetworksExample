#' ---
#' title: "Gene Of Interest First Degree Neighbour"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(here)
  library(igraph)
})

#' * Data
edgelist <- read.table(here("~/Git/GeneNetworksExample/backbone-edgelist.txt"),
                       header=FALSE,as.is=TRUE)

#' # Network
#' We create a graf
graf <- graph.edgelist(as.matrix(edgelist))

#' We can check the structure of the graf
#' 
#' There are 68 subgrafs, one with almost all genes, and the rest involving mainly 2 genes each
#' (look at the $no and $csize list elements)
clusters(graf)$csize

#' We can look at the neighborhood of RB1 and TP53
goi <- c("TP53","RB1")

#' Extract the first degree neighbours of these genes from the network
subgraf <- make_ego_graph(graf,1,
                          get.vertex.attribute(graf,"name") %in% goi)

#' very few genes
barplot(table(sapply(lapply(subgraf,clusters),"[[","csize")),
        las=2,main="Gene of interest cluster size",
        ylab="occurence",xlab="csize")

#' combine all these networks together
fdn <- Reduce("%u%",subgraf)

#' Look at how many clusters we get and how many nodes are involved
clusters(fdn)

#' Look at a static plot
plot(fdn)

#' Extend to second degree neighbours
subgraf <- make_ego_graph(graf,2,
                          get.vertex.attribute(graf,"name") %in% goi)

#' Look at a static plot for TP53
plot(subgraf[[1]])

#' Look at a static plot for RB1
plot(subgraf[[2]])

#' Let's export the RB1 data for visualisation
write_graph(subgraf[[2]],format = "graphml",file="RB1SecondDegreeNeighbour.graphml")

#' # Session Info
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```

