####Packages####
library(igraph)
library(zoo)
library(CTDquerier)
library(CINNA)

####Network####
#Construction of the protein-protein interactions network from string database
setwd("path to data folder")
network_actions <-read.table(file ="9606.protein.actions.v10.5.txt",header = T,sep = "\t",stringsAsFactors = F )

#annotation_string permit to obtain  to convert Ensembl name to protein names.
annotation<-read.table(file ="annotation_string.txt",header = T,sep = "\t",stringsAsFactors = F,quote="" )
annotation<-annotation[,1:2]

#replace Ensembl names by official protein names
network_actions<-merge(x = network_actions, y = annotation,by.x="item_id_a",by.y="protein_external_id")
network_actions<-network_actions[,-which(colnames(network_actions)=="item_id_a")]
colnames(network_actions)[which(colnames(network_actions)=="preferred_name")]<-"item_id_a"
network_actions<-merge(x = network_actions, y = annotation,by.x="item_id_b",by.y="protein_external_id")
network_actions<-network_actions[,-which(colnames(network_actions)=="item_id_b")]
colnames(network_actions)[which(colnames(network_actions)=="preferred_name")]<-"item_id_b"

#Select network with chosen parameters
network_actions_directed<-network_actions[network_actions$is_directional=="t",
                                          -which(colnames(network_actions)=="is_directional")]
network_actions_directed_acting<-network_actions_directed[network_actions_directed$a_is_acting=="t",
                                                          -which(colnames(network_actions_directed)=="a_is_acting")]
network_actions_directed_acting_900<-network_actions_directed_acting[which(network_actions_directed_acting$score>=900),
                                                                     -which(colnames(network_actions_directed_acting)=="score")]
network_actions_directed_acting_900<-network_actions_directed_acting_900[,c("item_id_a","item_id_b")]

#remove redundant relations and self loops
network_actions_directed_acting_900<-unique(network_actions_directed_acting_900)
sum(network_actions_directed_acting_900$item_id_a==network_actions_directed_acting_900$item_id_b)
proteins_actions_directed_acting_900<-unique(c(network_actions_directed_acting_900$item_id_a,network_actions_directed_acting_900$item_id_b))

####input predicted AD associated genes####
setwd("path to result folder")
DE_proteins<-read.csv('Predicted_AD_associated_genes.csv')
DE_proteins<-DE_proteins[,2]

####Contexctualization of the network####
network<-network_actions_directed_acting_900[network_actions_directed_acting_900$item_id_a%in%DE_proteins & network_actions_directed_acting_900$item_id_b%in%DE_proteins,]
network<-unique(network)

#The network is transform as an igraph object
graph<-graph_from_data_frame(network)

####Cliques calculation and select the largest clique#######
clique<-max_cliques(graph = graph,min=3)
proteins_clustered<- unique(names(unlist(clique)))
Fnetwork<- network[network$item_id_a%in%proteins_clustered & network$item_id_b%in%proteins_clustered,]
Fnetwork<-as.matrix(Fnetwork)
Fnetwork<-unique(Fnetwork)
g <- graph_from_edgelist(Fnetwork, directed=TRUE)
c <- components(g)
V(g)$group <- c$membership
BigComp = which.max(c$csize)
Main_g = induced_subgraph(g, which(V(g)$group == BigComp))
pr_cent<-proper_centralities(Main_g)
result<-calculate_centralities(Main_g, include = pr_cent[29])
result<-as.data.frame(do.call(cbind, result))

###output the result###
write.csv(result, file='AD_gene_centrality-scores.csv')

