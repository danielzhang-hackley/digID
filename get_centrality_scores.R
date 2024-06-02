####Packages####
library(igraph)
library(zoo)
library(CTDquerier)
library(CINNA)

#Construction of the protein-protein interactions network from String database
setwd('path to data files')
network_actions <-read.table(file ="9606.protein.actions.v10.5.txt",header = T,sep = "\t",stringsAsFactors = F )

#convert Ensembl name to protein names.
annotation<-read.table(file ="annotation_string.txt",header = T,sep = "\t",stringsAsFactors = F,quote="" )
annotation<-annotation[,1:2]

#replace Ensembl names by regular protein names
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

####input predicted AD genes####
setwd('path to the predicted AD gene file')
DE_proteins<-read.csv('./predicted_AD_genes.csv')
DE_proteins<-DE_proteins[,1]

#add relations with DE_proteins under a score of 900 to the network
network<-network_actions_directed_acting_900[network_actions_directed_acting_900$item_id_a%in%DE_proteins & network_actions_directed_acting_900$item_id_b%in%DE_proteins,]
network<-unique(network)

#transform the network to an igraph object
graph<-graph_from_data_frame(network)

####Cliques calculation and select the largest clique#######
clique<-max_cliques(graph = graph,min=3)
proteins_clustered<- unique(names(unlist(clique)))
Fnetwork<- network[network$item_id_a%in%proteins_clustered & network$item_id_b%in%proteins_clustered,]
Fnetwork<-as.matrix(Fnetwork)
Fnetwork<-unique(Fnetwork)

## calculate centrality scores
g <- graph_from_edgelist(Fnetwork, directed=TRUE)
c <- components(g)
V(g)$group <- c$membership
BigComp = which.max(c$csize)
Main_g = induced_subgraph(g, which(V(g)$group == BigComp))
pr_cent<-proper_centralities(Main_g)
result<-calculate_centralities(Main_g, include = pr_cent[29])
result<-as.data.frame(do.call(cbind, result))
result=cbind(rownames(result),result[,1])
colnames(result)=c('gene_name', 'centrality_score')
result=as.data.frame(result)
result = result[order(result$centrality_score, decreasing=TRUE),]

# output the result
write.csv(result, file='path to result file',row.names=FALSE)

