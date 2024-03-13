#' Create a dynamic network representing the relationship between two groups of variables
#'
#' @param data A data frame containing the relationship between the two groups to be represented in the network 
#' @param g1 Name of the column containing the labels of the first group that will be used to create the network
#' @param g2 Name of the column containing the labels of the second group that will be used to create the network
#' @param col1 Color of the nodes that will represent the first group represented in the network. The default value is aquamarine
#' @param col2 Color of the nodes that will represent the second group represented in the network The default value is red
#' @param edge_col Color of the edges that will connect the nodes in the network. The default value is gray
#' @param remove_label If is required to omit the labels for some of the groups, this argument receives the column name informed  the g1 or g2 arguments. The default value is NULL
#' @param node_size A vector with the node sizes to represent g1 and g2. The defaul values are 15 and 40, respectively
#' @param font_size The size of the font of the labels of each node (The default value is 45)
#' @param edge_width The width of the edges connecting the nodes in the network
#' @details This function returns a dynamic network, using visNetwork, representing the connection between two groups. For example, the output from the find_genes_qtls_around_markers() function can be used here to represent the connections between markers and QTLs. Another option is to combine the data frames with both gene and QTL annotation around markers to reprsent the connections between genes and QTLs.     
#' @return A dynamic network representing the connection between two groups. 
#' @importFrom igraph graph.empty
#' @importFrom igraph graph.data.frame
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom visNetwork toVisNetworkData
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visOptions
#' @importFrom visNetwork visNodes
#' @importFrom visNetwork visEdges
#' @importFrom visNetwork visIgraphLayout
#' @importFrom visNetwork visPhysics
#' @name NetVis
#' @export

NetVis<-function(data,g1,g2,col1="aquamarine", col2="red",edge_col="gray",remove_label=NULL,node_size=c(15,40),font_size=45, edge_width=1){
  
  data<-data[!is.na(data[,g1]),]
  
  data<-data[!is.na(data[,g2]),]

  g <- igraph::graph.empty()
    
  g <- igraph::graph.data.frame(data[,c(g1,g2)], directed = FALSE)
    
  igraph::V(g)$color <- ifelse(names(V(g))%in%data[,g1], col1,col2)
    
  igraph::E(g)$weight<-edge_width

  visGraph <- toVisNetworkData(g)
  
  edges<-visGraph$edges
  
  edges$color<-edge_col
  
  nodes<-visGraph$nodes
  
  if(!is.null(remove_label)){
  nodes[which(nodes$label %in% data[,remove_label]),"label"] <- ""
  }
  

  nodes$size<-ifelse(nodes$label%in%data[,g1],node_size[1],node_size[2])
  
  visNetwork(nodes = nodes, edges = edges) %>%
    visOptions(highlightNearest = T, nodesIdSelection = TRUE, clickToUse =T) %>%
    visNodes(font = list(size = font_size)) %>%
    visEdges (width = edge_width) %>%
    visIgraphLayout() %>%
    visPhysics(solver = "forceAtlas2Based",
               forceAtlas2Based = list(gravitationalConstant = -100))
  
}