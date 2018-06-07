function plotweb(A,s,tl)
    R"""
    set.seed(1);
    library(igraph)
    library(RColorBrewer)
    pal <- brewer.pal(3,"Set1")
    fw_g <- graph.adjacency($(A'));
    trophic <- as.numeric($tl);
    coords <- cbind(runif(vcount(fw_g)),trophic);
    plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='#6495ED30',vertex.label=NA,vertex.frame.color='black', vertex.color=grey(($s+1)/2))
    """
end

