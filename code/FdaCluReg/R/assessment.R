## assessment standard and information utilization function
ari=function(comm1,comm2){return(igraph::compare(comm1,comm2,method = "adjusted.rand"))}
ri=function(comm1,comm2){return(igraph::compare(comm1,comm2,method = "rand"))}
nmi=function(comm1,comm2){return(igraph::compare(comm1,comm2,method = "nmi"))}
#find mode
find_mode=function(x){return(names(sort(table(x),decreasing = T))[1])}#mode function
