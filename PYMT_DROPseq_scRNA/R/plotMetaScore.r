# requires Seurat object epidat
plotMetaScore<-function(genelist){
    require(Seurat) # for playing with object
    require(ggplot2) # just to make sure can plot
    require(matlab) # for the jet color palette
    require(dplyr) # for the gene contribution transformation
    message(table(genelist %in% rownames(epidat@data))[2]," out of ",length(genelist)," genes entered were used to generate score\n") # just shows how many are actually contributing to the score
    hm<-rowSums(epidat@scale.data[rownames(epidat@scale.data) %in% genelist,]) #
    message("top contributing genes (by percentage) contributing to signature")
    print((hm*100/sum(hm)) %>% sort(decreasing=T) %>% signif(2) %>% head(10) )
    print(ggplot(data.frame(cbind(epidat@dr$tsne@cell.embeddings),metascore=colSums(epidat@scale.data[rownames(epidat@scale.data) %in% genelist,])),aes(x=tSNE_1,y=tSNE_2,color=metascore)) + geom_point(alpha=.5) + scale_colour_gradientn(colours = jet.colors(7))) #generates the object on the fly and plots it
}