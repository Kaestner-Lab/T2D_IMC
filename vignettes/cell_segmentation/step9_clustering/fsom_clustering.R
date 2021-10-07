rm(list = ls())
library(ggplot2)
library(reshape2)
library(FlowSOM)
home_dir <- "/Users/wuminghui/MYDATA/Data0803/"
setwd(home_dir)
marker_dir <- "Cluster_marker.csv"
input_dir <- "normalized/data_norm.rds"
output_dir<- "clustering/"
dir.create(output_dir)
output_filename <- "clustered"

mydata <- readRDS(input_dir)
marker_df <- read.csv(marker_dir)
marker_reclus <- marker_df$marker[1:27]
mydata_dm <- data.matrix(mydata)
flowSOM.res <- ReadInput(mydata_dm, transform = F, scale = F)
# building the 15x15 (225 clusters) self-organizing map
fsom <- BuildSOM(flowSOM.res, colsToUse = marker_reclus, xdim = 15, ydim = 15, rlen = 10)
nmc <- 40
# combined into 40 groups through consensus clustering 
metaClustering <-  metaClustering_consensus(fsom$map$codes, k=nmc)
metaClustering_perCell <- GetMetaclusters(fsom, metaClustering)
metafsom_df <- cbind(mydata, metaClustering_perCell)
saveRDS(metafsom_df, paste0(output_dir,output_filename,".rds"))

# keep a record of number of cells in each cluster for each image
fre <- table(metafsom_df$Image_name, metafsom_df$metaClustering_perCell)
sum <- apply(fre,2,function(x){sum(x)})
fre_df <- rbind(fre,sum)
write.csv(fre_df, paste0(output_dir,output_filename, "frequ.csv"))
## ====== clust end ======= ##

# plot boxplot showing distribution of mean-intensity per marker for each cluster
# this plot is used to guide the cluster annotation 
## plot the boxplot 
marker <- c("HLA.ABC","C.peptide","Nestin","Glucagon","pan.Keratin","CD11b","CD44","PDX.1",
            "CD45","CD56","beta.Actin","CD4","NKX6.1","CD68","Somatostatin","CD20","CD8",
            "CD99","CA2","NFkb","GnzB","Ki67","CD57", "CD31","CD14", "Foxp3","p16","CD3",
            "pS6","CD45RO","HLA.DR","PP","GHRL")
data_FlowSOM <- data.frame(metafsom_df)
melt_data <- melt(data_FlowSOM, measure.vars = marker,
                  variable.name = "target", value.name = "value")
melt_data$metaClustering_perCell <- as.factor(melt_data$metaClustering_perCell)
## plot the boxplot of each marker 
ggplot(melt_data, aes( x = metaClustering_perCell, y = value, group = metaClustering_perCell)) +
  geom_boxplot( aes(color=metaClustering_perCell), alpha = 0.3, size = 0.3, show.legend = F)+
  facet_wrap(~ target, ncol = 4, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
save_name <- paste0(output_dir, output_filename, "boxreclumetaMar.png")
ggsave(save_name, width = 1000, height = 1000, units = 'mm')
## boxplot of each cluster
ggplot(melt_data, aes( x = target, y = value, group = target )) +
  geom_boxplot( aes(color=target), alpha = 0.3, size = 0.3, show.legend = F)+
  facet_wrap(~ metaClustering_perCell, ncol = 4, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
save_name <- paste0(output_dir, output_filename, "boxreclumetaClu.png")
ggsave(save_name, width = 1000, height = 1000, units = 'mm')
## ========== ##

# Another validation method is to reconstruct the image by projecting the cells based on their x and y coordinate and see the spatial distribution of cells (and potentially compare with the raw images). 
# select cells in an image
cluster_6009h2 <- subset(metafsom_df, metafsom_df$Image_name=="NPOD6444_Tail_ROI2")
# subset cells in a specific cluster
cluster_sel_6009h2 <- subset(cluster_6009h2, cluster_6009h2$metaClustering_perCell %in% c(4))
png("cell_spatial_dist.png")
plot(1:1100, 1:1100, ylim = rev(range(1:1100)),type = "n",
     main = paste(unique(cluster_sel_6009h2$metaClustering_perCell),
                  unique(cluster_6009h2$Image_name), nrow(cluster_sel_6009h2)))
text(x = cluster_sel_6009h2$Location_Center_X, y = cluster_sel_6009h2$Location_Center_Y, 
     label = "o",col = "red", cex = 0.8)
dev.off()
