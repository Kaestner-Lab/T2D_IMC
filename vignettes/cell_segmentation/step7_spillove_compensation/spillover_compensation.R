library("CATALYST")
library("dplyr")
library("SingleCellExperiment")
library("data.table")
library("LSD")

comp_datimg <- function(datimg, sm, method='nnls'){
  img_mat = as.matrix(datimg)
  sce <- prepData(img_mat %>%flowCore::flowFrame())
  sce <- compCytof(sce, sm, method = method, overwrite = T)
  return(data.frame(t(counts(sce))))
}

# read in the metal ion vs. gene table
tag_gene_list <- read.csv("imc_pannel.csv")
print(head(tag_gene_list))
sm_table <- read.csv("Isotope_Purity_Matrix.csv")
print(head(sm_table))
# read in file
image_6009 <- read.csv("NPOD6009_Head_ROI1.csv")
print(head(image_6009))
# check resulting dimension
dim(image_6009)
## metal as sm_table rowname
rownames(sm_table) <- sm_table[,1]
# grab the two columns of interest
tag_gene_list<- tag_gene_list[,3:4]
# create a new dataframe (selecting antibody columns)
image_6009_comp <- image_6009[,6:ncol(image_6009)]
rownames(image_6009_comp) <- image_6009[,1]
image_6009 <- image_6009[,6:ncol(image_6009)]
# change column names to their corresponding metal tags
colnames(image_6009_comp) <- tag_gene_list$Metal.Tag2[unlist(lapply(colnames(image_6009_comp), function(x){which(x == tag_gene_list$Target)}))]
#image_6009_comp <- image_6009_comp[,!is.na(colnames(image_6009_comp))]
# compensate 
sm_table[is.na(sm_table)] <- 0
sm_table <- sm_table[,-1]/100
image_6009_comp_2 <- comp_datimg(image_6009_comp,sm_table[colnames(image_6009_comp),colnames(image_6009_comp)])
# reset column names to gene names
colnames(image_6009_comp_2) <- colnames(image_6009)
# save output
write.csv(image_6009_comp_2, "~/Desktop/image_6009_comp.csv")


# heatscatter plot 
heatscatter(image_6009_comp$Nd145, image_6009_comp$Nd146)
heatscatter(image_6009_comp_2$C.peptide, image_6009_comp_2$Nestin)

# check with previous compensation 
data_control_full <- readRDS("~/Desktop/kaestner_lab/results/2020-03-30_cp_downstream_analysis_all/data_control_full.rds")
# compare raw data
summary(data_control_full$raw.data$CD57[1:8000])
summary(image_6009$CD57[1:8000])

# compare compensated data
summary(data_control_full$comp.data$CD57[1:8000])
summary(image_6009_comp$CD57[1:8000])

# View compensation result 
require(ggplot2)
par(mfrow=c(2,2))
# CD14
p<-list(ggplot(mapping = aes(x =1:nrow(image_6009), y= image_6009$C.peptide)) + geom_point(),
#CPEPTIDE
ggplot(mapping = aes(x =1:nrow(image_6009), y= image_6009$Nestin)) + geom_point(),
# after compensation
ggplot(mapping = aes(x =1:nrow(image_6009), y= image_6009_comp$C.peptide)) + geom_point(),
ggplot(mapping = aes(x =1:nrow(image_6009_comp), y= image_6009_comp$Nestin)) + geom_point())
gridExtra::grid.arrange(grobs=p)

