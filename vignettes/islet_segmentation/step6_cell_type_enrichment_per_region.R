## the coordination less than 1 will be error, set them as 1 will be good
rm(list = ls())
Sys.time()
library(EBImage)
library(stringr)
home_dir <- "/Users/wuminghui/MYDATA/"
setwd(home_dir)
mask_dir <- "Cellprofiler_output/CD99"
cell_dir <- "Data0108/anno/annot_0108.rds"
output_dir <- "Data0108/Imu_islet_dis_0108/"
dir.create(output_dir)
output_filename <- "imuToislet_0108"
mask_all <- list.files(mask_dir,full.names = T)
mydata <- readRDS(cell_dir)
imu <- grep('CD4|CD8|Macrophage|MemoryT',mydata$CellType5)
imu_df <- mydata[imu,]
dis_df <- c()
for (i in 1:length(mask_all)) { 
  mask_i <- readImage(mask_all[i])
  dmap_i <- distmap(x=1- mask_i, metric="euclidean")
  mask_name <- str_match(mask_all[i], "CD99/s*(.*?)s*_CD99")
  img_name <- mask_name[,2]
  img_imu <- subset(imu_df, imu_df$Image_name == img_name)
  DistanceToIslet <- c()
  for (loc_i in 1:length(img_imu$CellID)) {
    dis_i <- dmap_i[img_imu$Location_Center_X[loc_i], img_imu$Location_Center_Y[loc_i]]
    DistanceToIslet <- rbind(DistanceToIslet, dis_i)
  }
  dis_ls <- cbind(img_imu, DistanceToIslet)
  dis_df <- rbind(dis_df, dis_ls)
  print(img_name)
  print(dim(dis_df))
  print(i)
}
Sys.time()
saveRDS(dis_df,paste0(output_dir, output_filename, ".rds"))

