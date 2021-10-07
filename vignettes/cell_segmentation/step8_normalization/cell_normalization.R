## filter data transformation
rm(list = ls())
library(ggplot2)
library(reshape2)
library(fsbrain)

## The function subtracts all values with the the intensity value with maximum density (assuming data is log-tranformed values) 
density_fun <- function(x){
  den <- density(x,n=2^16)
  y_peak <- which.max(den$y)
  x_peak <- den$x[y_peak]
  print(x_peak)
  x <- x - x_peak
}

home_dir <- "/Users/wuminghui/MYDATA/Data0803/"
setwd(home_dir)
output_dir<- "normalized/"
output_filename <- "data_norm"
dir.create(output_dir)

## filter the cell which the area are less than or equal to 25
mydata <- readRDS('All_comp1031.rds')
mydata <- subset(mydata, clean_df$Cell_Area > 25)

## Log transformation 
log_data <- log(mydata[, 3:75]+1, base = 2)
log_df <- cbind(mydata[, 1:2], mydata[,76:89], log_data)
## subtract the peak of each image ##
images <- unique(mydata$Image_name)
norm_df <- c()
# perform normalization on each image individually
for (imag_i in images) {
  print(imag_i)
  imag_df <- subset(log_df,log_df$Image_name == imag_i)
  norm_intens <- apply(imag_df[, 17:89], 2, function(x){density_fun(x)})
  norm_data <- cbind(imag_df[,1:16],norm_intens)
  norm_df <- rbind(norm_df, norm_data)
}
# set all negative values to zero
norm_intens[norm_intens < 0] <- 0
saveRDS(norm_df, paste0(output_dir, output_filename,".rds"))

## plot data distribution after normalization  ##
clip_df <-norm_df
great_scale <- melt(data = clip_df,
                    id.vars = colnames(clip_df[,-c(54:87)]),
                    variable.name = "target", value.name = "intensity")
ggplot(great_scale, aes(x = intensity,  color = Image_name)) +
  geom_density(alpha = 0.4,size=0.5, show.legend = F) +
  facet_wrap(~ target, nrow = 4, scales = "free") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1))
ggsave(paste0(output_dir, output_filename,".png"), width = 900, height = 500, units = 'mm')
# ==== END ==== #

