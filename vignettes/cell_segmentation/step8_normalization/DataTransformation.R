## After compensation,this script is used for the data transformation before clustering
## 1.log-transformation: all data add 1 and do base 2 logarithm.
## Add 1 is good for visualization, I tried base exp, but base 2 is better to amply the dim marker.
## 2.subtracting the lowest mode: subtract the max-density intensity.  
## For my markers, the distribution shape of same markers in different images look like similar.
## All of them(except the pS6 and CD11b, they are no-special stating in some tissues, I don't take them in
## statistical analysis, but I also transform them.) have the obvious one peak, then I assume it as background.
## I refer the zero-centered normalization, I subtract the the max-density intensity. 
## Other: For my data, some channels(such as PP and Ghrelin) have 0 positive cell in some donors, I tried 99th percentile, 
## z-score, 0-1 normalization, they caused many false positive cells after FlowSOM algorithm clustering.
rm(list = ls())
library(ggplot2)
library(reshape2)
## The function is subtracting the max-density intensity
density_fun <- function(x){
  cut <- subset(x, x > 0)
  y_peak <- which.max(density(cut)$y)
  x_peak <- density(cut)$x[y_peak]
  print(x_peak)
  x <- x - x_peak
}
home_dir <- "/Users/wuminghui/Documents/Minghui_Method/"
setwd(home_dir)
input_dir <- "Testdata.rds"
output_dir<- "data_ad1base2/"
dir.create(output_dir)
mydata <- readRDS(input_dir)
## Log transformation 
log_data <- log(mydata[, 3:38]+1, base = 2)
log_df <- cbind(mydata[, 1:2], log_data)
saveRDS(log_df, paste0(output_dir, "All_cpad1log2.rds"))
## plot log-transformation data distribution ##
great_scale <- melt(data = log_df,
                    id.vars = colnames(log_df[,-c(3:38)]),
                    variable.name = "target",
                    value.name = "intensity")
ggplot(great_scale, aes(x = intensity,  color = Image_name)) +
  geom_density(alpha = 0.4,size=0.5, show.legend = F)+
  facet_wrap(~ target, nrow = 4, scales = "free")+ theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1))
ggsave(paste0(output_dir, "All_cpad1log2.png"), width = 900, height = 500, units = 'mm')

## subtract the peak of each image ##
images <- unique(mydata$Image_name)
norm_df <- c()
for (imag_i in images) {
  print(imag_i)
  imag_df <- subset(log_df,log_df$Image_name == imag_i)
  norm_intens <- apply(imag_df[, 3:38], 2, function(x){density_fun(x)})
  norm_data <- cbind(imag_df[,1:2], norm_intens)
  norm_df <- rbind(norm_df, norm_data)
}
saveRDS(norm_df, paste0(output_dir, "All_cpad1log2den1014.rds"))

##  data distribution  ##
great_scale <- melt(data = norm_df,
                    id.vars = colnames(norm_df[,-c(3:38)]),
                    variable.name = "target",
                    value.name = "intensity")
ggplot(great_scale, aes(x = intensity,  color = Image_name)) +
  geom_density(alpha = 0.4,size=0.5, show.legend = F)+
  facet_wrap(~ target, nrow = 4, scales = "free")+ theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1))
ggsave(paste0(output_dir, "All_cpad1log2den1014.png"), width = 900, height = 500, units = 'mm')
# ==== END ==== #







