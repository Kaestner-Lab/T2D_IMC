# Removing streaks or hot pixels in IMC images 
# Using the first two noise/background channel for points detection
# Then for rest of the channels, check if these points are greater than the mean of local neighborhood(5x5), if yes, remove.

# To run the script, please see the end of the file. 

# `tiff` the image data.
# n=5 is the size of the mask.  The value has to be odd so that there is a
# center row to drop.  The flashes we have seen are in the horizontal
# orientation so we drop a row (not column.)

# `apply_to` sets an intensity threshold.  Pixels over this threshold are
# considered to be in a flash.  If the value is negative, it is treated as a
# 'negative' quantile and is used to set the threshold value from the tiff data.
# For example, the default value of -0.1 means that the threshold will be set at
# the 90% quantile of pixel brightness in the image.

# `min_rel_change` this is the fold change needed to trigger replacement.  The
# purpose is to allow `apply_to` to be set liberally, but to only modify pixels
# that are 'clearly' flashes.


# ---
# free_mask_x and free_mask_y are relative column and row values from define the
# mask relative to the center pixel.  For example for n=5 the rows will have
# values in -2, -1, 1, and 2 (once the center row is dropped) and columns will
# be -2,-1,0,1, and 2.
# to speed things up we define a list of the flash pixel coordinates, oi.  We do
# not need to apply the filter everywhere so this speeds up the calculation
# by skipping useless work.

# as we iterate through the flash pixel coordinates we 'bind' the mask to the
# pixel's location and store in bound_mask.
# using the bound mask we calculate the local mean.

# Do we update the pixel?  If it was zero we do (this is an extra non-intuitive
# goal)  If the mean is 'much' less than the pixel, then we replace it in place.
# once we have iterated over the pixels, we return the modified tiff.

tiff.deflash_pt <- function ( tiff, n=5, apply_to=-0.02, min_rel_change=2.0 ) {
  if ( n %% 2 == 0) n <- n + 1
  half_n <- (n-1)/2
  offset_i <- 1:n - half_n - 1
  free_mask_x <- rep( offset_i, n-1 )
  free_mask_y <- rep( offset_i, each=n )
  free_mask_y <- free_mask_y[ free_mask_y != 0 ]
  # free_mask   <- cbind( free_mask_y, free_mask_x ) # backwards to match row,col indexing
  inv_min_rel_change <- 1/min_rel_change
  rv_tiff     <- tiff
  # convert quantile cutoff to absolute
  if (apply_to < 0) {
    apply_to <- as.numeric(quantile(as.numeric(tiff), 1+apply_to))
  }
  if(apply_to==0){
    apply_to<-1
  }
  # get a data.frame of the pixel coordinates for the flash pixels.
  oi   <- (data.frame( which( tiff >= apply_to, arr.ind=TRUE ))
           %>% dplyr::filter( row >= half_n, col >= half_n )
           %>% dplyr::filter( row <= nrow(tiff) - half_n, col <= ncol(tiff) - half_n))
  oi_n <- nrow(oi)
  # work through the flash pixels.
  print(paste0("going through ",oi_n))
  n_changed <- 0 
  # points to be changed
  oi_tc <- c()
  if(oi_n>=1){
    for (i in 1:oi_n) {
      #print(i)
      bound_mask <- cbind( free_mask_y + oi$row[i], free_mask_x + oi$col[i] ) 
      old_value  <- tiff[oi$row[i],oi$col[i]]
      new_value  <- mean( tiff[bound_mask] ) + 1
      if ( old_value == 0 ) {
        update_px <- TRUE
      } else if ( new_value/old_value < inv_min_rel_change ) {
        update_px <- TRUE
        n_changed<-n_changed+1
        oi_tc<-c(oi_tc,i)
      } else {
        update_px <- FALSE
      }
      #if (update_px) rv_tiff[ oi$row[i], oi$col[i] ] <- new_value
    }
    print(length(oi_tc))
    return(oi[oi_tc,])
  }
}

write_multi_tiff<-function(data,filename,bits_per_sample,width,height){
  write_tif(result_tiff_obj[[1]][,,1,],paste0(output_directory,i),
            bits_per_sample = 16,overwrite = TRUE)
  tiff(filename = paste0(output_directory,i),
       width = round(width/1000*0.4,digits = 1), height = round(height/1000*0.4,digits = 1), units = "inch", pointsize = 12,
       compression = "none", bit)
}

tiff_deflash_multistack<- function(tiff, n=5, apply_to=-0.02, min_rel_change=5.0 ){
  # first two noise channels, the minimal change requirement is much less, only 2 fold-change required
  points_1 <- tiff.deflash_pt(tiff[,,,1],n=n, apply_to=apply_to, min_rel_change=2.0 )
  points_2 <- tiff.deflash_pt(tiff[,,,2],n=n, apply_to=apply_to, min_rel_change=2.0 )
  oi_all <- intersect(points_1,points_2)
  print(paste0("Number of common signals detected in channel 1 and 2: ", nrow(oi_all)))
  if ( n %% 2 == 0) n <- n + 1
  half_n <- (n-1)/2
  offset_i <- 1:n - half_n - 1
  free_mask_x <- rep( offset_i, n-1 )
  free_mask_y <- rep( offset_i, each=n )
  free_mask_y <- free_mask_y[ free_mask_y != 0 ]
  inv_min_rel_change <- 1/min_rel_change
  
  oi_n <- nrow(oi_all)
  # work through the flash pixels.
  all_changed <- 0 
  res_tiff <- tiff
  if(oi_n>=1){
    print(paste0("removing streaks: ", oi_n, " points"))
    for(j in 1:dim(tiff)[4]){
      n_changed<-0
      tiff_i<-tiff[,,,j]
      rv_tiff     <- tiff_i
      for (i in 1:oi_n) {
        #print(i)
        bound_mask <- cbind( free_mask_y + oi_all$row[i], free_mask_x + oi_all$col[i] ) 
        old_value  <- tiff_i[oi_all$row[i],oi_all$col[i]]
        new_value  <- mean( tiff_i[bound_mask] )
        if ( old_value == 0 ) {
          update_px <- TRUE
        } else if ( new_value/old_value < inv_min_rel_change ) {
          update_px <- TRUE
          n_changed <- n_changed+1
          all_changed <- all_changed+1
        } else {
          update_px <- FALSE
        }
        if (update_px) rv_tiff[ oi_all$row[i], oi_all$col[i] ] <- new_value
      }
      res_tiff[,,,j]<-round(rv_tiff,digits = 0)
      print(paste0("In channel ", j, ": ",n_changed," points changed."))
    }
    print(paste0("All channels: ",all_changed," points changed."))
    }
    return(list(res_tiff,(all_changed>=5)))
}

deflash_folder_2<- function(input_directory,output_directory,n=5, apply_to=-0.02, min_rel_change=5.0){
  # copy the whole directory into the output folder
  file.copy(input_directory,output_directory ,recursive = T,overwrite = T)
  # get all tiff files that need to be processed
  files <- list.files(output_directory,recursive = T,pattern = ".tiff")
  log<-c()
  for(i in files){
    print(paste0("processing ", i))
    flash_tiff <- read_tif(paste0(output_directory,i))
    
    result_tiff_obj <- tiff_deflash_multistack(flash_tiff, n ,apply_to, min_rel_change)
    write_tif(result_tiff_obj[[1]][,,1,],paste0(output_directory,i),
              bits_per_sample = 16,overwrite = TRUE)
    tiff(filename = "Rplot%03d.tiff",
         width = 480, height = 480, units = "px", pointsize = 12,
         compression = c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"),
         bg = "white", res = NA,  ...,
         type = c("cairo", "Xlib", "quartz"), antialias)
    if(result_tiff_obj[[2]]){
      log<-c(log, i)
    }
  }
  return(log)
}

#==== Run here=====
# input folder should contain a list of .tiff files
input_directory<-"~/Desktop/NPOD6259_Body/"
# output folder location, will mimic the folder structure from input folder
output_directory<-"~/Desktop/tmp/"

log <- deflash_folder_2(input_directory,output_directory)
# print the name of the images that have been changed 
print(paste0("Images with streaks: ",log))
