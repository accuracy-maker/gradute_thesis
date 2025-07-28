# Load Dataset and Maxpooling
# Author: Haitao Gao
# Date: 11st July 2025

library(parallel)
library(evd) 
################ Load Dataset #####################
# The dataset contains 4 parts
# 1. coord.grid: the spatial coordinates of the data
# 2. data: the index and rainfall values (list of list)
# 3. data.fit.sum: the extremes which are above 99.95% of mean value
# 4. data.fit_max: the extreme whcih are above 99.95% of max value

load("../data/Florida_new2.RData")


############## Convert extremes into three sets ##############
# 1. set 1: intersection set of data.fit.sum and data.fit.max; label 0
# 2. set 2: data.fit.sum - set 1; label 1
# 3. set 3: data.fit.max - set 1; label 2

create_fitting_data <- function(data, threshold = 0.9995){
  # transform to unit-pareto
  data_pareto <- mclapply(data,function(x){list(x[[1]],qgpd(x[[2]],1,1,1))})
  
  # Compute sum and max vectors
  data_mean <- unlist(lapply(data_pareto, function(x) mean(x[[2]])))
  data_max <- unlist(lapply(data_pareto, function(x) max(x[[2]])))
  
  # threshold filtering
  data_fit_max <- data_pareto[data_max > quantile(data_max, threshold)]
  data_fit_mean <- data_pareto[data_mean > quantile(data_mean, threshold)]
  
  # extract the index sets
  get_idx_str <- function(x) paste(sort(x[[1]]), collapse = "-")
  mean_keys <- sapply(data_fit_mean, get_idx_str)
  max_keys  <- sapply(data_fit_max,  get_idx_str)
  
  # find the intersections and differences
  intersect_keys <- intersect(mean_keys, max_keys)
  only_mean_keys <- setdiff(mean_keys, intersect_keys)
  only_max_keys  <- setdiff(max_keys, intersect_keys)
  
  # label data
  label_data <- function(data_list, keys, label){
    result <- list()
    for(i in seq_along(data_list)){
      if(get_idx_str(data_list[[i]]) %in% keys){
        result[[length(result)+1]] <- list(data_list[[i]][[1]],
                                           data_list[[i]][[2]],
                                           label)
      }
    }
    return (result)
  }
  
  data_0 <- label_data(data_fit_mean, intersect_keys, 0)
  data_1 <- label_data(data_fit_mean, only_mean_keys, 1)
  data_2 <- label_data(data_fit_max,  only_max_keys,  2)
  
  # combine the dataset
  combined_data <- c(data_0, data_1, data_2)
  
  # print the set size
  cat("Set sizes:\n")
  cat("Intersection (label 0):", length(data_0), "\n")
  cat("Only mean (label 1):   ", length(data_1), "\n")
  cat("Only max (label 2):    ", length(data_2), "\n")
  
  return(combined_data)
}

# get the new dataset
# new_dataset <- create_fitting_data(data_down_raw, threshold = 0.999)

