library(ggplot2)
library(factoextra)
library(NbClust)
library(cluster)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)

#for running from IDE
args <- c("/Users/adm515/Desktop/v46.2/v46.2_pca_sorted.csv", TRUE)

#when calling, make sure that group_data is a df of samples within a particular group_ID
calc_gap_stat_opt_clust <- function(group_data) {
  max_clus <- length(group_data[,1]) - 1
  if (max_clus == 0) {
    max_clus <- 1
  }
  gap_stat <- clusGap(group_data, FUN = kmeans, nstart = 25, K.max = max_clus, B = 100)
  opt_clus <- maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"])
  return(opt_clus)
}

run_kmeans <- function(group_data, num_clus) {
  km.curr <- kmeans(group_data, centers = num_clus, nstart = 25)
  return(km.curr)
}

detect_outliers <- function(group_data, kmeans_obj) {
  cluster_centers <- kmeans_obj$centers[kmeans_obj$cluster, ]
  cluster_dists <- sqrt(rowSums(group_data - cluster_centers)^2)
  
  cluster_means <- ave(cluster_dists, kmeans_obj$cluster, FUN = function(x) mean(x, na.rm = TRUE))
  cluster_SDs <- ave(cluster_dists, kmeans_obj$cluster, FUN = function(x) sd(x, na.rm = TRUE))
  
  cluster_z_score <- (cluster_dists - cluster_means)/cluster_SDs
  
  
  cluster_outlier_2 <- ifelse(abs(cluster_z_score) > 2, T, F)
  cluster_outlier_1 <- ifelse(abs(cluster_z_score) > 1, T, F)
  cluster_outlier_0.75 <- ifelse(abs(cluster_z_score) > 0.75, T, F)
  cluster_outlier_0.5 <- ifelse(abs(cluster_z_score) > 0.5, T, F)
  
  return(cbind.data.frame(kmeans_obj$cluster, cluster_z_score, cluster_outlier_2, cluster_outlier_1, cluster_outlier_0.75, cluster_outlier_0.5))
}

if (length(args) < 1) {
  stop("USAGE: [PCA_source_table.csv] [Invert? TRUE/FALSE]")
}

pca_df <- read.csv(args[1], header=TRUE)
row.names(pca_df) <- pca_df[, 1]

#figure out how to pass in bool from CLI
if (args[2] == "TRUE") {
  pca_df$PC_1 <- pca_df$PC_1 * -1
  pca_df$PC_2 <- pca_df$PC_2 * -1
}

set.seed(2021)

group_ids <- unique(pca_df$Group_ID)

results <- data.frame()

for (ID in group_ids) {
  group_df <- pca_df[pca_df$Group_ID == ID, 2:3]
  group_opt_clust <- 1
  
  km.group <- run_kmeans(group_data = group_df, num_clus = group_opt_clust)
  
  clustering_result_df <- detect_outliers(group_data = group_df, kmeans_obj = km.group)
  
  cbind.data.frame(clustering_result_df, rep(ID,nrow(group_df)))
  results <- rbind.data.frame(results, clustering_result_df)
  if (any(is.na(clustering_result_df$cluster_z_score)) == FALSE) {
    pdf(file = paste("/Users/adm515/Desktop/v46.2/outlier_PDF/",ID,".pdf",sep=""), width = 6, height= 4)
    p <- ggplot(group_df,aes(PC_1, PC_2, color=factor(clustering_result_df$cluster_outlier_0.5))) + labs(color='Outlier?') + xlim(-.1, .1) + ylim(-.1, .1) + geom_point()
    print(p)
    dev.off()
    }
}
