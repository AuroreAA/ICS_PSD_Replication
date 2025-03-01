

# Libraries --------------------------------------------------------------

library(ICSOutlier)
library(ICSClust)
library(geigen)
library(MASS)
library(moments)
library(rrcov)
library(ggthemes)


# Functions
source("code/utilities.R")
# To generate the mixture of gaussian distributions
mixture_sim2 <- function (pct_clusters = c(0.5, 0.5), n = 500, p = 10, 
                          mean_lst = list(), sd_lst = list()){
  
  
  if (sum(pct_clusters) != 1) {
    stop("the sum of the groups is not equal to 1!")
  }
  n_groups = floor(n * pct_clusters)
  if (sum(n_groups) != n) {
    n_groups[length(pct_clusters)] <- rev(n - cumsum(n_groups)[1:(length(pct_clusters) - 
                                                                    1)])[1]
  }
  if (sum(n_groups) != n) {
    warning(paste("the total number of observation is not equal to n", 
                  paste(round(pct_clusters, 2), collapse = " - ")))
  }
  
  X_list <- lapply(1:length(pct_clusters), function(i) {
    n_gr <- n_groups[i]
    if (n > 0) {
      data.frame(mvtnorm::rmvnorm(n = n_gr, mean = mean_lst[[i]], 
                                  sigma = sd_lst[[i]]))
    }
  })
  data.frame(cluster = rep(paste0("Group", 1:length(pct_clusters)), 
                           n_groups), do.call(rbind, (X_list)))
}



colors <- ggthemes::colorblind_pal()(4)[c(3, 2)]
set.seed(20240725)



# Data: 2 groups collinear  --------------------------------------------------------------------
p <- 4
n <- 1000
mean_lst <- list(rep(0,p), rep(0,p))

W1 = diag(rep(1,p/2))
W2 = diag(c(2,rep(0,(p/2)-1)))

Sigma1 = rbind(cbind(W1,matrix(0,nrow(W1),ncol(W2))), matrix(0,p-ncol(W1),p))
Sigma2 = rbind(cbind(W1,matrix(0,nrow(W1),ncol(W2))), cbind(matrix(0,nrow(W2),p-ncol(W2)), W2))

sd_lst <- list(Sigma1, Sigma2)
df <- mixture_sim2(pct_clusters = c(0.98, 0.02), n = n, p = p, 
                   mean_lst = mean_lst, sd_lst = sd_lst)

clusters <- df[,1]
X_ini <- as.matrix(df[,-1])

plot_ini <- component_plot(data.frame(X_ini), clusters = clusters, 
                           colors = colors)
plot_ini

# Transformation
G <- toeplitz(p:1)/p

X <- X_ini %*% G
plot_ini <- component_plot(data.frame(X), clusters = clusters, 
                           colors = colors)
plot_ini


file_plot = "figures/illustrative_example/data_G.%s"
pdf(sprintf(file_plot, "pdf"), width = 6.5, height = 5.75)
print(plot_ini)
dev.off()


# V2 is the theoretically empirical scatter matrix
eps = 0.02
V2 = rbind(cbind(W1,matrix(0,nrow(W1),ncol(W2))), cbind(matrix(0,nrow(W2),p-ncol(W2)), eps*W2))
# V1 is the theoretically robust scatter matrix
V1 = Sigma1
X.V1 = sqrt(V1)
X.V2 = sqrt(V2)

V1 = G %*%V1%*%t(G)
V2 = G %*%V2%*%t(G)

X.V1 = X.V1%*%t(G)
X.V2 = X.V2%*%t(G)


# not working
try(ics_res <- ICS(X, S1 = V1, S2 = V2,  algorithm = "QR"))
try(ics_res <- ICS(X, S1 = V1, S2 = V2, algorithm = "standard"))





# Pseudo inverse - Ginv ----------------------------------------------------------
### COV-COV4 -----
ics_res <- ICS_ginv(X, S1 = V1, S2 = V2, 
                    algorithm = "standard",
                    tolerance = .Machine$double.eps)
val <- sprintf(fmt = "%0.02e", ics_res$gen_kurtosis)
val
# colnames(ics_res$scores) <- paste("GINV",  paste0("~rho", "==", val),
#                                   colnames(ics_res$scores), sep = "-")
colnames(ics_res$scores) <- paste("GINV", colnames(ics_res$scores), sep = "-")
ics_res$gen_kurtosis
plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )+
  scale_x_continuous(
    # labels =   ~ sprintf(fmt = "%0.01e", .),
    guide  = guide_axis(check.overlap = T)) + 
  scale_y_continuous(
    # labels =   ~ sprintf(fmt = "%0.01e", .),
    guide  = guide_axis(check.overlap = T))


# save plot to pdf
file_plot = "figures/illustrative_example/data_G_%s.%s"
pdf(sprintf(file_plot,"GINV_COVCOV4", "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()




# GSVD --------------------------------------------------------------------

ICS_X_V <- function(x, V){
  
  location <- rep(0, ncol(V))
  out <- list(location = location, scatter = V, label = "")
  class(out) <- "ICS_scatter"
  out
}

ics_res <- ICS_generalized(X, S1 = ICS_X_V , S2 = ICS_X_V,
                           S1_args = list(V = X.V1), S2_args = list(V = X.V2),
                           algorithm = c("GSVD"), center = FALSE,
                           na.action = na.fail) 
colnames(ics_res$scores) <- paste("GSVD", colnames(ics_res$scores), sep = "-")
ics_res$gen_kurtosis_all
plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )


# save plot to pdf
file_plot = "figures/illustrative_example/data_G_%s.%s"
pdf(sprintf(file_plot,"GSVD_COVCOV4", "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()

ics_res <- ICS_generalized(X, S1 = cov4mean_ginv_K, S2 = cov_K, 
                           algorithm = c("GSVD"), center = FALSE, 
                           fix_signs = c("scores", "W"), na.action = na.fail) 
colnames(ics_res$scores) <- paste("GSVD", colnames(ics_res$scores), sep = "-")
ics_res$gen_kurtosis_all
plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )
# save plot to pdf
file_plot = "figures/illustrative_example/data_G_%s.%s"
pdf(sprintf(file_plot,"GSVD_COVCOV4_sym", "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()






# GEP ---------------------------------------------------------------------
# ics_res <- ICS_generalized(X, S1 = ICS_cov, S2 = cov4mean_ginv, 
#                            algorithm = c("GEP"), center = FALSE, 
#                            fix_signs = c("scores", "W"), na.action = na.fail) 
# 
# select_plot(ics_res)
# component_plot(ics_res, clusters = clusters, colors = colors)
# gen_kurtosis(ics_res)
# ics_res$gen_kurtosis_all