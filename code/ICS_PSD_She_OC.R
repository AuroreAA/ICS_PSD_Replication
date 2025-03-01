

# Libraries --------------------------------------------------------------

library(ICSOutlier)
library(ICSClust)
library(geigen)
library(MASS)
library(moments)
library(rrcov)
library(ggthemes)


set.seed(20240727)
source("code/utilities.R")

colors <- ggthemes::colorblind_pal()(4)[c(3, 2)]

# Data --------------------------------------------------------------------
n = 100
p = 5
r = 3
mu = 0
n_outlier <- 4 # 4, 10, 16
L <- 3.5 # 4.5
# row type outliers
S <- L * matrix(c(rep(1, n_outlier*(p-r)), rep(0,(n-n_outlier)*(p-r))),
                nrow = n, ncol = p-r, byrow = TRUE)

#Forming singular values 
d <- c(100, 60, 2)
d <- c(1000, 400, 200)

# Generate uncorrelated data of low rank
#Forming orthogonal left-singular vectors
U <- matrix(rnorm(n*r), n, r)
decomp <- qr(U)
U <- qr.Q(decomp)

crossprod(U)

#Forming orthogonal right-singular vectors
V <- matrix(rnorm(p*r), p, r)
decomp <- qr(V)
V <- qr.Q(decomp)
crossprod(V)
tcrossprod(V)
Proj_V <- V %*% t(V)


Proj_OC <- diag(rep(1,p)) - Proj_V
#V_OC <- diag(rep(1,p)) - V%*%solve(t(V)%*%V)%*%
decomp <- eigen(Proj_OC)
V_OC <- decomp$vectors[,1:(p-r)]
# check
Proj_OC - V_OC%*%t(V_OC)

#Forming matrix x
Xr <- U %*% diag(d) %*% t(V)
X_OC <- S %*% t(V_OC)
X <- Xr + X_OC 
#X <- Xr + X_OC + E
clusters <- rep(paste0("Group", 1:2), 
                c(n_outlier, n-n_outlier))


plot_ics <- component_plot(data.frame(X), clusters = clusters, 
                           colors = colors)
# save plot to pdf
file_plot = "figures/OC_5/data.%s"
pdf(sprintf(file_plot,  "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()

plot_ics <- component_plot(data.frame(X%*%Proj_OC), clusters = clusters)
# save plot to pdf
file_plot = "figures/OC_5/data_proj_OC.%s"
pdf(sprintf(file_plot,  "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()


plot_ics <- component_plot(data.frame(X_OC), clusters = clusters)
# save plot to pdf
file_plot = "figures/OC_5/data_X_OC.%s"
pdf(sprintf(file_plot,  "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()


# Reduction dimension -----------------------------------------------------
## SVD ---------------------------------------------------------------------

### Rang = 4 -----
SVD <- svd(X)
val <- SVD$d
val
rank <- sum(val > max(sqrt(.Machine$double.eps) * val[1]))
rank
rank <- sum(val > max(dim(X)) * .Machine$double.eps * val[1])
rank

loadings <- SVD$v[,1:rank]
X_red <- as.matrix(X) %*% loadings
dim(X_red)

# working
ics_res <- ICS(X_red)
colnames(ics_res$scores) <- paste("DR", colnames(ics_res$scores), sep = "-")
plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors)

# save plot to pdf
file_plot = "figures/OC_5/data_%s_%d.%s"
pdf(sprintf(file_plot,  "COVCOV4_DR", rank, "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()

try(ics_res <- ICS(X_red, algorithm = "standard", S1 = ICS_mcd_raw, S2 = ICS_cov))
try(ics_res <- ICS(X_red, algorithm = "standard", S1 = ICS_cov, S2 = ICS_mcd_raw))



### Rank var -----
cumsum(val^2)/sum(val^2)*100
threshold <- 95
rank <- which(cumsum(val^2)/sum(val^2)*100 > threshold)[1]
rank
loadings <- SVD$v[,1:rank]
X_red <- as.matrix(X) %*% loadings
dim(X_red)

ics_res <- ICS(X_red)
colnames(ics_res$scores) <- paste("DR", colnames(ics_res$scores), sep = "-")
plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors)

# save plot to pdf
file_plot = "figures/OC_5/data_%s_%d.%s"
pdf(sprintf(file_plot,  "COVCOV4_DR", rank, "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()



ics_res <- ICS(X_red, algorithm = "standard", 
               S1 = ICS_mcd_raw, S2 = ICS_cov)
colnames(ics_res$scores) <- paste("DR", colnames(ics_res$scores), sep = "-")
plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors)

# save plot to pdf
file_plot = "figures/OC_5/data_%s_%d.%s"
pdf(sprintf(file_plot,  "MCDCOV_DR", rank, "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()

# Regularized ---------------------------------------------------------------------
## RMCD --------------------------------------------------------------------
RMCD <- CovMrcd(X)
ics_res <- ICS(X, S1 = RMCD@cov, S2 = cov(X), algorithm = "standard")
colnames(ics_res$scores) <- paste("MRCD", colnames(ics_res$scores), sep = "-")

plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors)+
  scale_x_continuous(
    # labels =   ~ sprintf(fmt = "%0.01e", .),
    guide  = guide_axis(check.overlap = T)) + 
  scale_y_continuous(
    # labels =   ~ sprintf(fmt = "%0.01e", .),
    guide  = guide_axis(check.overlap = T))

# save plot to pdf
file_plot = "figures/OC_5/data_%s.%s"
pdf(sprintf(file_plot,  "MRCDCOV", "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()




# Pseudo inverse - Ginv ----------------------------------------------------------
ics_res <- ICS_ginv(X, S1 = ICS_cov, S2 = ICS_cov4, 
                        algorithm = "whiten",
                        tolerance = .Machine$double.eps)
colnames(ics_res$scores) <- paste("GINV", colnames(ics_res$scores), sep = "-")

plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors)+
  scale_x_continuous(
    # labels =   ~ sprintf(fmt = "%0.01e", .),
    guide  = guide_axis(check.overlap = T)) + 
  scale_y_continuous(
    # labels =   ~ sprintf(fmt = "%0.01e", .),
    guide  = guide_axis(check.overlap = T))

# save plot to pdf
file_plot = "figures/OC_5/data_%s.%s"
pdf(sprintf(file_plot,  "COVCOV4_GINV",  "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()


try(ics_res <- ICS_ginv(X, S1 = ICS_mcd_raw, S2 = ICS_cov, 
                    algorithm = "whiten",
                    tolerance = .Machine$double.eps))


# GSVD --------------------------------------------------------------------
ics_res <- ICS_generalized(X, S1 = cov_K, S2 = cov4mean_ginv_K, 
                           S2_args = list(tolerance = .Machine$double.eps),
                algorithm = c("GSVD"), center = FALSE, 
                fix_signs = c("scores", "W"), na.action = na.fail) 
colnames(ics_res$scores) <- paste("GSVD", colnames(ics_res$scores), sep = "-")

plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors)

# save plot to pdf
file_plot = "figures/OC_5/data_%s.%s"
pdf(sprintf(file_plot,  "COVCOV4_GSVD",  "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()



ics_res <- ICS_generalized(X, S1 = cov4mean_ginv_K, S2 = cov_K, 
                           algorithm = c("GSVD"), center = FALSE, 
                           fix_signs = c("scores", "W"), na.action = na.fail) 
colnames(ics_res$scores) <- paste("GSVD", colnames(ics_res$scores), sep = "-")

plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors)

# save plot to pdf
file_plot = "figures/OC_5/data_%s.%s"
pdf(sprintf(file_plot,  "COVCOV4_GSVD_sym",  "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()




