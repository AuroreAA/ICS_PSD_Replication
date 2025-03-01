

# Libraries --------------------------------------------------------------

library(ICSOutlier)
library(ICSClust)
library(geigen)
library(MASS)
library(moments)
library(rrcov)
library(ggplot2)
library(dplyr)


source("code/utilities.R")

colors <- ggthemes::colorblind_pal()(4)[c(3, 2)]
# Data  --------------------------------------------------------------------

data(HTP3)
X <- as.matrix(HTP3)
outliers <- c(32)
clusters <- rep("normal", nrow(X))
clusters[outliers] <- "outlier"
boxplot(X, horizontal = TRUE)
dim(X)

#X <- scale(X, center = TRUE, scale = FALSE)

# dimension ok

#working
ics_res <- ICS(X, algorithm = "QR")

df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,1]^2)

# Create the y-axis label based on the combination of scatter matrices and
# the number of selected components

k_details <- 1
y_axis <- bquote("ICSQR - ICSD"[ ~ k == .(k_details) ] ^2 )


plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top")
plot_ics


# save plot to pdf
file_plot = "figures/HTP3/data_%s.%s"
pdf(sprintf(file_plot,  "COVCOV4_QR", "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()



# Reduction dimension ------------------------------------------------------
## SVD ---------------------------------------------------------------------
SVD <- svd(X)
val <- SVD$d
val
rank <- sum(val > max(sqrt(.Machine$double.eps) * val[1]))
rank
# rank <- sum(val > max(dim(X)) * .Machine$double.eps * val[1])
# rank



### Rank 23 -----
loadings <- SVD$v[,1:rank]
X_red <- as.matrix(X) %*% loadings
dim(X_red)

# working
ics_res <- ICS(X_red, algorithm = "standard")

df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = rowSums(ics_res$scores[,1:2]^2))
k=1:2
k_details <- paste(length(k), ifelse(length(k)==1, k,
                                    paste(range(k), collapse = ":")),
                   sep = ", ")
y_axis <- bquote("DR rank = 23 - ICSD"[ ~ k == .(k_details) ] ^2 )

plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top")


# save plot to pdf
file_plot = "figures/HTP3/data_%s_%d.%s"
pdf(sprintf(file_plot,  "COVCOV4_DR", rank, "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()

ics_res <- ICS(X_red, S1= ICS_cov4, S2 = ICS_cov, algorithm = "standard")
df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,23]^2)
k_details <- 23
y_axis <- bquote("DR rank = 23 - ICSD"[ ~ k == .(k_details) ] ^2 )
plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top")


# save plot to pdf
file_plot = "figures/HTP3/data_%s_%d.%s"
pdf(sprintf(file_plot,  "COVCOV4_DR_sym", rank, "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()


try(ics_res <- ICS(X_red, algorithm = "standard", 
                   S1 = ICS_mcd_raw, S2 = ICS_cov))


ics_res <- ICS(X_red, algorithm = "whiten", 
                   S1 = ICS_cov, S2 = ICS_mcd_raw)

df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,23]^2)
k_details <- 23
y_axis <- bquote("DR rank = 23 - ICSD"[ ~ k == .(k_details) ] ^2 )
plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top")


# save plot to pdf
file_plot = "figures/HTP3/data_%s_%d.%s"
pdf(sprintf(file_plot,  "COVMCD_DR", rank, "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()

### Rank var -----
cumsum(val^2)/sum(val^2)*100
threshold <- 99
rank <- which(cumsum(val^2)/sum(val^2)*100 > threshold)[1]
rank
loadings <- SVD$v[,1:rank]
X_red <- as.matrix(X) %*% loadings
dim(X_red)

# working
ics_res <- ICS(X_red, algorithm = "standard")


df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,1]^2)
k_details <- 1
y_axis <- bquote("DR rank = 3 - ICSD"[ ~ k == .(k_details) ] ^2 )
plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top")


# save plot to pdf
file_plot = "figures/HTP3/data_%s_%d.%s"
pdf(sprintf(file_plot,  "COVCOV4_DR", rank, "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()


ics_res <- ICS(X_red, algorithm = "standard", 
                   S1 = ICS_mcd_raw, S2 = ICS_cov)

df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,1]^2)

plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top")


# save plot to pdf
file_plot = "figures/HTP3/data_%s_%d.%s"
pdf(sprintf(file_plot,  "MCDCOV_DR", rank, "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()




# Regularized ---------------------------------------------------------------------

## COV-COV4 ----------------------------------------------------------------
V2 <- cov4mean_ginv(X)
# solve(V2)
# ics_res <- ICS(X, S1 = cov, S2 = cov4mean_ginv,
#                S2_args = list(tolerance = sqrt(.Machine$double.eps))) 
# ics_res <- ICS(X, S1 = V2, S2 = cov(X)) 


## RMCD --------------------------------------------------------------------
RMCD <- CovMrcd(X)
ics_res <- ICS(X, S1 = RMCD@cov, S2 = cov(X), algorithm = "standard")

df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,1]^2)
y_axis <- bquote("MRCD - ICSD"[ ~ k == .(k_details) ] ^2 )
plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top")
# save plot to pdf
file_plot = "figures/HTP3/data_%s.%s"
pdf(sprintf(file_plot,  "MRCD", "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()



k_details <- 33
y_axis <- bquote("ICSD"[ ~ k == .(k_details) ] ^2 )
try(ics_res <- ICS(X, S1 = cov(X), S2 = RMCD@cov, algorithm = "standard"))




# Pseudo inverse - Ginv ----------------------------------------------------------

try(ics_res <- ICS_ginv(X, S1 = ICS_cov, S2 = ICS_cov4, 
                     algorithm = "whiten",
                     tolerance = .Machine$double.eps))



try(ics_res <- ICS_ginv(X, S1 = ICS_cov4, S2 = ICS_cov, 
                    algorithm = "whiten",
                    tolerance = .Machine$double.eps))
# on a des résultats mais c'est juste pour des valeurs propres toutes petites.

try(ics_res <- ICS_ginv(X, S1 = ICS_mcd_raw, S2 = ICS_cov, 
                    algorithm = "whiten",
                    tolerance = .Machine$double.eps))



# GSVD --------------------------------------------------------------------
ics_res <- ICS_generalized(X, S1 = cov_K, S2 = cov4mean_ginv_K, 
                           S2_args = list(tolerance = .Machine$double.eps),
                algorithm = c("GSVD"), center = FALSE, 
                fix_signs = c("scores", "W"), na.action = na.fail) 
k_details <- 1
y_axis <- bquote("GSVD - ICSD"[ ~ k == .(k_details) ] ^2 )

df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,1]^2)

plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top")


# save plot to pdf
file_plot = "figures/HTP3/data_%s.%s"
pdf(sprintf(file_plot,  "COVCOV4_GSVD", "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()


ics_res <- ICS_generalized(X, S1 = cov4mean_ginv_K, S2 = cov_K, 
                           algorithm = c("GSVD"), center = FALSE, 
                           fix_signs = c("scores", "W"), na.action = na.fail) 

df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,33]^2)

k_details <- 33
y_axis <- bquote("GSVD - ICSD"[ ~ k == .(k_details) ] ^2 )

plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top")


# save plot to pdf
file_plot = "figures/HTP3/data_%s.%s"
pdf(sprintf(file_plot,  "COVCOV4_GSVD_sym", "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()


# GEP ---------------------------------------------------------------------
# ics_res <- ICS_generalized(X, S1 = ICS_cov, S2 = cov4mean_ginv, 
#                            algorithm = c("GEP"), center = FALSE, 
#                            fix_signs = c("scores", "W"), na.action = na.fail) 
# 
# select_plot(ics_res)
# component_plot(ics_res, clusters = clusters)
# gen_kurtosis(ics_res)
# ics_res$gen_kurtosis_all
