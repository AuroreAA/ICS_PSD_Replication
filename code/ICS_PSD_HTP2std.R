

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

# Data: 2 groups collinear  --------------------------------------------------------------------


data(HTP2)
X <- as.matrix(HTP2)
outliers <- c(28)
clusters <- rep("normal", nrow(X))
clusters[outliers] <- "outliers"
boxplot(X, horizontal = TRUE)
dim(X)

# standardization of the data
X <- scale(X, center = TRUE, scale = TRUE)

# not working
try(ics_res <- ICS(X, algorithm = "QR"))
try(ics_res <- ICS(X, algorithm = "standard", S1 = ICS_mcd_raw, S2 = ICS_cov))




# Reduction dimension -----------------------------------------------------
## SVD ---------------------------------------------------------------------

### X ----
SVD <- svd(X)
val <- SVD$d
val
rank <- sum(val > max(sqrt(.Machine$double.eps) * val[1]))
rank
rank <- sum(val > max(dim(X)) * .Machine$double.eps * val[1])
rank
cumsum(val^2)/sum(val^2)*100
threshold <- 98
rank <- which(cumsum(val^2)/sum(val^2)*100 > threshold)[1]
rank


rank <- sum(val > max(sqrt(.Machine$double.eps) * val[1]))
rank
loadings <- SVD$v[,1:rank]
X_red <- as.matrix(X) %*% loadings
dim(X_red)

# working
ics_res <- ICS(X_red, algorithm = "standard")


df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,1]^2)
k=1
k_details <- paste(length(k), ifelse(length(k)==1, k,
                                     paste(range(k), collapse = ":")),
                   sep = ", ")
y_axis <- bquote("DR rank = 141, ICSD"[ ~ k == .(k_details) ] ^2 )

plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top")+ 
  scale_x_continuous(
    # labels =   ~ sprintf(fmt = "%0.01e", .),
    guide  = guide_axis(check.overlap = T)) + 
  scale_y_continuous(
    # labels =   ~ sprintf(fmt = "%0.01e", .),
    guide  = guide_axis(check.overlap = T))


# save plot to pdf
file_plot = "figures/HTP2std/data_%s_%d.%s"
pdf(sprintf(file_plot,  "COVCOV4_DR", rank, "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()





try(ics_res <- ICS(X_red, algorithm = "standard", S1 = ICS_mcd_raw, S2 = ICS_cov))
try(ics_res <- ICS(X_red, algorithm = "standard", S1 = ICS_cov, S2 = ICS_mcd_raw))





# Regularized ---------------------------------------------------------------------

## RMCD --------------------------------------------------------------------
RMCD <- CovMrcd(X)
ics_res <- ICS(X, S1 = RMCD@cov, S2 = cov(X), algorithm = "standard")
component_plot(ics_res, clusters = clusters)

df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,1]^2)
k=1
k_details <- paste(length(k), ifelse(length(k)==1, k,
                                     paste(range(k), collapse = ":")),
                   sep = ", ")
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
file_plot = "figures/HTP2std/data_%s.%s"
pdf(sprintf(file_plot,  "MRCDCOV",  "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()



# Pseudo inverse - Ginv ----------------------------------------------------------

try(ics_res <- ICS_ginv(X, S1 = ICS_cov, S2 = ICS_cov4, 
                     algorithm = "whiten",
                     tolerance = .Machine$double.eps))



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
file_plot = "figures/HTP2std/data_%s.%s"
pdf(sprintf(file_plot,  "COVCOV4_GSVD", "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()


ics_res <- ICS_generalized(X, S1 = cov4mean_ginv_K, S2 = cov_K, 
                           algorithm = c("GSVD"), center = FALSE, 
                           fix_signs = c("scores", "W"), na.action = na.fail) 
df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,141]^2)

k_details <- 141
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
file_plot = "figures/HTP2std/data_%s.%s"
pdf(sprintf(file_plot,  "COVCOV4_GSVD_sym", "pdf"), width = 6.5, height = 5.75)
print(plot_ics)
dev.off()


