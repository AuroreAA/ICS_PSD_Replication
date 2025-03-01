

# Libraries --------------------------------------------------------------

library(ICSOutlier)
library(ICSClust)
library(geigen)
library(MASS)
library(moments)
library(rrcov)
library(ggthemes)



source("code/utilities.R")

colors <- ggthemes::colorblind_pal()(4)[c(3, 2)]
set.seed(20240725)
# Data: 2 groups collinear  --------------------------------------------------------------------

df <- mixture_sim(pct_clusters = c(0.5, 0.5), n = 1000, p = 3, delta = 10)
df <- data.frame(df, X4 = df[,2]-3*df[,3], X5 = df[,3]+5*df[,4])
clusters <- df[,1]
X_ini <- as.matrix(df[,-1])



for(transformation in c("", "std")){
  if (transformation == ""){
    X <- X_ini
  }else if (transformation == "std"){
    X <- scale(X_ini, center = TRUE, scale = TRUE)
  }
  
  plot_ini <- component_plot(data.frame(X), clusters = clusters, colors = colors)
  
  # save plot to pdf
  file_plot = "figures/collinear/data_%s.%s"
  pdf(sprintf(file_plot, transformation, "pdf"), width = 6.5, height = 5.75)
  print(plot_ini)
  dev.off()
  
  
  # not working
  # ics_res <- ICS(X, algorithm = "QR")
  # ics_res <- ICS(X, algorithm = "standard", S1 = ICS_mcd_raw, S2 = ICS_cov)
  # ics_res <- ICS(X, algorithm = "standard")
  # ics_res <- ICS(X, algorithm = "whiten")
  
  
  # Reduction dimension -----------------------------------------------------
  
  
  ## SVD ---------------------------------------------------------------------
  SVD <- svd(X)
  val <- SVD$d
  val
  rank <- sum(val > max(sqrt(.Machine$double.eps) * val[1]))
  rank
  rank <- sum(val > max(dim(X)) * .Machine$double.eps * val[1])
  rank
  # cumsum(val^2)/sum(val^2)*100
  # threshold <- 98
  # rank <- which(cumsum(val^2)/sum(val^2)*100 > threshold)[1]
  # rank
  
  loadings <- SVD$v[,1:rank]
  X_red <- as.matrix(X) %*% loadings
  dim(X_red)
  
  
  
  
  ### COV-COV4 ----
  ics_res <- ICS(X_red, algorithm = "standard")
  colnames(ics_res$scores) <- paste("DR", colnames(ics_res$scores), sep = "-")
  #select_plot(ics_res)
  plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )
  
  # save plot to pdf
  file_plot = "figures/collinear/data_%s_%s.%s"
  pdf(sprintf(file_plot, transformation, "dim_reduc_COVCOV4", "pdf"), width = 6.5, height = 5.75)
  print(plot_ics)
  dev.off()
  
  
  ics_res <- ICS(X_red, S1 = ICS_cov4, S2 = ICS_cov, algorithm = "standard")
  colnames(ics_res$scores) <- paste("DR", colnames(ics_res$scores), sep = "-")
  
  #select_plot(ics_res)
  plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )
  
  # save plot to pdf
  file_plot = "figures/collinear/data_%s_%s.%s"
  pdf(sprintf(file_plot, transformation, "dim_reduc_COVCOV4_sym", "pdf"), width = 6.5, height = 5.75)
  print(plot_ics)
  dev.off()
  
  
  
  
  ### MCD-COV ----
  
  ics_res <- ICS(X_red, algorithm = "standard", S1 = ICS_mcd_raw, S2 = ICS_cov)
  colnames(ics_res$scores) <- paste("DR", colnames(ics_res$scores), sep = "-")
  
  plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )
  # save plot to pdf
  file_plot = "figures/collinear/data_%s_%s.%s"
  pdf(sprintf(file_plot, transformation, "dim_reduc_MCDCOV", "pdf"), width = 6.5, height = 5.75)
  print(plot_ics)
  dev.off()
  
  
  ics_res <- ICS(X_red, algorithm = "standard", S1 = ICS_cov, S2 = ICS_mcd_raw)
  colnames(ics_res$scores) <- paste("DR", colnames(ics_res$scores), sep = "-")
  
  plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )
  # save plot to pdf
  file_plot = "figures/collinear/data_%s_%s.%s"
  pdf(sprintf(file_plot, transformation, "dim_reduc_MCDCOV_sym", "pdf"), width = 6.5, height = 5.75)
  print(plot_ics)
  dev.off()
  
  
  
  
  # Regularized ---------------------------------------------------------------------
  
  ## COV-COV4 ----------------------------------------------------------------
  # K_V2 <- cov4mean_ginv_K(X)
  # V2 <- crossprod(K_V2)
  # solve(V2)
  # ics_res <- ICS(X, S1 = cov(X), S2 = V2) 
  # ics_res <- ICS(X, S1 = V2, S2 = cov(X)) 
  # select_plot(ics_res)
  # component_plot(ics_res, clusters = clusters, colors = colors)
  
  
  
  ## RMCD --------------------------------------------------------------------
  RMCD <- CovMrcd(X)
  try({ics_res <- ICS(X, S1 = RMCD@cov, S2 = cov(X), algorithm = "standard")
  colnames(ics_res$scores) <- paste("MRCD", colnames(ics_res$scores), sep = "-")
  
  plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )+
    scale_x_continuous(
      # labels =   ~ sprintf(fmt = "%0.01e", .),
      guide  = guide_axis(check.overlap = T)) + 
    scale_y_continuous(
      # labels =   ~ sprintf(fmt = "%0.01e", .),
      guide  = guide_axis(check.overlap = T))
  
  # save plot to pdf
  file_plot = "figures/collinear/data_%s_%s.%s"
  pdf(sprintf(file_plot, transformation, "regul_RMCDCOV", "pdf"), width = 6.5, height = 5.75)
  print(plot_ics)
  dev.off()})
  
  # not working
  try({
    ics_res <- ICS(X, S1 = cov(X), S2 = RMCD@cov, algorithm = "standard")
    colnames(ics_res$scores) <- paste("MRCD", colnames(ics_res$scores), sep = "-")
    
    plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )+
      scale_x_continuous(
        # labels =   ~ sprintf(fmt = "%0.01e", .),
        guide  = guide_axis(check.overlap = T)) + 
      scale_y_continuous(
        # labels =   ~ sprintf(fmt = "%0.01e", .),
        guide  = guide_axis(check.overlap = T))
    
    # save plot to pdf
    file_plot = "figures/collinear/data_%s_%s.%s"
    pdf(sprintf(file_plot, transformation, "regul_RMCDCOV_sym", "pdf"), width = 6.5, height = 5.75)
    print(plot_ics)
    dev.off()
  })
  
  
  # Pseudo inverse - Ginv ----------------------------------------------------------
  ### COV-COV4 -----
  try({ics_res <- ICS_ginv(X, S1 = ICS_cov, S2 = ICS_cov4, 
                           algorithm = "whiten",
                           tolerance = .Machine$double.eps)
  colnames(ics_res$scores) <- paste("GINV", colnames(ics_res$scores), sep = "-")
  
  
  
  plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors ) +
    scale_x_continuous(
      # labels =   ~ sprintf(fmt = "%0.01e", .),
      guide  = guide_axis(check.overlap = T)) + 
    scale_y_continuous(
      # labels =   ~ sprintf(fmt = "%0.01e", .),
      guide  = guide_axis(check.overlap = T))
  
  
  
  # save plot to pdf
  file_plot = "figures/collinear/data_%s_%s.%s"
  pdf(sprintf(file_plot, transformation, "GINV_COVCOV4", "pdf"), width = 6.5, height = 5.75)
  print(plot_ics)
  dev.off()
  })
  
  try({
    ics_res <- ICS_ginv(X, S1 = ICS_cov4, S2 = ICS_cov, 
                        algorithm = "whiten",
                        tolerance = .Machine$double.eps)# save plot to pdf
    colnames(ics_res$scores) <- paste("GINV", colnames(ics_res$scores), sep = "-")
    plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )+
      scale_x_continuous(
        # labels =   ~ sprintf(fmt = "%0.01e", .),
        guide  = guide_axis(check.overlap = T)) + 
      scale_y_continuous(
        # labels =   ~ sprintf(fmt = "%0.01e", .),
        guide  = guide_axis(check.overlap = T))
    
    file_plot = "figures/collinear/data_%s_%s.%s"
    pdf(sprintf(file_plot, transformation, "GINV_COVCOV4_sym", "pdf"), width = 6.5, height = 5.75)
    print(plot_ics)
    dev.off()
  })
  
  
  
  ### MCD-COV ----
  try(ics_res <- ICS_ginv(X, S1 = ICS_mcd_raw, S2 = ICS_cov, 
                          algorithm = "whiten",
                          tolerance = .Machine$double.eps))
  
  
  try(ics_res <- ICS_ginv(X, S1 = ICS_cov, S2 = ICS_mcd_raw, 
                          algorithm = "whiten",
                          tolerance = .Machine$double.eps))
  
  
  # GSVD --------------------------------------------------------------------
  ics_res <- ICS_generalized(X, S1 = cov_K, S2 = cov4mean_ginv_K, 
                             algorithm = c("GSVD"), center = FALSE, 
                             fix_signs = c("scores", "W"), na.action = na.fail) 
  colnames(ics_res$scores) <- paste("GSVD", colnames(ics_res$scores), sep = "-")
  
  plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )
  # save plot to pdf
  file_plot = "figures/collinear/data_%s_%s.%s"
  pdf(sprintf(file_plot, transformation, "GSVD_COVCOV4", "pdf"), width = 6.5, height = 5.75)
  print(plot_ics)
  dev.off()
  
  ics_res <- ICS_generalized(X, S1 = cov4mean_ginv_K, S2 = cov_K, 
                             algorithm = c("GSVD"), center = FALSE, 
                             fix_signs = c("scores", "W"), na.action = na.fail) 
  colnames(ics_res$scores) <- paste("GSVD", colnames(ics_res$scores), sep = "-")
  plot_ics <- component_plot(ics_res, clusters = clusters, colors = colors )
  # save plot to pdf
  file_plot = "figures/collinear/data_%s_%s.%s"
  pdf(sprintf(file_plot, transformation, "GSVD_COVCOV4_sym", "pdf"), width = 6.5, height = 5.75)
  print(plot_ics)
  dev.off()
  
  
}




# GEP ---------------------------------------------------------------------
# ics_res <- ICS_generalized(X, S1 = ICS_cov, S2 = cov4mean_ginv, 
#                            algorithm = c("GEP"), center = FALSE, 
#                            fix_signs = c("scores", "W"), na.action = na.fail) 
# 
# select_plot(ics_res)
# component_plot(ics_res, clusters = clusters, colors = colors)
# gen_kurtosis(ics_res)
# ics_res$gen_kurtosis_all