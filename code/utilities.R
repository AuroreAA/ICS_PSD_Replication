svd_qr <- function(X,  ref = "previous", tol_eps = 10^(-8)){
  # Here there is a singularity issue. One solution is to first reduce the 
  # dimension. To ensure higher numerical stability of the subsequent methods
  # we suggest to permute the data and to use the QR decomposition instead of
  # the regular SVD decomposition.
  Xt <- X
  # Normalization by the mean 
  Xt.c <- sweep(X, 2, colMeans(X), "-")
  
  
  n <- nrow(Xt.c)
  p <- ncol(Xt.c)
  
  # Permutation by rows 
  # decreasing by infinity norm:  absolute maximum
  norm_inf <- apply(Xt.c, 1, function(x) max(abs(x)))
  order_rows <- order(norm_inf, decreasing = TRUE)
  Xt_row_per <- Xt.c[order_rows,]
  
  # QR decomposition of Xt with column pivoting from LAPACK 
  qr_Xt <- qr(1/sqrt(nrow(Xt.c)-1)*Xt_row_per, LAPACK = TRUE)
  
  # Estimation of rank q 
  # R is nxp, but with only zero for rows > p
  # the diag of R is already in decreasing order and is a good approximation
  # of the rank of X.c. To decide on which singular values are zero we use
  # a relative criteria based on previous values.
  # R should be pxp
  R <- qr.R(qr_Xt)
  r_all <- abs(diag(R))
  r_ratios <- r_all[2:length(r_all)]/r_all[1:(length(r_all)-1)]
  q <- which(r_ratios < max(dim(Xt.c)) *.Machine$double.eps)[1]
  q <- ifelse(is.na(q), length(r_all), q)
  q
  
  if (ref =="first"){
    r_ratios <-  r_all/r_all[1]
    q <- sum(r_ratios > max(n,p) *tol_eps)
    q
  }else{
    r_ratios <- r_all[2:length(r_all)]/r_all[1:(length(r_all)-1)]
    q <- which(r_ratios < max(n,p) *tol_eps)[1]
    q
  }
  print(q)
  
  # Q should be nxp but we are only interested in nxq
  Q1 <- qr.Q(qr_Xt)[,1:q]
  
  # QR decomposition of Rt 
  R_q <- R[1:q, ]
  qr_R <- qr(t(R_q), LAPACK = TRUE)
  Tau <- qr.Q(qr_R)[1:q, ]
  Omega1 <- qr.R(qr_R)[1:q, 1:q]
  
  # New X tilde 
  # permutation matrices
  # permutation of rows
  Pi2 <- data.frame(model.matrix(~ . -1, data = data.frame(row=as.character(order_rows))))
  Pi2 <- Pi2[,order(as.numeric(substr(colnames(Pi2), start = 4, stop = nchar(colnames(Pi2)))))]
  colnames(Pi2) <- rownames(Xt)
  
  # permutation of cols
  Pi3 <- data.frame(model.matrix(~ . -1, data = data.frame(col=as.character( qr_R$pivot))))
  Pi3 <- t(Pi3[,order(as.numeric(substr(colnames(Pi3), start = 4, stop = nchar(colnames(Pi3)))))])
  
  X_tilde <- sqrt(nrow(Xt)-1)* Tau %*% t(Pi3) %*% t(Q1)
  
  Xt_tilde <- t(Pi2) %*% t(X_tilde)
  
  return(Xt_tilde)
  
}



ICS_ginv <- function (X, S1 = ICS_cov, S2 = ICS_cov4, S1_args = list(), S2_args = list(), 
                      algorithm = c("whiten", "standard"), center = FALSE, 
                      fix_signs = c("scores", "W"), na.action = na.fail,
                      tolerance = .Machine$double.eps) 
{
  X <- na.action(X)
  X <- as.matrix(X)
  p <- ncol(X)
  if (p < 2L) 
    stop("'X' must be at least bivariate")
  algorithm <- match.arg(algorithm)
  whiten <- ifelse(algorithm == "whiten", TRUE, FALSE)
  S1_label <- deparse(substitute(S1))
  S2_label <- deparse(substitute(S2))
  
  if (isTRUE(whiten)) {
    if (!is.function(S2)) {
      warning("whitening requires 'S2' to be a function; ", 
              "proceeding without whitening")
      whiten <- FALSE
      algorithm <- "standard"
    }
  }
  center <- isTRUE(center)
  fix_signs <- match.arg(fix_signs)
  if (is.function(S1)) {
    S1_X <- ICS:::get_scatter(X, fun = S1, args = S1_args, label = S1_label)
  }
  else {
    S1_X <- ICS:::to_ICS_scatter.matrix(S1, label = S1_label)
    if (length(S1_args) > 0L) {
      warning("'S1' is not a function; ignoring additional arguments")
      S1_args <- list()
    }
  }
  S1_label <- S1_X$label
  
  W1 <- mat_sqrt_ginv(S1_X$scatter, inverse = TRUE, tolerance = tolerance)
  
  if (whiten) {
    Y <- X %*% W1
    S2_Y <- ICS:::get_scatter(Y, fun = S2, args = S2_args, 
                              label = S2_label)
    S2_label <- S2_Y$label
    missing_T2_X <- TRUE
  }else {
    if (is.function(S2)) {
      S2_X <- ICS:::get_scatter(X, fun = S2, args = S2_args, 
                                label = S2_label)
    }
    else {
      S2_X <- ICS:::to_ICS_scatter.matrix(S2, label = S2_label)
      if (length(S2_args) > 0L) {
        warning("'S2' is not a function; ignoring additional arguments")
        S2_args <- list()
      }
    }
    S2_label <- S2_X$label
    S2_Y <- ICS:::to_ICS_scatter.matrix(W1 %*% S2_X$scatter %*% W1, 
                                        label = S2_label)
    T2_X <- S2_X$location
    missing_T2_X <- is.null(T2_X)
  }
  if (center) {
    T1_X <- S1_X$location
    if (is.null(T1_X)) {
      warning("location component in 'S1' required for centering the data; ", 
              "proceeding without centering")
      center <- FALSE
    }
    else X_centered <- sweep(X, 2L, T1_X, "-")
  }
  
  S2_Y_eigen <- eigen(S2_Y$scatter, symmetric = TRUE)
  gen_kurtosis <- S2_Y_eigen$values
  if (!all(is.finite(gen_kurtosis))) {
    warning("some generalized kurtosis values are non-finite")
  }
  W <- crossprod(S2_Y_eigen$vectors, W1)
  if (fix_signs == "scores") {
    if (center && !(missing_T2_X || isTRUE(all.equal(T1_X, 
                                                     T2_X)))) {
      T1_Z <- T1_X %*% W
      T2_Z <- T2_X %*% W
      gen_skewness <- as.vector(T1_Z - T2_Z)
    }
    else {
      Z <- if (center) 
        tcrossprod(X_centered, W)
      else tcrossprod(X, W)
      gen_skewness <- colMeans(Z) - apply(Z, 2L, median)
    }
    skewness_signs <- ifelse(gen_skewness >= 0, 1, -1)
    gen_skewness <- skewness_signs * gen_skewness
    W_final <- sweep(W, 1L, skewness_signs, "*")
  }
  else {
    row_signs <- apply(W, 1L, .sign.max)
    row_norms <- sqrt(rowSums(W^2))
    W_final <- sweep(W, 1L, row_norms * row_signs, "/")
    gen_skewness <- NULL
  }
  if (center) 
    Z_final <- tcrossprod(X_centered, W_final)
  else Z_final <- tcrossprod(X, W_final)
  IC_names <- paste("IC", seq_len(p), sep = ".")
  names(gen_kurtosis) <- IC_names
  dimnames(W_final) <- list(IC_names, colnames(X))
  dimnames(Z_final) <- list(rownames(X), IC_names)
  if (!is.null(gen_skewness)) 
    names(gen_skewness) <- IC_names
  res <- list(gen_kurtosis = gen_kurtosis, W = W_final, scores = Z_final, 
              gen_skewness = gen_skewness, S1_label = S1_label, S2_label = S2_label, 
              S1_args = S1_args, S2_args = S2_args, algorithm = algorithm, 
              center = center, fix_signs = fix_signs)
  class(res) <- "ICS"
  res
}



ICS_generalized <- function (X, S1 = cov_K, S2 = cov4mean_ginv_K, S1_args = list(), S2_args = list(), 
                             algorithm = c("GSVD", "GEP"), center = FALSE, 
                             fix_signs = c("scores", "W"), na.action = na.fail) 
{
  X <- na.action(X)
  X <- as.matrix(X)
  p <- ncol(X)
  if (p < 2L) 
    stop("'X' must be at least bivariate")
  algorithm <- match.arg(algorithm)
  GSVD <- ifelse(algorithm == "GSVD", TRUE, FALSE)
  GEP <- ifelse(algorithm == "GEP", TRUE, FALSE)
  S1_label <- deparse(substitute(S1))
  S2_label <- deparse(substitute(S2))
  
  center <- isTRUE(center)
  fix_signs <- match.arg(fix_signs)
  if (is.function(S1)) {
    S1_X <- ICS:::get_scatter(X, fun = S1, args = S1_args, label = S1_label)
  }
  else {
    S1_X <- ICS:::to_ICS_scatter(S1, label = S1_label)
    if (length(S1_args) > 0L) {
      warning("'S1' is not a function; ignoring additional arguments")
      S1_args <- list()
    }
  }
  S1_label <- S1_X$label
  T1_X <- S1_X$location
  if (is.null(T1_X)) 
    T1_X <- colMeans(X)
  if (is.function(S2)) {
    S2_X <- ICS:::get_scatter(X, fun = S2, args = S2_args, 
                              label = S2_label)
  }
  else {
    S2_X <- ICS:::to_ICS_scatter(S2, label = S2_label)
    if (length(S2_args) > 0L) {
      warning("'S2' is not a function; ignoring additional arguments")
      S2_args <- list()
    }
  }
  S2_label <- S2_X$label
  missing_T2_X <- TRUE
  
  if (GSVD){

    # if (!all(dim(X)==dim(S1_X$scatter))) {
    #   stop("'S1' does not return a nxp matrix")
    # }
    # if (!all(dim(X)==dim(S2_X$scatter))) {
    #   stop("'S2' does not return a nxp matrix")
    # }
    # GSVD: B^-1A
    out_gsvd <- gsvd(A = S2_X$scatter, B = S1_X$scatter)
    gen_kurtosis_all <- data.frame(alpha = out_gsvd$alpha, 
                                   beta = out_gsvd$beta,
                                   rho = out_gsvd$alpha^2/out_gsvd$beta^2)
    
    # generalized eigenvectors
    R <- gsvd.R(out_gsvd)
    R0 <- gsvd.oR(out_gsvd)
    s <- ncol(R0) - ncol(R)
    Rinv_all <- matrix(0, ncol(out_gsvd$Q), ncol(out_gsvd$Q))
    diag(Rinv_all) <- 1
    Rinv_all[(s+1):ncol(out_gsvd$Q), (s+1):ncol(out_gsvd$Q)] <- solve(R)
    W <- out_gsvd$Q %*% Rinv_all
    
    # Warning: the generalized eigenvectors are not ordered
    # the first ones are associated to eigenvalues = NaN
    
    ind_nan <- 0:sum(is.na(gen_kurtosis_all$rho))
    if (sum(ind_nan) > 0){
      W <- W[,-ind_nan]
    }
    
    val_order <- order(gen_kurtosis_all$rho, decreasing = TRUE, na.last = NA)
    gen_kurtosis <- gen_kurtosis_all$rho[val_order]
    gen_kurtosis_all <- gen_kurtosis_all[unique(c(val_order,length(val_order)+ind_nan)),]
    W <- t(W[,val_order])
  }else if (GEP){
    
    out_GEP <- geigen(A = S2_X$scatter, B  = S1_X$scatter, symmetric = FALSE)
    gen_kurtosis_all <- out_GEP$values
    
    ind_real_rs <- which(Im(out_GEP$values)==0)
    val <- Re(out_GEP$values)[ind_real_rs]
    val_order <- order(val, decreasing = TRUE)
    gen_kurtosis <- val[val_order]
    
    W <- Re(out_GEP$vectors[,ind_real_rs])
    W <- t(W[,val_order])
  }
  
  if (center) {
    T1_X <- S1_X$location
    if (is.null(T1_X)) {
      warning("location component in 'S1' required for centering the data; ", 
              "proceeding without centering")
      center <- FALSE
    }
    else X_centered <- sweep(X, 2L, T1_X, "-")
  }
  if (!all(is.finite(gen_kurtosis))) {
    warning("some generalized kurtosis values are non-finite")
  }
  
  if (fix_signs == "scores") {
    if (center && !(missing_T2_X || isTRUE(all.equal(T1_X, 
                                                     T2_X)))) {
      T1_Z <- T1_X %*% W
      T2_Z <- T2_X %*% W
      gen_skewness <- as.vector(T1_Z - T2_Z)
    }
    else {
      Z <- if (center) 
        tcrossprod(X_centered, W)
      else tcrossprod(X, W)
      gen_skewness <- colMeans(Z) - apply(Z, 2L, median)
    }
    skewness_signs <- ifelse(gen_skewness >= 0, 1, -1)
    gen_skewness <- skewness_signs * gen_skewness
    W_final <- sweep(W, 1L, skewness_signs, "*")
  }
  else {
    row_signs <- apply(W, 1L, .sign.max)
    row_norms <- sqrt(rowSums(W^2))
    W_final <- sweep(W, 1L, row_norms * row_signs, "/")
    gen_skewness <- NULL
  }
  if (center) 
    Z_final <- tcrossprod(X_centered, W_final)
  else Z_final <- tcrossprod(X, W_final)
  IC_names <- paste("IC", seq_len(ncol(Z_final)), sep = ".")
  names(gen_kurtosis) <- IC_names
  dimnames(W_final) <- list(IC_names, colnames(X))
  dimnames(Z_final) <- list(rownames(X), IC_names)
  if (!is.null(gen_skewness)) 
    names(gen_skewness) <- IC_names
  res <- list(gen_kurtosis = gen_kurtosis, W = W_final, scores = Z_final, 
              gen_skewness = gen_skewness, S1_label = S1_label, S2_label = S2_label, 
              S1_args = S1_args, S2_args = S2_args, algorithm = algorithm, 
              center = center, fix_signs = fix_signs, 
              gen_kurtosis_all = gen_kurtosis_all
  )
  class(res) <- "ICS"
  res
}

mat_sqrt_ginv <- function (A, inverse = TRUE, tolerance = .Machine$double.eps) 
{
  power <- if (inverse) 
    -0.5
  else 0.5
  eigen_A <- eigen(A, symmetric = TRUE)
  r1 <- which(eigen_A$values > max(tolerance *eigen_A$values[1]))
  eigen_A$vectors[,r1] %*% tcrossprod(diag(eigen_A$values^power)[r1,r1], 
                                      eigen_A$vectors[,r1])
}


cov_K <- function(X) {1/sqrt(nrow(X)-1)*X}

cov4mean_ginv_K <- function (X, tolerance = sqrt(.Machine$double.eps)) 
{
  X <- as.matrix(X)
  p <- dim(X)[2]
  n <- dim(X)[1]
  eigen.cov <-  eigen(cov(X), symmetric = TRUE)
  r1 <- which(eigen.cov$values > max(tolerance))
  cov.inv.demi <- eigen.cov$vectors[,r1] %*% tcrossprod(diag(eigen.cov$values[r1]^(-0.5)), 
                                                        eigen.cov$vectors[,r1])
  data.centered <- sweep(X, 2, colMeans(X), "-")
  radius <- sqrt(rowSums((data.centered %*% cov.inv.demi)^2))
  y <- radius * data.centered
  V <- (1/sqrt((n * (p + 2)))) * y
  return(V)
}

cov4mean_ginv <- function (X, tolerance = sqrt(.Machine$double.eps)){
  K <- cov4mean_ginv_K(X, tolerance = tolerance)
  crossprod(K)
}

