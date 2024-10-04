cibersort <-
function (sig_matrix, mixture_file, perm = 0, QN = TRUE) 
{
    if (is.character(sig_matrix)) {
        X <- read.delim(sig_matrix, header = T, sep = "\t", row.names = 1, 
            check.names = F)
        X <- data.matrix(X)
    }
    else {
        X <- sig_matrix
    }
    if (is.character(mixture_file)) {
        Y <- read.delim(mixture_file, header = T, sep = "\t", 
            row.names = 1, check.names = F)
        Y <- data.matrix(Y)
    }
    else {
        Y <- mixture_file
    }
    X <- X[order(rownames(X)), ]
    Y <- Y[order(rownames(Y)), ]
    P <- perm
    if (max(Y) < 50) {
        Y <- 2^Y
    }
    if (QN == TRUE) {
        tmpc <- colnames(Y)
        tmpr <- rownames(Y)
        Y <- normalize.quantiles(Y)
        colnames(Y) <- tmpc
        rownames(Y) <- tmpr
    }
    Xgns <- row.names(X)
    Ygns <- row.names(Y)
    YintX <- Ygns %in% Xgns
    Y <- Y[YintX, ]
    XintY <- Xgns %in% row.names(Y)
    X <- X[XintY, ]
    X <- (X - mean(X))/sd(as.vector(X))
    if (P > 0) {
        nulldist <- sort(doPerm(P, X, Y)$dist)
    }
    header <- c("Mixture", colnames(X), "P-value", "Correlation", 
        "RMSE")
    output <- matrix()
    itor <- 1
    mixtures <- dim(Y)[2]
    pval <- 9999
    while (itor <= mixtures) {
        y <- Y[, itor]
        y <- (y - mean(y))/sd(y)
        result <- tryCatch({
          CoreAlg(X, y)
        }, error = function(e) {
          cat("Error in sample:", colnames(Y)[itor], ":", e$message, "\n")
          return(list(w = rep(NA, ncol(X)), mix_rmse = NA, mix_r = NA))  # Handle error gracefully
        })
        w <- result$w
        mix_r <- result$mix_r
        mix_rmse <- result$mix_rmse
        if (P > 0) {
            pval <- 1 - (which.min(abs(nulldist - mix_r))/length(nulldist))
        }
        out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
        if (itor == 1) {
            output <- out
        }
        else {
            output <- rbind(output, out)
        }
        itor <- itor + 1
    }
    obj <- rbind(header, output)
    obj <- obj[, -1]
    obj <- obj[-1, ]
    obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
    rownames(obj) <- colnames(Y)
    colnames(obj) <- c(colnames(X), "P-value", "Correlation", 
        "RMSE")
    obj
}
