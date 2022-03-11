#' Eliminating the Unmeasured Confounders and Estimating
#' causal effect and lagged causal effect
#'
#' @param data an optional data frame containing the variables in the model.
#' @param x the variable name of exposure at time t2

#' @param y1 the variable name of outcome at time t1
#' @param y2 the variable name of outcome at time t2
#' @param y3 the variable name of outcome at time t3
#' @param centre the variable name of centre
#'





ncor <- function(data, centre, y1, y2, y3, x) {

    coef_ncor1 <- NULL
    se_ncor1 <- NULL
    coef_ncor2 <- NULL
    coef_ncor3 <- NULL
    centrename <- unique(data[, centre])
    len <- length(centrename)
    for (i in 1:len) {

        data_i_sample <- data[data[, centre] == centrename[i], ]
        data_i <- data_i_sample[, c(y1, y2, y3, x)]
        colnames(data_i) <- c("y1", "y2", "y3", "x")

        model_ncor1 <- lm(y3 ~ x, data = data_i)
        coef_ncor1[i] <- summary(model_ncor1)$coefficients[2, 1]
        se_ncor1[i] <- summary(model_ncor1)$coefficients[2, 2]

        model_ncor2 <- lm(y2 ~ x, data = data_i)
        coef_ncor2[i] <- summary(model_ncor2)$coefficients[2, 1]

        model_ncor3 <- lm(y1 ~ x, data = data_i)
        coef_ncor3[i] <- summary(model_ncor3)$coefficients[2, 1]

    }

    nn_centre <- length(centrename)
    we_matrix <- matrix(rep(0, nn_centre * nn_centre),
                        nrow = nn_centre, ncol = nn_centre)
    cen_data <- matrix(rep(1, nn_centre * nn_centre),
                       nrow = nn_centre, ncol = nn_centre)
    cen_data[upper.tri(cen_data)] <- 0
    cen_data[lower.tri(cen_data)] <- 0
    for (ppp in 1:nn_centre) {
        for (qqq in 1:nn_centre) {
            we_matrix[ppp, qqq] <- (1 / se_ncor1[ppp]) * (1 / se_ncor1[qqq]) *
                cen_data[ppp, qqq]
        }
    }

    xo_xo <- data.frame(1, coef_ncor2, coef_ncor3)
    xo_xo <- as.matrix(xo_xo)
    omega <- we_matrix
    beta1 <- solve(t(xo_xo) %*% omega %*% xo_xo) %*%
      t(xo_xo) %*% omega %*% coef_ncor1

    causal_ncor <- beta1[1]
    direct_ncor <- beta1[2]

    all_coef <- data.frame(coef_ncor1, coef_ncor2, coef_ncor3)
    causal_ncor_boot <- NULL

    direct_ncor_boot <- NULL

    for (boo in 1:1000) {
        hanghao <- sample(1:nn_centre, nn_centre, replace = T)
        we_matrix_boot <- matrix(rep(0, nn_centre * nn_centre),
                                 nrow = nn_centre, ncol = nn_centre)
        for (ppp in 1:nn_centre) {
            we_matrix_boot[ppp, ppp] <- we_matrix[hanghao[ppp], hanghao[ppp]]
        }
        data_boot <- all_coef[hanghao, ]

        xo_xo_boot <- data.frame(1, data_boot[, 2], data_boot[, 3])
        xo_xo_boot <- as.matrix(xo_xo_boot)
        omega_boot <- we_matrix_boot
        beta1_boot <- try(solve(t(xo_xo_boot) %*% omega_boot %*% xo_xo_boot) %*%
            t(xo_xo_boot) %*% omega_boot %*% data_boot[, 1], silent = T)
        if ("try-error" %in% class(beta1_boot)) {
            (next)()
        }
        causal_ncor_boot[boo] <- beta1_boot[1]
        direct_ncor_boot[boo] <- beta1_boot[2]
    }

    se_causal_ncor <- sd(causal_ncor_boot, na.rm = T)
    se_direct_ncor <- sd(direct_ncor_boot, na.rm = T)


    causal_ncor_25 <- quantile(causal_ncor_boot, c(0.025, 0.975), na.rm = T)[1]
    causal_ncor_975 <- quantile(causal_ncor_boot, c(0.025, 0.975), na.rm = T)[2]
    direct_ncor_25 <- quantile(direct_ncor_boot, c(0.025, 0.975), na.rm = T)[1]
    direct_ncor_975 <- quantile(direct_ncor_boot, c(0.025, 0.975), na.rm = T)[2]

    allresult <- c(causal_ncor, se_causal_ncor, causal_ncor_25, causal_ncor_975,
        direct_ncor, se_direct_ncor, direct_ncor_25, direct_ncor_975)

    return(allresult)

}
