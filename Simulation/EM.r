EM <- function(data, starting_values = NULL, tol = 1e-5, max_iter = 50) {
    start_time <- Sys.time()
    
    if (is.null(starting_values)) {
        avg <- mean(data)
        var <- var(data)
        starting_values <- numeric(2)
        starting_values[1] <- avg * avg / (var - avg) # 初始化 r
        starting_values[2] <- avg / var # 初始化 p
    }

    n <- length(data)
    params <- starting_values
    converged <- FALSE
    inner_iter <- 0

    for (iter in 1:max_iter) {
        r <- params[1]
        p <- params[2]
        beta <- p / (1 - p)

        # E步: 计算隐变量lambda的期望
        z <- (r + data) / (beta + 1)
        log_z <- digamma(r + data) - log(beta + 1)

        # M步: 更新参数
        beta_new <- n * r / sum(z)

        r_func <- function(a) {
            log(a) - digamma(a) - log(mean(z)) + mean(log_z)
        }
        res <- tryCatch(
            {
                uniroot(r_func, c(0.01, 100), tol = .Machine$double.eps)
            },
            error = function(e) {
                warning("Failed to converge in uniroot. Stopping iterations.")
                break
            }
        )
        r_new <- res$root
        inner_iter <- inner_iter + res$iter

        p_new <- beta_new / (1 + beta_new)

        params_new <- c(r_new, p_new)
        if (sum(abs(params_new - params)) < tol) {
            params <- params_new
            converged <- TRUE
            break
        }

        # 检查新参数是否在合理范围内
        if (params_new[1] <= 0 || params_new[2] <= 0 || params_new[2] >= 1) {
            warning("Parameters out of bounds. Iteration stopped.")
            break
        }

        params <- params_new
    }

    end_time <- Sys.time() # 结束时间
    duration <- difftime(end_time, start_time, units = "secs") # 计算时间

    return(list(params = params, iterations = iter + inner_iter, converged = converged, duration = duration))
}