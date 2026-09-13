## Is the frequentist clmm over-correcting?
## The concern: tau-hat = 3.36 (logit SD) may be inflated by quasi-separation
## (studies contributing many effect sizes all of one sign push their random
## intercept toward +/- infinity), and a too-large tau over-attenuates the
## marginal probabilities.
## Test: refit as a Bayesian cumulative model with priors that regularise tau,
## and compare the marginal P(positive) with the frequentist estimates.
suppressMessages({library(brms); library(posterior)})
options(mc.cores = 4)

d <- read.csv("raw.csv")
d$Vote <- factor(tolower(trimws(d$Vote_counting)),
                  levels = c("negative","neutral","ambiguous","positive"), ordered = TRUE)
d$study <- factor(d$biblio_internal_id)
d$Driver_change <- factor(d$Driver_change)
d <- droplevels(subset(d, Driver_change != "Coastal-erosion"))   # single study

gh <- function(n = 60) { i <- 1:(n-1); a <- sqrt(i); J <- diag(0, n)
  J[cbind(i,i+1)] <- a; J[cbind(i+1,i)] <- a; e <- eigen(J, symmetric = TRUE)
  list(x = rev(e$values), w = rev(e$vectors[1,]^2)) }
GH <- gh(60)

marginal_ppos <- function(draws, lev) {
  th3 <- draws$`b_Intercept[3]`; sdu <- draws$sd_study__Intercept
  bn <- grep("^b_Driver_change", names(draws), value = TRUE)   # in level order
  stopifnot(length(bn) == length(lev) - 1)
  B <- cbind(0, as.matrix(as.data.frame(draws)[, bn, drop = FALSE]))
  colnames(B) <- lev
  out <- sapply(lev, function(l)
    vapply(seq_along(th3), function(s)
      sum(GH$w * (1 - plogis(th3[s] - B[s, l] - GH$x * sdu[s]))), 0))
  out
}

fit_one <- function(tag, prior_sd) {
  f <- brm(Vote ~ Driver_change + (1 | study), data = d,
           family = cumulative("logit"),
           prior = c(prior_string(prior_sd, class = "sd"),
                     prior(normal(0, 2.5), class = "b")),
           chains = 4, iter = 3000, warmup = 1000, seed = 1,
           control = list(adapt_delta = 0.95), refresh = 0, silent = 2)
  dr <- as_draws_df(f)
  lev <- levels(d$Driver_change)
  P <- marginal_ppos(dr, lev)
  q <- apply(P, 2, quantile, c(0.5, 0.025, 0.975))
  cat("\n=====", tag, "=====\n")
  cat("prior on sd:", prior_sd, "\n")
  cat("posterior sd(study): median", round(median(dr$sd_study__Intercept), 2),
      " 95% CrI [", round(quantile(dr$sd_study__Intercept, .025), 2), ",",
      round(quantile(dr$sd_study__Intercept, .975), 2), "]\n")
  o <- data.frame(driver = lev,
                  P_pos = sprintf("%.2f [%.2f-%.2f]", q[1,], q[2,], q[3,]),
                  `P(<0.5)` = sprintf("%.2f", colMeans(P < 0.5)),
                  verdict = ifelse(q[2,] > .5, "POS", ifelse(q[3,] < .5, "NEG", "ns")),
                  check.names = FALSE, stringsAsFactors = FALSE)
  print(o, row.names = FALSE)
  invisible(cbind(prior = tag, o))
}

a <- fit_one("weakly informative  (half-t(3,0,2.5))", "student_t(3, 0, 2.5)")
b <- fit_one("regularising        (half-normal(0,1))", "normal(0, 1)")
c_ <- fit_one("strongly regularising (exponential(1))", "exponential(1)")
write.csv(rbind(a, b, c_), "bayes_drivers.csv", row.names = FALSE)
