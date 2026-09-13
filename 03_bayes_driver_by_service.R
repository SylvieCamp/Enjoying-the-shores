## Driver x CES blocks: the frequentist clmm is degenerate where the reported
## outcomes are close to separated (tau^2 = 104 for management, 29 for pollution).
## We therefore estimate these blocks in a Bayesian framework with a regularising
## half-normal(0,1) prior on the random-effect SD, which is defined in that regime.
suppressMessages({library(brms); library(posterior)}); options(mc.cores = 4)

d <- read.csv("raw.csv")
d$Vote <- factor(tolower(trimws(d$Vote_counting)),
                 levels = c("negative","neutral","ambiguous","positive"), ordered = TRUE)
d$study <- factor(d$biblio_internal_id)

gh <- function(n=60){i<-1:(n-1);a<-sqrt(i);J<-diag(0,n);J[cbind(i,i+1)]<-a;J[cbind(i+1,i)]<-a
  e<-eigen(J,symmetric=TRUE); list(x=rev(e$values), w=rev(e$vectors[1,]^2))}
GH <- gh(60)

block <- function(driver) {
  D <- droplevels(subset(d, Driver_change == driver))
  k <- tapply(D$study, D$Service, function(z) length(unique(z)))
  D <- droplevels(D[D$Service %in% names(k)[k >= 2], ])
  D$Service <- factor(D$Service); D$Vote <- droplevels(D$Vote)
  f <- brm(Vote ~ Service + (1 | study), data = D, family = cumulative("logit"),
           prior = c(prior(normal(0,1), class = "sd"), prior(normal(0,2.5), class = "b")),
           chains = 4, iter = 3000, warmup = 1000, seed = 1,
           control = list(adapt_delta = 0.97), refresh = 0, silent = 2)
  dr <- as_draws_df(f); lev <- levels(D$Service)
  nth <- nlevels(D$Vote) - 1
  th <- dr[[paste0("b_Intercept[", nth, "]")]]; sdu <- dr$sd_study__Intercept
  bn <- grep("^b_Service", names(dr), value = TRUE)
  B <- cbind(0, as.matrix(as.data.frame(dr)[, bn, drop = FALSE])); colnames(B) <- lev
  P <- sapply(lev, function(l) vapply(seq_along(th), function(s)
        sum(GH$w * (1 - plogis(th[s] - B[s,l] - GH$x * sdu[s]))), 0))
  q <- apply(P, 2, quantile, c(.5,.025,.975))
  cat("\n---", driver, "--- posterior sd(study) median",
      round(median(sdu),2), "CrI [", round(quantile(sdu,.025),2), ",",
      round(quantile(sdu,.975),2), "]\n")
  o <- data.frame(driver = driver, level = lev,
    studies = as.integer(tapply(D$study, D$Service, function(z) length(unique(z)))[lev]),
    k = as.integer(table(D$Service)[lev]),
    pBayes = q[1,], loBayes = q[2,], hiBayes = q[3,],
    pPos = colMeans(P > 0.5),
    verdict = ifelse(q[2,] > .5, "positive", ifelse(q[3,] < .5, "negative", "not distinguishable")),
    sd_median = round(median(sdu),2), stringsAsFactors = FALSE)
  print(o[, c("level","studies","k","pBayes","loBayes","hiBayes","pPos","verdict")],
        row.names = FALSE, digits = 3)
  o
}

## single-cell drivers: only Recreation has >= 2 studies, so an intercept-only model
cell <- function(driver, service) {
  D <- droplevels(subset(d, Driver_change == driver & Service == service)); D$Vote <- droplevels(D$Vote)
  if (nlevels(D$Vote) >= 3) {
    f <- brm(Vote ~ 1 + (1 | study), data = D, family = cumulative("logit"),
             prior = prior(normal(0,1), class = "sd"),
             chains = 4, iter = 3000, warmup = 1000, seed = 1,
             control = list(adapt_delta = 0.97), refresh = 0, silent = 2)
    dr <- as_draws_df(f); nth <- nlevels(D$Vote) - 1
    th <- dr[[paste0("b_Intercept[", nth, "]")]]; sdu <- dr$sd_study__Intercept
    P <- vapply(seq_along(th), function(s) sum(GH$w * (1 - plogis(th[s] - GH$x * sdu[s]))), 0)
  } else {   # two outcome levels only: the cumulative model reduces to a logistic one
    D$pos <- as.integer(D$Vote == "positive")
    f <- brm(pos ~ 1 + (1 | study), data = D, family = bernoulli("logit"),
             prior = prior(normal(0,1), class = "sd"),
             chains = 4, iter = 3000, warmup = 1000, seed = 1,
             control = list(adapt_delta = 0.97), refresh = 0, silent = 2)
    dr <- as_draws_df(f); b <- dr$b_Intercept; sdu <- dr$sd_study__Intercept
    P <- vapply(seq_along(b), function(s) sum(GH$w * plogis(b[s] + GH$x * sdu[s])), 0)
  }
  q <- quantile(P, c(.5,.025,.975))
  cat("\n---", driver, "/", service, "--- sd median", round(median(sdu),2), "\n")
  data.frame(driver = driver, level = service,
    studies = length(unique(D$study)), k = nrow(D),
    pBayes = q[1], loBayes = q[2], hiBayes = q[3], pPos = mean(P > 0.5),
    verdict = ifelse(q[2] > .5, "positive", ifelse(q[3] < .5, "negative", "not distinguishable")),
    sd_median = round(median(sdu),2), stringsAsFactors = FALSE)
}
out <- rbind(block("Management"), block("Pollution"), block("LSUC"),
             cell("Climate-change", "Recreation"), cell("Direct-exploitation", "Recreation"))
write.csv(out, "bayes_blocks.csv", row.names = FALSE)
