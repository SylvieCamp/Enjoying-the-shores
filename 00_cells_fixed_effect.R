## Cell-level estimates for Figure 2a: one intercept-only model per driver x service
## cell, fitted where the cell has at least two independent studies.
suppressMessages({library(ordinal); library(MASS)}); set.seed(20260909)
d <- read.csv("raw.csv")
d$Vote <- factor(tolower(trimws(d$Vote_counting)),
                 levels = c("negative","neutral","ambiguous","positive"), ordered = TRUE)
d$study <- factor(d$biblio_internal_id)

gh <- function(n=60){i<-1:(n-1);a<-sqrt(i);J<-diag(0,n);J[cbind(i,i+1)]<-a;J[cbind(i+1,i)]<-a
  e<-eigen(J,symmetric=TRUE); list(x=rev(e$values), w=rev(e$vectors[1,]^2))}
GH <- gh(60)
pm <- function(th, s) sum(GH$w * (1 - plogis(th - GH$x * s)))

out <- list()
for (dr in unique(d$Driver_change)) for (sv in unique(d$Service)) {
  z <- droplevels(subset(d, Driver_change == dr & Service == sv))
  if (nrow(z) == 0) next
  ns <- length(unique(z$study)); z$Vote <- droplevels(z$Vote)
  row <- data.frame(driver = dr, service = sv, k = nrow(z), studies = ns,
                    pA = NA_real_, loA = NA_real_, hiA = NA_real_,
                    pB = NA_real_, loB = NA_real_, hiB = NA_real_)
  if (nlevels(z$Vote) >= 2) {
    m <- try(clm(Vote ~ 1, data = z, Hess = TRUE), silent = TRUE)
    if (!inherits(m, "try-error")) {
      nth <- nlevels(z$Vote) - 1; p <- m$alpha
      V <- tryCatch(vcov(m), error = function(e) NULL)
      row$pA <- pm(p[nth], 0)
      if (!is.null(V) && all(is.finite(V))) {
        S <- MASS::mvrnorm(4000, p, V)
        q <- quantile(1 - plogis(S[, nth]), c(.025, .975))
        row$loA <- q[1]; row$hiA <- q[2]
      }
    }
    if (ns >= 2) {
      mm <- try(clmm(Vote ~ 1 + (1 | study), data = z, Hess = TRUE), silent = TRUE)
      if (!inherits(mm, "try-error")) {
        nth <- nlevels(z$Vote) - 1; p <- mm$alpha; s <- as.numeric(mm$ST$study)
        V <- tryCatch(vcov(mm)[names(p), names(p), drop = FALSE], error = function(e) NULL)
        row$pB <- pm(p[nth], s)
        if (!is.null(V) && all(is.finite(V)) && s < 8) {
          S <- MASS::mvrnorm(4000, p, V)
          vals <- vapply(seq_len(nrow(S)), function(i) pm(S[i, nth], s), 0)
          q <- quantile(vals, c(.025, .975)); row$loB <- q[1]; row$hiB <- q[2]
        }
      }
    }
  }
  out[[length(out) + 1]] <- row
}
res <- do.call(rbind, out)
write.csv(res, "cells_final.csv", row.names = FALSE)
cat("cells:", nrow(res), "| with clmm estimate:", sum(!is.na(res$pB)), "\n")
