################################################################################
##  Enjoying the shores — master analysis on the definitive data file
##  Input : Raw_data_14-02-24.xlsx  (exported to raw.csv)
##  Output: verification against published results, Table S11, Table S12,
##          and the data behind Figures 1 and 2.
################################################################################
suppressMessages({library(ordinal); library(metafor); library(MASS)})
set.seed(20260909)

d <- read.csv("raw.csv", check.names = FALSE)
d$Vote  <- factor(tolower(trimws(d$Vote_counting)),
                  levels = c("negative","neutral","ambiguous","positive"), ordered = TRUE)
d$study <- factor(d$biblio_internal_id)
d$yi    <- suppressWarnings(as.numeric(as.character(d$Effect_size)))

cat("### CORPUS\n")
cat("effect sizes:", nrow(d), "| studies:", nlevels(d$study), "\n")
print(table(d$Vote))
cat("validity mean:", round(mean(d$validity), 3), "\n")
cat("effect sizes per study: median", median(table(d$study)),
    "| max", max(table(d$study)), "\n")
cat("dispersion reported (Data_with_SD):\n"); print(table(d$Data_with_SD, useNA = "ifany"))
cat("\nmetric families:\n"); print(table(d$Value_type, useNA = "ifany"))

################################ helpers ######################################
gh <- function(n = 60) { i <- 1:(n-1); a <- sqrt(i); J <- diag(0, n)
  J[cbind(i,i+1)] <- a; J[cbind(i+1,i)] <- a; e <- eigen(J, symmetric = TRUE)
  list(x = rev(e$values), w = rev(e$vectors[1,]^2)) }
GH <- gh(60)
pm  <- function(th, eta, s) vapply(eta, function(e)
        sum(GH$w * (1 - plogis(th - e - GH$x * s))), 0)
sim <- function(mod, s, nl, nth, n = 5000) {
  p <- c(mod$alpha, mod$beta)
  V <- tryCatch(vcov(mod)[names(p), names(p), drop = FALSE], error = function(e) NULL)
  est <- pm(p[nth], c(0, p[-(1:nth)]), s)
  if (is.null(V) || any(!is.finite(V))) return(list(e = est, ci = matrix(NA_real_, 2, nl)))
  S <- t(apply(MASS::mvrnorm(n, p, V), 1, function(r) pm(r[nth], c(0, r[-(1:nth)]), s)))
  if (nl == 1) S <- matrix(S, ncol = 1)
  list(e = est, ci = apply(S, 2, quantile, c(.025, .975))) }
fmt <- function(o) ifelse(is.na(o$ci[1,]), sprintf("%.2f (n.e.)", o$e),
        sprintf("%.2f (%.2f-%.2f)", o$e, o$ci[1,], o$ci[2,]))
vd  <- function(o) ifelse(is.na(o$ci[1,]), "-", ifelse(o$ci[1,] > .5, "positive",
        ifelse(o$ci[2,] < .5, "negative", "not distinguishable")))

################################ ordinal ######################################
ordinal_block <- function(tag, group, D) {
  D[[group]] <- factor(D[[group]]); D$Vote <- droplevels(D$Vote); D <- droplevels(D)
  lev <- levels(D[[group]]); nl <- length(lev); nth <- nlevels(D$Vote) - 1
  mA <- clm(as.formula(paste("Vote ~", group)), data = D, Hess = TRUE)
  A  <- sim(mA, 0, nl, nth)
  k  <- tapply(D$study, D[[group]], function(z) length(unique(z)))
  Dr <- droplevels(D[D[[group]] %in% names(k)[k >= 2], ]); Dr$Vote <- droplevels(Dr$Vote)
  levr <- levels(Dr[[group]]); nthr <- nlevels(Dr$Vote) - 1
  mAr <- clm(as.formula(paste("Vote ~", group)), data = Dr, Hess = TRUE)   # same data as mB
  mB <- clmm(as.formula(paste("Vote ~", group, "+ (1|study)")), data = Dr, Hess = TRUE)
  sdu <- as.numeric(mB$ST$study); B <- sim(mB, sdu, length(levr), nthr)
  agg <- do.call(rbind, lapply(split(Dr, list(Dr$study, Dr[[group]]), drop = TRUE),
    function(z) { t <- table(z$Vote); data.frame(g = z[[group]][1],
      Vote = factor(names(t)[which.max(t)], levels = levels(D$Vote), ordered = TRUE)) }))
  agg$g <- factor(agg$g, levels = levr); agg$Vote <- droplevels(agg$Vote)
  mC <- clm(Vote ~ g, data = agg, Hess = TRUE)
  C  <- sim(mC, 0, length(levr), nlevels(agg$Vote) - 1)
  m <- match(lev, levr)
  data.frame(model = tag, level = lev,
    studies = as.integer(k[lev]), k = as.integer(table(D[[group]])[lev]),
    pA = A$e, loA = A$ci[1,], hiA = A$ci[2,],
    pB = sapply(m, function(i) if (is.na(i)) NA else B$e[i]),
    loB = sapply(m, function(i) if (is.na(i)) NA else B$ci[1,i]),
    hiB = sapply(m, function(i) if (is.na(i)) NA else B$ci[2,i]),
    A = fmt(A),
    B = sapply(m, function(i) if (is.na(i)) "-" else fmt(B)[i]),
    C = sapply(m, function(i) if (is.na(i)) "-" else fmt(C)[i]),
    verdictA = vd(A),
    verdictB = sapply(m, function(i) if (is.na(i)) "single study" else vd(B)[i]),
    tau2 = round(sdu^2, 1), lrt = format.pval(anova(mAr, mB)$`Pr(>Chisq)`[2], digits = 3),
    stringsAsFactors = FALSE) }

S11 <- rbind(
  ordinal_block("All CES, by driver",  "Driver_change", d),
  ordinal_block("All drivers, by CES", "Service",       d),
  ordinal_block("Pollution, by CES",   "Service", subset(d, Driver_change == "Pollution")),
  ordinal_block("Management, by CES",  "Service", subset(d, Driver_change == "Management")),
  ordinal_block("LSUC, by CES",        "Service", subset(d, Driver_change == "LSUC")))
write.csv(S11, "S11_final.csv", row.names = FALSE)
cat("\n### TABLE S11 (published spec A vs study-clustered spec B)\n")
print(S11[, c("model","level","studies","k","A","B","C","verdictA","verdictB")], row.names = FALSE)

## counts for the figures
cnt <- rbind(
  data.frame(grouping = "driver",  level = as.character(d$Driver_change), Vote = d$Vote),
  data.frame(grouping = "service", level = as.character(d$Service),       Vote = d$Vote))
write.csv(as.data.frame(table(cnt$grouping, cnt$level, cnt$Vote)), "counts_final.csv", row.names = FALSE)

################################ meta-analysis ################################
dm <- d[!is.na(d$yi) & !is.na(d$Value_type), ]
families <- c("% annual change", "% Total change", "-5 to 5; value of change",
              "Annual change", "total change")
families <- intersect(families, unique(dm$Value_type))
schemes <- list(
  used         = function(y) rep(max(0, 1, abs(y * 0.10)), length(y)),
  proportional = function(y) pmax(0.001, abs(y * 0.10)),
  half         = function(y) rep(max(0, 1, abs(y * 0.10)) * 0.5, length(y)),
  double       = function(y) rep(max(0, 1, abs(y * 0.10)) * 2.0, length(y)))

res <- list()
for (fam in c("% annual change", "% Total change", "-5 to 5; value of change")) {
  D <- dm[dm$Value_type == fam, ]
  cells <- subset(as.data.frame(table(D$Service, D$Driver_change, D$measurement_prediction)),
                  Freq > 1)
  for (r in seq_len(nrow(cells))) {
    z <- D[D$Service == as.character(cells$Var1[r]) &
           D$Driver_change == as.character(cells$Var2[r]) &
           D$measurement_prediction == as.character(cells$Var3[r]), ]
    row <- list(family = fam, service = as.character(cells$Var1[r]),
                driver = as.character(cells$Var2[r]), type = as.character(cells$Var3[r]),
                k = nrow(z), studies = length(unique(z$study)))
    for (s in names(schemes)) {
      z$vi <- schemes[[s]](z$yi); z$precision <- as.numeric(z$validity)
      m <- try(suppressWarnings(rma.mv(yi, vi, W = precision, random = ~1 | study,
                                       data = z, method = "REML")), silent = TRUE)
      if (inherits(m, "try-error")) { row[[paste0(s, "_sig")]] <- "err" } else {
        row[[paste0(s, "_sig")]] <- ifelse(m$pval < 0.05,
          ifelse(m$beta > 0, "POS*", "NEG*"), "ns")
        if (s == "used") { row$beta <- as.numeric(m$beta); row$se <- m$se
          row$ci.lb <- m$ci.lb; row$ci.ub <- m$ci.ub; row$pval <- m$pval } } }
    res[[length(res) + 1]] <- as.data.frame(row, stringsAsFactors = FALSE) } }
S12 <- do.call(rbind, res)
S12$stable <- ifelse(S12$used_sig == S12$proportional_sig &
                     S12$used_sig == S12$half_sig & S12$used_sig == S12$double_sig,
                     "stable", "changes")
write.csv(S12, "S12_final.csv", row.names = FALSE)
cat("\n### TABLE S12 / FIGURE 2b\n")
print(S12[, c("family","service","driver","type","k","studies","beta","pval",
              "used_sig","proportional_sig","half_sig","double_sig","stable")],
      row.names = FALSE, digits = 4)
cat("\nmodel fits:", nrow(S12),
    "| distinct CES-driver combinations:", nrow(unique(S12[, c("service","driver")])),
    "| involving pollution:", sum(!duplicated(S12[, c("service","driver")]) & S12$driver == "Pollution"),
    "| significant:", sum(grepl("\\*", S12$used_sig)), "\n")
