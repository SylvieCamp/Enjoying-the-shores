suppressMessages({library(ordinal); library(MASS)}); set.seed(1)
d <- read.csv("raw.csv"); d$study <- factor(d$biblio_internal_id)
d$Vote <- factor(tolower(trimws(d$Vote_counting)), levels=c("negative","neutral","ambiguous","positive"), ordered=TRUE)
d3 <- droplevels(subset(d, Vote != "ambiguous")); d3$Vote <- droplevels(d3$Vote)
gh <- function(n=60){i<-1:(n-1);a<-sqrt(i);J<-diag(0,n);J[cbind(i,i+1)]<-a;J[cbind(i+1,i)]<-a;e<-eigen(J,symmetric=TRUE);list(x=rev(e$values),w=rev(e$vectors[1,]^2))}
GH<-gh(60); pm<-function(th,eta,s) vapply(eta,function(e) sum(GH$w*(1-plogis(th-e-GH$x*s))),0)
for (D in list(d, d3)) {
  D$Driver_change <- factor(D$Driver_change); k <- tapply(D$study, D$Driver_change, function(z) length(unique(z)))
  D <- droplevels(D[D$Driver_change %in% names(k)[k>=2],]); nth <- nlevels(D$Vote)-1
  m <- clmm(Vote ~ Driver_change + (1|study), data=D, Hess=TRUE); p <- c(m$alpha,m$beta); s <- as.numeric(m$ST$study)
  est <- pm(p[nth], c(0,p[-(1:nth)]), s); names(est) <- levels(D$Driver_change)
  cat(sprintf("%-20s", ifelse(nlevels(D$Vote)==4,"with ambiguous","ambiguous excluded")), paste(sprintf("%s=%.2f", names(est), est), collapse="  "), "\n")
}
