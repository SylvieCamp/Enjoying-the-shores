## Figure 2 — original design, corrected values.
## a) one boxed block per driver: services, counts, diverging proportion bars and
##    the probability of a positive effect on the shared bottom axis.
## b) mean effects from the meta-analytical models, printed by metric family and
##    by measurement / prediction, (original layout).
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)
                  library(grid); library(png)})

GREEN <- "#8CCAA0"; RED <- "#EB8C8D"; NEU <- "#D2D2D2"; AMB <- "#8C8C8C"
BOX <- "#EDEDED"; INK <- "#1a1a1a"
GTXT <- "#2E7D4F"; RTXT <- "#C0392B"
BORDER <- c("Climate-change"="#C1622A","Direct-exploitation"="#2E5FA3",
            "LSUC"="#E0AE28","Pollution"="#3E7A45","Coastal-erosion"="#6E6E6E",
            "Management"="#A64BA6")
TITLE <- c("Climate-change"="Climate change","Direct-exploitation"="Direct Exploitation",
           "LSUC"="Land/ sea use change","Pollution"="Pollution",
           "Coastal-erosion"="Coastal erosion","Management"="Management")
ORD <- names(BORDER)

d <- read.csv("raw.csv")
d$Vote <- factor(tolower(trimws(d$Vote_counting)),
                 levels = c("negative","neutral","ambiguous","positive"))
cells <- read.csv("cells_final.csv", stringsAsFactors = FALSE) %>% select(driver, service, k, studies, pA, loA, hiA)
bayes <- read.csv("bayes_blocks.csv", stringsAsFactors = FALSE) %>%
  transmute(driver, service = level, pB = pBayes, loB = loBayes, hiB = hiBayes)
cells <- full_join(cells, bayes, by = c("driver","service"))
s12 <- read.csv("S12_final.csv", stringsAsFactors = FALSE)

cnt <- d %>% count(driver = Driver_change, service = Service, Vote, name = "n") %>%
  tidyr::pivot_wider(names_from = Vote, values_from = n, values_fill = 0)
for (v in c("negative","neutral","ambiguous","positive"))
  if (!v %in% names(cnt)) cnt[[v]] <- 0
cnt <- cnt %>% mutate(tot = negative + neutral + ambiguous + positive,
                      p_neg = negative/tot, p_neu = neutral/tot,
                      p_amb = ambiguous/tot, p_pos = positive/tot)
X <- left_join(cnt, cells, by = c("driver","service")) %>%
  mutate(lab = ifelse(service == "Cultural", "Non-specified CES", service),
         icon = paste0("icons/a_", service, ".png"))

tx <- function(p) 2*p - 1
FAM <- c("% Total change","% annual change","-5 to 5; value of change")
XPOS <- c(0.14, 0.28, 0.47, 0.61, 0.80, 0.94)   # M,P for each of the three families

block <- function(dr, axis_bottom) {
  x <- X %>% filter(driver == dr) %>% arrange(is.na(pB), pB, desc(p_neg))
  n <- nrow(x); x$y <- rev(seq_len(n))
  seg <- rbind(
    data.frame(y=x$y, xmin=0, xmax=x$p_pos, f="positive"),
    data.frame(y=x$y, xmin=-x$p_neu, xmax=0, f="neutral"),
    data.frame(y=x$y, xmin=-(x$p_neu+x$p_amb), xmax=-x$p_neu, f="ambiguous"),
    data.frame(y=x$y, xmin=-(x$p_neu+x$p_amb+x$p_neg), xmax=-(x$p_neu+x$p_amb), f="negative"))
  seg$f <- factor(seg$f, levels=c("negative","ambiguous","neutral","positive"))

  left <- ggplot(x) +
    annotate("rect", xmin=.62, xmax=1.0, ymin=.4, ymax=n+.6, fill=BOX, colour=NA) +
    geom_text(aes(.60, y, label = lab), hjust=1, size=2.85, colour=INK) +
    geom_text(aes(.74, y, label = studies), size=2.8, colour=INK) +
    geom_text(aes(.92, y, label = k),       size=2.8, colour=INK) +
    scale_x_continuous(limits=c(.02,1.0), expand=c(0,0)) +
    scale_y_continuous(limits=c(.4,n+.6), expand=c(0,0)) +
    theme_void() + theme(plot.margin=margin(2,0,2,2))

  mid <- ggplot() +
    geom_rect(data=seg, aes(xmin=xmin, xmax=xmax, ymin=y-.34, ymax=y+.34, fill=f)) +
    geom_vline(xintercept=0, colour="black", linewidth=.45) +
    geom_segment(data=x, aes(x=tx(loB), xend=tx(hiB), y=y, yend=y),
                 colour="black", linewidth=.45, na.rm=TRUE) +
    geom_point(data=x, aes(tx(pB), y), shape=23, size=1.75, stroke=.45,
               colour="black", fill="black", na.rm=TRUE) +
    scale_fill_manual(values=c(negative=RED, ambiguous=AMB, neutral=NEU,
                               positive=GREEN), guide="none") +
    scale_x_continuous(limits=c(-1.02,1.02), expand=c(0,0),
      breaks=c(-1,-.5,0,.5,1), labels=c("0%","25%","50%","75%","100%")) +
    scale_y_continuous(limits=c(.4,n+.6), expand=c(0,0)) +
    labs(x=NULL, y=NULL) + theme_minimal(base_size=9) +
    theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_line(colour="#e6e6e6", linewidth=.3),
          axis.text.y=element_blank(),
          axis.text.x=if (axis_bottom) element_text(colour=INK, size=8) else element_blank(),
          plot.margin=margin(2,4,2,2),
          plot.background=element_rect(fill="white", colour=NA),
          panel.background=element_rect(fill="white", colour=NA))

  ic <- ggplot(x) + scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
    scale_y_continuous(limits=c(.4,n+.6), expand=c(0,0)) + theme_void() +
    theme(plot.margin=margin(2,0,2,2))
  for (i in seq_len(n)) {
    img <- png::readPNG(x$icon[i])
    ic <- ic + annotation_custom(rasterGrob(img, interpolate=TRUE),
                                 xmin=.10, xmax=.90,
                                 ymin=x$y[i]-.42, ymax=x$y[i]+.42)
  }

  e <- s12 %>% filter(driver == dr) %>%
    mutate(col = ifelse(beta > 0, GTXT, RTXT), sig = pval < 0.05,
           fi = match(family, FAM), mi = ifelse(type == "measurement", 1, 2),
           xp = XPOS[(fi - 1) * 2 + mi],
           yp = match(service, x$service),
           yp = ifelse(is.na(yp), 1, x$y[yp]))
  right <- ggplot(e) +
    geom_text(aes(xp, yp, label = sprintf("%.1f", beta), colour = I(col),
                  fontface = ifelse(sig, "bold", "plain")), size = 2.75) +
    scale_x_continuous(limits=c(0,1.04), expand=c(0,0)) +
    scale_y_continuous(limits=c(.4,n+.6), expand=c(0,0)) +
    theme_void() + theme(plot.margin=margin(2,2,2,2))

  row <- left + ic + mid + right +
    plot_layout(widths = c(1.35, .30, 3.1, 1.55), nrow = 1) +
    plot_annotation(theme = theme(
      plot.background = element_rect(fill = NA, colour = BORDER[dr], linewidth = .9),
      plot.margin = margin(3, 3, 3, 3)))
  wrap_elements(full = row)
}

hdr <- ggplot() + xlim(0,1) + ylim(0,1) + theme_void() +
  annotate("text", x=.30, y=.62, label="a)", size=3.4, fontface="bold", colour=INK) +
  annotate("text", x=.55, y=.62, label="Proportion of effect-sizes (bars)", size=3.2, colour=INK) +
  annotate("text", x=.80, y=.62, label="b)   Mean effect", size=3.2, fontface="bold", colour=INK) +
  annotate("text", x=.185, y=.16, label="Cultural services   Articles   Data",
           size=2.7, colour=INK)

## panel-b occupies the last 1.55 of widths c(1.35,.30,3.1,1.55) = 6.30 total
BX0 <- (1.35 + .30 + 3.1) / 6.30
BW  <- 1.55 / 6.30
ftr <- ggplot() + xlim(0,1) + ylim(0,1) + theme_void() +
  annotate("text", x=(1.35 + .30 + 3.1/2)/6.30, y=.84,
           label="Probability of having a positive effect (diamonds)", size=3.2, colour=INK) +
  annotate("segment", x=BX0 + BW*(XPOS-0.06), xend=BX0 + BW*(XPOS+0.06),
           y=.80, yend=.80, colour="#555555", linewidth=.3) +
  annotate("text", x=BX0 + BW*XPOS, y=.68, label=rep(c("M","P"),3), size=2.6, colour=INK) +
  annotate("text", x=BX0 + BW*c(.21,.54,.87), y=.52,
           label=c("% total\nchange","% annual\nchange","Semi-\nquantitative"),
           size=2.4, colour=INK, lineheight=.95, vjust=1) +
  annotate("text", x=BX0 + BW*.54, y=.08,
           label="M: measurement; P: prediction", size=2.3, colour="#555555")

titles <- lapply(ORD, function(dr)
  ggplot() + xlim(0,1) + ylim(0,1) + theme_void() +
    annotate("text", x=.02, y=.5, hjust=0, label=TITLE[dr], size=3.3,
             fontface="bold", colour=BORDER[dr]))

blocks <- list()
for (i in seq_along(ORD)) {
  blocks[[length(blocks)+1]] <- titles[[i]]
  blocks[[length(blocks)+1]] <- block(ORD[i], i == length(ORD))
}
hts <- as.numeric(rbind(rep(.55, length(ORD)),
                        sapply(ORD, function(dr) sum(X$driver == dr))))

fig <- wrap_plots(c(list(hdr), blocks, list(ftr)), ncol = 1,
                  heights = c(1.0, hts, 2.4)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("Figure2_final.png", fig, width = 9.6, height = 11.7, dpi = 340, bg = "white")
ggsave("Figure2_final.pdf", fig, width = 9.6, height = 11.7, bg = "white")
cat("written Figure2_final.png / .pdf\n")
