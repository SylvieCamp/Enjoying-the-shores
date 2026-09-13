## Figure 1 — original design, corrected values.
## Diverging stacked proportion bars (top axis) with the probability of a positive
## effect overlaid as diamonds on the bottom axis, (original layout).
## COMPARE = FALSE -> Figure 1 (primary estimates only).
## COMPARE = TRUE  -> Figure S2 (primary estimates beside the fixed-effect ones).
if (!exists("COMPARE")) COMPARE <- FALSE
OUT <- if (COMPARE) "FigureS2_final" else "Figure1_final"
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)
                  library(grid); library(png)})

GREEN <- "#8CCAA0"; RED <- "#EB8C8D"; NEU <- "#D2D2D2"; AMB <- "#8C8C8C"
BOX   <- "#EDEDED"; INK <- "#1a1a1a"

d <- read.csv("raw.csv")
d$Vote <- factor(tolower(trimws(d$Vote_counting)),
                 levels = c("negative","neutral","ambiguous","positive"))
d$study <- factor(d$biblio_internal_id)
s11 <- read.csv("S11_final.csv", stringsAsFactors = FALSE)
bb  <- read.csv("bayes_blocks.csv", stringsAsFactors = FALSE)

ORD_A <- c("Spiritual","Existence","Recreation","Cultural","Heritage","Aesthetic",
           "Symbolic","Scientific","Religion","Bequest","Cognitive")
ORD_B <- c("Climate-change","Direct-exploitation","LSUC","Pollution",
           "Coastal-erosion","Management")
LAB_A <- c(Cultural = "Non-specified CES")
LAB_B <- c("Climate-change"="Climate change","Direct-exploitation"="Direct exploitation",
           "LSUC"="Land/sea use change","Coastal-erosion"="Coastal erosion")

build <- function(col, ord, model_tag, labmap, prefix) {
  cnt <- d %>% count(level = .data[[col]], Vote, name = "n") %>%
    tidyr::pivot_wider(names_from = Vote, values_from = n, values_fill = 0)
  for (v in c("negative","neutral","ambiguous","positive"))
    if (!v %in% names(cnt)) cnt[[v]] <- 0
  cnt <- cnt %>% mutate(tot = negative + neutral + ambiguous + positive,
                        p_neg = negative/tot, p_neu = neutral/tot,
                        p_amb = ambiguous/tot, p_pos = positive/tot)
  est <- s11 %>% filter(model == model_tag) %>%
    select(level, studies, k, pA, loA, hiA, pB, loB, hiB)
  x <- left_join(cnt, est, by = "level") %>%
    filter(level %in% ord) %>%
    mutate(level = factor(level, levels = ord),
           lab = ifelse(as.character(level) %in% names(labmap),
                        labmap[as.character(level)], as.character(level)),
           icon = paste0("icons/", prefix, "_", as.character(level), ".png"))
  arrange(x, level)
}
A <- build("Service", ORD_A, "All drivers, by CES", LAB_A, "a")
B <- build("Driver_change", ORD_B, "All CES, by driver", LAB_B, "b")

## ------------------------------------------------------------------ bar panel
bars <- function(x, toplab, botlab) {
  n <- nrow(x); x$y <- rev(seq_len(n))
  seg <- rbind(
    data.frame(y = x$y, xmin = 0,                     xmax = x$p_pos,                          f = "positive"),
    data.frame(y = x$y, xmin = -x$p_neu,              xmax = 0,                                f = "neutral"),
    data.frame(y = x$y, xmin = -(x$p_neu + x$p_amb),  xmax = -x$p_neu,                         f = "ambiguous"),
    data.frame(y = x$y, xmin = -(x$p_neu + x$p_amb + x$p_neg), xmax = -(x$p_neu + x$p_amb),    f = "negative"))
  seg$f <- factor(seg$f, levels = c("negative","ambiguous","neutral","positive"))
  tx <- function(p) 2 * p - 1                       # probability -> bar-axis units
  OFF <- if (COMPARE) .17 else 0
  ggplot() +
    geom_rect(data = seg, aes(xmin = xmin, xmax = xmax, ymin = y - .34, ymax = y + .34,
                              fill = f), colour = NA) +
    geom_vline(xintercept = 0, colour = "black", linewidth = .5) +
    { if (COMPARE) list(
      geom_segment(data = x, aes(x = tx(loA), xend = tx(hiA), y = y - .17, yend = y - .17),
                   colour = "#7d7d7d", linewidth = .35, na.rm = TRUE),
      geom_point(data = x, aes(tx(pA), y - .17), shape = 23, size = 1.7, stroke = .55,
                 colour = "#5a5a5a", fill = "white", na.rm = TRUE)) } +
    geom_segment(data = x, aes(x = tx(loB), xend = tx(hiB), y = y + OFF, yend = y + OFF),
                 colour = "black", linewidth = .5, na.rm = TRUE) +
    geom_point(data = x, aes(tx(pB), y + OFF), shape = 23, size = 2.1, stroke = .5,
               colour = "black", fill = "black", na.rm = TRUE) +
    scale_fill_manual(values = c(negative = RED, ambiguous = AMB,
                                 neutral = NEU, positive = GREEN), guide = "none") +
    scale_x_continuous(limits = c(-1.02, 1.02), expand = c(0, 0),
      breaks = c(-1,-.5,0,.5,1), labels = c("100%","50%","0%","50%","100%"),
      position = "top",
      sec.axis = sec_axis(~ (. + 1)/2, breaks = c(0,.25,.5,.75,1),
                          labels = c("0%","25%","50%","75%","100%"), name = botlab)) +
    scale_y_continuous(limits = c(.4, n + 1.0), expand = c(0, 0)) +
    labs(x = toplab, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#e3e3e3", linewidth = .3),
          axis.text.y = element_blank(),
          axis.text.x = element_text(colour = INK, size = 9),
          axis.title.x = element_text(colour = INK, size = 9.5),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(4, 6, 4, 2))
}

## ------------------------------------------------------------- label + counts
side <- function(x, title, header) {
  n <- nrow(x); x$y <- rev(seq_len(n))
  g <- ggplot(x) +
    annotate("rect", xmin = .60, xmax = 1.00, ymin = .4, ymax = n + .55,
             fill = BOX, colour = NA) +
    geom_text(aes(.585, y, label = lab), hjust = 1, size = 3.2, colour = INK) +
    geom_text(aes(.79, y, label = studies), size = 3.1, colour = INK) +
    geom_text(aes(.93, y, label = k),       size = 3.1, colour = INK) +
    scale_x_continuous(limits = c(.10, 1.00), expand = c(0, 0)) +
    scale_y_continuous(limits = c(.4, n + 1.0), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(4, 0, 4, 4),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(colour = INK, face = "bold", size = 10.5,
                                    hjust = 0, margin = margin(b = 2)))
  for (i in seq_len(n)) {
    img <- png::readPNG(x$icon[i])
    g <- g + annotation_custom(rasterGrob(img, interpolate = TRUE),
                               xmin = .612, xmax = .688,
                               ymin = x$y[i] - .40, ymax = x$y[i] + .40)
  }
  if (header)
    g <- g + annotate("text", x = c(.79, .93), y = n + .95,
                      label = c("Number of\nArticles", "of\nData"),
                      size = 2.7, colour = INK, vjust = 1, lineheight = .95) +
             coord_cartesian(clip = "off")
  g + ggtitle(title)
}

pa <- side(A, "a) Cultural Services", TRUE) + bars(A, "Proportion of effect-sizes (bars)", NULL) +
      plot_layout(widths = c(1, 1.55))
pb <- side(B, "b) Drivers", FALSE) + bars(B, NULL, "Probability of having a positive effect (diamonds)") +
      plot_layout(widths = c(1, 1.55))

leg <- ggplot() + xlim(0, 10) + ylim(0, 1) +
  annotate("rect", xmin = c(0.15,1.35,2.55,3.85), xmax = c(0.55,1.75,2.95,4.25),
           ymin = .40, ymax = .68, fill = c(RED, AMB, NEU, GREEN)) +
  annotate("text", x = c(0.65,1.85,3.05,4.35), y = .54,
           label = c("negative","ambiguous","neutral","positive"),
           hjust = 0, size = 2.9, colour = INK) +
  { if (COMPARE) list(
    annotate("point", x = 5.75, y = .54, shape = 23, size = 1.7, stroke = .55,
             colour = "#5a5a5a", fill = "white"),
    annotate("segment", x = 5.45, xend = 6.05, y = .54, yend = .54,
             colour = "#7d7d7d", linewidth = .35),
    annotate("text", x = 6.15, y = .54, hjust = 0, size = 2.9, colour = INK,
             label = "effect sizes treated as independent")) } +
  annotate("segment", x = 5.45, xend = 6.05, y = if (COMPARE) .18 else .54, yend = if (COMPARE) .18 else .54,
           colour = "black", linewidth = .5) +
  annotate("point", x = 5.75, y = if (COMPARE) .18 else .54, shape = 23, size = 2.1, colour = "black",
           fill = "black") +
  annotate("text", x = 6.15, y = if (COMPARE) .18 else .54, hjust = 0, size = 2.9, colour = INK,
           label = if (COMPARE) "study-level random intercept (primary analysis)" else "probability of a positive effect, 95% confidence interval") +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA),
                       plot.margin = margin(0, 4, 2, 4))

fig <- pa / pb / leg +
  plot_layout(heights = c(11, 6.2, 1.5)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(paste0(OUT, ".png"), fig, width = 9.6, height = 7.6, dpi = 400, bg = "white")
ggsave(paste0(OUT, ".pdf"), fig, width = 9.6, height = 7.6, bg = "white")
cat("written", OUT, "\n")
