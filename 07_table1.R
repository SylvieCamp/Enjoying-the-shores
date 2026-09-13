## Table 1 — composition of the evidence base (counts only; no model).
suppressMessages(library(dplyr))
d <- read.csv("raw.csv"); d$v <- tolower(trimws(d$Vote_counting))
row <- function(x, label) data.frame(
  label, studies = n_distinct(x$biblio_internal_id), k = nrow(x),
  measured = sum(x$measurement_prediction == "measurement"),
  predicted = sum(x$measurement_prediction == "prediction"),
  positive = sum(x$v == "positive"), neutral = sum(x$v == "neutral"),
  ambiguous = sum(x$v == "ambiguous"), negative = sum(x$v == "negative"))
DRV <- c("Climate-change","Coastal-erosion","Direct-exploitation","LSUC","Management","Pollution")
CES <- c("Aesthetic","Bequest","Cognitive","Existence","Heritage","Cultural","Recreation",
         "Religion","Scientific","Spiritual","Symbolic")
tab <- bind_rows(
  bind_rows(lapply(DRV, function(k) row(d[d$Driver_change == k, ], k))),
  bind_rows(lapply(CES, function(k) row(d[d$Service == k, ], k))),
  row(d, "All"))
write.csv(tab, "table1_counts.csv", row.names = FALSE); print(tab)
