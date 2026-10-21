# Internal helpers shared across the package

# Fast AUC (c-statistic) via the Wilcoxon rank-sum statistic, with ties
# counted as 1/2. Based on the implementation by Ben Van Calster. Counts are
# held as doubles so n1 * n2 cannot overflow; returns NA if either outcome
# group is empty.
.fast_auc <- function(p, y) {
  x1 <- p[y == 1]; n1 <- as.numeric(length(x1))
  x2 <- p[y == 0]; n2 <- as.numeric(length(x2))
  if (n1 == 0 || n2 == 0) return(NA_real_)
  r  <- rank(c(x1, x2))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}
