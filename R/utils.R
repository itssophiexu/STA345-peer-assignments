# utils.R — Shared helper functions

if (!requireNamespace("tidyverse", quietly = TRUE)) stop("tidyverse required")
if (!exists(".utils_loaded", envir = globalenv())) {
  library(tidyverse)
  assign(".utils_loaded", TRUE, envir = globalenv())
}

# Clamp a value between lo and hi
clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))

# Truncated geometric sampler
# Returns n IID samples from TruncGeom(q, imax)
sample_trunc_geom <- function(n, q, imax = 16) {
  probs <- q * (1 - q)^(0:(imax - 1))
  probs <- probs / sum(probs)
  sample.int(imax, size = n, replace = TRUE, prob = probs)
}

# Region bracket pairings (standard NCAA bracket order)
REGION_SEEDS <- matrix(c(
  1, 16,  8, 9,  5, 12,  4, 13,
  6, 11,  3, 14, 7, 10,  2, 15
), ncol = 2, byrow = TRUE)

# ESPN scoring: points per round
ESPN_POINTS <- c(10L, 20L, 40L, 80L, 160L, 320L)
