# generators.R — Bracket generation algorithms
# All generators from Ludden et al. 2020 + upset injection baseline.
# Any probability model's matrix can be plugged into any generator.

source("R/utils.R")

# =============================================================================
# CORE SIMULATION FUNCTIONS
# =============================================================================

sim_game <- function(seed_a, seed_b, rnd, prob_mat) {
  p_a <- prob_mat[seed_a, seed_b, rnd]
  if (runif(1) < p_a) seed_a else seed_b
}

simulate_region <- function(prob_mat) {
  r1 <- integer(8)
  for (i in 1:8) {
    r1[i] <- sim_game(REGION_SEEDS[i, 1], REGION_SEEDS[i, 2], 1L, prob_mat)
  }
  r2 <- integer(4)
  for (i in 1:4) r2[i] <- sim_game(r1[2*i - 1], r1[2*i], 2L, prob_mat)
  r3 <- integer(2)
  r3[1] <- sim_game(r2[1], r2[2], 3L, prob_mat)
  r3[2] <- sim_game(r2[3], r2[4], 3L, prob_mat)
  r4 <- sim_game(r3[1], r3[2], 4L, prob_mat)
  list(r1 = r1, r2 = r2, r3 = r3, r4 = r4)
}

simulate_bracket <- function(prob_mat) {
  reg <- lapply(1:4, function(i) simulate_region(prob_mat))
  f4 <- c(reg[[1]]$r4, reg[[2]]$r4, reg[[3]]$r4, reg[[4]]$r4)
  f4w <- c(sim_game(f4[1], f4[2], 5L, prob_mat),
           sim_game(f4[3], f4[4], 5L, prob_mat))
  champ <- sim_game(f4w[1], f4w[2], 6L, prob_mat)
  list(regions = reg, f4_seeds = f4, f4_winners = f4w, champion = champ)
}

# =============================================================================
# GENERATORS — all take a prob_mat and optional sampling params
# =============================================================================

gen_R64 <- function(prob_mat) simulate_bracket(prob_mat)

gen_PF <- function() {
  reg <- lapply(1:4, function(i) {
    r1 <- pmin(REGION_SEEDS[, 1], REGION_SEEDS[, 2])
    r2 <- integer(4)
    for (j in 1:4) r2[j] <- min(r1[2*j - 1], r1[2*j])
    r3 <- c(min(r2[1], r2[2]), min(r2[3], r2[4]))
    r4 <- min(r3[1], r3[2])
    list(r1 = r1, r2 = r2, r3 = r3, r4 = r4)
  })
  f4 <- sapply(reg, function(r) r$r4)
  f4w <- c(min(f4[1], f4[2]), min(f4[3], f4[4]))
  list(regions = reg, f4_seeds = f4, f4_winners = f4w, champion = min(f4w))
}

gen_E8 <- function(prob_mat, q_e8) {
  reg <- lapply(1:4, function(i) {
    r <- simulate_region(prob_mat)
    e8 <- sample_trunc_geom(2, q_e8)
    r$r3 <- e8
    r$r4 <- sim_game(e8[1], e8[2], 4L, prob_mat)
    r
  })
  f4 <- sapply(reg, function(r) r$r4)
  f4w <- c(sim_game(f4[1], f4[2], 5L, prob_mat),
           sim_game(f4[3], f4[4], 5L, prob_mat))
  champ <- sim_game(f4w[1], f4w[2], 6L, prob_mat)
  list(regions = reg, f4_seeds = f4, f4_winners = f4w, champion = champ)
}

gen_F4A <- function(prob_mat, q_f4) {
  f4 <- sample_trunc_geom(4, q_f4)
  reg <- lapply(1:4, function(i) {
    r <- simulate_region(prob_mat)
    r$r4 <- f4[i]
    r
  })
  f4w <- c(sim_game(f4[1], f4[2], 5L, prob_mat),
           sim_game(f4[3], f4[4], 5L, prob_mat))
  champ <- sim_game(f4w[1], f4w[2], 6L, prob_mat)
  list(regions = reg, f4_seeds = f4, f4_winners = f4w, champion = champ)
}

gen_F4B <- function(prob_mat, q_top, q_bot) {
  reg <- lapply(1:4, function(i) {
    r <- simulate_region(prob_mat)
    top_seed <- sample_trunc_geom(1, q_top, imax = 16)
    bot_seed <- sample_trunc_geom(1, q_bot, imax = 16)
    r$r3 <- c(top_seed, bot_seed)
    r$r4 <- sim_game(top_seed, bot_seed, 4L, prob_mat)
    r
  })
  f4 <- sapply(reg, function(r) r$r4)
  f4w <- c(sim_game(f4[1], f4[2], 5L, prob_mat),
           sim_game(f4[3], f4[4], 5L, prob_mat))
  champ <- sim_game(f4w[1], f4w[2], 6L, prob_mat)
  list(regions = reg, f4_seeds = f4, f4_winners = f4w, champion = champ)
}

gen_NCG <- function(prob_mat, q_ncg) {
  ncg <- sample_trunc_geom(2, q_ncg)
  reg <- lapply(1:4, function(i) simulate_region(prob_mat))
  f4 <- sapply(reg, function(r) r$r4)
  f4w <- ncg
  champ <- sim_game(ncg[1], ncg[2], 6L, prob_mat)
  list(regions = reg, f4_seeds = f4, f4_winners = f4w, champion = champ)
}

gen_UpsetInject <- function(upset_rates_by_round) {
  bracket <- gen_PF()
  for (reg_idx in 1:4) {
    for (i in 1:8) {
      if (runif(1) < upset_rates_by_round[1]) {
        bracket$regions[[reg_idx]]$r1[i] <- REGION_SEEDS[i, 2]
      }
    }
    r1 <- bracket$regions[[reg_idx]]$r1
    for (i in 1:4) {
      candidates <- c(r1[2*i - 1], r1[2*i])
      if (runif(1) < upset_rates_by_round[2]) {
        bracket$regions[[reg_idx]]$r2[i] <- max(candidates)
      } else {
        bracket$regions[[reg_idx]]$r2[i] <- min(candidates)
      }
    }
    r2 <- bracket$regions[[reg_idx]]$r2
    for (i in 1:2) {
      idx <- if (i == 1) c(1, 2) else c(3, 4)
      candidates <- r2[idx]
      if (runif(1) < upset_rates_by_round[3]) {
        bracket$regions[[reg_idx]]$r3[i] <- max(candidates)
      } else {
        bracket$regions[[reg_idx]]$r3[i] <- min(candidates)
      }
    }
    r3 <- bracket$regions[[reg_idx]]$r3
    if (runif(1) < upset_rates_by_round[4]) {
      bracket$regions[[reg_idx]]$r4 <- max(r3)
    } else {
      bracket$regions[[reg_idx]]$r4 <- min(r3)
    }
  }
  f4 <- sapply(bracket$regions, function(r) r$r4)
  bracket$f4_seeds <- f4
  for (i in 1:2) {
    pair <- f4[c(2*i - 1, 2*i)]
    if (runif(1) < upset_rates_by_round[5]) {
      bracket$f4_winners[i] <- max(pair)
    } else {
      bracket$f4_winners[i] <- min(pair)
    }
  }
  if (runif(1) < upset_rates_by_round[6]) {
    bracket$champion <- max(bracket$f4_winners)
  } else {
    bracket$champion <- min(bracket$f4_winners)
  }
  bracket
}
