# models.R — Probability models for NCAA bracket generation
#
# 1. Power Model (Ludden et al. 2020) — the paper's original, used as baseline
# 2. AdjEM Model — our improvement: logistic regression on team efficiency data,
#    using region-specific team AdjEM values for each tournament year

source("R/utils.R")

# =============================================================================
# SHARED HELPERS
# =============================================================================

compute_alpha <- function(p_bar, s1, s2) {
  if (s1 == s2) return(0)
  if (p_bar >= 1) return(2)
  if (p_bar <= 0) return(-2)
  logit_p <- log(p_bar / (1 - p_bar))
  alpha <- logit_p / (log(1/s1) - log(1/s2))
  clamp(alpha, -2, 2)
}

power_model_prob <- function(alpha, s1, s2) {
  if (s1 == s2) return(0.5)
  s2^alpha / (s1^alpha + s2^alpha)
}

build_prob_matrix_from_alphas <- function(r1_alphas, round_alphas) {
  mat <- array(0.5, dim = c(16, 16, 6))

  round_alpha_lookup <- rep(1, 6)
  for (i in seq_len(nrow(round_alphas))) {
    round_alpha_lookup[round_alphas$round[i]] <- round_alphas$alpha_bar[i]
  }

  default_r1_alpha <- if (1 %in% round_alphas$round) {
    round_alphas %>% filter(round == 1) %>% pull(alpha_bar)
  } else { 1.0 }
  if (length(default_r1_alpha) == 0) default_r1_alpha <- 1.0

  for (s1 in 1:15) {
    for (s2 in (s1+1):16) {
      row <- r1_alphas %>% filter(.data$s1 == !!s1, .data$s2 == !!s2)
      alpha <- if (nrow(row) > 0) row$alpha[1] else default_r1_alpha
      p <- power_model_prob(alpha, s1, s2)
      mat[s1, s2, 1] <- p
      mat[s2, s1, 1] <- 1 - p

      for (rnd in 2:6) {
        p <- power_model_prob(round_alpha_lookup[rnd], s1, s2)
        mat[s1, s2, rnd] <- p
        mat[s2, s1, rnd] <- 1 - p
      }
    }
  }
  mat
}

# =============================================================================
# 1. POWER MODEL — Original (Ludden et al. 2020)
#    Round 1: matchup-specific alphas derived from observed win rates.
#    Rounds 2-6: one pooled alpha per round.
#    As we show in the critique: R1 predictions = observed frequencies exactly.
# =============================================================================

build_power_model <- function(data, train_years) {
  train <- data %>%
    filter(year %in% train_years) %>%
    mutate(
      s1 = pmin(winning_team_seed, losing_team_seed),
      s2 = pmax(winning_team_seed, losing_team_seed),
      s1_won = (winning_team_seed == s1)
    )

  r1 <- train %>%
    filter(round == 1) %>%
    group_by(s1, s2) %>%
    summarise(p_bar = mean(s1_won), n = n(), .groups = "drop") %>%
    rowwise() %>%
    mutate(alpha = compute_alpha(p_bar, s1, s2)) %>%
    ungroup()

  later <- train %>%
    filter(round >= 2, s1 != s2) %>%
    group_by(round, s1, s2) %>%
    summarise(p_bar = mean(s1_won), n = n(), .groups = "drop") %>%
    rowwise() %>%
    mutate(alpha = compute_alpha(p_bar, s1, s2)) %>%
    ungroup()

  round_alphas <- later %>%
    group_by(round) %>%
    summarise(alpha_bar = weighted.mean(alpha, w = n), .groups = "drop")

  list(r1_alphas = r1, round_alphas = round_alphas, name = "Power Model")
}

build_power_prob_matrix <- function(model) {
  build_prob_matrix_from_alphas(model$r1_alphas, model$round_alphas)
}

# =============================================================================
# 2. AdjEM MODEL — Region-specific team efficiency logistic regression
#
#    Uses KenPom-style Adjusted Efficiency Margin (AdjEM = AdjOE - AdjDE)
#    for the ACTUAL teams in each region of a given tournament year.
#
#    Model: P(s1 wins) ~ em_diff * round + log(s2/s1) + seed_diff
#
#    Key difference from Power Model:
#    - Power treats all 1-seeds identically across all years and regions
#    - AdjEM knows that UConn (AdjEM=31) as a 1-seed is different from
#      Purdue (AdjEM=24) as a 1-seed, AND puts them in separate regions
#
#    This is what breaks through the seed-only performance ceiling.
# =============================================================================

# Train the AdjEM logistic model on historical games
train_adjem_model <- function(data, train_years) {
  train <- data %>%
    filter(year %in% train_years, !is.na(win_AdjEM), !is.na(lose_AdjEM)) %>%
    mutate(
      s1 = pmin(winning_team_seed, losing_team_seed),
      s2 = pmax(winning_team_seed, losing_team_seed),
      s1_won = as.integer(winning_team_seed == s1),
      em_diff = ifelse(s1_won == 1, win_AdjEM - lose_AdjEM, lose_AdjEM - win_AdjEM),
      log_seed_ratio = log(s2 / s1),
      seed_diff = s2 - s1
    ) %>%
    filter(s1 != s2)

  glm(s1_won ~ em_diff * round + log_seed_ratio + seed_diff,
      data = train, family = binomial())
}

# Build a probability matrix for one region using its actual team AdjEM values
build_region_prob_matrix <- function(fit, region_teams) {
  mat <- array(0.5, dim = c(16, 16, 6))
  for (s1 in 1:15) {
    for (s2 in (s1+1):16) {
      e1 <- region_teams$AdjEM[region_teams$Seed == s1]
      e2 <- region_teams$AdjEM[region_teams$Seed == s2]
      if (length(e1) == 0 || length(e2) == 0 || is.na(e1[1]) || is.na(e2[1])) next
      for (rnd in 1:6) {
        nd <- data.frame(em_diff = e1[1] - e2[1], round = rnd,
                         log_seed_ratio = log(s2 / s1), seed_diff = s2 - s1)
        p <- predict(fit, newdata = nd, type = "response")
        mat[s1, s2, rnd] <- clamp(p, 0.01, 0.99)
        mat[s2, s1, rnd] <- 1 - mat[s1, s2, rnd]
      }
    }
  }
  mat
}

# Get the 4 tournament regions and their teams from metrics data
get_tournament_regions <- function(metrics, yr) {
  tourney_regions <- c("East", "West", "Midwest", "South", "Southeast", "Southwest")

  teams <- metrics %>%
    filter(Season == yr, Region %in% tourney_regions, !is.na(Seed)) %>%
    mutate(Seed = as.integer(Seed), AdjEM = AdjOE - AdjDE) %>%
    filter(!is.na(Seed)) %>%
    select(team = `Mapped ESPN Team Name`, Seed, Region, AdjEM)

  regions <- teams %>% distinct(Region) %>% pull(Region)
  lapply(regions[1:min(4, length(regions))], function(reg) {
    teams %>% filter(Region == reg)
  })
}

# Build all 4 region-specific probability matrices for a tournament year
build_adjem_region_matrices <- function(data, metrics, train_years, test_year) {
  fit <- train_adjem_model(data, train_years)
  region_teams <- get_tournament_regions(metrics, test_year)
  lapply(region_teams, function(rt) build_region_prob_matrix(fit, rt))
}

# Region-aware bracket simulation (each region uses its own probability matrix)
simulate_bracket_regionaware <- function(region_mats) {
  reg <- lapply(1:4, function(i) simulate_region(region_mats[[i]]))
  f4 <- c(reg[[1]]$r4, reg[[2]]$r4, reg[[3]]$r4, reg[[4]]$r4)
  # F4 games: blend matrices of the two competing regions
  f4_mat12 <- (region_mats[[1]] + region_mats[[2]]) / 2
  f4_mat34 <- (region_mats[[3]] + region_mats[[4]]) / 2
  f4w1 <- sim_game(f4[1], f4[2], 5L, f4_mat12)
  f4w2 <- sim_game(f4[3], f4[4], 5L, f4_mat34)
  # Championship: blend all 4 region matrices
  champ_mat <- (region_mats[[1]] + region_mats[[2]] + region_mats[[3]] + region_mats[[4]]) / 4
  champ <- sim_game(f4w1, f4w2, 6L, champ_mat)
  list(regions = reg, f4_seeds = f4, f4_winners = c(f4w1, f4w2), champion = champ)
}
