# scoring.R — ESPN scoring system & evaluation metrics

source("R/utils.R")

# =============================================================================
# ESPN SCORING
# =============================================================================

precompute_actual <- function(data, target_year) {
  actual <- data %>% filter(year == target_year)
  round_of_map <- c(64, 32, 16, 8, 4, 2)
  lapply(1:6, function(rnd) {
    sort(actual %>% filter(round_of == round_of_map[rnd]) %>% pull(winning_team_seed))
  })
}

count_matches_fast <- function(sim_seeds, actual_seeds) {
  sim_t <- tabulate(sim_seeds, nbins = 16)
  act_t <- tabulate(actual_seeds, nbins = 16)
  sum(pmin(sim_t, act_t))
}

score_bracket_fast <- function(bracket, actual_rounds) {
  score <- 0L
  for (rnd in 1:3) {
    field <- c("r1", "r2", "r3")[rnd]
    sim <- unlist(lapply(bracket$regions, function(r) r[[field]]))
    score <- score + count_matches_fast(sim, actual_rounds[[rnd]]) * ESPN_POINTS[rnd]
  }
  sim_r4 <- sapply(bracket$regions, function(r) r$r4)
  score <- score + count_matches_fast(sim_r4, actual_rounds[[4]]) * ESPN_POINTS[4]
  score <- score + count_matches_fast(bracket$f4_winners, actual_rounds[[5]]) * ESPN_POINTS[5]
  score <- score + count_matches_fast(bracket$champion, actual_rounds[[6]]) * ESPN_POINTS[6]
  score
}

# =============================================================================
# HELPERS
# =============================================================================

compute_upset_rates <- function(data, train_years) {
  train <- data %>% filter(year %in% train_years)
  rates <- train %>%
    mutate(upset = (winning_team_seed > losing_team_seed)) %>%
    group_by(round) %>%
    summarise(upset_rate = mean(upset), .groups = "drop") %>%
    arrange(round)
  rates$upset_rate
}

compute_sampling_params <- function(data, train_years) {
  train_data <- data %>% filter(year %in% train_years)
  q_f4  <- 1 / mean(train_data %>% filter(round == 4) %>% pull(winning_team_seed))
  q_e8  <- 1 / mean(train_data %>% filter(round == 3) %>% pull(winning_team_seed))
  q_ncg <- 1 / mean(train_data %>% filter(round == 5) %>% pull(winning_team_seed))
  s16_winners <- train_data %>% filter(round == 3) %>% pull(winning_team_seed)
  q_top <- 1 / mean(pmax(s16_winners[s16_winners <= 5], 1))
  q_bot <- 1 / mean(pmax(s16_winners[s16_winners >= 2], 1))
  list(q_f4 = q_f4, q_e8 = q_e8, q_ncg = q_ncg, q_top = q_top, q_bot = q_bot)
}

# =============================================================================
# EXPERIMENT RUNNER
# =============================================================================

# Run experiment: test multiple models × generators for a single year
run_experiment <- function(data, test_year, N = 1000, R = 10,
                           model_names = c("power", "recency", "mov", "smoothed"),
                           gen_names = c("R64", "E8", "F4A", "NCG", "PickFav")) {
  cat("Year:", test_year, "\n")

  train_years <- setdiff(1985:(test_year - 1), 2020)
  train_years <- train_years[train_years %in% unique(data$year)]

  # Build all probability matrices
  prob_mats <- list()
  for (mn in model_names) {
    mm <- build_model_and_matrix(mn, data, train_years, target_year = test_year)
    prob_mats[[mn]] <- mm$prob_mat
  }

  actual <- precompute_actual(data, test_year)
  sp <- compute_sampling_params(data, train_years)
  upset_rates <- compute_upset_rates(data, train_years)

  # For each model × generator combo, run N*R brackets
  results <- list()

  for (mn in model_names) {
    pm <- prob_mats[[mn]]

    # Define generators for this probability matrix
    gens <- list(
      R64     = function() gen_R64(pm),
      E8      = function() gen_E8(pm, sp$q_e8),
      F4A     = function() gen_F4A(pm, sp$q_f4),
      F4B     = function() gen_F4B(pm, sp$q_top, sp$q_bot),
      NCG     = function() gen_NCG(pm, sp$q_ncg),
      PickFav = function() gen_PF(),
      UpsetInject = function() gen_UpsetInject(upset_rates)
    )
    gens <- gens[intersect(gen_names, names(gens))]

    for (gn in names(gens)) {
      t0 <- Sys.time()
      gen_fn <- gens[[gn]]

      rep_max_scores <- numeric(R)
      all_scores <- numeric(N * R)

      for (rep in 1:R) {
        offset <- (rep - 1) * N
        best <- 0L
        for (i in 1:N) {
          b <- gen_fn()
          s <- score_bracket_fast(b, actual)
          all_scores[offset + i] <- s
          if (s > best) best <- s
        }
        rep_max_scores[rep] <- best
      }

      elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
      cat("  ", mn, "+", gn, ":", elapsed, "sec\n")

      results[[paste(mn, gn, sep = "_")]] <- tibble(
        model          = mn,
        generator      = gn,
        year           = test_year,
        max_score_mean = mean(rep_max_scores),
        max_score_sd   = sd(rep_max_scores),
        mean_score     = mean(all_scores),
        median_score   = median(all_scores),
        scores         = list(all_scores),
        max_scores     = list(rep_max_scores)
      )
    }
  }

  bind_rows(results)
}
