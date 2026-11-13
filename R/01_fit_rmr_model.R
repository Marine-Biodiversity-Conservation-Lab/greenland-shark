# 01 - Interspecific resting metabolic rate model (panel B)
#
# Phylogenetic Bayesian regression of whole-animal resting metabolic power on
# body mass and measurement temperature across 26 shark species:
#
#   ln W_i = b0 + E * x_T,i + b_M * ln M_i + e,    e ~ MVN(0, sigma^2 V(lambda))
#
#   x_T = 1/(k T_ref) - 1/(k T_i)    Arrhenius (Boltzmann-Arrhenius) temperature,
#                                    k = 8.617333e-5 eV/K, T_ref = median
#                                    measurement temperature (294.69 K, 21.54 C).
#                                    Its coefficient E is an apparent activation
#                                    energy in eV.
#   V(lambda) = lambda * A + (1 - lambda) * I, A the phylogenetic correlation
#                                    matrix, lambda fixed at its profile-likelihood
#                                    maximum, 0.816.
#
# Lambda is FIXED rather than estimated: with one observation per species the
# phylogenetic and residual variances are not separately identifiable, so it was
# chosen by profile likelihood over a grid (support interval 0.381-0.980;
# lambda = 0 rejected, dAIC +9.18). The profile is not repeated here.
#
# Writes
#   outputs/fits/m2_rmr_lambda0816.rds   brms fit (git-ignored, ~MB, regenerable)
#   outputs/m2_posterior_draws.csv       fixed-effect draws, used by 02_figure.R
#   outputs/m2_summary.csv               posterior medians and 95% CrI

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ape)
  library(brms)
})

K_BOLTZ <- 8.617333e-5
LAMBDA  <- 0.816
SEED <- 20260803; ITER <- 4000; WARMUP <- 1000; CHAINS <- 4
options(mc.cores = min(CHAINS, max(1, parallel::detectCores() - 1)))

d <- read.csv(here("data", "rmr_interspecific_26spp.csv")) |>
  mutate(temp_k   = temp_c + 273.15,
         ln_watts = log(metabolic_power_w),
         ln_mass  = log(mass_kg))
T_REF <- median(d$temp_k)                     # 294.69 K
d <- d |> mutate(x_T = 1 / (K_BOLTZ * T_REF) - 1 / (K_BOLTZ * temp_k))

tree <- read.nexus(here("data", "phylogeny_26spp.nex"))
stopifnot(setequal(tree$tip.label, d$binomial))
A <- vcv.phylo(tree, corr = TRUE)
V <- LAMBDA * A + (1 - LAMBDA) * diag(nrow(A))
dimnames(V) <- dimnames(A)
V <- V[d$binomial, d$binomial]

dir.create(here("outputs", "fits"), showWarnings = FALSE, recursive = TRUE)

# Priors centred on metabolic theory (mass exponent 0.8, E 0.65 eV), weakly
# informative; the posterior is robust to them (prior sensitivity in the
# accompanying paper's supplement).
m2 <- brm(bf(ln_watts ~ x_T + ln_mass + fcor(V)),
          data = d, data2 = list(V = V), family = gaussian(),
          prior = c(prior(normal(0.8, 0.2),  class = "b", coef = "ln_mass"),
                    prior(normal(0.65, 0.3), class = "b", coef = "x_T"),
                    prior(student_t(3, 0, 1), class = "sigma")),
          chains = CHAINS, iter = ITER, warmup = WARMUP, seed = SEED,
          backend = "rstan", refresh = 0,
          file = here("outputs", "fits", "m2_rmr_lambda0816"),
          file_refit = "on_change")

draws <- posterior::as_draws_df(m2) |>
  as_tibble() |>
  select(.chain, .iteration, .draw, b_Intercept, b_x_T, b_ln_mass, sigma)
write.csv(draws, here("outputs", "m2_posterior_draws.csv"), row.names = FALSE)

summ <- draws |>
  pivot_longer(c(b_Intercept, b_x_T, b_ln_mass, sigma), names_to = "term") |>
  group_by(term) |>
  summarise(median = median(value), lo95 = quantile(value, 0.025),
            hi95 = quantile(value, 0.975), .groups = "drop") |>
  mutate(T_ref_K = T_REF, lambda = LAMBDA, n_species = nrow(d))
write.csv(summ, here("outputs", "m2_summary.csv"), row.names = FALSE)

print(summ, width = Inf)
np <- nuts_params(m2)
cat(sprintf("divergences %d | max Rhat %.4f\n",
            sum(subset(np, Parameter == "divergent__")$Value),
            max(rhat(m2), na.rm = TRUE)))
