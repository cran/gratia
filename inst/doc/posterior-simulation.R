## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-packages------------------------------------------------------------
pkgs <- c(
  "mgcv", "gratia", "dplyr", "tidyr", "ggplot2", "ggdist",
  "distributional", "tibble", "withr", "patchwork"
)
vapply(pkgs, library, logical(1L), logical.return = TRUE, character.only = TRUE)

## ----ss-data------------------------------------------------------------------
ss_df <- data_sim("eg1", seed = 42)
m_ss <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = ss_df, method = "REML")

## ----ss-sm-coefs--------------------------------------------------------------
s_x0 <- get_smooth(m_ss, "s(x0)")
smooth_coef_indices(s_x0)

## ----sample-smooths-m-ss------------------------------------------------------
sm_samp <- smooth_samples(m_ss, select = "s(x0)", n_vals = 100, n = 100,
  seed = 21)

## ----draw-smooth-samples-x0---------------------------------------------------
sm_samp |>
  draw(alpha = 0.3)

## ----plot-smooth-ci-and-posterior-smooths-------------------------------------
# evaluate the fitted smooth over x0 and add on a credible interval
sm_est <- smooth_estimates(m_ss, select = "s(x0)") |>
  add_confint()
# plot the smooth, credible interval, and posterior smooths
sm_est |>
  ggplot(aes(x = x0)) +
  geom_lineribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
    orientation = "vertical", fill = "#56B4E9", alpha = 0.5
  ) +
  geom_line(
    data = sm_samp,
    aes(y = .value, group = .draw), alpha = 0.2
  ) +
  geom_line(aes(y = .estimate), linewidth = 1, colour = "#E69F00") +
  labs(y = smooth_label(s_x0))

## ----luo-wabha-function-sim-plot----------------------------------------------
f <- function(x) {
  sin(2 * ((4 * x) - 2)) + (2 * exp(-256 * (x - 0.5)^2))
}
df <- data_sim("lwf6", dist = "normal", scale = 0.3, seed = 2)
plt <- df |>
  ggplot(aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_function(fun = f)
plt

## ----luo-wabha-fit-gam--------------------------------------------------------
m <- gam(y ~ s(x, k = 25, bs = "ad"), data = df, method = "REML")

## ----luo-wabha-new-data-------------------------------------------------------
new_df <- data_slice(m, x = evenly(x, lower = 0, upper = 1, n = 200)) |>
  mutate(.row = row_number())

## ----luo-wabha-fitted---------------------------------------------------------
fv <- fitted_values(m, data = new_df)

## ----luo-wabha-posterior------------------------------------------------------
fs <- fitted_samples(m, data = new_df, n = 10, seed = 4) |>
  left_join(new_df |> select(.row, x), by = join_by(.row == .row))

## ----luo-wabha-plot-samples---------------------------------------------------
plt +
  geom_ribbon(data = fv, aes(y = .fitted, ymin = .lower_ci, ymax = .upper_ci),
    fill = "red", alpha = 0.3) +
  geom_line(data = fs, aes(group = .draw, x = x, y = .fitted),
    colour = "yellow", alpha = 0.4)

## -----------------------------------------------------------------------------
df <- data_sim("gwf2", n = 400, scale = 1, dist = "normal", seed = 8)

## ----simulated-data-gaussian-gam----------------------------------------------
df |>
  ggplot(aes(x = x, y = y)) +
  geom_point() +
  geom_function(fun = gw_f2, colour = "#0072B2", linewidth = 1.5)

## ----fit-gaussian-gam---------------------------------------------------------
m <- gam(y ~ s(x), data = df, method = "REML", family = gaussian())

## ----predict-new-x-gaussian-gam-----------------------------------------------
mu <- predict(m, newdata = data.frame(x = 0.5))
mu

## ----scale-parameter-gaussian-gam---------------------------------------------
sigma <- m$scale
sigma

## -----------------------------------------------------------------------------
df |>
  ggplot(aes(x = x, y = y)) +
  stat_halfeye(aes(ydist = dist_normal(mean = mu, sd = sigma)),
    x = 0.5, scale = 0.2, slab_fill = "#E69F00", slab_alpha = 0.7
  ) +
  geom_point() +
  geom_function(fun = gw_f2, colour = "#0072B2", linewidth = 1.5) +
  geom_point(x = 0.5, y = mu, colour = "red")

## ----expected-value-uncertainty-gaussian-gam----------------------------------
fitted_values(m, data = data.frame(x = 0.5))

## ----data-slice-gaussian-gam--------------------------------------------------
ds <- data_slice(m, x = evenly(x, n = 200)) |>
  mutate(.row = row_number())

## ----fitted-values-gaussian-gam-----------------------------------------------
fv <- fitted_values(m, data = ds)

## ----posterior-samples-gaussian-gam-------------------------------------------
ps <- posterior_samples(m, n = 10000, data = ds, seed = 24,
  unconditional = TRUE) |>
  left_join(ds, by = join_by(.row == .row))
ps

## ----quantile-fun-------------------------------------------------------------
quantile_fun <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
  tibble::tibble(
    .value = quantile(x, probs = probs, ...),
    .q = probs * 100
  )
}

## ----process-posterior-samples-gaussian-gam-----------------------------------
p_int <- ps |>
  group_by(.row) |>
  reframe(quantile_fun(.response)) |>
  pivot_wider(
    id_cols = .row, names_from = .q, values_from = .value,
    names_prefix = ".q"
  ) |>
  left_join(ds, by = join_by(.row == .row))
p_int

## ----plot-prediction-intervals-gaussian-gam-----------------------------------
fv |>
  ggplot(aes(x = x, y = .fitted)) +
  # summarise the posterior samples
  geom_hex(
    data = ps, aes(x = x, y = .response, fill = after_stat(count)),
    bins = 50, alpha = 0.7
  ) +
  # add the lower and upper prediction intervals
  geom_line(data = p_int, aes(y = .q2.5), colour = "#56B4E9",
    linewidth = 1.5) +
  geom_line(data = p_int, aes(y = .q97.5), colour = "#56B4E9",
    linewidth = 1.5) +
  # add the lower and upper credible intervals
  geom_line(aes(y = .lower_ci), colour = "#56B4E9", linewidth = 1) +
  geom_line(aes(y = .upper_ci), colour = "#56B4E9", linewidth = 1) +
  # add the fitted model
  geom_line() +
  # add the observed data
  geom_point(data = df, aes(x = x, y = y)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(legend.position = "none") +
  labs(y = "Response")

## ----mh-gam-data--------------------------------------------------------------
ga_fail <- function(seed) {
  df <- tibble(y = c(
    rep(0, 89), 1, 0, 1, 0, 0, 1, rep(0, 13), 1, 0, 0, 1,
    rep(0, 10), 1, 0, 0, 1, 1, 0, 1, rep(0, 4), 1, rep(0, 3),
    1, rep(0, 3), 1, rep(0, 10), 1, rep(0, 4), 1, 0, 1, 0, 0,
    rep(1, 4), 0, rep(1, 5), rep(0, 4), 1, 1, rep(0, 46)
  )) |>
    mutate(
      x = withr::with_seed(
        seed,
        sort(c(0:10 * 5, rnorm(length(y) - 11) * 20 + 100))
      ),
      .row = row_number()
    ) |>
    relocate(.row, .before = 1L)
  df
}

## ----plot-ga-fail-data--------------------------------------------------------
df <- ga_fail(3)

df |>
  ggplot(aes(x = x, y = y)) +
  geom_point()

## ----fit-ga-fail-gam----------------------------------------------------------
m_logit <- gam(y ~ s(x, k = 15), data = df, method = "REML", family = binomial)

## -----------------------------------------------------------------------------
fs_ga <- fitted_samples(m_logit, n = 2000, seed = 2)
fs_mh <- fitted_samples(m_logit,
  n = 2000, seed = 2, method = "mh", thin = 2,
  rw_scale = 0.4
)

## ----ga-fail-intervals--------------------------------------------------------
int_ga <- fs_ga |>
  group_by(.row) |>
  median_qi(.width = c(0.5, 0.8, 0.95)) |>
  left_join(df, by = join_by(.row == .row))

int_mh <- fs_mh |>
  group_by(.row) |>
  median_qi(.width = c(0.5, 0.8, 0.95)) |>
  left_join(df, by = join_by(.row == .row))

## ----ga-fail-plot-plus-intervals----------------------------------------------
plt_ga <- df |>
  ggplot(aes(x = x, y = y)) +
  geom_point() +
  geom_lineribbon(
    data = int_ga,
    aes(x = x, y = .fitted, ymin = .lower, ymax = .upper)
  ) +
  scale_fill_brewer() +
  labs(title = "Gaussian approximation")

plt_mh <- df |>
  ggplot(aes(x = x, y = y)) +
  geom_point() +
  geom_lineribbon(
    data = int_mh,
    aes(x = x, y = .fitted, ymin = .lower, ymax = .upper)
  ) +
  scale_fill_brewer() +
  labs(title = "Metropolis Hastings sampler")

plt_ga + plt_mh + plot_layout(guides = "collect")

