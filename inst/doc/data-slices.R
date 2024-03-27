## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-packages, results = "hide"------------------------------------------
library("mgcv")
library("gratia")
library("dplyr")
library("ggplot2")
library("forcats")
library("datasets")

## ----load-co2-----------------------------------------------------------------
## data load and prep
data(CO2, package = "datasets")
plant <- CO2 |>
  as_tibble() |>
  rename(plant = Plant, type = Type, treatment = Treatment) |>
  mutate(plant = factor(plant, ordered = FALSE))

## ----plot-plant-data----------------------------------------------------------
plant_ylab <- expression(CO[2] ~ uptake ~ (mu * mol ~ m^{-3})) # nolint: brace_linter.
plant_xlab <- expression(CO[2] ~ concentration ~ (mL ~ L^{-1})) # nolint: brace_linter.

plant |>
  ggplot(aes(x = conc, y = uptake, group = plant, colour = treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~type) +
  labs(y = plant_ylab, x = plant_xlab, colour = "Treatment")

## ----fit-plant-gam------------------------------------------------------------
plant <- plant |>
  mutate(tt = fct_cross(treatment, type))
m_plant <- gam(uptake ~ treatment * type + s(conc, by = tt, k = 6) +
    s(plant, bs = "re"),
  data = plant, method = "REML", familly = Gamma(link = "log")
)
overview(m_plant)

## ----draw-plant-model---------------------------------------------------------
draw(m_plant, residuals = TRUE, scales = "fixed")

## ----plant-ds-1---------------------------------------------------------------
ds1 <- data_slice(m_plant,
  conc = evenly(conc, n = 100),
  type = level(type, "Quebec"), treatment = level(treatment, "chilled")
)
ds1

## ----plant-ds-1a--------------------------------------------------------------
ds1 <- data_slice(m_plant,
  conc = evenly(conc, n = 100),
  treatment = level(treatment, "chilled"), type = level(type, "Quebec"),
  tt = level(tt, "chilled:Quebec")
)
ds1

## ----plant-fitted-1-----------------------------------------------------------
fv1 <- fitted_values(m_plant, data = ds1, scale = "response", exclude = "s(plant)")
fv1

## ----plot-plant-fitted-1------------------------------------------------------
fv1 |>
  ggplot(aes(x = conc, y = .fitted)) +
  geom_point(
    data = plant |>
      filter(type == "Quebec", treatment == "chilled"),
    mapping = aes(y = uptake),
    alpha = 0.8, colour = "steelblue"
  ) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line() +
  labs(
    x = plant_xlab, y = plant_ylab,
    title = expression(Estimated ~ CO[2] ~ uptake),
    subtitle = "Chilled plants of the Quebec type"
  )

## ----plant-data-slice-2-------------------------------------------------------
ds2 <- data_slice(m_plant,
  conc = evenly(conc, n = 100),
  treatment = evenly(treatment), type = level(type, "Mississippi")
) |>
  mutate(tt = fct_cross(treatment, type, keep_empty = TRUE))
ds2

## ----draw-plant-fitted-values-2-----------------------------------------------
fitted_values(m_plant,
  data = ds2, scale = "response",
  exclude = "s(plant)"
) |>
  ggplot(aes(x = conc, y = .fitted, group = treatment)) +
  geom_point(
    data = plant |> filter(type == "Mississippi"),
    mapping = aes(y = uptake, colour = treatment),
    alpha = 0.8
  ) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = treatment),
    alpha = 0.2
  ) +
  geom_line(aes(colour = treatment)) +
  labs(
    x = plant_xlab, y = plant_ylab,
    title = expression(Estimated ~ CO[2] ~ uptake),
    subtitle = "Comparison of treatment in plants of the Mississippi type",
    colour = "Treatment", fill = "Treatment"
  )

## ----data-slice-args----------------------------------------------------------
args(gratia:::data_slice.gam)

## ----data-sim-and-fit---------------------------------------------------------
# simulate data from the bivariate surface
df <- data_sim("eg2", n = 1000, scale = 0.25, seed = 2)

# fit the GAM
m_biv <- gam(y ~ te(x, z), data = df, method = "REML")

## ----data-slice-biv-1---------------------------------------------------------
ds3 <- data_slice(m_biv,
  x = evenly(x, n = 100),
  z = quantile(z, probs = 0.25)
)

z_val <- with(ds3, round(quantile(z, probs = 0.25), 2))
ylab <- bquote(hat(f)(x, .(z_val)))

## ----smooth-estimates-bivar---------------------------------------------------
sm <- smooth_estimates(m_biv, smooth = "te(x,z)", data = ds3) |>
  add_confint()
sm

## ----plot-smooth-estimates-bivar----------------------------------------------
sm |>
  ggplot(aes(x = x, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line() +
  labs(
    title = "Evaluation of smooth te(x,z) at fixed z",
    y = ylab
  )

## ----data-slice-biv-2---------------------------------------------------------
ds4 <- data_slice(m_biv,
  x = evenly(x, n = 100),
  z = round(quantile(z, probs = c(0.25, 0.5, 0.75)), 2)
)

sm <- smooth_estimates(m_biv, smooth = "te(x,z)", data = ds4) |>
  add_confint() |>
  mutate(fz = factor(z))

sm |>
  ggplot(aes(x = x, y = .estimate, colour = fz, group = fz)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = fz, colour = NULL),
    alpha = 0.2
  ) +
  geom_line() +
  labs(
    title = "Evaluation of smooth te(x,z) at fixed z",
    y = expression(hat(f)(x, z)), colour = "z", fill = "z"
  )

## ----data-slice-biv-1-fitted-values-------------------------------------------
fitted_values(m_biv, data = ds3) |> # default is response scale, not link
  ggplot(aes(x = x, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line() +
  labs(
    title = "Fitted values from model",
    y = expression(hat(y))
  )

## ----data-slice-biv-2-fitted-values-------------------------------------------
fitted_values(m_biv, data = ds4) |>
  mutate(fz = factor(z)) |>
  ggplot(aes(x = x, y = .fitted, colour = fz, group = fz)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = fz, colour = NULL),
    alpha = 0.2
  ) +
  geom_line() +
  labs(
    title = "Fitted values from model",
    y = expression(hat(y)), colour = "z", fill = "z"
  )

