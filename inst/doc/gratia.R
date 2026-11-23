## ----setup-knitr, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library("gratia")
load_mgcv()

## ----data-sim-----------------------------------------------------------------
df <- data_sim("eg1", seed = 42)
df

## ----fit-gam------------------------------------------------------------------
m <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = df, method = "REML")
summary(m)

## ----draw-gam-----------------------------------------------------------------
draw(m)

## ----theme-smooths------------------------------------------------------------
library("ggplot2")
draw(m) & theme(
  geom.gratia.smooth = element_geom(colour = "forestgreen", linewidth = 1.5),
  geom.gratia.interval = element_geom(fill = "forestgreen")
)

## ----smooth-estimates---------------------------------------------------------
sm <- smooth_estimates(m)
sm

## ----smooths------------------------------------------------------------------
smooths(m)

## ----smooth-estimates-x2------------------------------------------------------
sm <- smooth_estimates(m, select = "s(x2)")
sm

## ----ggplot-smooth------------------------------------------------------------
sm |>
  add_confint() |>
  ggplot(aes(y = .estimate, x = x2)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
    alpha = 0.2, fill = "forestgreen"
  ) +
  geom_line(colour = "forestgreen", linewidth = 1.5) +
  labs(
    y = "Partial effect",
    title = expression("Partial effect of" ~ f(x[2])),
    x = expression(x[2])
  )

## ----appraise-----------------------------------------------------------------
appraise(m)

## ----draw-partial-residuals---------------------------------------------------
draw(m, residuals = TRUE)

