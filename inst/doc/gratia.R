## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library("gratia")
library("mgcv")

## ----data-sim-----------------------------------------------------------------
df <- data_sim("eg1", seed = 42)
df

## ----fit-gam------------------------------------------------------------------
m <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = df, method = "REML")
summary(m)

## ----draw-gam-----------------------------------------------------------------
draw(m)

## ----smooth-estimates---------------------------------------------------------
sm <- smooth_estimates(m)
sm

## ----smooths------------------------------------------------------------------
smooths(m)

## ----smooth-estimates-x2------------------------------------------------------
sm <- smooth_estimates(m, smooth = "s(x2)")
sm

## ----ggplot-smooth------------------------------------------------------------
library("ggplot2")
library("dplyr")
sm %>%
  add_confint() %>%
  ggplot(aes(y = est, x = x2)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
                alpha = 0.2, fill = "forestgreen") +
    geom_line(colour = "forestgreen", size = 1.5) +
    labs(y = "Partial effect",
         title = expression("Partial effect of" ~ f(x[2])),
         x = expression(x[2]))

## ----appraise-----------------------------------------------------------------
appraise(m)

## ----draw-partial-residuals---------------------------------------------------
draw(m, residuals = TRUE)

