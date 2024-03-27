## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library("gratia")
library("mgcv")
library("ggplot2")
library("dplyr")
library("patchwork")

# simulate data
n <- 400
eg1 <- data_sim("eg1", n = n, seed = 1)

# fit model
m <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3),
  data = eg1, method = "REML"
)

## ----draw-eg1-1---------------------------------------------------------------
p <- draw(m)
p

## ----draw-eg1-wrong-theme-----------------------------------------------------
p + theme_bw()

## ----draw-eg1-right-theme-----------------------------------------------------
p & theme_bw()

## ----draw-multiple-patchworks, error = TRUE, fig.show = "hide"----------------
p1 <- draw(m, select = "s(x0)")
p2 <- draw(m, select = "s(x1)")
p3 <- draw(m, select = "s(x2)")
p1 + p2 + p3

## ----draw-multiple-patchworks-layout------------------------------------------
p1 + p2 + p3 + plot_layout(ncol = 3)

## ----draw-multiple-patchworks-draw--------------------------------------------
draw(m, select = c("s(x0)", "s(x1)", "s(x2)"), ncol = 3)

## ----custom-plot-eval-and-add-confint-----------------------------------------
# evaluate the smooths
sm <- smooth_estimates(m) %>%
  add_confint()
sm

## ----custom-plot-add-partial-resids-------------------------------------------
# add partial residuals to data
eg1 <- eg1 %>%
  add_partial_residuals(m)

## ----custom-plot-names--------------------------------------------------------
names(eg1)

## ----custom-plot-sx2----------------------------------------------------------
p_sx2 <- sm %>%
  filter(.smooth == "s(x2)") %>%
  ggplot() +
  geom_rug(aes(x = x2),
    data = eg1,
    sides = "b", length = grid::unit(0.02, "npc")
  ) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = x2),
    alpha = 0.2
  ) +
  geom_point(aes(x = x2, y = `s(x2)`),
    data = eg1, cex = 1.5, colour = "steelblue3"
  ) +
  geom_line(aes(x = x2, y = .estimate), lwd = 1.2) +
  labs(y = "Partial effect", title = "s(x2)")
p_sx2

## ----custom-plot-other-smooths, echo = FALSE----------------------------------
p_sx0 <- sm %>%
  filter(.smooth == "s(x0)") %>%
  ggplot() +
  geom_rug(aes(x = x0),
    data = eg1,
    sides = "b", length = grid::unit(0.02, "npc")
  ) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = x0),
    alpha = 0.2
  ) +
  geom_point(aes(x = x0, y = `s(x0)`),
    data = eg1, cex = 1.5, colour = "steelblue3"
  ) +
  geom_line(aes(x = x0, y = .estimate), lwd = 1.2) +
  labs(y = "Partial effect", title = "s(x0)")

p_sx1 <- sm %>%
  filter(.smooth == "s(x1)") %>%
  ggplot() +
  geom_rug(aes(x = x1),
    data = eg1,
    sides = "b", length = grid::unit(0.02, "npc")
  ) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = x1),
    alpha = 0.2
  ) +
  geom_point(aes(x = x1, y = `s(x1)`),
    data = eg1, cex = 1.5, colour = "steelblue3"
  ) +
  geom_line(aes(x = x1, y = .estimate), lwd = 1.2) +
  labs(y = "Partial effect", title = "s(x1)")

p_sx3 <- sm %>%
  filter(.smooth == "s(x3)") %>%
  ggplot() +
  geom_rug(aes(x = x3),
    data = eg1,
    sides = "b", length = grid::unit(0.02, "npc")
  ) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = x3),
    alpha = 0.2
  ) +
  geom_point(aes(x = x3, y = `s(x3)`),
    data = eg1, cex = 1.5, colour = "steelblue3"
  ) +
  geom_line(aes(x = x3, y = .estimate), lwd = 1.2) +
  labs(y = "Partial effect", title = "s(x3)")

## ----custom-plot-final--------------------------------------------------------
p_sx0 + p_sx1 + p_sx2 + p_sx3 + plot_layout(ncol = 2)

## ----add-factor-to-data-------------------------------------------------------
set.seed(12)
eg1 <- eg1 %>%
  mutate(fac = sample(letters[1:4], nrow(.), replace = TRUE))

## ----plot-partial-resids-with-colours-by-fac----------------------------------
plt <- sm %>%
  filter(.smooth == "s(x2)") %>%
  ggplot() +
  geom_rug(aes(x = x2),
    data = eg1,
    sides = "b", length = grid::unit(0.02, "npc")
  ) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = x2),
    alpha = 0.2
  ) +
  geom_line(aes(x = x2, y = .estimate), lwd = 1.2) +
  labs(y = "Partial effect", title = "s(x2)")

plt +
  geom_point(
    aes(
      x = x2, y = `s(x2)`,
      colour = fac
    ), # <-- map fac to colour aesthetic
    data = eg1, cex = 1.5
  )

## ----plot-partial-resids-with-colours-by-x1-----------------------------------
plt +
  geom_point(
    aes(
      x = x2, y = `s(x2)`,
      colour = x1, size = x1
    ), # <-- map fac to colour aesthetic
    data = eg1, alpha = 0.3
  ) + # <-- deleted cex!!
  scale_colour_viridis_c(option = "plasma")

