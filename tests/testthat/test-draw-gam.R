# Test draw.gam()
# A lot of the current tests for this are over in test-draw-methods.R

test_that("draw.gam works for m_1_smooth", {
  expect_silent(plt <- draw(m_1_smooth, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_1_smooth", plt)
})

test_that("draw.gam works for m_1_smooth with angled labels", {
  expect_silent(plt <- draw(m_1_smooth, angle = 45, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_1_smooth rotated", plt)
})

test_that("draw.gam works for m_1_smooth with partial residuals", {
  expect_silent(plt <- draw(m_1_smooth,
    residuals = TRUE,
    resid_col = "hotpink", rug = FALSE
  ))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_1_smooth hotpink residuals", plt)
})

test_that("draw.gam works for m_gam", {
  expect_silent(plt <- draw(m_gam, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_gam", plt)
})

test_that("draw.gam works for m_gam rotated", {
  expect_silent(plt <- draw(m_gam, rug = FALSE, angle = 45))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_gam rotated", plt)
})

test_that("draw.gam works for m_gamm", {
  skip_on_ci() # testing without as moved to mac os x
  expect_silent(plt <- draw(m_gamm, rug = FALSE))
  expect_doppelganger("draw_gam m_gamm", plt)
})

test_that("draw.gam works for m_gamm4", {
  skip_on_ci() # testing without as moved to mac os x
  skip_on_cran()
  expect_silent(plt <- draw(m_gamm4, rug = FALSE))
  expect_doppelganger("draw_gam m_gamm4", plt)
})

test_that("draw.gam works for m_bam", {
  expect_silent(plt <- draw(m_bam, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_bam", plt)
})

test_that("draw.gam works for m_gaulss", {
  expect_silent(plt <- draw(m_gaulss, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_gaulss", plt)
})

test_that("draw.gam works for m_scat", {
  expect_silent(plt <- draw(m_scat, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_scat", plt)
})

test_that("draw.gam works for m_gamgcv", {
  expect_silent(plt <- draw(m_gamgcv, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_gamgcv", plt)
})

test_that("draw.gam works for rm1", {
  expect_silent(plt <- draw(rm1, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam rm1", plt)
})

test_that("draw.gam works for rm2", {
  expect_silent(plt <- draw(rm2, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam rm2", plt)
})

test_that("draw.gam works for dlnm_m", {
  # skip_on_ci() # testing without as moved to mac os x
  expect_silent(plt <- draw(dlnm_m, rug = FALSE))

  skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam dlnm_m", plt)
})

test_that("draw.gam works for m_ar1", {
  expect_silent(plt <- draw(m_ar1, rug = FALSE, n = 50))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_ar1", plt)
})

test_that("draw.gam works for m_ar1_by", {
  expect_silent(plt <- draw(m_ar1_by, rug = FALSE, n = 50))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_ar1_by", plt)
})

test_that("draw.gam works for m_2_fac", {
  expect_silent(plt <- draw(m_2_fac, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_2_fac", plt)
})

test_that("draw.gam works for m_para_sm", {
  expect_silent(plt <- draw(m_para_sm, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_para_sm", plt)
})

test_that("draw.gam works for m_2_fac", {
  expect_message(
    plt <- draw(m_2_fac,
      parametric = TRUE, rug = FALSE,
      data = df_2_fac, envir = teardown_env()
    ),
    "Interaction terms are not currently supported."
  )
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_2_fac parametric", plt)
})

test_that("draw.gam works for m_para_sm", {
  expect_message(
    plt <- draw(m_para_sm,
      parametric = TRUE, rug = FALSE,
      data = df_2_fac, envir = teardown_env()
    ),
    "Interaction terms are not currently supported."
  )
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_para_sm parametric", plt)
})

test_that("draw.gam with constant works for GAM with 1 smooth", {
  constant <- coef(m_1_smooth)[1]
  expect_silent(plt <- draw(m_1_smooth, constant = constant, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_1_smooth constant", plt)
})

test_that("draw.gam with constant works for GAM with multiple smooths", {
  constant <- coef(m_gam)[1]
  expect_silent(plt <- draw(m_gam, constant = constant, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_gam constant", plt)
})

test_that("draw.gam with fun works for GAM with 1 smooth", {
  expect_silent(plt <- draw(m_1_smooth, fun = exp, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_1_smooth fun", plt)
})

test_that("draw.gam with fun works for GAM with multiple smooths", {
  expect_silent(plt <- draw(m_gam, fun = exp, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_gam fun", plt)
})

test_that("draw.gam with constant + fun works for GAM with 1 smooth", {
  constant <- coef(m_gam)[1]
  expect_silent(plt <- draw(m_1_smooth,
    constant = constant, fun = exp,
    rug = FALSE
  ))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_1_smooth constant fun", plt)
})

test_that("draw.gam with constant + fun works for GAM with multiple smooths", {
  constant <- coef(m_gam)[1]
  expect_silent(plt <- draw(m_gam,
    constant = constant, fun = exp,
    rug = FALSE
  ))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam m_gam constant fun", plt)
})

test_that("draw.gam works for a trivariate smooth", {
  skip_on_os(os = "win")
  # skip_on_os(os = "mac")
  expect_silent(plt <- draw(su_m_trivar_te, n = 50))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam trivar te", plt)
})

test_that("draw.gam works for a quadvariate smooth", {
  skip_on_os(os = "win")
  # skip_on_os(os = "mac")
  expect_silent(plt <- draw(su_m_quadvar_te, n = 50))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam quadvar te", plt)
})

# test_that("draw.gam works for a trivariate smooth rotated labels", {
#     skip_on_os(os = "win")
#     skip_on_os(os = "mac")
#     expect_silent(plt <- draw(su_m_trivar_te, angle = 45))
#     expect_doppelganger("draw_gam trivar te rotated labels", plt)
# })

# test_that("draw.gam works for a quadvariate smooth rotated labels", {
#     skip_on_os(os = "win")
#     skip_on_os(os = "mac")
#     expect_silent(plt <- draw(su_m_quadvar_te, angle = 45))
#     expect_doppelganger("draw_gam quadvar te rotated labels", plt)
# })

test_that("draw.gam issues message for parametric only model", {
  expect_message(
    plt <- draw(m_only_para),
    "Unable to draw any of the model terms."
  )
})

test_that("draw.gam works for a parametric only model", {
  # rotated
  expect_message(
    plt <- draw(m_only_para,
      parametric = TRUE, angle = 90,
      rug = FALSE, data = df_2_fac, envir = teardown_env()
    ),
    "Interaction terms are not currently supported."
  )
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("draw_gam parametric only model rotated labels", plt)
})

test_that("draw.gam works with grouped by option", {
  expect_silent(plt1 <- draw(su_m_factor_by, grouped_by = TRUE, rug = FALSE))
  expect_silent(plt2 <- draw(m_ordered_by, grouped_by = TRUE, rug = FALSE))
  # skip_on_ci() # testing without as moved to mac os x
  expect_doppelganger("grouped by gam", plt1)
  expect_doppelganger("ordered grouped by gam", plt1)
})

test_that("draw gam works with univar tensor products #260", {
  expect_silent(plt_uni_te <- m_univar_te |> draw(rug = FALSE))
  expect_silent(plt_uni_ti <- m_univar_ti |> draw(rug = FALSE))
  expect_silent(plt_uni_t2 <- m_univar_t2 |> draw(rug = FALSE))

  skip_on_cran()

  expect_doppelganger("draw.gam univariate te smooth", plt_uni_te)
  expect_doppelganger("draw.gam univariate ti smooth", plt_uni_ti)
  expect_doppelganger("draw.gam univariate t2 smooth", plt_uni_t2)
})
