# Test lp_matrix() & methods

test_that("lp_matrix works for a GAM", {
  expect_silent(xp <- lp_matrix(m_gam))
  expect_s3_class(xp, "lp_matrix")
  expect_s3_class(xp, "matrix")
  expect_silent(ds <- data_slice(m_gam,
    x2 = evenly(x2, n = 50),
    data = su_eg1, envir = teardown_env()
  ))
  expect_silent(xp <- lp_matrix(m_gam, data = ds))
  expect_s3_class(xp, "lp_matrix")
  expect_s3_class(xp, "matrix")
})

test_that("print() method returns output invisibly", {
  # skip_on_ci() # testing without as moved to mac os x

  xp <- lp_matrix(m_gam)
  expect_output(ret <- withVisible(print(xp)))
  expect_false(ret$visible)
})

test_that("print() output is as expected", {
  # skip_on_ci() # testing without as moved to mac os x
  skip_on_cran()

  expect_snapshot({
    print(lp_matrix(m_gam))
  })
})
