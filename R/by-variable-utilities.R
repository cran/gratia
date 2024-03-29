#' Tests for by variable smooths
#'
#' Functions to check if a smooth is a by-variable one and to test of the type
#' of by-variable smooth is a factor-smooth or a continous-smooth interaction.
#'
#' @param smooth an object of class `"mgcv.smooth"`
#'
#' @return A logical vector.
#'
#' @author Gavin L. Simpson
#'
#' @export
#' @rdname is_by_smooth
`is_by_smooth` <- function(smooth) {
  check_is_mgcv_smooth(smooth)
  is_factor_by_smooth(smooth) | is_continuous_by_smooth(smooth)
}

#' @export
#' @rdname is_by_smooth
`is_factor_by_smooth` <- function(smooth) {
  check_is_mgcv_smooth(smooth)
  by.level <- smooth[["by.level"]]
  !is.null(by.level)
}

#' @export
#' @rdname is_by_smooth
`is_continuous_by_smooth` <- function(smooth) {
  check_is_mgcv_smooth(smooth)
  by.level <- by_level(smooth)
  by.var <- by_variable(smooth)
  !(is.null(by.level) & by.var == "NA")
}

#' @export
#' @rdname is_by_smooth
`by_variable` <- function(smooth) {
  check_is_mgcv_smooth(smooth)
  as.character(smooth[["by"]])
}

#' @export
#' @rdname is_by_smooth
`by_level` <- function(smooth) {
  check_is_mgcv_smooth(smooth)
  smooth[["by.level"]]
}


#' @importFrom tibble add_column
`add_missing_by_info_to_smooth` <- function(smooth) {
  smooth <- add_column(smooth,
    by_variable = factor(rep(NA_character_, nrow(smooth))),
    .after = 1L
  )
}

## Adds data to an object in a consistent way
##
## Variables in `data` are selected using `vars`. Then they are renamed
## `.x1`, `.x2` etc.
##
#' @importFrom dplyr bind_cols
#' @importFrom rlang !! :=
add_smooth_var_data <- function(x, vars, data) {
  sm_data <- data[vars]
  names(sm_data) <- paste0(".x", seq_len(NCOL(sm_data)))
  sm_data <- bind_cols(x, sm_data)
  sm_data
}

## Adds information about factor by variables to smooths
##
#' @importFrom tibble add_column
#' @importFrom rlang !! :=
`add_by_data` <- function(x, n = NULL, by_name, by_data, before = 1L) {
  ## n is number of observations to add
  ## by_name should be the name of the factor variable, which should
  ## be in by_data, a data frame/tibble of data used to create `x`
  if (is.null(n)) {
    n <- NROW(x)
  }

  # check if by_name is "NA" ! yes, that's how Simon stores it
  if (by_name == "NA") {
    by_name <- NA_character_
  }

  # add on the by variable info
  x <- add_by_var_column(x, by_var = by_name, n = n)

  if (!is.na(by_name)) { # no longer "NA" remember
    x <- add_column(x, !!(by_name) := by_data[[by_name]], .after = before)
  }
  x
}

# !! newer, used in several eval_smooth methods
#' @importFrom tibble add_column
`add_by_var_column` <- function(object, by_var, n = NULL) {
  if (is.null(n)) {
    n <- NROW(object)
  }
  add_column(object, .by = rep(by_var, times = n), .after = 1L)
}

# !! newer, used in several eval_smooth methods
#' @importFrom tibble add_column
`add_smooth_type_column` <- function(object, sm_type, n = NULL) {
  if (is.null(n)) {
    n <- NROW(object)
  }
  add_column(object, .type = rep(sm_type, times = n), .after = 1L)
}
