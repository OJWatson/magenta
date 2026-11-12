#' @noRd
match_clean <- function(a, b, quiet = TRUE) {
  a <- gsub("[[:punct:][:space:]]", "", tolower(stringi::stri_trans_general(a, "latin-ascii")))
  b <- gsub("[[:punct:][:space:]]", "", tolower(stringi::stri_trans_general(b, "latin-ascii")))
  ret <- match(a, b)
  if (sum(is.na(ret) > 0)) {
    dists <- stringdist::seq_distmatrix(lapply(a, utf8ToInt), lapply(b, utf8ToInt))
    ret[is.na(ret)] <- apply(dists[which(is.na(ret)), , drop = FALSE], 1, which.min)
    if (!quiet) {
      print(unique(cbind(a, b[ret])))
    }
  }
  return(ret)
}


#---
#' Match admin region
#'
#' @param country Character for country within which admin2 is in.
#'   Default = NULL
#' @param admin Character for admin region. Some fuzzy logic will be used to
#'   match. If not provided then no seasonality is introduced. Default = NULL
#' @param quiet Whether function should be quiet. If FALSE (default) the
#'   returned country and admin are printed.
#'
#' @return Row number of the matching region in
#'   \code{magenta::admin_units_seasonal}
#'
#' @keywords internal
admin_match <- function(admin = NULL, country = NULL, quiet = FALSE) {
  ads <- magenta::admin_units_seasonal

  # intiialise admin match as no match
  admin_matches <- 0

  if (!is.null(admin)) {
    # if there is no country given then search for the admin unit
    if (is.null(country)) {
      # find exact match
      admin_matches <- grep(paste("^", admin, "$", sep = ""),
        ads$admin1,
        ignore.case = TRUE
      )

      # if exact does not match try closest match
      if (length(admin_matches) != 1) {
        admin_matches <- match_clean(admin, ads$admin1)
      }

      # if we do have a country though find that match first and then find admin
    } else {
      # first find an exact match
      country_matches <- grep(paste("^", country, "$", sep = ""),
        ads$country,
        ignore.case = TRUE
      )

      if (length(unique(ads$country[country_matches])) == 1) {
        chosen.country <- unique(ads$country[country_matches])
      } else {
        # if exact does not match try closest
        country_matches <- match_clean(country, ads$country)
        chosen.country <- unique(ads$country[country_matches])
      }

      # find exact match
      admin_sub_matches <- grep(paste("^", admin, "$", sep = ""),
        ads$admin1[country_matches],
        ignore.case = TRUE
      )

      # if exact does not match try closest dist
      if (length(admin_sub_matches) != 1) {
        admin_sub_matches <- match_clean(admin, ads$admin1[country_matches])
      }
      admin_matches <- country_matches[admin_sub_matches]
    }

    message(
      "Requested: ", admin, ", ", country,
      "\nReturned: ", ads$admin1[admin_matches], ", ", ads$country[admin_matches]
    )
  }

  return(admin_matches)
}


#---
#' Create seasonal theta return
#'
#' @param country Character for country within which admin2 is in.
#'   Default = NULL
#' @param admin Character for admin region. Some fuzzy logic will be used
#'   to match. If not provided then no seasonality is introduced. Default = NULL
#'
#' @keywords internal
seasonal_profile <- function(admin = NULL, country = NULL) {
  admin_matches <- admin_match(admin, country)
  admin_units_seasonal <- magenta::admin_units_seasonal

  if (admin_matches != 0) {
    ssa0 <- admin_units_seasonal$a0[admin_matches]
    ssa1 <- admin_units_seasonal$a1[admin_matches]
    ssa2 <- admin_units_seasonal$a2[admin_matches]
    ssa3 <- admin_units_seasonal$a3[admin_matches]
    ssb1 <- admin_units_seasonal$b1[admin_matches]
    ssb2 <- admin_units_seasonal$b2[admin_matches]
    ssb3 <- admin_units_seasonal$b3[admin_matches]
    theta_c <- admin_units_seasonal$theta_c[admin_matches]

    t <- 1:365
    y <- sapply(t, function(x) {
      max((ssa0 +
        ssa1 * cos(2 * pi * x / 365) +
        ssa2 * cos(2 * 2 * pi * x / 365) +
        ssa3 * cos(3 * 2 * pi * x / 365) +
        ssb1 * sin(2 * pi * x / 365) +
        ssb2 * sin(2 * 2 * pi * x / 365) +
        ssb3 * sin(3 * 2 * pi * x / 365)) / theta_c, 0.001)
    })
  } else {
    ssa0 <- ssa1 <- ssa2 <- ssa3 <- ssb1 <- ssb2 <- ssb3 <- theta_c <- 0
    y <- rep(1, 365)
  }

  return(y)
}

#---
#' Plot country seasonal profiles
#'
#' @param country Character for country
#'
#' @keywords internal
country_seasonal_profiles <- function(country) {
  ads <- magenta::admin_units_seasonal
  if (!is.element(country, unique(ads$country))) {
    stop(paste0(
      "Country is not one of ",
      paste(unique(ads$country), collapse = ", ")
    ))
  }
  admins <- as.character(ads$admin1[which(ads$country == country)])

  if (length(admins) > 16) {
    par(mfrow = c(5, 5))
  } else if (length(admins) > 9) {
    par(mfrow = c(4, 4))
  } else if (length(admins) > 4) {
    par(mfrow = c(3, 3))
  } else {
    par(mfrow = c(2, 2))
  }

  max_ses <- max(sapply(admins, function(x) {
    max(seasonal_profile(x, country))
  }))

  ret <- sapply(admins, function(x) {
    invisible(plot(seasonal_profile(x, country),
      main = x,
      ylim = c(0, max_ses + 0.2),
      ylab = "Seasonal",
      xlab = "Day"
    ))
  })
}
