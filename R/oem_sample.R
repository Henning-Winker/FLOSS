# Sampling: generate observed data from operating models.

#' Generate observed data
#'
#' Samples catch, effort, CPUE, survey index, length frequencies and
#' biological samples from one or several operating models, for the chosen
#' iterations, and returns them with the truth (answer key).
#'
#' Each data stream of each species uses its own reproducible random stream,
#' so changing e.g. the biological-sample settings does not change the CPUE.
#'
#' @param om a `floss_om` or `floss_oms` (several species).
#' @param oem [oem_settings()].
#' @param iters iterations to sample (one dataset per iteration).
#' @param seed base seed.
#' @return an object of class `floss_data`: a list of data frames `catch`,
#'   `effort`, `cpue`, `survey`, `lfd`, `bio`, `lhpars`, `truth`, `refpts`,
#'   plus `settings` (fishery, species and sampling specs).
#' @export
sample_data <- function(om, oem, iters = 1:10, seed = 123) {
  oms <- if (inherits(om, "floss_om")) structure(list(om), names = om$species$id,
                                                 class = "floss_oms") else om
  nsim <- dims(oms[[1]]$stock)$iter
  if (any(iters > nsim)) stop("iters must be <= ", nsim)
  per_sp <- lapply(oms, function(o) {
    tag <- function(x) seed_for(seed, paste(o$species$id, x))
    list(
      catch  = with_seed(tag("catch"), sample_catch(o, oem, iters)),
      cpue   = with_seed(tag("cpue"), sample_cpue(o, oem, iters)),
      survey = with_seed(tag("survey"), sample_survey(o, oem, iters)),
      lfd    = with_seed(tag("lfd"), sample_lfd(o, oem, iters)),
      bio    = with_seed(tag("bio"), sample_bio(o, oem, iters)),
      lhpars = lh_table(o$species),
      truth  = truth_table(o, iters),
      refpts = refpts_table(o))
  })
  bind <- function(what) {
    out <- do.call(rbind, lapply(per_sp, `[[`, what))
    rownames(out) <- NULL
    out
  }
  fy <- oms[[1]]$fishery
  structure(list(
    catch = bind("catch"),
    effort = with_seed(seed_for(seed, "effort"),
                       sample_effort(oms[[1]], oem, iters)),
    cpue = bind("cpue"), survey = bind("survey"), lfd = bind("lfd"),
    bio = bind("bio"), lhpars = bind("lhpars"), truth = bind("truth"),
    refpts = bind("refpts"),
    settings = list(fishery = fy, species = lapply(oms, `[[`, "species"),
                    oem = oem, iters = iters, seed = seed,
                    treatment = oms[[1]]$treatment,
                    artefacts = planted_artefacts(oem, fy, oms[[1]]$treatment),
                    nsim = nsim, floss_version =
                      as.character(utils::packageVersion("FLOSS")))),
    class = "floss_data")
}

#' @exportS3Method base::print
print.floss_data <- function(x, ...) {
  cat("FLOSS data:", length(x$settings$species), "species,",
      length(x$settings$fishery$fleets), "fleets,", length(x$settings$iters),
      "iterations\n")
  n <- vapply(x[c("catch", "effort", "cpue", "survey", "lfd", "bio")], NROW, 1)
  print(n)
  invisible(x)
}

# ---- individual streams -----------------------------------------------------

#' @describeIn sample_data observed catch by fleet (t).
#' @export
sample_catch <- function(om, oem, iters = 1:10) {
  cf <- catch_fleet(om)
  do.call(rbind, lapply(names(cf), function(f) {
    tru <- iter(cf[[f]], iters)
    sd <- cv2sd(oem$catch$cv[[f]])
    bias <- rep(oem$catch$bias[[f]], dims(tru)$year)
    if (!is.null(oem$catch$bias_from))           # bias only from a given year
      bias[an(dimnames(tru)$year) < oem$catch$bias_from[[f]]] <- 0
    obs <- tru * (1 - bias)
    obs[] <- c(obs) * stats::rlnorm(length(obs), -sd^2 / 2, sd)
    long_df(obs, tru, om$species$species, f, "catch_t", "catch_true")
  }))
}

#' @describeIn sample_data observed nominal effort by fleet (relative).
#' @export
sample_effort <- function(om, oem, iters = 1:10) {
  do.call(rbind, lapply(names(om$effort), function(f) {
    tru <- iter(om$effort[[f]], iters)
    sd <- cv2sd(oem$effort$cv[[f]])
    obs <- tru
    obs[] <- c(tru) * stats::rlnorm(length(tru), -sd^2 / 2, sd)
    df <- long_df(obs, tru, NA, f, "effort", "effort_true")
    df$species <- NULL
    df
  }))
}

#' @describeIn sample_data CPUE by fleet with AR1 lognormal residuals,
#'   scaled to mean 1 over the sampled years; `index_true` is the true
#'   vulnerable biomass (without effort creep).
#' @export
sample_cpue <- function(om, oem, iters = 1:10) {
  stk <- iter(om$stock, iters)
  do.call(rbind, lapply(om$fishery$fleets, function(f) {
    yrs <- year_window(stk, oem$years[[f]])
    vb <- vuln_biomass(stk, c(om$sel[[f]]$sel_a), 0.5)[, yrs]
    for (p in fleet_periods(om, f)[-1]) {          # selectivity changes
      yy <- intersect(yrs, ac(p$years))
      if (length(yy)) vb[, yy] <- vuln_biomass(stk, sel_a_for(om, p$sel), 0.5)[, yy]
    }
    beta <- oem$cpue$beta[[f]] %||% 1
    signal <- vb^beta
    if (isTRUE(oem$cpue$creep)) signal <- signal * iter(om$efficiency[[f]], iters)[, yrs]
    qs <- oem$cpue$q_step[[f]]
    if (!is.null(qs)) {                          # artefact: catchability step
      yy <- yrs[an(yrs) >= qs[1]]
      signal[, yy] <- signal[, yy] * qs[2]
    }
    obs <- ar1_obs(signal, oem$cpue$cv[[f]], oem$cpue$rho[[f]])
    # index_true = true vulnerable biomass (no creep), for the answer key
    rel_index_df(obs, vb, om$species$species, f)
  }))
}

#' @describeIn sample_data survey biomass index (relative to its mean).
#' @export
sample_survey <- function(om, oem, iters = 1:10) {
  sv <- oem$survey
  if (is.null(sv)) return(NULL)
  stk <- iter(om$stock, iters)
  yrs <- year_window(stk, sv$years)
  g <- floss_gears(om, oem)$survey
  n <- FLRef::survey_n_gear(stk, g, years = yrs, timing = sv$timing)
  b <- quantSums(n * stock.wt(stk)[, yrs])
  signal <- b
  if (!is.null(sv$q_step)) {                     # artefact: catchability step
    yy <- yrs[an(yrs) >= sv$q_step[1]]
    signal[, yy] <- signal[, yy] * sv$q_step[2]
  }
  obs <- ar1_obs(signal, sv$cv, sv$rho)
  rel_index_df(obs, b, om$species$species, "survey")
}

#' @describeIn sample_data two-stage length samples by fleet and survey
#'   ([FLRef::lfd.sim()]); `len_mm` is the lower bin limit.
#' @export
sample_lfd <- function(om, oem, iters = 1:10) {
  stk <- iter(om$stock, iters)
  gears <- floss_gears(om, oem)
  cn <- catch_fleet(om, "n")
  out <- lapply(om$fishery$fleets, function(f) {
    yrs <- year_window(stk, oem$years[[f]])
    do.call(rbind, lapply(fleet_periods(om, f), function(p) {
      yy <- intersect(yrs, ac(p$years))
      if (!length(yy)) return(NULL)
      g <- if (identical(p$sel, om$species$sel[[f]])) gears[[f]] else
        gear_for(om, oem, f, p$sel)
      x <- iter(cn[[f]], iters)[, yy]
      lfd_df(FLRef::lfd.sim(x, g, scale = FALSE), om$species$species, f)
    }))
  })
  sv <- oem$survey
  if (!is.null(sv) && isTRUE(sv$lfd)) {
    yrs <- year_window(stk, sv$years)
    n <- FLRef::survey_n_gear(stk, gears$survey, years = yrs, timing = sv$timing)
    out <- c(out, list(lfd_df(FLRef::lfd.sim(n, gears$survey, scale = FALSE),
                              om$species$species, "survey")))
  }
  do.call(rbind, out)
}

#' @describeIn sample_data biological samples: one row per fish.
#' @export
sample_bio <- function(om, oem, iters = 1:10) {
  stk <- iter(om$stock, iters)
  sp <- om$species
  gears <- floss_gears(om, oem)
  cn <- catch_fleet(om, "n")
  bio <- oem$bio
  L50 <- sp$lh$lm
  L95 <- sp$lh$lm * bio$mat_L95
  sdw <- cv2sd(bio$wt_cv)
  src <- c(stats::setNames(om$fishery$fleets, om$fishery$fleets),
           if (!is.null(oem$survey) && oem$survey$bio_n > 0) c(survey = "survey"))
  rows <- lapply(names(src), function(s) {
    if (s == "survey") {
      yrs <- year_window(stk, oem$survey$years)
      x <- FLRef::survey_n_gear(stk, gears$survey, years = yrs,
                                timing = oem$survey$timing)
      nfish <- oem$survey$bio_n
    } else {
      yrs <- year_window(stk, oem$years[[s]])
      x <- iter(cn[[s]], iters)[, yrs]
      nfish <- bio$n[[s]]
    }
    if (nfish <= 0) return(NULL)
    key_for <- function(y) {
      if (s == "survey") return(gears$survey$condALK)
      for (p in fleet_periods(om, s)[-1])
        if (an(y) %in% p$years) return(gear_for(om, oem, s, p$sel)$condALK)
      gears[[s]]$condALK
    }
    K <- gears[[s]]$condALK
    lens <- as.numeric(colnames(K))
    bw <- gears[[s]]$grid$bin
    ages <- as.numeric(rownames(K))
    keys <- lapply(yrs, key_for)
    names(keys) <- yrs
    do.call(rbind, lapply(seq_along(iters), function(j) {
      do.call(rbind, lapply(yrs, function(y) {
        K <- keys[[y]]
        p <- c(x[, y, , , , j])
        if (sum(p) <= 0) return(NULL)
        a <- sample(ages, nfish, replace = TRUE, prob = p)
        l <- vapply(a, function(ai) {
          k <- sample(seq_along(lens), 1, prob = K[match(ai, ages), ])
          lens[k] + stats::runif(1) * bw
        }, 0)
        w <- sp$lh$a * l^sp$lh$b * stats::rlnorm(nfish, -sdw^2 / 2, sdw)
        pm <- 1 / (1 + exp(-log(19) * (l - L50) / (L95 - L50)))
        # observed age: true age x (1 + ageing bias) + reading error
        age_obs <- pmax(0, round(a * (1 + (bio$age_bias %||% 0)) +
                                   stats::rnorm(nfish, 0, bio$age_cv * pmax(a, 1))))
        if (identical(bio$age_scheme, "stratified")) {
          # bins on the reported (mm-rounded) length, as an observer would
          bin_id <- floor(round(l * 10) / (bw * 10))
          keep <- stats::ave(seq_along(l), bin_id, FUN = function(i)
            seq_along(i) %in% sample(seq_along(i), min(length(i), bio$n_per_bin)))
          age_obs[keep == 0] <- NA
        } else {
          age_obs[stats::runif(nfish) > bio$p_aged] <- NA
        }
        data.frame(species = sp$species, iter = iters[j], year = as.integer(y),
                   source = s, fish_id = seq_len(nfish),
                   len_mm = round(l * 10), age = age_obs,
                   wt_g = round(w, 1), mature = stats::rbinom(nfish, 1, pm),
                   age_true = a)
      }))
    }))
  })
  do.call(rbind, rows)
}

# ---- truth ------------------------------------------------------------------

truth_table <- function(om, iters) {
  s <- iter(om$stock, iters)
  rp <- om$stock@refpts
  df <- function(x, nm) {
    d <- as.data.frame(x)[, c("year", "iter", "data")]
    names(d)[3] <- nm
    d
  }
  out <- Reduce(function(a, b) merge(a, b, by = c("year", "iter")), list(
    df(ssb(s), "ssb"), df(rec(s), "rec"), df(fbar(s), "fbar"),
    df(catch(s), "catch"), df(ssb(s) / c(rp["Bmsy"]), "B_Bmsy"),
    df(fbar(s) / c(rp["Fmsy"]), "F_Fmsy")))
  out$iter <- iters[as.integer(as.character(out$iter))]
  cbind(species = om$species$species, out[order(out$iter, out$year), ],
        row.names = NULL)
}

refpts_table <- function(om) {
  rp <- om$stock@refpts
  data.frame(species = om$species$species, t(stats::setNames(c(rp), rownames(rp))),
             target_msy = om$species$msy, check.names = FALSE)
}

lh_table <- function(sp) {
  lh <- sp$lh
  data.frame(species = sp$species, linf_cm = lh$linf, k = lh$k, t0 = lh$t0,
             a = lh$a, b = lh$b, lm_cm = lh$lm, tmax = lh$tmax, M = lh$M,
             steepness = lh$s)
}

# ---- helpers ----------------------------------------------------------------

#' FLRef gear objects for an operating model
#'
#' Builds [FLRef::build_gear()] objects (length selectivity, capture-
#' conditioned length key, sample sizes) for every fleet and the survey, on
#' the length grid of the sampling settings.
#' @param om a `floss_om`.
#' @param oem [oem_settings()].
#' @return named list of gears.
#' @export
floss_gears <- function(om, oem) {
  sp <- om$species
  lh <- FLPar(linf = sp$lh$linf, k = sp$lh$k, t0 = sp$lh$t0, a = sp$lh$a,
              b = sp$lh$b)
  bin <- oem$lfd$bin %||% max(1, round(sp$lh$linf / 40))
  cfg <- lapply(om$fishery$fleets, function(f)
    c(unclass(sp$sel[[f]]), f_mult = 1, ess_age = oem$lfd$ess_age[[f]],
      ess_len = oem$lfd$ess_len[[f]]))
  names(cfg) <- om$fishery$fleets
  sv <- oem$survey
  if (!is.null(sv)) {
    s <- sv$sel
    rel <- intersect(names(s), c("L50", "L95", "Lpeak", "sd", "sd_left", "sd_right"))
    s[rel] <- lapply(s[rel], function(v) v * sp$lh$linf)
    cfg$survey <- c(s, obs = "survey", ess_age = sv$ess_age,
                    ess_len = sv$ess_len, timing = sv$timing)
  }
  FLRef::build_gear(cfg, lhpar = lh, age = 0:sp$lh$tmax, timing = 0.5,
                    lmin = bin, lmax_mult = 1.2, bin = bin, cv = oem$lfd$cv)
}

cv2sd <- function(cv) sqrt(log(1 + cv^2))

# gear for one fleet with a different selectivity (e.g. after a gear switch)
gear_for <- function(om, oem, fleet, sel) {
  om2 <- om
  om2$species$sel[[fleet]] <- sel
  floss_gears(om2, oem)[[fleet]]
}

# selectivity-at-age (ages 0:tmax) of a selectivity spec for the OM species
sel_a_for <- function(om, sel) {
  c(FLRef::sel_la(species_lhpar(om$species), amin = 0, amax = om$species$lh$tmax,
                  type = sel$type, L50 = sel$L50, L95 = sel$L95, Lpeak = sel$Lpeak,
                  sd = sel$sd, sd_left = sel$sd_left, sd_right = sel$sd_right)$sel_a)
}

year_window <- function(stk, yrs) {
  have <- an(dimnames(stk)$year)
  ac(have[have >= yrs[1] & have <= yrs[2]])
}

vuln_biomass <- function(stk, sel_a, timing = 0.5) {
  z <- harvest(stk) + m(stk)
  quantSums(stock.n(stk) * exp(-z * timing) * sel_a * stock.wt(stk))
}

ar1_obs <- function(x, cv, rho) {
  sd <- cv2sd(cv)
  ny <- dim(x)[2]
  ni <- dim(x)[6]
  e <- matrix(0, ny, ni)
  e[1, ] <- stats::rnorm(ni, 0, sd)
  if (ny > 1) for (i in 2:ny)
    e[i, ] <- rho * e[i - 1, ] + sqrt(1 - rho^2) * stats::rnorm(ni, 0, sd)
  out <- x
  out[] <- c(x) * exp(c(e) - sd^2 / 2)
  out
}

long_df <- function(obs, tru, species, source, nm, nm_true) {
  d <- as.data.frame(obs)
  t <- as.data.frame(tru)
  out <- data.frame(species = species, iter = as.integer(as.character(d$iter)),
                    year = as.integer(as.character(d$year)), fleet = source)
  out[[nm]] <- d$data
  out[[nm_true]] <- t$data
  out
}

rel_index_df <- function(obs, tru, species, source) {
  d <- as.data.frame(obs)
  t <- as.data.frame(tru)
  it <- as.character(d$iter)
  out <- data.frame(species = species, iter = as.integer(it),
                    year = as.integer(as.character(d$year)), source = source,
                    index = d$data / stats::ave(d$data, it),
                    index_true = t$data / stats::ave(t$data, it))
  out
}

lfd_df <- function(x, species, source) {
  d <- as.data.frame(x)
  d <- d[d$data > 0, ]
  data.frame(species = species, iter = as.integer(as.character(d$iter)),
             year = as.integer(as.character(d$year)), source = source,
             len_mm = round(as.numeric(as.character(d$len)) * 10),
             n = as.integer(round(d$data)))
}
