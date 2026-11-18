# Growth-type groups (growth morphs) alongside the age-structured operating model
#
# Selectivity is length-based: within an age, fast growers reach the selected lengths sooner and are removed sooner.
# An age-structured OM with selectivity at age (the age-length key averaged over lengths) removes all fish of an age
# at the same rate, so survivors keep the unexploited length distribution. With `gtg`, build_om() tracks growth
# morphs (asymptotic lengths Linf_g = Linf (1 + cv z_g), recruits split by normal weights; k and M shared) through
# each cohort's fleet F history, and runs the age-structured OM with the resulting effective F, weight and maturity
# at age. Numbers, catch and SSB at age of the FLStock match the morph dynamics, so the OM stays a standard FLR
# object (FLBRP, FLRef, mse). The length samples are drawn from the morphs in each fleet's catch at age.
# The morph composition of a cohort depends only on its F history, not on recruitment, which is what makes this
# separation exact up to the plus group (cohorts mixed per recruit).

#' Settings for growth-type groups in the operating model
#'
#' @param n number of growth morphs.
#' @param maxsd range of the morphs in standard deviations of Linf.
#' @param cv total coefficient of variation of Linf among fish (length-at-age variability).
#' @param cv_res CV of length within a morph; default `cv / sqrt(n)`. The between-morph CV is
#'   `sqrt(cv^2 - cv_res^2)`, so the total spread equals `cv` (split after M. Taylor's FLR platoon code).
#' @param nsub number of sub-annual steps of the morph dynamics (default 12, monthly): within the year each fleet
#'   selects on the length each morph has in that step, so fishing follows growth through the year (as in a
#'   continuous-time length-based model); the FLStock stays annual. Seasons of the sampler group these steps
#'   (12, 4 or 1 seasons with nsub = 12).
#' @return a list for `build_om(gtg = )` / `build_oms(gtg = )`, with `cv_between` added.
#' @export
gtg_settings <- function(n = 15, maxsd = 2.5, cv = 0.1, cv_res = NULL, nsub = 12) {
  cv_res <- cv_res %||% (cv / sqrt(n))
  if (cv_res >= cv) stop("cv_res must be smaller than cv (the total CV of Linf)")
  if (nsub < 1) stop("nsub must be >= 1")
  list(n = n, maxsd = maxsd, cv = cv, cv_res = cv_res, cv_between = sqrt(cv^2 - cv_res^2), nsub = as.integer(nsub))
}

# selectivity at length of a species' fleet spec (FLRef shapes), scaled to a maximum of 1
sel_len_spec <- function(s, L, Lmax) {
  f <- function(L) switch(s$type,
    logistic = 1 / (1 + exp(-log(19) * (L - s$L50) / (s$L95 - s$L50))),
    normal   = exp(-0.5 * ((L - s$Lpeak) / s$sd)^2),
    dnormal  = ifelse(L <= s$Lpeak, exp(-0.5 * ((L - s$Lpeak) / s$sd_left)^2),
                      exp(-0.5 * ((L - s$Lpeak) / s$sd_right)^2)),
    stop("unknown selectivity type ", s$type))
  f(L) / max(f(seq(0.1, Lmax, by = 0.1)))
}

# morph structure: lengths, weights, maturity and selectivity at age by morph (morph x age matrices)
gtg_setup <- function(sp, fleets, gtg, wscale = 1) {
  lh <- sp$lh
  ages <- 0:lh$tmax
  gtg <- gtg_settings(gtg$n, gtg$maxsd, gtg$cv, gtg$cv_res, gtg$nsub %||% 12)
  z <- seq(-gtg$maxsd, gtg$maxsd, length.out = gtg$n)
  recP <- stats::dnorm(z); recP <- recP / sum(recP)
  z <- z / sqrt(sum(recP * z^2))                           # grid SD exactly 1 (truncation at maxsd)
  linf_g <- lh$linf * (1 + gtg$cv_between * z)
  len <- function(t) pmax(outer(linf_g, 1 - exp(-lh$k * (t - lh$t0))), 1e-6)
  vb <- function(t) lh$linf * (1 - exp(-lh$k * (t - lh$t0)))
  L50 <- vb(lh$tm); L95 <- vb(lh$tm + 1)
  L_beg <- len(ages); L_mid <- len(ages + 0.5)
  Lmax <- 1.5 * max(linf_g)
  nsub <- gtg$nsub; tsub <- (seq_len(nsub) - 0.5) / nsub    # sub-annual steps (middle of each)
  L_sub <- array(vapply(tsub, function(t) len(ages + t), L_beg), c(dim(L_beg), nsub))
  sel_sub <- lapply(sp$sel[fleets], function(s) array(sel_len_spec(s, L_sub, Lmax), dim(L_sub)))
  list(ages = ages, recP = recP, linf_g = linf_g, len = len, M = lh$M, cv_res = gtg$cv_res,
       L50 = L50, L95 = L95, fleets = fleets,
       w_beg = wscale * lh$a * L_beg^lh$b, w_mid = wscale * lh$a * L_mid^lh$b,
       w_sub = wscale * lh$a * L_sub^lh$b, nsub = nsub, tsub = tsub, L_sub = L_sub,
       mat_beg = cbind(0, (1 / (1 + exp(-log(19) * (L_beg - L50) / (L95 - L50))))[, -1, drop = FALSE]),  # age 0 immature (FLasher ignores age-0 SSB)
       sel_sub = sel_sub,                                    # morph x age x sub-step, at the length in that step
       sel = lapply(sel_sub, function(x) apply(x, c(1, 2), mean)),   # mean over the year: annual Z and survival
       sel_fun = lapply(sp$sel[fleets], function(s) function(L) sel_len_spec(s, L, Lmax)),
       settings = gtg)
}

# catch numbers by fleet, morph, age and sub-annual step for numbers l (morph x age) at the start of the year:
# within the year, each fleet selects on the length of each morph in that step. Survival over the year equals
# exp(-Z) with Z from the mean selectivity (G$sel), as in gtg_track().
gtg_catch_sub <- function(G, l, Ff) {
  fl <- names(Ff); dt <- 1 / G$nsub
  N <- l
  out <- lapply(fl, function(f) array(0, c(dim(l), G$nsub))); names(out) <- fl
  for (s in seq_len(G$nsub)) {
    Fs <- lapply(fl, function(f) G$sel_sub[[f]][, , s] * Ff[[f]])
    Zs <- G$M + Reduce(`+`, Fs)
    D <- N * (1 - exp(-Zs * dt)) / Zs
    for (i in seq_along(fl)) out[[i]][, , s] <- D * Fs[[i]]
    N <- N * exp(-Zs * dt)
  }
  out
}

# per-recruit numbers by morph and age at constant fleet F (named vector by fleet)
gtg_pr <- function(G, Ff) {
  Z <- G$M + Reduce(`+`, lapply(names(Ff), function(f) G$sel[[f]] * Ff[[f]]))
  na <- length(G$ages)
  l <- matrix(0, length(G$recP), na); l[, 1] <- G$recP
  for (a in 2:na) l[, a] <- l[, a - 1] * exp(-Z[, a - 1])
  l[, na] <- l[, na] / (1 - exp(-Z[, na]))
  list(l = l, Z = Z)
}

# equilibrium SPR at constant fleet F (spawning at the start of the year, maturity at length)
gtg_spr <- function(G, Ff) {
  f0 <- stats::setNames(rep(0, length(Ff)), names(Ff))
  sum(gtg_pr(G, Ff)$l * G$w_beg * G$mat_beg) / sum(gtg_pr(G, f0)$l * G$w_beg * G$mat_beg)
}

# track the morph composition of each cohort through the fleet F history of one iteration
# Fy: matrix years x fleets (apical F). Per recruit; year 1 starts from the equilibrium at its F.
gtg_track <- function(G, Fy) {
  ny <- nrow(Fy); na <- length(G$ages); ng <- length(G$recP); fl <- colnames(Fy)
  l <- Z <- array(0, c(ng, na, ny))
  for (y in seq_len(ny)) Z[, , y] <- G$M + Reduce(`+`, lapply(fl, function(f) G$sel[[f]] * Fy[y, f]))
  l[, , 1] <- gtg_pr(G, Fy[1, ])$l
  for (y in seq_len(ny - 1)) {
    l[, 1, y + 1] <- G$recP
    l[, 2:na, y + 1] <- l[, 1:(na - 1), y] * exp(-Z[, 1:(na - 1), y])
    l[, na, y + 1] <- l[, na, y + 1] + l[, na, y] * exp(-Z[, na, y])
  }
  list(l = l, Z = Z)
}

# aggregates by age and year for one iteration: effective F (survival-equivalent), weight and maturity at age,
# catch numbers, catch shares and catch weight by fleet
gtg_aggregate <- function(G, Fy, tr = gtg_track(G, Fy), keep_sub = FALSE) {
  l <- tr$l; Z <- tr$Z; fl <- colnames(Fy); ny <- nrow(Fy)
  N <- apply(l, c(2, 3), sum)
  surv <- apply(l * exp(-Z), c(2, 3), sum)
  Feff <- pmax(-log(surv / N) - G$M, 0)
  wb <- array(G$w_beg, dim(l)); mb <- array(G$mat_beg, dim(l))
  wt <- apply(l * wb, c(2, 3), sum) / N
  mat <- apply(l * wb * mb, c(2, 3), sum) / apply(l * wb, c(2, 3), sum)
  csub <- lapply(fl, function(f) array(0, c(dim(l), G$nsub))); names(csub) <- fl
  cn <- lapply(fl, function(f) array(0, dim(l))); names(cn) <- fl
  cw <- lapply(fl, function(f) matrix(0, dim(l)[2], ny)); names(cw) <- fl
  for (y in seq_len(ny)) {
    cy <- gtg_catch_sub(G, l[, , y], stats::setNames(Fy[y, ], fl))
    for (f in fl) {
      if (keep_sub) csub[[f]][, , y, ] <- cy[[f]]
      cn[[f]][, , y] <- apply(cy[[f]], c(1, 2), sum)
      cw[[f]][, y] <- colSums(apply(cy[[f]] * G$w_sub, c(1, 2), sum))
    }
  }
  c_f <- lapply(cn, function(x) apply(x, c(2, 3), sum))
  ctot <- Reduce(`+`, c_f)
  list(Feff = Feff, wt = wt, mat = mat, share = lapply(c_f, function(x) ifelse(ctot > 0, x / ctot, 0)),
       cwt = lapply(fl, function(f) ifelse(c_f[[f]] > 0, cw[[f]] / c_f[[f]], 0)) |> stats::setNames(fl),
       cwt_tot = ifelse(ctot > 0, Reduce(`+`, cw) / ctot, 0), cn = cn, csub = if (keep_sub) csub)
}

# within-morph length distributions in each sub-annual step, weighted by the gear's selectivity at length and
# normalised within each morph: list over steps of lists over ages of morph x length matrices. Within a morph,
# length is normal around the morph's VB length in that step (SD cv_res x length).
gtg_len_dist <- function(G, lens, bin, fleet) {
  sl <- G$sel_fun[[fleet]](lens + bin / 2)
  lapply(seq_len(G$nsub), function(s) {
    Lm <- G$L_sub[, , s]; sd <- pmax(G$cv_res * Lm, 1e-3)
    lapply(seq_along(G$ages), function(a) {
      P <- stats::pnorm(outer(Lm[, a], lens + bin, "-") / -sd[, a]) -
           stats::pnorm(outer(Lm[, a], lens, "-") / -sd[, a])
      P <- sweep(P, 2, sl, `*`)
      P / pmax(rowSums(P), 1e-300)
    })
  })
}

# seasonal catch fractions at age (age x season) and age-length keys (list over seasons of age x length) of one
# fleet in one year, from its catch by morph, age and sub-annual step (csub_y: morph x age x step). Seasons group
# the sub-annual steps (nsub = 12 and 12, 4 or 1 seasons: months, quarters, year).
# Building catch at length from morph numbers follows genLF() in M. Taylor's FLR platoon code.
gtg_season_keys <- function(G, csub_y, S, ld) {
  na <- length(G$ages); nl <- ncol(ld[[1]][[1]])
  if (S <= G$nsub) {                                       # seasons group the steps
    grp <- split(seq_len(G$nsub), ceiling(seq_len(G$nsub) * S / G$nsub)); share <- rep(1, S)
  } else {                                                 # more seasons than steps: each season in its step
    st <- ceiling(season_mid(S) * G$nsub); grp <- as.list(st); share <- 1 / tabulate(st, G$nsub)[st]
  }
  ca <- apply(csub_y, 2, sum)
  fr <- sapply(seq_along(grp), function(q) share[q] * apply(csub_y[, , grp[[q]], drop = FALSE], 2, sum))
  fr <- matrix(fr, na); fr <- fr / ifelse(ca > 0, ca, 1)
  keys <- lapply(grp, function(ss) {
    K <- matrix(0, na, nl)
    for (a in seq_len(na)) for (s in ss) K[a, ] <- K[a, ] + as.numeric(csub_y[, a, s] %*% ld[[s]][[a]])
    K / ifelse(rowSums(K) > 0, rowSums(K), 1)
  })
  list(fr = fr, keys = keys)
}

# equilibrium per recruit at constant fleet F: spawning biomass, yield (weight) and effective F at age
gtg_eq <- function(G, Ff) {
  pr <- gtg_pr(G, Ff); l <- pr$l; Z <- pr$Z
  cs <- gtg_catch_sub(G, l, Ff)
  ypr <- sum(vapply(cs, function(x) sum(x * G$w_sub), 0))
  N <- colSums(l); surv <- colSums(l * exp(-Z))
  list(phi = sum(l * G$w_beg * G$mat_beg), ypr = ypr, Feff = pmax(-log(surv / N) - G$M, 0))
}

# reference points of the growth-morph equilibrium for a fleet F pattern p (apical F by fleet), with the OM's
# Beverton-Holt SRR (rec = a S / (b + S)): weight, maturity and selectivity at age change with F as the morph
# composition changes, which FLBRP on a fixed (fished) pattern does not capture. F on the fbar scale of the FLStock
# (effective F averaged over the fbar ages).
gtg_refpts <- function(G, a, b, p, fbar_ages, spr_ref = 0.4) {
  eqx <- function(x) gtg_eq(G, x * p)
  fb <- function(e) mean(e$Feff[match(fbar_ages, G$ages)])
  phi0 <- eqx(0)$phi
  ssb <- function(e) max(a * e$phi - b, 0)
  Y <- function(x) { e <- eqx(x); s <- ssb(e); if (s <= 0) 0 else s / e$phi * e$ypr }
  xup <- 1; while (ssb(eqx(xup)) > 0 && xup < 1e3) xup <- xup * 2
  xmsy <- stats::optimize(Y, c(0, xup), maximum = TRUE)$maximum
  em <- eqx(xmsy)
  xr <- stats::uniroot(function(x) eqx(x)$phi / phi0 - spr_ref, c(0, xup))$root
  er <- eqx(xr)
  c(Fmsy = fb(em), Bmsy = ssb(em), MSY = Y(xmsy), B0 = a * phi0 - b,
    Fspr = fb(er), Bspr = ssb(er), SPRmsy = em$phi / phi0, Fmult_msy = xmsy, Fmult_spr = xr)
}

#' Reference points with growth morphs from steepness and unfished recruitment
#'
#' Fmsy, Bmsy, MSY and B0 of the growth-morph equilibrium, for a Beverton-Holt stock-recruitment relationship
#' given by steepness `h` and unfished recruitment `R0`. Weight, maturity and effective selectivity at age
#' change with F as fast growers are fished out, which equilibrium methods on a fixed selectivity at age miss
#' (see [build_om()] with `gtg`).
#'
#' The relationship `R = a S / (b + S)` follows from `h` and `R0` with the unfished spawning biomass per
#' recruit `phi0` of the morphs: `B0 = R0 phi0`, `a = 4 h R0 / (5 h - 1)`, `b = B0 (1 - h) / (5 h - 1)`.
#' F is moved along the fleet pattern `p` (apical F by fleet, scaled to sum to 1): `F_f = x p_f`.
#' @param x a species spec ([load_species()]) or a `floss_om` built with `gtg`.
#' @param h steepness. Default: the species' `s`; for an OM, its own SRR when `h` and `R0` are both `NULL`.
#' @param R0 unfished recruitment (default 1, so biomass is per recruit). For an OM with `h` given and
#'   `R0 = NULL`, the OM's R0.
#' @param pattern apical F by fleet (named vector; any scale). Default: the species' `relq` for a spec; the
#'   median fleet F of the last three years for an OM.
#' @param gtg morph settings ([gtg_settings()]), used for a species spec.
#' @param fbar_ages ages over which effective F is averaged for Fmsy (default 1 to tmax - 1, as in FLOSS OMs;
#'   for an OM, its fbar range).
#' @param spr_ref SPR of the SPR-based proxy (default 0.4).
#' @return an `FLPar` with `Fmsy` (on the fbar scale, effective F), `Bmsy`, `MSY`, `B0`, `BmsyB0`, `SPRmsy`,
#'   `Fspr` and `Bspr` (at `spr_ref`), `h`, `R0`, and the apical F of each fleet at MSY (`Fmsy_<fleet>`).
#'   For a species spec, weights are in kg (FLOSS specs give g), so biomass is in tonnes when `R0` is in
#'   thousands.
#' @examples
#' \dontrun{
#' sp <- load_species("siganus_sutor")
#' gtg_brp(sp, h = 0.7, R0 = 1000)
#' gtg_brp(sp, h = 0.7, pattern = c(traps = 1, gillnet = 0, linefish = 0))   # traps only
#' }
#' @export
gtg_brp <- function(x, h = NULL, R0 = NULL, pattern = NULL, gtg = gtg_settings(), fbar_ages = NULL,
                    spr_ref = 0.4) {
  a <- b <- NULL                                          # set directly only for an OM's own SRR
  if (inherits(x, "floss_om")) {
    if (is.null(x$gtg)) stop("the OM was built without growth morphs (gtg)")
    G <- x$gtg$G
    if (is.null(pattern)) {
      ly <- utils::tail(dimnames(x$stock)$year, 3)
      pattern <- vapply(G$fleets, function(f) stats::median(c(apply(x$fapic[[f]][, ly], 6, mean))), 0)
    }
    rng <- range(x$stock)
    fbar_ages <- fbar_ages %||% (rng[["minfbar"]]:rng[["maxfbar"]])
    a0 <- an(params(x$sr)["a"]); b0 <- an(params(x$sr)["b"])
    if (is.null(h) && is.null(R0)) { a <- a0; b <- b0 }
    else {
      p0 <- stats::setNames(rep(0, length(G$fleets)), G$fleets)
      phi0 <- gtg_eq(G, p0)$phi
      R0 <- R0 %||% ((a0 * phi0 - b0) / phi0)
      h <- h %||% { B0 <- a0 * phi0 - b0; 0.2 * (b0 + B0) / (b0 + 0.2 * B0) }
    }
  } else {
    sp <- x
    fl <- names(pattern %||% sp$relq)
    G <- gtg_setup(sp, fl, do.call(gtg_settings, gtg[intersect(names(gtg), c("n", "maxsd", "cv", "cv_res", "nsub"))]),
                   wscale = 1 / 1000)
    pattern <- pattern %||% unlist(sp$relq)[fl]
    fbar_ages <- fbar_ages %||% (1:(sp$lh$tmax - 1))
    h <- h %||% sp$lh$s
    R0 <- R0 %||% 1
  }
  pattern <- pattern[G$fleets]
  if (any(!is.finite(pattern)) || sum(pattern) <= 0) stop("'pattern' needs positive F for at least one fleet")
  p <- pattern / sum(pattern)
  if (is.null(a)) {
    if (h <= 0.2 || h > 1) stop("steepness must be in (0.2, 1]")
    phi0 <- gtg_eq(G, p * 0)$phi
    B0 <- R0 * phi0
    a <- 4 * h * R0 / (5 * h - 1)
    b <- B0 * (1 - h) / (5 * h - 1)
  }
  rp <- gtg_refpts(G, a, b, p, fbar_ages = fbar_ages, spr_ref = spr_ref)
  B0 <- rp[["B0"]]
  h_out <- 0.2 * (b + B0) / (b + 0.2 * B0)
  R0_out <- a * B0 / (b + B0)
  out <- c(rp[c("Fmsy", "Bmsy", "MSY", "B0")], BmsyB0 = rp[["Bmsy"]] / B0, rp[c("SPRmsy", "Fspr", "Bspr")],
           h = h_out, R0 = R0_out, stats::setNames(rp[["Fmult_msy"]] * p, paste0("Fmsy_", names(p))))
  do.call(FLPar, as.list(out))
}
