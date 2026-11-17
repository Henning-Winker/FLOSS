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
#' @return a list for `build_om(gtg = )` / `build_oms(gtg = )`, with `cv_between` added.
#' @export
gtg_settings <- function(n = 15, maxsd = 2.5, cv = 0.1, cv_res = NULL) {
  cv_res <- cv_res %||% (cv / sqrt(n))
  if (cv_res >= cv) stop("cv_res must be smaller than cv (the total CV of Linf)")
  list(n = n, maxsd = maxsd, cv = cv, cv_res = cv_res, cv_between = sqrt(cv^2 - cv_res^2))
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
  gtg <- gtg_settings(gtg$n, gtg$maxsd, gtg$cv, gtg$cv_res)
  z <- seq(-gtg$maxsd, gtg$maxsd, length.out = gtg$n)
  recP <- stats::dnorm(z); recP <- recP / sum(recP)
  z <- z / sqrt(sum(recP * z^2))                           # grid SD exactly 1 (truncation at maxsd)
  linf_g <- lh$linf * (1 + gtg$cv_between * z)
  len <- function(t) pmax(outer(linf_g, 1 - exp(-lh$k * (t - lh$t0))), 1e-6)
  vb <- function(t) lh$linf * (1 - exp(-lh$k * (t - lh$t0)))
  L50 <- vb(lh$tm); L95 <- vb(lh$tm + 1)
  L_beg <- len(ages); L_mid <- len(ages + 0.5)
  Lmax <- 1.5 * max(linf_g)
  list(ages = ages, recP = recP, linf_g = linf_g, len = len, M = lh$M, cv_res = gtg$cv_res,
       L50 = L50, L95 = L95, fleets = fleets,
       w_beg = wscale * lh$a * L_beg^lh$b, w_mid = wscale * lh$a * L_mid^lh$b,
       mat_beg = cbind(0, (1 / (1 + exp(-log(19) * (L_beg - L50) / (L95 - L50))))[, -1, drop = FALSE]),  # age 0 immature (FLasher ignores age-0 SSB)
       sel = lapply(sp$sel[fleets], function(s) matrix(sel_len_spec(s, L_mid, Lmax), nrow(L_mid))),
       sel_fun = lapply(sp$sel[fleets], function(s) function(L) sel_len_spec(s, L, Lmax)),
       settings = gtg)
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
gtg_aggregate <- function(G, Fy, tr = gtg_track(G, Fy)) {
  l <- tr$l; Z <- tr$Z; fl <- colnames(Fy); ny <- nrow(Fy)
  N <- apply(l, c(2, 3), sum)
  surv <- apply(l * exp(-Z), c(2, 3), sum)
  Feff <- pmax(-log(surv / N) - G$M, 0)
  wb <- array(G$w_beg, dim(l)); mb <- array(G$mat_beg, dim(l)); wm <- array(G$w_mid, dim(l))
  wt <- apply(l * wb, c(2, 3), sum) / N
  mat <- apply(l * wb * mb, c(2, 3), sum) / apply(l * wb, c(2, 3), sum)
  cn <- lapply(fl, function(f) {
    Ff <- aperm(array(Fy[, f], c(ny, dim(l)[1:2])), c(2, 3, 1))
    l * array(G$sel[[f]], dim(l)) * Ff / Z * (1 - exp(-Z))
  })
  names(cn) <- fl
  c_f <- lapply(cn, function(x) apply(x, c(2, 3), sum))
  cw_f <- lapply(cn, function(x) apply(x * wm, c(2, 3), sum))
  ctot <- Reduce(`+`, c_f)
  list(Feff = Feff, wt = wt, mat = mat, share = lapply(c_f, function(x) ifelse(ctot > 0, x / ctot, 0)),
       cwt = lapply(fl, function(f) ifelse(c_f[[f]] > 0, cw_f[[f]] / c_f[[f]], 0)) |> stats::setNames(fl),
       cwt_tot = ifelse(ctot > 0, Reduce(`+`, cw_f) / ctot, 0), cn = cn)
}

# age-length keys of one fleet's catch (age x length x year x time), from the morphs in the catch at age.
# Within a morph, length at time t in the year is normal around the morph's VB length (SD cv_res x length), and
# the gear selects on that length; cn_f is the morph catch (morph x age x year) from gtg_aggregate().
# Building catch at length from morph numbers follows genLF() in M. Taylor's FLR platoon code.
gtg_keys <- function(G, cn_f, tt, lens, bin, fleet) {
  na <- length(G$ages); ny <- dim(cn_f)[3]; nl <- length(lens)
  sl <- G$sel_fun[[fleet]](lens + bin / 2)
  out <- array(0, c(na, nl, ny, length(tt)))
  for (k in seq_along(tt)) {
    Lm <- G$len(G$ages + tt[k])                            # morph x age
    sd <- pmax(G$cv_res * Lm, 1e-3)
    for (a in seq_len(na)) {
      P <- stats::pnorm(outer(Lm[, a], lens + bin, "-") / -sd[, a]) -
           stats::pnorm(outer(Lm[, a], lens, "-") / -sd[, a])       # morph x length
      K <- crossprod(matrix(cn_f[, a, ], ncol = ny), P)             # year x length
      K <- sweep(K, 2, sl, `*`)
      s <- rowSums(K)
      K <- K / ifelse(s > 0, s, 1)
      out[a, , , k] <- t(K)
    }
  }
  out
}

gtg_key <- function(G, cn_f, y, t, lens, bin, fleet)
  gtg_keys(G, cn_f[, , y, drop = FALSE], t, lens, bin, fleet)[, , 1, 1]

# equilibrium per recruit at constant fleet F: spawning biomass, yield (weight) and effective F at age
gtg_eq <- function(G, Ff) {
  pr <- gtg_pr(G, Ff); l <- pr$l; Z <- pr$Z
  ypr <- sum(vapply(names(Ff), function(f)
    sum(l * G$sel[[f]] * Ff[[f]] / Z * (1 - exp(-Z)) * G$w_mid), 0))
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
    Fspr = fb(er), Bspr = ssb(er), SPRmsy = em$phi / phi0)
}
