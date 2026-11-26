# Sex-structured operating model with growth morphs (growth-type groups)
#
# om_sex(gtg = ) tracks growth morphs of each sex (asymptotic lengths Linf_s (1 + cv z_g), k, M, weight and
# maturity of that sex) through each iteration's apical F history. Fleets select at length, with the apical F at
# length shared by the sexes (as in stock_sex()), so within an age fast growers are removed sooner (Rosa Lee) and
# survivors no longer keep the unexploited length distribution. The two-sex age-structured population is then
# projected with the resulting effective F, weight and maturity at age of each sex: recruitment from female
# spawning biomass (Beverton-Holt, deviations shared by the sexes), recruits split equally between the sexes.
# Reference points come from the joint morph equilibrium (female spawning biomass, yield of both sexes), length
# samples from the morphs in each sex's catch. The morph helpers of gtg.R (gtg_pr, gtg_track, gtg_aggregate,
# gtg_catch_sub, gtg_eq, gtg_len_dist) are used unchanged on the structure of each sex.

# morph structure of one sex of a sex-structured stock (same fields as gtg_setup()); fleet selectivity at length
# scaled by relq and by the joint apical value at length, so the apical F at length is shared by the sexes
gtg_setup_sex1 <- function(l, sel, relq, ages, gtg, female, mx) {
  z <- seq(-gtg$maxsd, gtg$maxsd, length.out = gtg$n)
  recP <- stats::dnorm(z); recP <- recP / sum(recP)
  z <- z / sqrt(sum(recP * z^2))
  linf_g <- l$linf * (1 + gtg$cv_between * z)
  len <- function(t) pmax(outer(linf_g, 1 - exp(-l$k * (t - l$t0))), 1e-6)
  L_beg <- len(ages)
  nsub <- gtg$nsub; tsub <- (seq_len(nsub) - 0.5) / nsub
  L_sub <- array(vapply(tsub, function(t) len(ages + t), L_beg), c(dim(L_beg), nsub))
  fl <- names(sel)
  sel_fun <- lapply(fl, function(f) { s <- sel[[f]]; q <- relq[[f]]; function(L) q * sel_len_floss(s, L) / mx })
  names(sel_fun) <- fl
  sel_sub <- lapply(sel_fun, function(fn) array(fn(L_sub), dim(L_sub)))
  matL <- function(L) 1 / (1 + exp(-log(19) * (L - l$L50) / (l$L95 - l$L50)))
  mat_beg <- if (female) cbind(0, matL(L_beg)[, -1, drop = FALSE]) else L_beg * 0
  list(ages = ages, recP = recP, linf_g = linf_g, len = len, M = l$M, cv_res = gtg$cv_res,
       L50 = l$L50, L95 = l$L95, fleets = fl,
       w_beg = l$a * L_beg^l$b, w_sub = l$a * L_sub^l$b, nsub = nsub, tsub = tsub, L_sub = L_sub,
       mat_beg = mat_beg, sel_sub = sel_sub, sel = lapply(sel_sub, function(x) apply(x, c(1, 2), mean)),
       sel_fun = sel_fun, settings = gtg)
}

#' Growth morphs of a sex-structured stock
#'
#' Morph structures of both sexes for [om_sex()] with `gtg`: lengths, weights, maturity (females only) and fleet
#' selectivity at length by morph, age and sub-annual step. Fleet selectivity is scaled by `relq` and by the
#' maximum of the joint selectivity at length, so the apical F at length is the same for both sexes (the F of
#' [stock_sex()] and of length-based assessments such as FLicc).
#' @param stk a sex-structured stock from [stock_sex()].
#' @param gtg morph settings ([gtg_settings()]) or a list of its arguments; `TRUE` for the defaults.
#' @return a list with the morph structures `F` and `M` (class `floss_gtg_sex`).
#' @export
gtg_sex <- function(stk, gtg = gtg_settings()) {
  if (isTRUE(gtg)) gtg <- list()
  if (!all(c("cv_between", "nsub") %in% names(gtg)))
    gtg <- do.call(gtg_settings, gtg[intersect(names(gtg), c("n", "maxsd", "cv", "cv_res", "nsub"))])
  lh <- attr(stk, "lh"); sel <- attr(stk, "sel")
  if (is.null(lh) || is.null(sel)) stop("stk must come from stock_sex() (attributes 'lh' and 'sel')")
  relq <- attr(stk, "relq") %||% stats::setNames(rep(1, length(sel)), names(sel))
  ages <- as.numeric(dimnames(stk)$age)
  Lg <- seq(0.1, 1.5 * max(lh$F$linf, lh$M$linf) * (1 + 3 * gtg$cv), by = 0.05)
  mx <- max(Reduce(`+`, lapply(names(sel), function(f) relq[[f]] * sel_len_floss(sel[[f]], Lg))))
  structure(list(F = gtg_setup_sex1(lh$F, sel, relq, ages, gtg, TRUE, mx),
                 M = gtg_setup_sex1(lh$M, sel, relq, ages, gtg, FALSE, mx), settings = gtg),
            class = "floss_gtg_sex")
}

# fleet F vector at apical F x (the apical F at length is shared by the fleets' scaled selectivities)
ff_sex <- function(G, x) stats::setNames(rep(x, length(G$fleets)), G$fleets)

#' Sex-structured equilibrium with growth morphs
#'
#' Per-recruit and equilibrium quantities of a sex-structured stock with growth morphs at apical F (recruits
#' split equally between the sexes): female spawning biomass per recruit, SPR, yield per recruit of each sex and
#' in total, and with `sr` recruitment, spawning biomass and yield.
#' @param Gs morph structures from [gtg_sex()].
#' @param sr optional stock-recruitment model ([sr_sex()]).
#' @param F apical F values.
#' @return a data frame, as [eq_sex()].
#' @export
eq_sex_gtg <- function(Gs, sr = NULL, F = seq(0, 2, by = 0.01)) {
  sb0 <- 0.5 * gtg_eq(Gs$F, ff_sex(Gs$F, 0))$phi
  out <- do.call(rbind, lapply(F, function(f) {
    eF <- gtg_eq(Gs$F, ff_sex(Gs$F, f)); eM <- gtg_eq(Gs$M, ff_sex(Gs$M, f))
    sbpr <- 0.5 * eF$phi
    data.frame(F = f, sbpr = sbpr, spr = sbpr / sb0, ypr = 0.5 * (eF$ypr + eM$ypr),
               ypr_F = 0.5 * eF$ypr, ypr_M = 0.5 * eM$ypr)
  }))
  if (!is.null(sr)) {
    p <- params(sr); a <- an(p["a"]); b <- an(p["b"])
    out$rec <- pmax(a - b / out$sbpr, 0)
    out$ssb <- out$rec * out$sbpr
    out$yield <- out$rec * out$ypr
    out$yield_F <- out$rec * out$ypr_F; out$yield_M <- out$rec * out$ypr_M
  }
  out
}

#' Reference points of a sex-structured stock with growth morphs
#'
#' Fmsy, Bmsy, MSY, B0, R0, Fcrash, SPR at MSY, Bmsy/B0 and SPR-based F from the joint morph equilibrium
#' (female spawning biomass, yield of both sexes), on the apical-F scale of [om_sex()]. Weight, maturity and
#' effective selectivity at age change with F as fast growers are fished out, which [refpts_sex()] on a fixed
#' selectivity at age does not capture.
#' @param Gs morph structures from [gtg_sex()].
#' @param sr stock-recruitment model ([sr_sex()], with `gtg` so that B0 matches the morphs).
#' @param spr SPR levels for the SPR-based F.
#' @return an `FLPar`, as [refpts_sex()].
#' @export
refpts_sex_gtg <- function(Gs, sr, spr = c(0.3, 0.4)) {
  e <- function(f, s = sr) eq_sex_gtg(Gs, s, F = f)
  e0 <- e(0)
  fhi <- 1; while (e(fhi)$ssb > 0 && fhi < 50) fhi <- fhi * 2
  fcrash <- stats::uniroot(function(f) e(f)$ssb - 1e-6 * e0$ssb, c(0, fhi))$root
  fmsy <- stats::optimize(function(f) e(f)$yield, c(0, fcrash), maximum = TRUE, tol = 1e-8)$maximum
  em <- e(fmsy)
  fs <- sapply(spr, function(p) stats::uniroot(function(f) e(f, NULL)$spr - p, c(0, 50))$root)
  out <- c(Fmsy = fmsy, Bmsy = em$ssb, MSY = em$yield, B0 = e0$ssb, R0 = e0$rec,
           Fcrash = fcrash, SPRmsy = em$spr, Bmsy_B0 = em$ssb / e0$ssb)
  out <- c(out, stats::setNames(fs, paste0("Fspr", round(100 * spr))))
  do.call(FLPar, as.list(out))
}

# two-sex projection with growth morphs for om_sex(gtg = ): returns the FLStock (units F, M) and the morph tracks
om_sex_gtg_project <- function(stk, sr, Gs, fq, deviances, years) {
  ages <- Gs$F$ages; na <- length(ages); ny <- length(years); nsim <- dim(fq)[6]
  a <- an(params(sr)["a"]); b <- an(params(sr)["b"])
  sx <- c("F", "M")
  arr <- function() array(0, c(na, ny, 2, nsim), dimnames = list(age = ages, year = years, unit = sx, iter = seq_len(nsim)))
  N <- Fe <- W <- Mt <- CW <- CR <- arr()     # CR: morph catch numbers per number at the start of the year
  ng <- length(Gs$F$recP)
  trk <- lapply(sx, function(s) array(0, c(ng, na, ny, nsim)))
  names(trk) <- sx
  rec <- matrix(0, ny, nsim)
  for (i in seq_len(nsim)) {
    fy <- c(fq[, , , , , i])
    for (s in sx) {
      G <- Gs[[s]]
      Fy <- matrix(fy, ny, length(G$fleets), dimnames = list(NULL, G$fleets))
      tr <- gtg_track(G, Fy); ag <- gtg_aggregate(G, Fy, tr)
      Fe[, , s, i] <- ag$Feff; W[, , s, i] <- ag$wt; CW[, , s, i] <- ag$cwt_tot
      Mt[, , s, i] <- if (s == "F") ag$mat else 0
      trk[[s]][, , , i] <- tr$l
      npr <- apply(tr$l, c(2, 3), sum)
      cpr <- Reduce(`+`, lapply(ag$cn, function(x) apply(x, c(2, 3), sum)))
      CR[, , s, i] <- ifelse(npr > 0, cpr / npr, 0)
    }
    # first year: equilibrium at its F
    e1 <- eq_sex_gtg(Gs, sr, F = fy[1])
    rec[1, i] <- e1$rec
    for (s in sx) N[, 1, s, i] <- 0.5 * e1$rec * colSums(gtg_pr(Gs[[s]], ff_sex(Gs[[s]], fy[1]))$l)
    for (y in seq_len(ny)) {
      if (y > 1) {
        ssb <- sum(N[, y, "F", i] * W[, y, "F", i] * Mt[, y, "F", i])
        rec[y, i] <- (a * ssb / (b + ssb)) * c(deviances[, ac(years[y]), "F", , , i])
        N[1, y, , i] <- 0.5 * rec[y, i]
      }
      if (y < ny) for (s in sx) {
        z <- Fe[, y, s, i] + Gs[[s]]$M
        n <- N[, y, s, i] * exp(-z)
        N[2:na, y + 1, s, i] <- n[1:(na - 1)]
        N[na, y + 1, s, i] <- N[na, y + 1, s, i] + n[na]
      }
    }
  }
  Mq <- arr(); for (s in sx) Mq[, , s, ] <- Gs[[s]]$M
  # catch numbers from the morph catch (fleets select on each morph's length within the year), so catch biomass
  # (with the morph catch weight) matches the morph equilibrium; Baranov with the annual effective F would not
  CN <- N * CR
  dmn <- list(age = ages, year = years, unit = sx, season = "all", area = "unique", iter = seq_len(nsim))
  q <- function(x, u = "NA") FLQuant(array(x, c(na, ny, 2, 1, 1, nsim)), dimnames = dmn, units = u)
  st <- FLStock(stock.n = q(N), catch.n = q(CN), landings.n = q(CN), discards.n = q(0 * CN),
                m = q(Mq), mat = q(Mt), stock.wt = q(W), catch.wt = q(CW), landings.wt = q(CW),
                discards.wt = q(CW), harvest = q(Fe, "f"))
  harvest.spwn(st)[] <- 0; m.spwn(st)[] <- 0
  range(st, c("minfbar", "maxfbar")) <- range(stk)[c("minfbar", "maxfbar")]
  landings(st) <- computeLandings(st); discards(st) <- computeDiscards(st); catch(st) <- computeCatch(st)
  stock(st) <- computeStock(st)
  list(stock = st, tracks = trk, rec = rec)
}

# morph numbers (morph x age) of one sex, year and iteration: the projected numbers at age split by the morph
# composition of each cohort (exact for the plus group as well)
gtg_sex_numbers <- function(om, s, y, i) {
  y <- ac(y); yi <- match(y, ac(om$years))
  l <- om$gtg$tracks[[s]][, , yi, i]
  n <- c(stock.n(om$stock)[, y, s, , , i])
  sweep(l, 2, ifelse(colSums(l) > 0, n / colSums(l), 0), `*`)
}

# expected catch at length (fleet) or survey numbers at length of one sex with growth morphs
gtg_sex_len <- function(om, s, y, i, lens, fleet = NULL, survey_sel = NULL, timing = 0.5) {
  G <- om$gtg$G[[s]]; bin <- diff(lens[1:2]); y <- ac(y); Nm <- gtg_sex_numbers(om, s, y, i)
  if (!is.null(fleet)) {
    fy <- c(om$fapic[, y, , , , i])
    cy <- gtg_catch_sub(G, Nm, ff_sex(G, fy))[[fleet]]
    ld <- gtg_len_dist(G, lens, bin, fleet)
    out <- numeric(length(lens))
    for (st in seq_len(G$nsub)) for (a in seq_along(G$ages)) out <- out + as.numeric(cy[, a, st] %*% ld[[st]][[a]])
    return(out)
  }
  fy <- c(om$fapic[, y, , , , i])
  Z <- G$M + Reduce(`+`, lapply(G$fleets, function(f) G$sel[[f]] * fy))
  Nt <- Nm * exp(-timing * Z)
  Lm <- G$len(G$ages + timing); sd <- pmax(G$cv_res * Lm, 1e-3)
  sl <- sel_len_floss(survey_sel, lens + bin / 2)
  out <- numeric(length(lens))
  for (a in seq_along(G$ages)) {
    P <- stats::pnorm(outer(Lm[, a], lens + bin, "-") / -sd[, a]) - stats::pnorm(outer(Lm[, a], lens, "-") / -sd[, a])
    out <- out + as.numeric(Nt[, a] %*% P)
  }
  out * sl
}
