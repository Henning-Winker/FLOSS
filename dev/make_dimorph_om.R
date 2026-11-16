# Builds data(dimorph_om): the sex-structured operating model used in the FLOSS guides,
# a dimorphic stock with a trawl fishery (unsexed samples) and a survey (sexed samples),
# with the life history and sampling of FLicc's data(dimorph).
# Run from the package root:  source("dev/make_dimorph_om.R")
library(FLOSS); library(FLasher)
set.seed(1234)
hist_years <- 1976:2025
nsim <- 10
lh <- lh_sex(list(linf = 60, k = 0.15, t0 = -0.3, a = 1e-5, b = 3, L50 = 35, M = 0.3),
             linf = 45, k = 0.25, L50 = 28)
sel <- list(Trawl = selectivity("dnormal", Lpeak = 34, sd_left = 6, sd_right = 20))
srv_sel <- selectivity("logistic", L50 = 12, L95 = 16)
stk <- stock_sex(lh, sel, tmax = 20, years = hist_years)
sr <- sr_sex(stk, s = 0.7, R0 = 1e6)
rp <- refpts_sex(stk, sr)
f_hist <- fapic_fudc(hist_years, nsim, fref = rp["Fmsy"], f0 = 0.1, fhi = 2.5, flo = 1.5,
                     sigmaF = 0.15, breaks = c(0.6, 0.8), spread = 0.2, seed = 1234)
devs <- rlnormar1(rho = 0.3, years = hist_years, n = nsim, meanlog = 0, sdlog = 0.5)
om <- om_sex(stk, sr, fapic = f_hist, deviances = devs)
lfd_years <- 2011:2025
lfd <- FLQuants(Trawl = lfd_sex(om, "Trawl", n = 400, years = lfd_years, iter = 1),
                Survey = lfd_sex(om, survey_sel = srv_sel, n = 250, sexed = TRUE,
                                 years = lfd_years, iter = 1))
lfd$Survey[, ac(2011:2015)] <- 0
dimorph_om <- list(om = om, lh = lh, sel = sel, survey_sel = srv_sel, lfd = lfd,
                   lhpar = lhpar_flicc_sex(lh), truth = om$truth[om$truth$iter == "1", ])
save(dimorph_om, file = "data/dimorph_om.rda", compress = "xz")
