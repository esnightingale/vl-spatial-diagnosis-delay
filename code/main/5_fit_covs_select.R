################################################################################
# Description: Preliminary variable selection on patient characteristics, with
# non-spatial IID effects. 
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))

figdir <- "figures/fit/univariate"
outdir <- "output/univariate"

covs_all <- list("age_s","sex","comorb","prv_tx","caste4_r","occ4_cat", "poss_acd",
                 "block_endm_2017", "IRS_2017","inc_2017_gt0",
                 "traveltime_s","traveltime_t_s", "rain")
#c("age_cat.12.25.","age_cat.25.42.","age_cat.42.95.") 
covs_pat1 <- list("age_s","sex","comorb","prv_tx","caste4_r","occ4_cat") 
covs_pat <- c(covs_pat1, list("poss_acd"))
covs_aware <- list("block_endm_2017", "IRS_2017","inc_2017_gt0")
covs_accessD <- list("traveltime_s", "rain")
covs_accessT <- list("traveltime_t_s", "rain") 

covs.multi.list <- list(Patient = covs_pat1, 
                        `Patient + detection` = covs_pat, 
                        Awareness = covs_aware, 
                        `Access: diag` = covs_accessD, 
                        `Access: trt` = covs_accessT)

mesh <- readRDS(here::here("data/analysis","mesh.rds")) 
# stk <- readRDS(here::here("data/analysis","stack.rds"))

# Original sf object for plotting vgm
dat <- read_data() %>% 
  mutate(y = delay)

# Specify likelihood family = NB or Poisson
family <- "poisson"

# ---------------------------------------------------------------------------- #
# Null fit

f <- as.formula(paste0("y ~ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01))"))

fit.null <- inla(f,
                  family = family,
                  # data = inla.stack.data(stk),
                  data = dat,
                  control.predictor = list(
                    compute = TRUE, link = 1), #,A = inla.stack.A(stk)),
                  control.compute = list(waic = TRUE, 
                                         config = TRUE,
                                         cpo = FALSE,
                                         return.marginals.predictor = TRUE),
                  control.fixed = list(mean = 0, 
                                       prec = 0.1, 
                                       mean.intercept = 0, 
                                       prec.intercept = 0.1),
                  verbose = TRUE)
fit.null <- list(f = f, fit = fit.null)

p_vgm_null <- plot_vgm(fit.null$fit$summary.random$id$mean, dat, title = "Fitted IID effects: Null model",ylim = c(0.5,1))
#   model     psill    range kappa
# 1   Nug 0.5942502  0.00000   0.0
# 2   Mat 0.2564780 30.99597   0.5
# model     psill    range kappa
# 1   Nug 0.6265227  0.00000   0.0
# 2   Mat 0.1938768 18.59652   0.5

# Estimated variogram not sensitive to choice of prior precision

saveRDS(fit.null, here::here(outdir, "fit_null.rds"))

# ---------------------------------------------------------------------------- #
# Univariate fits

fit_covs <- function(cov) {
  
  # Define formula
  f <- as.formula(paste0("y ~", 
                         paste0(cov, collapse = " + "), 
                         "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01))"))
  
  print(f)
  
  # Fit model
  fit <- init_inla(f, data = dat, family = family)$fit #data.stack = stk, 
  
  res <- list(f = f, fit = fit)
  return(res)
  
}

fits.uni <- plyr::llply(covs_all, fit_covs)
# saveRDS(fits.uni, here::here(outdir, "fits_covs_univar_nonspatial.rds"))

# ---------------------------------------------------------------------------- #
# Domain fits

fits.domain <- plyr::llply(covs.multi.list, fit_covs)

p_vgms_domain <- purrr::imap(fits.domain, function(x, nm) plot_vgm(x$fit$summary.random$id$mean, dat, 
                                                                  title = paste("Fitted IID effects:", nm),
                                                                  ylim = c(0.5,1)))

# model     psill    range kappa
# 1   Nug 0.7174815  0.00000   0.0
# 2   Mat 0.1601273 21.44282   0.5
# model     psill    range kappa
# 1   Nug 0.7151361  0.00000   0.0
# 2   Mat 0.1484297 20.78582   0.5
# model     psill    range kappa
# 1   Nug 0.6839024  0.00000   0.0
# 2   Mat 0.1858052 20.16757   0.5
# model     psill    range kappa
# 1   Nug 0.6562865  0.00000   0.0
# 2   Mat 0.1909908 18.69911   0.5
# model     psill   range kappa
# 1   Nug 0.6662924  0.0000   0.0
# 2   Mat 0.1902892 18.9436   0.5

# Compare travel time to diagnosis vs treatment
plyr::llply(fits.domain, function(x) summary(x$fit))
plyr::llply(fits.domain, function(x) x$fit$waic$waic) 
# $Patient
# [1] 28573.9
# 
# $`Patient + detection`
# [1] 28574.22
# 
# $Awareness
# [1] 28581.63
# 
# $`Access: diag`
# [1] 28580.84
# 
# $`Access: trt`
# [1] 28580.26

# Treatment facility travel time has a somewhat larger coefficient, but neither 
# significant on 95% CrI
  
# ---------------------------------------------------------------------------- #
# Full multivariate fit

covs_full <- c(covs_pat, covs_aware, covs_accessT) 

f <- as.formula(paste0("y ~ ", 
                       paste0(covs_full, collapse = " + "), 
                       "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01))"))

fit.full <- init_inla(f, data = dat, family = family)

# saveRDS(fit.full, here::here(outdir, "fit_covs_multivar_nonspatial.rds"))

p_vgm_full <- plot_vgm(fit.full$fit$summary.random$id$mean, dat, title = "Fitted IID effects: Full model", ylim = c(0.5,1))
# model     psill    range kappa
# 1   Nug 0.7232051  0.00000   0.0
# 2   Mat 0.1319103 20.76699   0.5

# ---------------------------------------------------------------------------- #
# Compare univariate versus multivariate coefficients

fits.all <- rlist::list.append(c(fits.uni,
                                 setNames(fits.domain[-1], NULL)),
                               fit.full)
# rlist::list.save(fits.all, here::here(outdir, "fits_covs_full_nonspatial.rds"))

c("Intercept",
  "Age (std)",
  "Sex: Male",
  "HIV positive",
  "Prv. VL/PKDL",
  "Marginalised caste (SC/ST)",
  "Occupation: Unskilled",
  "Skilled",
  "Salaried/self-employed",
  "Detection: ACD",
  "Block endemic (2017)",
  "Village IRS targeting (2017)",
  "Village incidence > 0 (2017)",
  "Travel time - diagnosis facility (std)",
  "Travel time - treatment facility (std)",
  "Diagnosis season: rain"
) -> axis.labs

Efxplot(plyr::llply(fits.all, function(x) x$fit),
        ModelGroups = c(rep("Univariate",length(fits.uni)),
                        rep("Domain multivariate",length(fits.domain)-1),
                        "Full multivariate"),
        VarNames = rev(axis.labs),
        Intercept = FALSE,
        Size = 1.5,
        # exp = FALSE
        ) +
  geom_vline(xintercept = c(3.5, 6.5), col = "lightgrey", lwd = 0.5) +
  scale_colour_viridis_d(option = "plasma", direction = -1, end = 0.8) +
  theme(
    # legend.position = "none",
    axis.text.y = element_text(face = rev(c("plain","bold","plain","bold",rep("plain",5),"bold","bold","plain","bold","plain","bold","plain"))),
    # legend.text = element_text(face = c("plain","bold","plain"))
  ) + 
  labs(title = "Estimated covariate effects: non-spatial models",
       colour = "")

ggsave(here::here(figdir, "covs_domain_efx_nonspatial.png"), height = 6, width = 8, units = "in", dpi = 300)

# Look at spatial correlation after including patient characteristics and detection
png(here::here(figdir, "vgms_domain_nonspatial.png"), height = 400, width = 2500)
gridExtra::grid.arrange(grobs = list.append(list.prepend(p_vgms_domain,p_vgm_null), 
                                            p_vgm_full), nrow = 1)
dev.off()

# ---------------------------------------------------------------------------- #
# Multivariate fit with selection

covs_select <- c("age_s","comorb", "poss_acd",
                 "block_endm_2017","inc_2017_gt0","traveltime_t_s") 

f <- as.formula(paste0("y ~ ", 
                       paste0(covs_select, collapse = " + "), 
                       "+ f(id, model = 'iid',
                             prior = 'pc.prec', 
                             param = c(1, 0.01))"))

fit.select <- init_inla(f, data = dat, family = family)

# saveRDS(fit.select, here::here(outdir, "fit_covs_multivar_select_nonspatial.rds"))

p_vgm_select <- plot_vgm(fit.select$fit$summary.random$id$mean, dat, title = "Fitted IID effects: Full model with selection", ylim = c(0.5,1))
# model     psill    range kappa
# 1   Nug 0.7225837  0.00000   0.0
# 2   Mat 0.1331138 19.94885   0.5

png(here::here(figdir, "vgms_select_nonspatial.png"), height = 1000, width = 500)
gridExtra::grid.arrange(p_vgm_null, p_vgm_select, nrow = 2)
dev.off()

# ---------------------------------------------------------------------------- #
# Compare univariate, multivariate and multivariate with selection

fits.all <- rlist::list.append(c(fits.uni, 
                                 setNames(fits.domain[-1], NULL)),
                               fit.select)
rlist::list.save(fits.all, here::here(outdir, "fits_covs_select_all_nonspatial.rds"))

c("Intercept",
  "Age (std)",
  "Sex: Male",
  "HIV positive",
  "Prv. VL/PKDL treatment",
  "Marginalised caste (SC/ST)",
  "Occupation: Unskilled",
  "Skilled",
  "Salaried/self-employed",
  "Detection: ACD",
  "Block endemic (2017)",
  "Village IRS targeting (2017)",
  "Village incidence > 0 (2017)",
  "Travel time - diagnosis facility (std)",
  "Travel time - treatment facility (std)",
  "Diagnosis season: rain"
) -> axis.labs


png(here::here(figdir, "covs_select_efx_nonspatial.png"), height = 8, width = 8, unit = "in", res = 320)
Efxplot(plyr::llply(fits.all, function(x) x$fit),
                   ModelGroups = c(rep("Univariate",length(fits.uni)), 
                                   rep("Domain multivariate",length(fits.domain)-1), 
                                   "Selected multivariate"),
                   VarNames = rev(axis.labs),
                   Intercept = FALSE, Size = 1.5) +
  geom_vline(xintercept = c(3.5, 6.5), col = "lightgrey", lwd = 1) +
  scale_colour_viridis_d(option = "plasma", end= 0.8, direction = -1) +
  scale_y_continuous(breaks = seq(0.6,1.8, by = 0.2)) +
  theme(legend.position = c(0.75,0.2),
        legend.title = element_blank(),
        legend.box.background = element_rect(color="white", size=2),
        axis.text.y = element_text(face = rev(c("plain","bold","plain","bold",rep("plain",5),"bold","bold","plain","bold","plain","bold","plain"))),
        # legend.text = element_text(face = c("plain","bold","plain"))
        ) + 
  labs(#title = "Estimated covariate effects: non-spatial models",
       # colour = "",
       caption = "Selected covariates from each domain are highlighted in bold.")
dev.off()

# ---------------------------------------------------------------------------- #
# Estimated coefficients

fits.all[[1]]$fit -> fit_age
fits.all[[3]]$fit -> fit_hiv
fits.all[[7]]$fit -> fit_acd

fits.all[[8]]$fit -> fit_end
fits.all[[10]]$fit -> fit_inc

fits.all[[12]]$fit -> fit_trt

fits.all[[14]]$fit -> fit_pat
fits.all[[15]]$fit -> fit_awa
fits.all[[17]]$fit -> fit_acc

fits.all[[18]]$fit -> fit_sel

lapply(list(fit_age, fit_hiv, fit_acd, 
            fit_end, fit_inc, 
            fit_trt,
            fit_pat, fit_awa, fit_acc,
            fit_sel), 
       get_coeffs)

# [[1]]
# effect       mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#   1 (Intercept) 19.1  0.296       18.5      18.9     19.1      19.3       19.6 
#   2 age_s        1.14 0.0178       1.11      1.13     1.14      1.15       1.18
# 
# [[2]]
# effect       mean    sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#   1 (Intercept) 18.7  0.298      18.1      18.5     18.7      18.9       19.3 
#   2 comorb1      1.54 0.126       1.31      1.46     1.54      1.63       1.81
# 
# [[3]]
# effect         mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#   1 (Intercept)  21.5   0.428      20.7      21.2     21.5      21.8       22.3  
#   2 poss_acdTRUE  0.738 0.0233      0.693     0.722    0.738     0.754      0.785
#
# [[4]]
# effect                mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#   1 (Intercept)         20.4   0.430      19.6      20.1     20.4      20.7       21.3  
#   2 block_endm_2017TRUE  0.854 0.0267      0.803     0.835    0.853     0.872      0.907
# 
# [[5]]
# effect             mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#   1 (Intercept)      20.7   0.480      19.8      20.4     20.7      21.0       21.7  
#   2 inc_2017_gt0TRUE  0.855 0.0267      0.804     0.836    0.854     0.873      0.908
# 
# [[6]]
# effect          mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#   1 (Intercept)    19.0  0.298      18.4       18.8     19.0      19.2       19.6 
#   2 traveltime_t_s  1.02 0.0159      0.991      1.01     1.02      1.03       1.05
# 
# [[7]]
# effect         mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#   1 (Intercept)  21.6   0.641      20.3      21.1     21.6      22.0       22.9  
# 2 age_s         1.14  0.0187      1.10      1.12     1.14      1.15       1.17 
# 3 sex1          0.947 0.0320      0.886     0.925    0.947     0.969      1.01 
# 4 comorb1       1.37  0.113       1.16      1.29     1.36      1.44       1.60 
# 5 prv_txTRUE    0.959 0.0527      0.859     0.922    0.957     0.993      1.07 
# 6 caste4_r1     1.06  0.0349      0.995     1.04     1.06      1.09       1.13 
# 7 occ4_cat1     0.982 0.0383      0.909     0.956    0.981     1.01       1.06 
# 8 occ4_cat2     1.03  0.0681      0.901     0.981    1.03      1.07       1.17 
# 9 occ4_cat3     0.992 0.0650      0.871     0.947    0.990     1.04       1.13 
# 10 poss_acdTRUE  0.749 0.0238      0.703     0.732    0.748     0.764      0.796
# 
# [[8]]
# effect                mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#   1 (Intercept)         21.7   0.722      20.3      21.2     21.7      22.2       23.1  
# 2 block_endm_2017TRUE  0.885 0.0290      0.829     0.865    0.884     0.904      0.943
# 3 IRS_2017TRUE         0.988 0.0409      0.910     0.959    0.987     1.01       1.07 
# 4 inc_2017_gt0TRUE     0.890 0.0319      0.829     0.868    0.889     0.911      0.954
# 
# [[9]]
# effect          mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#   1 (Intercept)    18.9  0.360      18.2       18.6     18.9      19.1       19.6 
# 2 traveltime_t_s  1.02 0.0159      0.990      1.01     1.02      1.03       1.05
# 3 rainTRUE        1.03 0.0342      0.961      1.00     1.03      1.05       1.09
# 
# [[10]]
# effect                mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#   1 (Intercept)         23.6   0.659      22.3      23.1     23.6      24.0       24.9  
# 2 age_s                1.13  0.0175      1.09      1.12     1.13      1.14       1.16 
# 3 comorb1              1.28  0.104       1.08      1.20     1.27      1.34       1.49 
# 4 poss_acdTRUE         0.748 0.0236      0.702     0.731    0.747     0.763      0.795
# 5 block_endm_2017TRUE  0.860 0.0280      0.806     0.840    0.859     0.878      0.916
# 6 inc_2017_gt0TRUE     0.937 0.0301      0.880     0.917    0.937     0.957      0.998
# 7 traveltime_t_s       1.02  0.0159      0.993     1.01     1.02      1.03       1.06 


################################################################################
################################################################################
