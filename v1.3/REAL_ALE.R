#' @title Estimating selection coefficients and testing their changes from ancient DNA data
#' @author Xiaoyang Dai, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.3
#' Phenotypes controlled by a single gene
#' Non-constant natural selection and non-constant demographic histories
#' Prior knowledge from modern samples (gene polymorphism)

#' Allele frequency data

#' Horse base coat colours (ASIP & MC1R) and white coat patterns (KIT13 & KIT16 & TRMP1)

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2019/HE2021-WFM-2L-DiffusApprox-PMMH1-MolEcol")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

#install.packages("emdbook")
library("emdbook")

# call R functions
source("./Code/Code v1.0/Code 1L/Code v1.3/RFUN_ALE.R")

################################################################################

#' ASIP

#' Grouped data of Wutke et al. (2016) from 21520 BC
load("./Data/REAL_GRP_COL.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
smp_gen
smp_gen * 8 + 2000
smp_siz
smp_cnt <- smp_ale_cnt_ASIP[2, ]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_COL_PMMHa_ASIP.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_GRP_COL_PMMHa_ASIP.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_COL_PMMHa_ASIP_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_COL_PMMHa_ASIP_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_COL_PMMHa_ASIP_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 12496 BC
load("./Data/REAL_RAW_COL.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0
pop_siz <- pop_siz[-(1:(smp_gen[4] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
smp_gen <- smp_gen[-(1:3)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:3)]
smp_cnt <- smp_ale_cnt_ASIP[2, -(1:3)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMHa_ASIP.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMHa_ASIP.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMHa_ASIP_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMHa_ASIP_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMHa_ASIP_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_RAW_COL.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0
pop_siz <- pop_siz[-(1:(smp_gen[4] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
smp_gen <- smp_gen[-(1:4)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:4)]
smp_cnt <- smp_ale_cnt_ASIP[2, -(1:4)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMH1a_ASIP.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMH1a_ASIP.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMH1a_ASIP_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMH1a_ASIP_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMH1a_ASIP_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' MC1R

#' Grouped data of Wutke et al. (2016) from 21520 BC
load("./Data/REAL_GRP_COL.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
smp_gen
smp_gen * 8 + 2000
smp_siz
smp_cnt <- smp_ale_cnt_MC1R[2, ]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_COL_PMMHa_MC1R.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_GRP_COL_PMMHa_MC1R.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_COL_PMMHa_MC1R_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_COL_PMMHa_MC1R_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_COL_PMMHa_MC1R_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 12496 BC
load("./Data/REAL_RAW_COL.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0
pop_siz <- pop_siz[-(1:(smp_gen[4] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
smp_gen <- smp_gen[-(1:3)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:3)]
smp_cnt <- smp_ale_cnt_MC1R[2, -(1:3)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMHa_MC1R.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMHa_MC1R.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMHa_MC1R_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMHa_MC1R_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMHa_MC1R_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_RAW_COL.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0
pop_siz <- pop_siz[-(1:(smp_gen[4] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
smp_gen <- smp_gen[-(1:4)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:4)]
smp_cnt <- smp_ale_cnt_MC1R[2, -(1:4)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMH1a_MC1R.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMH1a_MC1R.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMH1a_MC1R_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMH1a_MC1R_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_COL_PMMH1a_MC1R_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' KIT13

#' Grouped data of Wutke et al. (2016) from 3552 BC
load("./Data/REAL_GRP_PTN.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz <- pop_siz[-(1:(smp_gen[3] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round(400 - 2000) / 8 # 400 AD (the Middle Ages)
smp_gen <- smp_gen[-(1:2)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:2)]
smp_cnt <- smp_ale_cnt_KIT13[2, -(1:2)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_KIT13.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_KIT13.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_KIT13_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after the Middle Ages")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_KIT13_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient after the Middle Ages",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the Middle Ages")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after the Middle Ages")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_KIT13_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 3648 BC
load("./Data/REAL_RAW_PTN.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz <- pop_siz[-(1:(smp_gen[14] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round(400 - 2000) / 8 # 400 AD (the Middle Ages)
smp_gen <- smp_gen[-(1:13)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:13)]
smp_cnt <- smp_ale_cnt_KIT13[2, -(1:13)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_KIT13.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_KIT13.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_KIT13_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after the Middle Ages")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_KIT13_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient after the Middle Ages",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the Middle Ages")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after the Middle Ages")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_KIT13_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 3504 BC (Domestication 3500 BC)
load("./Data/REAL_RAW_PTN.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz <- pop_siz[-(1:(smp_gen[14] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round(400 - 2000) / 8 # 400 AD (the Middle Ages)
smp_gen <- smp_gen[-(1:18)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:18)]
smp_cnt <- smp_ale_cnt_KIT13[2, -(1:18)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_KIT13.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_KIT13.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_KIT13_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after the Middle Ages")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_KIT13_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient after the Middle Ages",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the Middle Ages")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after the Middle Ages")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_KIT13_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' KIT16

#' Grouped data of Wutke et al. (2016) from 3552 BC
load("./Data/REAL_GRP_PTN.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz <- pop_siz[-(1:(smp_gen[3] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round(400 - 2000) / 8 # 400 AD (the Middle Ages)
smp_gen <- smp_gen[-(1:2)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:2)]
smp_cnt <- smp_ale_cnt_KIT16[2, -(1:2)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_KIT16.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_KIT16.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_KIT16_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after the Middle Ages")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_KIT16_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient after the Middle Ages",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the Middle Ages")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after the Middle Ages")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_KIT16_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 3648 BC
load("./Data/REAL_RAW_PTN.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz <- pop_siz[-(1:(smp_gen[14] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round(400 - 2000) / 8 # 400 AD (the Middle Ages)
smp_gen <- smp_gen[-(1:13)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:13)]
smp_cnt <- smp_ale_cnt_KIT16[2, -(1:13)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_KIT16.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_KIT16.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_KIT16_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after the Middle Ages")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_KIT16_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient after the Middle Ages",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the Middle Ages")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after the Middle Ages")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_KIT16_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 3504 BC (Domestication 3500 BC)
load("./Data/REAL_RAW_PTN.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz <- pop_siz[-(1:(smp_gen[14] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round(400 - 2000) / 8 # 400 AD (the Middle Ages)
smp_gen <- smp_gen[-(1:18)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:18)]
smp_cnt <- smp_ale_cnt_KIT16[2, -(1:18)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_KIT16.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_KIT16.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_KIT16_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after the Middle Ages")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_KIT16_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the Middle Ages", ylab = "Selection coefficient after the Middle Ages",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the Middle Ages")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after the Middle Ages")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_KIT16_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

############################################################

#' TRMP1

#' Grouped data of Wutke et al. (2016) from 32848 BC
load("./Data/REAL_GRP_PTN_LP.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-1600 - 2000) / 8) # 1600 BC (the early Bronze Age)
smp_gen
smp_gen * 8 + 2000
smp_siz
smp_cnt <- smp_ale_cnt_TRPM1[2, ]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_TRPM1.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_TRPM1.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_TRPM1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the early Bronze Age")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after the early Bronze Age")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_TRPM1_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the early Bronze Age", ylab = "Selection coefficient after the early Bronze Age",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the early Bronze Age")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after the early Bronze Age")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_GRP_PTN_PMMHa_TRPM1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Raw data of Wutke et al. (2016) from 14496 BC
load("./Data/REAL_RAW_PTN_LP.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz <- pop_siz[-(1:(smp_gen[3] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round((-1600 - 2000) / 8) # 1600 BC (the early Bronze Age)
smp_gen <- smp_gen[-(1:2)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:2)]
smp_cnt <- smp_ale_cnt_TRPM1[2, -(1:2)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_TRPM1.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_TRPM1.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_TRPM1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the early Bronze Age")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after the early Bronze Age")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_TRPM1_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the early Bronze Age", ylab = "Selection coefficient after the early Bronze Age",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the early Bronze Age")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after the early Bronze Age")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMHa_TRPM1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

#######################################

#' Raw data of Wutke et al. (2016) from 9320 BC (Holocene 9700 BC)
load("./Data/REAL_RAW_PTN_LP.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz <- pop_siz[-(1:(smp_gen[3] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
smp_gen <- smp_gen[-(1:4)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:4)]
smp_cnt <- smp_ale_cnt_TRPM1[2, -(1:4)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_TRPM1.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_TRPM1.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_TRPM1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after domestication")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_TRPM1_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before domestication", ylab = "Selection coefficient after domestication",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before domestication")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after domestication")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH1a_TRPM1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

####################

#' Raw data of Wutke et al. (2016) from 3504 BC (Domestication 3500 BC)
load("./Data/REAL_RAW_PTN_LP.rda")

set.seed(1)

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz <- pop_siz[-(1:(smp_gen[3] - smp_gen[1]))]
ref_siz <- 1.6e+04
evt_gen <- round((-1600 - 2000) / 8) # 1600 BC (the early Bronze Age)
smp_gen <- smp_gen[-(1:17)]
smp_gen * 8 + 2000
smp_siz <- smp_siz[-(1:17)]
smp_cnt <- smp_ale_cnt_TRPM1[2, -(1:17)]
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

# system.time(sel_cof_chn <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))
#
# save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, sel_cof_chn,
#      file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH2a_TRPM1.rda")

load("./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH2a_TRPM1.rda")

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH2a_TRPM1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the early Bronze Age")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient after the early Bronze Age")
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

sel_cof_est <- rowMeans(sel_cof_chn)

sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH2a_TRPM1_Posterior_SelCoeff.pdf", width = 24, height = 12)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
      xlab = "Selection coefficient before the early Bronze Age", ylab = "Selection coefficient after the early Bronze Age",
      main = "Posterior for selection coefficient")
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
# HPDregionplot(as.mcmc(sel_cof_chn), vars = 1:2, n = grd_num, prob = 0.95, col = "blue", lwd = 2, add = TRUE)

hist(sel_cof_chn[1, ], breaks = seq(min(sel_cof_chn[1, ]), max(sel_cof_chn[1, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient before the early Bronze Age")
lines(density(sel_cof_chn[1, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[1], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1, 2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_chn[2, ], breaks = seq(min(sel_cof_chn[2, ]), max(sel_cof_chn[2, ]), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient after the early Bronze Age")
lines(density(sel_cof_chn[2, ]), lwd = 2, col = 'black')
abline(v = sel_cof_est[2], col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2, 2], col = 'blue', lty = 2, lwd = 2)
dev.off()

dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

dif_sel_est <- mean(dif_sel_chn)

dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/REAL v1.3/REAL_RAW_PTN_PMMH2a_TRPM1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
