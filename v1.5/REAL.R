#' @title Estimating selection coefficients and testing their changes from ancient DNA data II: extension for modelling sampling uncertainties
#' @author Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.1
#' Phenotypes controlled by a single gene
#' Non-constant natural selection and non-constant demographic histories
#' Prior knowledge from modern samples (gene polymorphism)
#' Joint estimation of the underlying trajectory of mutant allele frequencies and unknown alleles

#' Genotype frequency data

#' Horse base coat colours (ASIP & MC1R) and white coat patterns (KIT13 & KIT16 & TRMP1)

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2020/HE2021-WFM-2L-DiffusApprox-PMMH2-MolEcol")

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
source("./Code/Code v1.0/Code 1L/Code v1.1/RFUN.R")

################################################################################

#' ASIP

#' Grouped data of Wutke et al. (2016) from 12500 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- ASIP
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
max(raw_smp$age_mean[which(rowSums(raw_smp[, 5:9]) != 0)])
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, 5:9]) != 0)])), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_ASIP.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_ASIP.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_ASIP_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_ASIP_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_ASIP_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_ASIP_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

########################################

#' Raw data of Wutke et al. (2016) from 9700 BC (Holocene)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- ASIP
raw_smp <- raw_smp[which(raw_smp$age_mean <= 9700 + 2000), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_ASIP.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_ASIP.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_ASIP_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_ASIP_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_ASIP_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_ASIP_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

############################################################

#' MC1R

#' Grouped data of Wutke et al. (2016) from 4300 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- MC1R
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
max(raw_smp$age_mean[which(rowSums(raw_smp[, 5:9]) != 0)])
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, 5:9]) != 0)])), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_MC1R.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_MC1R.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_MC1R_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_MC1R_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_MC1R_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMHg_MC1R_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

########################################

#' Raw data of Wutke et al. (2016) from 9700 BC (Holocene)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- MC1R
raw_smp <- raw_smp[which(raw_smp$age_mean <= 9700 + 2000), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_MC1R.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_MC1R.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_MC1R_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before domestication")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_MC1R_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_MC1R_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_COL_PMMH1g_MC1R_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

############################################################

#' KIT13

#' Raw data of Wutke et al. (2016) from 3645 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT13
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
max(raw_smp$age_mean[which(rowSums(raw_smp[, 5:9]) != 0)])
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, 5:9]) != 0)])), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT13.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT13.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT13_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT13_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT13_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT13_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

########################################

#' Raw data of Wutke et al. (2016) from 3500 BC (Domestication)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT13
raw_smp <- raw_smp[which(raw_smp$age_mean <= 3500 + 2000), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT13.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT13.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT13_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT13_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT13_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT13_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

############################################################

#' KIT16

#' Raw data of Wutke et al. (2016) from 3500 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT16
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
max(raw_smp$age_mean[which(rowSums(raw_smp[, 5:9]) != 0)])
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, 5:9]) != 0)])), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT16.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT16.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT16_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT16_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT16_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_KIT16_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

########################################

#' Raw data of Wutke et al. (2016) from 3500 BC (Domestication 3500 BC)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT16
raw_smp <- raw_smp[which(raw_smp$age_mean <= 3500 + 2000), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT16.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT16.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT16_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the Middle Ages")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT16_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT16_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_KIT16_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

############################################################

#' TRPM1

#' Raw data of Wutke et al. (2016) from 14500 BC
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- TRPM1
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
max(raw_smp$age_mean[which(rowSums(raw_smp[, 5:9]) != 0)])
raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, 5:9]) != 0)])), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-1600 - 2000) / 8) # 1600 BC (the early Bronze Age)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_TRPM1.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_TRPM1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_TRPM1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the early Bronze Age")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the early Bronze Age")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_TRPM1_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_TRPM1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMHg_TRPM1_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

#######################################

#' Raw data of Wutke et al. (2016) from 9700 BC (Holocene)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- TRPM1
raw_smp <- raw_smp[which(raw_smp$age_mean <= 9700 + 2000), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_TRPM1.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_TRPM1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_TRPM1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the early Bronze Age")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the early Bronze Age")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_TRPM1_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_TRPM1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH1g_TRPM1_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

####################

#' Raw data of Wutke et al. (2016) from 3500 BC (Domestication)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- TRPM1
raw_smp <- raw_smp[which(raw_smp$age_mean <= 3500 + 2000), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1
pop_siz
ref_siz <- 1.6e+04
evt_gen <- round((-1600 - 2000) / 8) # 1600 BC (the early Bronze Age)
raw_smp <- raw_smp[, -(2:3)]
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH2g_TRPM1.rda")

load("./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH2g_TRPM1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
imp_smp_chn <- PMMH$imp_smp_chn

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH2g_TRPM1_Traceplot_SelCoeff.pdf", width = 12, height = 12)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the early Bronze Age")

plot(1:itn_num, sel_cof_chn[2, 1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient before the early Bronze Age")
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH2g_TRPM1_Posterior_SelCoeff.pdf", width = 24, height = 12)
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

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH2g_TRPM1_Posterior_SelChange.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(dif_sel_chn, breaks = seq(min(dif_sel_chn), max(dif_sel_chn), length.out = 50), freq = FALSE,
     xlab = "Change in selection coefficient",
     main = "Posterior for change in selection coefficient")
lines(density(dif_sel_chn), lwd = 2, col = 'black')
abline(v = dif_sel_est, col = 'black', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = dif_sel_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]

# MMSE estimate for selection coefficients and mutant allele frequencies
frq_pth_est <- colMeans(frq_pth_chn)
# grd_num <- 1e+03
# frq_pth_est <- rep(NA, length.out = dim(frq_pth_chn)[2])
# for (i in 1:dim(frq_pth_chn)[2]) {
#    frq_pth_pdf <- density(frq_pth_chn[, i], n = grd_num)
#    frq_pth_est[i] <- frq_pth_pdf$x[which(frq_pth_pdf$y == max(frq_pth_pdf$y))]
# }

# 95% HPD interval for selection coefficients and mutant allele frequencies
frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
for (i in 1:dim(frq_pth_chn)[2]) {
   frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
}

pdf(file = "./Output/Output v1.0/REAL v1.1/REAL_PTN_PMMH2g_TRPM1_Posterior_Traj.pdf", width = 12, height = 6)
par(mar = c(5.1, 5.1, 4.1, 1.1), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
smp_gen <- min(unique(raw_smp$age_mean)):max(unique(raw_smp$age_mean))
plot(0, type = 'n', xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Posterior for underlying trajectory of mutant allele")

for (i in 1:dim(frq_pth_chn)[1]) {
   lines(min(smp_gen):max(smp_gen), frq_pth_chn[i, ], col = 'grey', lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = 'black', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, ], col = 'blue', lty = 2, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, ], col = 'blue', lty = 2, lwd = 2)
dev.off()

# burn-in and thinning
imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

# posterior probability for genotypes
imp_smp_est <- matrix(NA, nrow = nrow(raw_smp), ncol = 3)
for (i in 1:nrow(imp_smp_est)) {
   imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
}

################################################################################
