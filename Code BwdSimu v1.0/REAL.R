#' @title Estimating selection coefficients and testing their changes from ancient DNA data II
#' @author Zhangyi He, Wenyang Lyu, Xiaoyang Dai, Mark Beaumont, Feng Yu

#' Phenotypes controlled by a single gene
#' Non-constant natural selection and non-constant demographic histories

#' Input: genotype likelihoods
#' Output: posteriors for the selection coefficient and the genotype frequency trajectories of the population

#' Horse base coat colours (ASIP & MC1R) and white spotting patterns (KIT13 & KIT16 & TRPM1)

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2021/





WL2023-PopGenet-StatAdv-BPE-WFM1-PMMHwGibbs-MolBiolEvol")

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
source("./Code/Code v1.0/REAL_RFUN.R")

################################################################################

#' ASIP

#' Raw data of Wutke et al. (2016)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- ASIP
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
# raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
# raw_smp[which(raw_smp[, 7] == 1), 4:5] <- 1 / 2
# raw_smp[which(raw_smp[, 8] == 1), 5:6] <- 1 / 2
# raw_smp[which(raw_smp[, 9] == 1), 4:6] <- 1 / 3
raw_smp[which(raw_smp[, 7] == 1), 4] <- 1 / 2
raw_smp[which(raw_smp[, 7] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 6] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 4] <- 1 / 4
raw_smp[which(raw_smp[, 9] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 6] <- 1 / 4
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
n_x <- round(2 * max(pop_siz) * 0.02) + 1
n_t <- round(n_x * n_x / 8 / min(pop_siz)) + 1
ptn_num <- c(n_t, n_x - 24, 26)
unif_grd <- TRUE
# coord_x <- calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)
mut_rat <- c(7.24e-09, 7.24e-09)
pcl_num <- 1e+03
itn_num <- 2e+01
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptBackPMMH(mut_rat, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/MANU_COL_ASIP_REAL1.rda")

########################################

#' Raw data of Wutke et al. (2016)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- ASIP
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
# raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
# raw_smp[which(raw_smp[, 7] == 1), 4:5] <- 1 / 2
# raw_smp[which(raw_smp[, 8] == 1), 5:6] <- 1 / 2
# raw_smp[which(raw_smp[, 9] == 1), 4:6] <- 1 / 3
raw_smp[which(raw_smp[, 7] == 1), 4] <- 1 / 2
raw_smp[which(raw_smp[, 7] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 6] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 4] <- 1 / 4
raw_smp[which(raw_smp[, 9] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 6] <- 1 / 4
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- max(raw_smp$age_mean)
n_x <- round(2 * max(pop_siz) * 0.02) + 1
n_t <- round(n_x * n_x / 8 / min(pop_siz)) + 1
ptn_num <- c(n_t, n_x - 24, 26)
unif_grd <- TRUE
coord_x <- calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)
mut_rat <- c(7.24e-09, 7.24e-09)
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptBackPMMH(mut_rat, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/MANU_COL_ASIP_REAL2.rda")

########################################

load("./Output/Output v1.0/MANU_COL_ASIP_REAL1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
sel_cof_est <- c(sel_cof_pdf$x[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[1]], sel_cof_pdf$y[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[2]])
frq_pth_est <- colMeans(frq_pth_chn)

#' Show the call rate and accuracy accoss different data qualities and call thresholds
model <- "WFM"
sel_cof <- sel_cof_est
dom_par
pop_siz
int_con <- frq_pth_est[1]
evt_gen
smp_lab <- rownames(raw_smp)
smp_gen <- raw_smp$age_mean
smp_qua1 <- (0:100) / 100
smp_qua2 <- (1:50) / 10
thr_val <- (1:5) * 2
sim_num <- 500

sim_cal_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val), sim_num))
sim_err_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val), sim_num))
avg_cal_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val)))
avg_err_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val)))
for (i in 1:length(smp_qua1)) {
  print(i)
  for (j in 1:length(smp_qua2)) {
    print(j)
    for (k in 1:length(thr_val)) {
      print(k)
      for (r in 1:sim_num) {
        sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, c(smp_qua1[i], smp_qua2[j]), thr_val[k])
        sim_cal_rat[i, j, k, r] <- sim_HMM_WFM$cal_rat
        sim_err_rat[i, j, k, r] <- sim_HMM_WFM$err_rat
      }
      avg_cal_rat[i, j, k] <- mean(sim_cal_rat[i, j, k, ])
      avg_err_rat[i, j, k] <- mean(sim_err_rat[i, j, k, ])
    }
  }
}

save(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua1, smp_qua2, thr_val, sim_num, sim_cal_rat, sim_err_rat, avg_cal_rat, avg_err_rat,
     file = "./Output/Output v1.0/MANU_COL_ASIP_SIMU_CALL.rda")

############################################################

#' MC1R

#' Raw data of Wutke et al. (2016)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- MC1R
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
# raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
# raw_smp[which(raw_smp[, 7] == 1), 4:5] <- 1 / 2
# raw_smp[which(raw_smp[, 8] == 1), 5:6] <- 1 / 2
# raw_smp[which(raw_smp[, 9] == 1), 4:6] <- 1 / 3
raw_smp[which(raw_smp[, 7] == 1), 4] <- 1 / 2
raw_smp[which(raw_smp[, 7] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 6] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 4] <- 1 / 4
raw_smp[which(raw_smp[, 9] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 6] <- 1 / 4
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
n_x <- round(2 * max(pop_siz) * 0.02) + 1
n_t <- round(n_x * n_x / 8 / min(pop_siz)) + 1
ptn_num <- c(n_t, n_x - 24, 26)
unif_grd <- TRUE
coord_x <- calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)
mut_rat <- c(7.24e-09, 7.24e-09)
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptBackPMMH(mut_rat, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/MANU_COL_MC1R_REAL1.rda")

########################################

#' Raw data of Wutke et al. (2016)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- MC1R
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
# raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
# raw_smp[which(raw_smp[, 7] == 1), 4:5] <- 1 / 2
# raw_smp[which(raw_smp[, 8] == 1), 5:6] <- 1 / 2
# raw_smp[which(raw_smp[, 9] == 1), 4:6] <- 1 / 3
raw_smp[which(raw_smp[, 7] == 1), 4] <- 1 / 2
raw_smp[which(raw_smp[, 7] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 6] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 4] <- 1 / 4
raw_smp[which(raw_smp[, 9] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 6] <- 1 / 4
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 0e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- max(raw_smp$age_mean)
raw_smp
n_x <- round(2 * max(pop_siz) * 0.02) + 1
n_t <- round(n_x * n_x / 8 / min(pop_siz)) + 1
ptn_num <- c(n_t, n_x - 24, 26)
unif_grd <- TRUE
coord_x <- calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)
mut_rat <- c(7.24e-09, 7.24e-09)
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptBackPMMH(mut_rat, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/MANU_COL_MC1R_REAL2.rda")

########################################

load("./Output/Output v1.0/MANU_COL_MC1R_REAL1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
sel_cof_est <- c(sel_cof_pdf$x[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[1]], sel_cof_pdf$y[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[2]])
frq_pth_est <- colMeans(frq_pth_chn)

#' Show the call rate and accuracy accoss different data qualities and call thresholds
model <- "WFM"
sel_cof <- sel_cof_est
dom_par
pop_siz
int_con <- frq_pth_est[1]
evt_gen
smp_lab <- rownames(raw_smp)
smp_gen <- raw_smp$age_mean
smp_qua1 <- (0:100) / 100
smp_qua2 <- (1:50) / 10
thr_val <- (1:5) * 2
sim_num <- 500

sim_cal_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val), sim_num))
sim_err_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val), sim_num))
avg_cal_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val)))
avg_err_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val)))
for (i in 1:length(smp_qua1)) {
  print(i)
  for (j in 1:length(smp_qua2)) {
    print(j)
    for (k in 1:length(thr_val)) {
      print(k)
      for (r in 1:sim_num) {
        sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, c(smp_qua1[i], smp_qua2[j]), thr_val[k])
        sim_cal_rat[i, j, k, r] <- sim_HMM_WFM$cal_rat
        sim_err_rat[i, j, k, r] <- sim_HMM_WFM$err_rat
      }
      avg_cal_rat[i, j, k] <- mean(sim_cal_rat[i, j, k, ])
      avg_err_rat[i, j, k] <- mean(sim_err_rat[i, j, k, ])
    }
  }
}

save(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua1, smp_qua2, thr_val, sim_num, sim_cal_rat, sim_err_rat, avg_cal_rat, avg_err_rat,
     file = "./Output/Output v1.0/MANU_COL_MC1R_SIMU_CALL.rda")

############################################################

#' KIT13

#' Raw data of Wutke et al. (2016)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT13
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
# raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
# raw_smp[which(raw_smp[, 7] == 1), 4:5] <- 1 / 2
# raw_smp[which(raw_smp[, 8] == 1), 5:6] <- 1 / 2
# raw_smp[which(raw_smp[, 9] == 1), 4:6] <- 1 / 3
raw_smp[which(raw_smp[, 7] == 1), 4] <- 1 / 2
raw_smp[which(raw_smp[, 7] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 6] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 4] <- 1 / 4
raw_smp[which(raw_smp[, 9] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 6] <- 1 / 4
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp
n_x <- round(2 * max(pop_siz) * 0.02) + 1
n_t <- round(n_x * n_x / 8 / min(pop_siz)) + 1
ptn_num <- c(n_t, n_x - 24, 26)
unif_grd <- TRUE
coord_x <- calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)
mut_rat <- c(7.24e-09, 7.24e-09)
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptBackPMMH(mut_rat, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/MANU_PTN_KIT13_REAL1.rda")

########################################

#' Raw data of Wutke et al. (2016)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT13
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
# raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
# raw_smp[which(raw_smp[, 7] == 1), 4:5] <- 1 / 2
# raw_smp[which(raw_smp[, 8] == 1), 5:6] <- 1 / 2
# raw_smp[which(raw_smp[, 9] == 1), 4:6] <- 1 / 3
raw_smp[which(raw_smp[, 7] == 1), 4] <- 1 / 2
raw_smp[which(raw_smp[, 7] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 6] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 4] <- 1 / 4
raw_smp[which(raw_smp[, 9] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 6] <- 1 / 4
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 1e+00
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- max(raw_smp$age_mean)
raw_smp
n_x <- round(2 * max(pop_siz) * 0.02) + 1
n_t <- round(n_x * n_x / 8 / min(pop_siz)) + 1
ptn_num <- c(n_t, n_x - 24, 26)
unif_grd <- TRUE
coord_x <- calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)
mut_rat <- c(7.24e-09, 7.24e-09)
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptBackPMMH(mut_rat, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/MANU_PTN_KIT13_REAL2.rda")

########################################

load("./Output/Output v1.0/MANU_PTN_KIT13_REAL1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
sel_cof_est <- c(sel_cof_pdf$x[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[1]], sel_cof_pdf$y[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[2]])
frq_pth_est <- colMeans(frq_pth_chn)

#' Show the call rate and accuracy accoss different data qualities and call thresholds
model <- "WFM"
sel_cof <- sel_cof_est
dom_par
pop_siz
int_con <- frq_pth_est[1]
evt_gen
smp_lab <- rownames(raw_smp)
smp_gen <- raw_smp$age_mean
smp_qua1 <- (0:100) / 100
smp_qua2 <- (1:50) / 10
thr_val <- (1:5) * 2
sim_num <- 500

sim_cal_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val), sim_num))
sim_err_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val), sim_num))
avg_cal_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val)))
avg_err_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val)))
for (i in 1:length(smp_qua1)) {
  print(i)
  for (j in 1:length(smp_qua2)) {
    print(j)
    for (k in 1:length(thr_val)) {
      print(k)
      for (r in 1:sim_num) {
        sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, c(smp_qua1[i], smp_qua2[j]), thr_val[k])
        sim_cal_rat[i, j, k, r] <- sim_HMM_WFM$cal_rat
        sim_err_rat[i, j, k, r] <- sim_HMM_WFM$err_rat
      }
      avg_cal_rat[i, j, k] <- mean(sim_cal_rat[i, j, k, ])
      avg_err_rat[i, j, k] <- mean(sim_err_rat[i, j, k, ])
    }
  }
}

save(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua1, smp_qua2, thr_val, sim_num, sim_cal_rat, sim_err_rat, avg_cal_rat, avg_err_rat,
     file = "./Output/Output v1.0/MANU_PTN_KIT13_SIMU_CALL.rda")

############################################################

#' KIT16

#' Raw data of Wutke et al. (2016)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT16
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
# raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
# raw_smp[which(raw_smp[, 7] == 1), 4:5] <- 1 / 2
# raw_smp[which(raw_smp[, 8] == 1), 5:6] <- 1 / 2
# raw_smp[which(raw_smp[, 9] == 1), 4:6] <- 1 / 3
raw_smp[which(raw_smp[, 7] == 1), 4] <- 1 / 2
raw_smp[which(raw_smp[, 7] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 6] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 4] <- 1 / 4
raw_smp[which(raw_smp[, 9] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 6] <- 1 / 4
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((400 - 2000) / 8) # 400 AD (the Middle Ages)
raw_smp
n_x <- round(2 * max(pop_siz) * 0.02) + 1
n_t <- round(n_x * n_x / 8 / min(pop_siz)) + 1
ptn_num <- c(n_t, n_x - 24, 26)
unif_grd <- TRUE
coord_x <- calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)
mut_rat <- c(7.24e-09, 7.24e-09)
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptBackPMMH(mut_rat, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/MANU_PTN_KIT16_REAL1.rda")

########################################

#' Raw data of Wutke et al. (2016)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- KIT16
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
# raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
# raw_smp[which(raw_smp[, 7] == 1), 4:5] <- 1 / 2
# raw_smp[which(raw_smp[, 8] == 1), 5:6] <- 1 / 2
# raw_smp[which(raw_smp[, 9] == 1), 4:6] <- 1 / 3
raw_smp[which(raw_smp[, 7] == 1), 4] <- 1 / 2
raw_smp[which(raw_smp[, 7] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 6] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 4] <- 1 / 4
raw_smp[which(raw_smp[, 9] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 6] <- 1 / 4
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- max(raw_smp$age_mean)
raw_smp
n_x <- round(2 * max(pop_siz) * 0.02) + 1
n_t <- round(n_x * n_x / 8 / min(pop_siz)) + 1
ptn_num <- c(n_t, n_x - 24, 26)
unif_grd <- TRUE
coord_x <- calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)
mut_rat <- c(7.24e-09, 7.24e-09)
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptPMMH(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/MANU_PTN_KIT16_REAL2.rda")

########################################

load("./Output/Output v1.0/MANU_PTN_KIT16_REAL1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
sel_cof_est <- c(sel_cof_pdf$x[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[1]], sel_cof_pdf$y[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[2]])
frq_pth_est <- colMeans(frq_pth_chn)

#' Show the call rate and accuracy accoss different data qualities and call thresholds
model <- "WFM"
sel_cof <- sel_cof_est
dom_par
pop_siz
int_con <- frq_pth_est[1]
evt_gen
smp_lab <- rownames(raw_smp)
smp_gen <- raw_smp$age_mean
smp_qua1 <- (0:100) / 100
smp_qua2 <- (1:50) / 10
thr_val <- (1:5) * 2
sim_num <- 500

sim_cal_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val), sim_num))
sim_err_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val), sim_num))
avg_cal_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val)))
avg_err_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val)))
for (i in 1:length(smp_qua1)) {
  print(i)
  for (j in 1:length(smp_qua2)) {
    print(j)
    for (k in 1:length(thr_val)) {
      print(k)
      for (r in 1:sim_num) {
        sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, c(smp_qua1[i], smp_qua2[j]), thr_val[k])
        sim_cal_rat[i, j, k, r] <- sim_HMM_WFM$cal_rat
        sim_err_rat[i, j, k, r] <- sim_HMM_WFM$err_rat
      }
      avg_cal_rat[i, j, k] <- mean(sim_cal_rat[i, j, k, ])
      avg_err_rat[i, j, k] <- mean(sim_err_rat[i, j, k, ])
    }
  }
}

save(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua1, smp_qua2, thr_val, sim_num, sim_cal_rat, sim_err_rat, avg_cal_rat, avg_err_rat,
     file = "./Output/Output v1.0/MANU_PTN_KIT16_SIMU_CALL.rda")

############################################################

#' TRPM1

#' Raw data of Wutke et al. (2016)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- TRPM1
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
# raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
# raw_smp[which(raw_smp[, 7] == 1), 4:5] <- 1 / 2
# raw_smp[which(raw_smp[, 8] == 1), 5:6] <- 1 / 2
# raw_smp[which(raw_smp[, 9] == 1), 4:6] <- 1 / 3
raw_smp[which(raw_smp[, 7] == 1), 4] <- 1 / 2
raw_smp[which(raw_smp[, 7] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 6] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 4] <- 1 / 4
raw_smp[which(raw_smp[, 9] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 6] <- 1 / 4
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- round((-3500 - 2000) / 8) # 3500 BC (domestication)
raw_smp
n_x <- round(2 * max(pop_siz) * 0.02) + 1
n_t <- round(n_x * n_x / 8 / min(pop_siz)) + 1
ptn_num <- c(n_t, n_x - 24, 26)
unif_grd <- TRUE
coord_x <- calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)
mut_rat <- c(7.24e-09, 7.24e-09)
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptBackPMMH(mut_rat, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/MANU_PTN_TRPM1_REAL1.rda")

########################################

#' Raw data of Wutke et al. (2016)
load("./Data/REAL.rda")

set.seed(1)
raw_smp <- TRPM1
raw_smp <- raw_smp[which(rowSums(raw_smp[, 4:9]) != 0), ]
int_gen <- -round(max(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
lst_gen <- -round(min(raw_smp$age_mean, raw_smp$age_lower, raw_smp$age_upper) / 8)
# raw_smp <- raw_smp[which(raw_smp$age_mean <= max(raw_smp$age_mean[which(rowSums(raw_smp[, c(5, 6, 8)]) != 0)])), ]
max(raw_smp$age_mean) - 2000
raw_smp$age_mean <- -round(raw_smp$age_mean / 8)
raw_smp$age_lower <- -round(raw_smp$age_lower / 8)
raw_smp$age_upper <- -round(raw_smp$age_upper / 8)
# raw_smp[which(raw_smp[, 7] == 1), 4:5] <- 1 / 2
# raw_smp[which(raw_smp[, 8] == 1), 5:6] <- 1 / 2
# raw_smp[which(raw_smp[, 9] == 1), 4:6] <- 1 / 3
raw_smp[which(raw_smp[, 7] == 1), 4] <- 1 / 2
raw_smp[which(raw_smp[, 7] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 8] == 1), 6] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 4] <- 1 / 4
raw_smp[which(raw_smp[, 9] == 1), 5] <- 1 / 2
raw_smp[which(raw_smp[, 9] == 1), 6] <- 1 / 4
raw_smp <- raw_smp[order(raw_smp$age_mean), ]
raw_smp <- raw_smp[, -c(2, 3, 7, 8, 9)]

sel_cof <- c(0e+00, 0e+00)
dom_par <- 5e-01
pop_siz <- pop_siz[min(raw_smp$age_mean - int_gen + 1):max(raw_smp$age_mean - int_gen + 1)]
ref_siz <- tail(pop_siz, n = 1)
evt_gen <- max(raw_smp$age_mean)
raw_smp
n_x <- round(2 * max(pop_siz) * 0.02) + 1
n_t <- round(n_x * n_x / 8 / min(pop_siz)) + 1
ptn_num <- c(n_t, n_x - 24, 26)
unif_grd <- TRUE
coord_x <- calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)
mut_rat <- c(7.24e-09, 7.24e-09)
pcl_num <- 1e+03
itn_num <- 2e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(PMMH <- cmprunAdaptBackPMMH(mut_rat, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto))

save(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, PMMH,
     file = "./Output/Output v1.0/MANU_PTN_TRPM1_REAL2.rda")

########################################

load("./Output/Output v1.0/MANU_PTN_TRPM1_REAL1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
frq_pth_chn <- PMMH$frq_pth_chn
brn_num <- 0.5 * dim(sel_cof_chn)[2] + 1
sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]
frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]
grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
sel_cof_est <- c(sel_cof_pdf$x[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[1]], sel_cof_pdf$y[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[2]])
frq_pth_est <- colMeans(frq_pth_chn)

#' Show the call rate and accuracy accoss different data qualities and call thresholds
model <- "WFM"
sel_cof <- sel_cof_est
dom_par
pop_siz
int_con <- frq_pth_est[1]
evt_gen
smp_lab <- rownames(raw_smp)
smp_gen <- raw_smp$age_mean
smp_qua1 <- (0:100) / 100
smp_qua2 <- (1:50) / 10
thr_val <- (1:5) * 2
sim_num <- 500

sim_cal_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val), sim_num))
sim_err_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val), sim_num))
avg_cal_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val)))
avg_err_rat <- array(NA, dim = c(length(smp_qua1), length(smp_qua2), length(thr_val)))
for (i in 1:length(smp_qua1)) {
  print(i)
  for (j in 1:length(smp_qua2)) {
    print(j)
    for (k in 1:length(thr_val)) {
      print(k)
      for (r in 1:sim_num) {
        sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, c(smp_qua1[i], smp_qua2[j]), thr_val[k])
        sim_cal_rat[i, j, k, r] <- sim_HMM_WFM$cal_rat
        sim_err_rat[i, j, k, r] <- sim_HMM_WFM$err_rat
      }
      avg_cal_rat[i, j, k] <- mean(sim_cal_rat[i, j, k, ])
      avg_err_rat[i, j, k] <- mean(sim_err_rat[i, j, k, ])
    }
  }
}

save(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua1, smp_qua2, thr_val, sim_num, sim_cal_rat, sim_err_rat, avg_cal_rat, avg_err_rat,
     file = "./Output/Output v1.0/MANU_PTN_TRPM1_SIMU_CALL.rda")

################################################################################
