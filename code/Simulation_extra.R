# Model-assisted calibration using
# Generalized entropy calibration in survey sampling

# Simulation setup

if (!interactive()) {
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
} else{
  args <- c(5)
}

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

# install.packages( "MatrixModels", type="win.binary" )

library(nleqslv)
library(PracTools)
suppressMessages(library(foreach))
suppressMessages(library(doRNG))
# library(caret)
library(tidyverse)
library(xtable)
library(GECal)
# library(mice)
# library(kableExtra)

load("../data/nhis.Rdata")

is.nmin = T # True for different sample size; False for different sampling fraction


set.seed(11)
SIMNUM = args[1]

if (!interactive()) {
  suppressMessages(library(doMC))
  dir.create(timenow0)
  setwd(timenow0)
  
  sink(timenow, append = TRUE)
  
  ## doMC backend registration
  nc <- parallel::detectCores(logical = TRUE)
  reserve <- 2L
  target <- max(1L, nc - reserve)
  want_workers <- 500L
  workers <- max(1L, min(want_workers, SIMNUM))
  
  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
  registerDoMC(workers)
  message(sprintf("doMC backend: detectCores=%d, using workers=%d", nc, getDoParWorkers()))
  
  # (Optional) Somke test
  .ok <- tryCatch({
    aa <- foreach(i = 1:getDoParWorkers(), .combine=c) %dopar% i
    length(aa) == getDoParWorkers()
  }, error=function(e){ message("Smoke test failed: ", conditionMessage(e)); FALSE })
  if (!.ok) stop("doMC Smoke test failed")
}else{
  suppressMessages(library(doParallel))
  
  # ## setup parallel backend to use many processors
  print(paste("detectCores =", detectCores()))
  cores = min(detectCores() - 3, 120)
  # cores = 2
  print(paste("cores =", cores))
  # cl <- makeCluster(cores, outfile = timenow) #not to overload your computer
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}
print(paste("SIMNUM =", SIMNUM))

print("y is Hemo")
# nhis$Hemo <- nhis$OralExam; print("y is OralExam")  # Variable of interest y is Smoking

# run_sim <- function(denom, getn = FALSE) { # is.nmin = FALSE
run_sim <- function(nh_min, getn = FALSE) { # is.nmin = TRUE
  if(is.nmin){
    denom = 1 
    if(!getn) print(paste("nh_min =", nh_min));    
  }else{
    nh_min = 5
    if(!getn) print(paste("denom =", denom)); 
  }
  # 
  nhis_sub <- nhis[sample(1:nrow(nhis), round(nrow(nhis) / denom)), ]
  N <- nrow(nhis_sub)
  tab1 <- table(nhis_sub$AgeGroup, nhis_sub$REGION1, nhis_sub$SEX)
  # tab1 <- ifelse(tab1 > 3 * nh_min, nh_min, round(tab1 / 3))
  tab1 = ifelse(tab1 > 9, nh_min, round(tab1 / 4))
  n <- sum(tab1)
  if(getn) return(n)
  theta <- mean(nhis_sub$Hemo)
  
  final_res <- foreach(simnum = 1:SIMNUM,
                       .packages = c("nleqslv","PracTools","GECal"),
                       .errorhandling = "pass") %dorng% {
                         theta = mean(nhis_sub$Hemo)
                         index = c()
                         pi_S = c()
                         pi = rep(0, N)
                         for(AgeGroup in unique(nhis_sub$AgeGroup)){
                           for(REGION1 in unique(nhis_sub$REGION1)){
                             for(SEX in unique(nhis_sub$SEX)){
                               idx_tmp = which(nhis_sub$AgeGroup == AgeGroup &
                                                 nhis_sub$REGION1 == REGION1 & 
                                                 nhis_sub$SEX == SEX)
                               n_h = tab1[AgeGroup,REGION1,SEX]
                               index = c(index, sample(idx_tmp, n_h, replace = FALSE))
                               pi_S = c(pi_S, rep(n_h / length(idx_tmp), n_h))
                               pi[idx_tmp] = n_h / length(idx_tmp)
                             }
                           }
                         }
                         # summary(pi)
                         
                         delta = as.integer(1:N %in% index)
                         Rmodel = glm(delta ~ AgeGroup + REGION1 + SEX, family = binomial,
                                      data = nhis_sub)
                         
                         # Rmodel = glm(delta ~ AgeGroup * REGION1 * SEX, family = binomial,
                         #              data = nhis_sub)
                         
                         # Rmodel = glm(delta ~ AgeGroup + REGION1 + SEX + HT + WT + Waist + Alcohol + SysBP + DiaBP+ FBS + Creatinine, family = binomial,
                         #              data = nhis_sub)
                         
                         # Rmodel = glm(delta ~ AgeGroup + SEX, family = binomial,
                         #              data = nhis_sub)
                         pihat0 = predict.glm(Rmodel, nhis_sub, type = "response")
                         
                         
                         findphi2 = function(phi, x0, Z, delta, ..., returnw = F){
                           pi_phi = drop(1 / (1 + exp(-x0 %*% phi)))
                           w_phi = ifelse(delta, (1 - pi_phi) / pi_phi, -1)
                           if(returnw){
                             return(drop(1 / pi_phi))
                           }else{
                             return(drop(w_phi %*% Z))
                           }
                         }
                         jacphi2 = function(phi, x0, Z, delta, ..., returnw = F){
                           pi_phi = drop(1 / (1 + exp(-x0 %*% phi)))
                           logit_phi = pi_phi / (1 - pi_phi)
                           return(-t(Z) %*% (x0 * ifelse(delta, 1 / logit_phi, 0) ))
                           # return(-t(x0) %*% (Z * ifelse(delta, 1 / logit_phi, logit_phi) ))
                         }
                         
                         xR = model.matrix(Rmodel)
                         nleqslv_res = nleqslv(Rmodel$coefficients, findphi2, jac = jacphi2, x0 = xR, Z = xR, delta = delta,
                                               control = list(maxit = 1e5, allowSingular = T), xscalm = "auto",
                                               method = "Newton")
                         
                         if(nleqslv_res$termcd != 1 & max(abs(findphi2(nleqslv_res$x, x0 = xR, Z = xR, delta = delta))) > 1e-5){
                           w_S = NA
                         }else{
                           w_S = findphi2(nleqslv_res$x, x0 = xR, Z = xR, delta = delta, returnw = T)
                         }
                         # drop(t(xR[Index_S,]) %*% w_S2[Index_S]); colSums(xR)
                         w_S2 = w_S
                         
                         # summary(lm(Smoking ~ ., data = nhis_sub))
                         
                         # Omodel = lm(Hemo ~ AgeGroup + SEX + HT + WT, data = nhis_sub)
                         
                         Omodel = lm(Hemo ~ AgeGroup + SEX + HT + WT + Waist + Alcohol + SysBP + DiaBP + FBS + Creatinine, data = nhis_sub)
                         yhat = predict.lm(Omodel, nhis_sub, type = "response")
                         
                         nhis_sub.samp <- nhis_sub[index, ]
                         # d_S       <- 1 / pi_S
                         y_S <- nhis_sub.samp$Hemo
                         
                         # omega2_S = matrix(0, nrow = length(index), ncol = length(index))
                         # for (cnt in 1:5) {
                         #   idx_strat2 = which(smho98.samp$stratum5 == cnt)
                         #   # omega2_S[idx_strat2, idx_strat2] <-
                         #   diag(omega2_S)[idx_strat2] <-
                         #     1 / pi_S[idx_strat2] * (1 / pi_S[idx_strat2] - 1)
                         # }
                         
                         theta_res = NULL
                         se_res = NULL
                         
                         # for(pimethod in c(0,1,3)){
                           for(pimethod in c(1,3)){                           
                           if(pimethod == 0){
                             pihat = pi
                           }else if(pimethod == 1){
                             pihat = rep(n / N, N)
                           }else if(pimethod == 2){
                             pihat = pihat0
                           }else if(pimethod == 3){
                             pihat = 1 / w_S2
                           }
                           d_S <- 1 / pihat[index]
                           
                           # #HT estimator
                           # const = numeric(0)
                           # calibration <- GEcalib(
                           #   ~ 0,
                           #   dweight = d_S,
                           #   data = nhis_sub.samp,
                           #   const = const,
                           #   entropy = 1,
                           #   method = "DS"
                           # )
                           # res_est = estimate(y_S ~ 1, calibration = calibration)$estimate
                           # 
                           # theta_res = c(theta_res, setNames(res_est[1] / N, "IPW"))
                           # se_res = c(se_res, setNames(res_est[2] / N, "IPW"))
                           
                           # Hajek estimator
                           # const = N
                           # calibration <- GEcalib(
                           #   ~ 1,
                           #   dweight = d_S,
                           #   data = nhis_sub.samp,
                           #   const = const,
                           #   entropy = 1,
                           #   method = "DS"
                           # )
                           # res_est = estimate(y_S ~ 1, calibration = calibration)$estimate
                           # 
                           # theta_res = c(theta_res, setNames(res_est[1] / N, "SIPW"))
                           # se_res = c(se_res, setNames(res_est[2] / N, "SIPW"))
                           
                           # SIPW estimator
                           theta_res = c(theta_res, setNames(sum((y_S) / pihat[index]) / sum(1 / pihat[index]), paste("IPW", ifelse(pimethod == 1, 2, pimethod), sep = ""))) # SIPW
                                                      
                           hhat = pihat * model.matrix(~1, data = nhis_sub)
                           kappa = solve(t(hhat[index,]) %*% (hhat[index,] * (1 / pihat[index] - 1) / pihat[index]),
                                         t(hhat[index,]) %*% ((y_S - yhat[index]) * (1 / pihat[index] - 1) / pihat[index]))
                           eta = yhat + drop(hhat %*% kappa)
                           eta[index] = eta[index] + (y_S - yhat[index] - drop(hhat %*% kappa)[index]) / pihat[index]

                           se_res = c(se_res, IPW = sqrt(var(eta) / N))
                           
                           # AIPW estimator
                           # theta_res = c(theta_res, AIPW = (sum(yhat) + sum((y_S - yhat[index]) / pihat[index])) / N) # AIPW
                           # 
                           # if(pimethod == 1){
                           #   hhat = pihat * model.matrix(Rmodel)
                           # }else if(pimethod == 2){
                           #   hhat = pihat / (1 - pihat) * model.matrix(Rmodel)
                           # }
                           # # hhat = pihat * model.matrix(Rmodel)
                           # kappa = solve(t(hhat[index,]) %*% (hhat[index,] * (1 / pihat[index] - 1) / pihat[index]),
                           #               t(hhat[index,]) %*% ((y_S - yhat[index]) * (1 / pihat[index] - 1) / pihat[index]))
                           # eta = yhat + drop(hhat %*% kappa)
                           # eta[index] = eta[index] + (y_S - yhat[index] - drop(hhat %*% kappa)[index]) / pihat[index]
                           # 
                           # se_res = c(se_res, AIPW = sqrt(var(eta) / N))
                           
                           vectmp <- as.character(formula(Omodel))[3]
                           # vectmp <- c("Y_IP")
                           
                           fortmp <- formula(paste("~", paste(vectmp, collapse = "+")))
                           fortmp2 <- formula(paste("~", paste(c(vectmp, "g(d_S)"), collapse = "+")))
                           
                           # Omodel_d = lm(reformulate(vectmp, response = "y_S"), weights = 1 / pi_S, data = smho98.samp)
                           # # Omodel_d = lm(reformulate(paste0("x.", 1:pcol), response = "y"), data = data_S)      
                           # 
                           # yhat_d = predict.lm(Omodel_d, smho98.samp, type = "response")
                           # e2 = (yhat_d - y_S)^2
                           # Vmodel_d = glm(reformulate(vectmp, response = "e2"),
                           #                weights = 1 / pi_S, data = smho98.samp)
                           # vhat = predict.glm(Vmodel_d, smho98, type = "response")
                           
                           for (entropy in list(-1, -1/2, 1)) {
                             if(entropy == "CE") pihat = ifelse(pihat > 0.6, 0.6, pihat) # To make CE convergent
                             d_S <- 1 / pihat[index]
                             # const = colSums(model.matrix(fortmp, nhis_sub))
                             # 
                             # calibration <- GEcalib(
                             #   fortmp,
                             #   dweight = d_S,
                             #   data = nhis_sub.samp,
                             #   const = const,
                             #   entropy = entropy,
                             #   method = "DS"
                             # )
                             # 
                             # res_est = estimate(y_S ~ 1, calibration = calibration)$estimate
                             # 
                             # theta_res = c(theta_res, setNames(res_est[1] / N, paste("DS", entropy, sep = "_"))) # DS
                             # se_res = c(se_res, setNames(res_est[2] / N, paste("DS", entropy, sep = "_")))
                             
                             if(pimethod == 1){
                               const = colSums(model.matrix(fortmp, nhis_sub))
                               calibration <- GEcalib(
                                 fortmp,
                                 dweight = d_S,
                                 data = nhis_sub.samp,
                                 const = const,
                                 entropy = entropy,
                                 method = "GEC0"
                               )
                             }else{
                               const = colSums(cbind(model.matrix(fortmp, nhis_sub), 
                                                     g(1 / pihat, entropy = entropy)))
                               calibration <- GEcalib(
                                 fortmp2,
                                 dweight = d_S,
                                 data = nhis_sub.samp,
                                 const = const,
                                 entropy = entropy,
                                 method = "GEC"
                               )
                             }
                             
                             calibration
                             
                             res_est = estimate(y_S ~ 1, calibration = calibration)$estimate
                             
                             if(entropy == -1){
                               entropy = "EL"
                             }else if(entropy == 0){
                               entropy = "ET"
                             }else if(entropy == -0.5){
                               entropy = "HD"
                             }else if(entropy == 1){
                               entropy = "SL"
                             }

                             
                             theta_res = c(theta_res, setNames(res_est[1] / N, paste(entropy, ifelse(pimethod == 1, 2, pimethod), sep = ""))) # DS
                             se_res = c(se_res, setNames(res_est[2] / N, paste(entropy, ifelse(pimethod == 1, 2, pimethod), sep = "")))
                           }
                         }
                         
                         CR_res = ifelse(abs(theta_res - theta) < qnorm(0.975) * se_res, 1, 0)
                         list(theta_res, se_res, CR_res)
                       }
  
  print(paste("# of failure:", sum(!sapply(lapply(final_res, function(x) x[[1]]), function(x) is.numeric(unlist(x)))))); 
  final_res0 = lapply(final_res, function(x) x[[1]])
  print(final_res0[!sapply(final_res0, function(x) is.numeric(unlist(x)))])
  
  final_res = final_res[sapply(final_res, length) == max(sapply(final_res, length))]
  
  final_res1 <- lapply(final_res, function(x) x[[1]])
  tmpres <- do.call(rbind, lapply(final_res1, c))
  # stop(tmpres)
  RMSE <- sqrt(colMeans((tmpres - theta)^2, na.rm = TRUE))
  return(RMSE)
}

if (interactive()) stopCluster(cl)

if(is.nmin){
  nhmin_vec <- c(1, 2, 3, 4, 5, 6, 7)
  # nhmin_vec <- c(5, 10)
  rmse_list <- lapply(nhmin_vec, run_sim)
}else{
  denom_vec <- c(1, 2, 4, 6.4, 10, 20, 40)
  denom_vec <- c(20, 40)
  rmse_list <- lapply(denom_vec, run_sim)
}

timenow2 = Sys.time()
print(timenow2 - timenow1)
if(!interactive()) save.image(paste(timenow0, ".RData", sep = ""))

rmse_df <- do.call(rbind, rmse_list)
rmse_df <- as.data.frame(rmse_df)
if(is.nmin){
  rmse_df$f <- sapply(nhmin_vec, run_sim, getn = TRUE)
}else{
  rmse_df$f <- sapply(denom_vec, run_sim, getn = TRUE)
}
rmse_df = rmse_df[,-8]


library(ggplot2)

rmse_df_long <- rmse_df %>%
  pivot_longer(-f, names_to = "Estimator", values_to = "RMSE")


p = ggplot(rmse_df_long, aes(x = f, y = RMSE, color = Estimator)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = round(rmse_df$f, ifelse(is.nmin, -2 , 2))) +
  labs(x = ifelse(is.nmin, "Sample size", "Sampling Fraction"),
       y = "RMSE") +  
  theme_minimal(base_size = 14)
p

if(!interactive()) ggsave(
  filename = "SamplingFraction.png",  # file name
  plot = p,                                   # which plot
  width = 8, height = 5, units = "in",        # size in inches
  dpi = 300                                   # high resolution
)
p
# ggplot(rmse_df_long, aes(x = f, y = RMSE, color = Estimator)) +
#   geom_line() + geom_point(size = 2) +
#   scale_x_continuous(
#     trans = "log2",
#     breaks = rmse_df$f,
#     labels = paste0(rmse_df$f)
#   ) +
#   labs(
#     x = "Sampling Fraction (n / N)",
#     y = "RMSE"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom")
