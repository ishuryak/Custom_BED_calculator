source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
OUT <- "C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/dd_grid_2026_10_18/"
R <- read.csv(paste0(OUT,"grid_rate_first.csv")); F <- R[R$band=="loose" & R$n_rate==13,]
F$ab <- round(F$alpha,2)
agg <- do.call(rbind, lapply(split(F, F$alpha), function(s) data.frame(alpha=s$alpha[1], di_loss=round(s$di_loss[1],2), ldd_min=round(min(s$lambda_dd),3), ldd_max=round(max(s$lambda_dd),3), onset_ref=paste(round(range(s$onset_ref),1),collapse="-"), n_onset_max=max(s$n_onset))))
print(agg[seq(1,nrow(agg),by=2),], row.names=FALSE)
write.csv(agg, paste0(OUT,"feasible_region_by_alpha.csv"), row.names=FALSE)
X <- read.csv("C:/Users/ishur/AppData/Local/Temp/dd_table.csv")[,c("id","name","d","n","T","ctrl_d","ctrl_n","ctrl_T","tcp","ctrl_tcp")]
ref <- default_ref(); K <- bed_lq(X$n,X$d,10)
pick <- function(a) { s <- F[abs(F$alpha-a)<1e-9,]; s$lambda_dd[which.min(abs(s$lambda_dd-median(s$lambda_dd)))] }
pts <- list(c(0.28, pick(0.28)), c(0.30, pick(0.30)), c(0.40, F$lambda_dd[F$alpha==0.4][which.max(F$n_onset[F$alpha==0.4]-F$onset_short[F$alpha==0.4]/100)]))
all <- NULL
for (pt in pts) { p <- default_params(); p$alpha <- pt[1]; p$lambda_dd <- pt[2]
  cal <- dd_calibration(ref,pt[1],10,0.2,pt[2],28); on <- cal$c_over_alpha*X$T/K; rate <- pt[2]*(1-exp(-pt[1]*K/X$T))
  o <- compute_all(X$n,X$d,X$T,p,ref); c0 <- compute_all(X$ctrl_n,X$ctrl_d,X$ctrl_T,p,ref)
  D <- data.frame(alpha=pt[1], lambda_dd=round(pt[2],3), onset_ref=round(cal$onset_ref,1), c_over_alpha=round(cal$c_over_alpha,2), id=X$id, name=X$name, d=X$d, n=X$n, T=X$T,
    onset=round(on,1), rep=on<X$T, tp=round(log(2)/rate,2), loss=round(rate/pt[1],2), di_loss=round(0.2/pt[1],2), bed_di=round(o$bed_di,1), bed_dd=round(o$bed_dd,1),
    dDI=round(o$bed_di-c0$bed_di,1), dDD=round(o$bed_dd-c0$bed_dd,1), dTCP=round(X$tcp-X$ctrl_tcp,1))
  cat(sprintf("\n== alpha %.2f lambda_dd %.3f onset_ref %.1f C/a %.1f DI loss %.2f ==\n", pt[1], pt[2], cal$onset_ref, cal$c_over_alpha, 0.2/pt[1])); print(D[,-(1:4)], row.names=FALSE); all <- rbind(all,D) }
write.csv(all, paste0(OUT,"candidate_points_regimens.csv"), row.names=FALSE)
