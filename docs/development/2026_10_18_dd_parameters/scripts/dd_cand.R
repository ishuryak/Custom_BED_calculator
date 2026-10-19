source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
OUT <- "C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/dd_grid_2026_10_18/"
X <- read.csv("C:/Users/ishur/AppData/Local/Temp/dd_table.csv")[,c("id","name","d","n","T","ctrl_d","ctrl_n","ctrl_T","tcp","ctrl_tcp")]
ref <- default_ref(); K <- bed_lq(X$n,X$d,10); all <- NULL
for (pt in list(c(0.30,0.375),c(0.40,0.394))) { p <- default_params(); p$alpha <- pt[1]; p$lambda_dd <- pt[2]
  cal <- dd_calibration(ref,pt[1],10,0.2,pt[2],28); on <- cal$c_over_alpha*X$T/K; rate <- pt[2]*(1-exp(-pt[1]*K/X$T))
  o <- compute_all(X$n,X$d,X$T,p,ref); c0 <- compute_all(X$ctrl_n,X$ctrl_d,X$ctrl_T,p,ref)
  D <- data.frame(alpha=pt[1], ldd=pt[2], onset_ref=round(cal$onset_ref,1), ca=round(cal$c_over_alpha,1), id=X$id, T=X$T, onset=round(on,1), tp=round(log(2)/rate,1), loss=round(rate/pt[1],2),
    tp_ok=log(2)/rate>=2 & log(2)/rate<=6, loss_ok=rate/pt[1]>=0.5 & rate/pt[1]<=1, on_ok=on>=18 | on>=X$T,
    bed_di=round(o$bed_di,1), bed_dd=round(o$bed_dd,1), dDI=round(o$bed_di-c0$bed_di,1), dDD=round(o$bed_dd-c0$bed_dd,1), dTCP=X$tcp-X$ctrl_tcp)
  print(D, row.names=FALSE); all <- rbind(all,D) }
write.csv(all, paste0(OUT,"candidate_points_regimens.csv"), row.names=FALSE)
# default for comparison of BED_DI with alpha change
p <- default_params(); cat("default ref DI", compute_all(35,2,47,p,ref)$bed_di, "\n")
