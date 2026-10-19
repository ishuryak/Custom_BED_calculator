source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
X <- read.csv("C:/Users/ishur/AppData/Local/Temp/dd_table.csv"); ref <- default_ref(); K <- bed_lq(X$n,X$d,10)
for (pt in list(c(0.35,0.375),c(0.35,0.38),c(0.30,0.38))) {
  cal <- dd_calibration(ref,pt[1],10,0.2,pt[2],28); on <- cal$c_over_alpha*X$T/K; rate <- pt[2]*(1-exp(-pt[1]*K/X$T)); tp <- log(2)/rate; loss <- rate/pt[1]; rep <- on<X$T
  bad <- rep & !(tp>=2 & tp<=6 & loss>=0.5 & loss<=1)
  p <- default_params(); p$alpha<-pt[1]; p$lambda_dd<-pt[2]
  a <- suggest_regimens(p, ref, min_gain=1, min_days_saved=7); b <- suggest_regimens(p, ref, min_gain=1, min_days_saved=7, early_cap=59)
  cat(sprintf("alpha %.2f ldd %.3f: onset_ref %.1f; Tp %.1f-%.1f; loss %.2f-%.2f; rate/loss implausible: %s; onsets<18: %s; DI loss %.2f; search strict %d, Fowler59 %d (DD-credited %d)\n",
    pt[1],pt[2],cal$onset_ref,min(tp[rep]),max(tp[rep]),min(loss[rep]),max(loss[rep]), paste(X$id[bad],collapse=","), paste(X$id[on<18 & rep],collapse=","), 0.2/pt[1], nrow(a), nrow(b), sum(b$d_dd>=1,na.rm=T)))
}
