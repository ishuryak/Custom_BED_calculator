source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
X <- read.csv("C:/Users/ishur/AppData/Local/Temp/dd_table.csv"); ref <- default_ref(); K <- bed_lq(X$n,X$d,10)
ldi <- log(2)/3
for (Tk in c(28,21)) for (l in seq(0.36,0.50,by=0.01)) { a <- 0.35
  cal <- dd_calibration(ref,a,10,ldi,l,Tk); if (cal$status!="ok") next; on <- cal$c_over_alpha*X$T/K; rate <- l*(1-exp(-a*K/X$T)); tp <- log(2)/rate; loss <- rate/a; rep <- on<X$T
  bad <- rep & !(tp>=2 & tp<=6 & loss>=0.5 & loss<=1)
  p <- default_params(); p$alpha<-a; p$lambda_di<-ldi; p$Tk<-Tk; p$lambda_dd<-l; b <- suggest_regimens(p, ref, min_gain=1, min_days_saved=7, early_cap=59)
  cat(sprintf("Tk %d ldd %.2f: onset_ref %.1f Tp %.1f-%.1f loss %.2f-%.2f fails:%s onsets<18:%s | Fowler59 %d\n", Tk, l, cal$onset_ref, min(tp[rep]),max(tp[rep]),min(loss[rep]),max(loss[rep]), paste(X$id[bad],collapse=""), paste(X$id[on<18&rep],collapse=""), nrow(b))) }
