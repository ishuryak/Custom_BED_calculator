source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
X <- read.csv("C:/Users/ishur/AppData/Local/Temp/dd_table.csv"); ref <- default_ref(); K <- bed_lq(X$n,X$d,10)
for (ldi in c(0.23,0.2)) { a <- 0.35; l <- 0.40
  cal <- dd_calibration(ref,a,10,ldi,l,28); on <- cal$c_over_alpha*X$T/K; rate <- l*(1-exp(-a*K/X$T)); tp <- log(2)/rate; loss <- rate/a; rep <- on<X$T
  bad <- rep & !(tp>=2 & tp<=6 & loss>=0.5 & loss<=1)
  p <- default_params(); p$alpha<-a; p$lambda_di<-ldi; p$lambda_dd<-l
  b <- suggest_regimens(p, ref, min_gain=1, min_days_saved=7, early_cap=59); s <- suggest_regimens(p, ref, min_gain=1, min_days_saved=7)
  cat(sprintf("ldi %.2f (DI Tp %.2f d, DI loss %.2f): onset_ref %.1f, C/a %.1f; DD Tp %.1f-%.1f, loss %.2f-%.2f, fails:%s; onsets<18:%s (range %.1f-%.1f); Fowler59 %d rows, strict %d\n",
    ldi, log(2)/ldi, ldi/a, cal$onset_ref, cal$c_over_alpha, min(tp[rep]),max(tp[rep]),min(loss[rep]),max(loss[rep]), paste(X$id[bad],collapse=""), paste(X$id[on<18&rep],collapse=""), min(on), max(on), nrow(b), nrow(s)))
  r <- compute_all(c(35,20),c(2,2.75),c(47,26),p,ref); cat("  ref BED DI/DD", round(r$bed_di[1],1), round(r$bed_dd[1],1), "| HYPNO test", round(r$bed_di[2],1), round(r$bed_dd[2],1), "\n") }
