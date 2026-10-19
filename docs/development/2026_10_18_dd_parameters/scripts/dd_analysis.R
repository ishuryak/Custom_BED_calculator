source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
p <- default_params(); ref <- default_ref()
K0 <- bed_lq(35,2,10); pk0 <- 1-exp(-p$alpha*K0/ref$T)
cat(sprintf("Reference: K0=%.1f, daily kill fraction pk0=%.3f\n", K0, pk0))
cal <- dd_calibration(ref,p$alpha,p$ab,p$lambda_di,p$lambda_dd,p$Tk)
cat(sprintf("Default lambda_dd=0.5: onset_ref=%.2f d, C/alpha=%.2f Gy, realized rate at ref=%.3f/d (Tp=%.2f d), loss %.2f Gy/d\n",
  cal$onset_ref, cal$c_over_alpha, 0.5*pk0, log(2)/(0.5*pk0), 0.5*pk0/p$alpha))
cat(sprintf("DI: rate 0.2/d Tp=%.2f d, loss %.2f Gy/d\n", log(2)/0.2, 0.2/0.2))
lmin <- p$lambda_di*smooth_max0(ref$T-p$Tk)/(pk0*ref$T)
cat(sprintf("Smallest feasible lambda_dd (onset at day ~0): %.3f\n", lmin))
cat(sprintf("Matched lambda_dd = lambda_di/pk0 = %.3f\n", p$lambda_di/pk0))
for (ld in c(0.3,0.4,0.5,p$lambda_di/pk0,1,2,5)) {
  c2 <- dd_calibration(ref,p$alpha,p$ab,p$lambda_di,ld,p$Tk)
  q <- p; q$lambda_dd <- ld
  s <- suggest_regimens(q, ref)
  s2 <- suggest_regimens(q, ref, criterion="dd")
  cat(sprintf("lambda_dd=%.3f: status=%s onset_ref=%.1f C/a=%.1f | either: %d rows (%d with dDD>=0) | dd-only: %d rows\n",
    ld, c2$status, c2$onset_ref, c2$c_over_alpha, nrow(s), sum(s$d_dd>=0,na.rm=TRUE), nrow(s2)))
}
cat("\nDefault suggestions with onsets under default and matched lambda_dd:\n")
s <- suggest_regimens(p, ref)
cm <- dd_calibration(ref,p$alpha,p$ab,p$lambda_di,p$lambda_di/pk0,p$Tk)
s$K <- bed_lq(s$m,s$d,10); s$pk <- 1-exp(-p$alpha*s$K/s$T)
s$onset_def <- cal$c_over_alpha*s$T/s$K; s$onset_match <- cm$c_over_alpha*s$T/s$K
s$rate_def <- 0.5*s$pk; s$rate_match <- (p$lambda_di/pk0)*s$pk
s$bed_dd_match <- bed_dd(s$m,s$d,s$T,p$alpha,10,p$lambda_di/pk0,cm$c_over_alpha)
print(round(s[,c("m","d","T","bed_di","bed_dd","bed_early","bed_late","onset_def","onset_match","rate_def","rate_match","bed_dd_match")],2))
r0dd <- bed_dd(35,2,47,p$alpha,10,p$lambda_di/pk0,cm$c_over_alpha); cat(sprintf("ref DD matched = %.2f (DI ref %.2f)\n", r0dd, bed_di(35,2,47,.2,10,.2,28)))
cat("\nBenchmarks / extremes: onset under default (C/a=", round(cal$c_over_alpha,1), ") and matched (C/a=", round(cm$c_over_alpha,1), ")\n")
X <- data.frame(m=c(20,15,10,5,5,33,30,25), d=c(2.75,3.1,4,6,8,2,2.2,2.4))
X$T <- overall_time(X$m); X$K <- bed_lq(X$m,X$d,10)
X$onset_def <- cal$c_over_alpha*X$T/X$K; X$onset_match <- cm$c_over_alpha*X$T/X$K
print(round(X,2))
cat("\nAcross the whole default search grid (5-34 fx, d up to 8): how many candidates have onset < 7 d, and how many of those pass both toxicity caps?\n")
g <- expand.grid(d=seq(2.05,8,by=0.05), m=5:34); g$T <- overall_time(g$m); g <- g[g$T<=46,]
a <- compute_all(g$m,g$d,g$T,p,ref); r0 <- compute_all(35,2,47,p,ref)
ok <- a$bed_early<=r0$bed_early+1e-9 & a$bed_late<=r0$bed_late+1e-9
for (CA in c(cal$c_over_alpha, cm$c_over_alpha)) { on <- CA*g$T/a$bed_di*0 + CA*g$T/bed_lq(g$m,g$d,10)
  cat(sprintf("C/a=%.1f: %d of %d grid points onset<7; of those passing toxicity caps: %d; min onset among cap-passing: %.1f d\n",
    CA, sum(on<7), nrow(g), sum(on<7 & ok), min(on[ok]))) }
cat("\nLiterature presets (lambda/alpha Gy/day): Withers 0.72; Fowler2010 a.35 Tp3 ->", round(log(2)/3/.35,3), "; Hartley a.3 Tp3 ->", round(log(2)/3/.3,3), "; calculator a.2 l.2 -> 1.0\n")
for (nm in list(list("fowler2010 Tk21",0.35,log(2)/3,21),list("fowler2010 Tk28",0.35,log(2)/3,28),list("hartley Tk22",0.3,log(2)/3,22),list("calc a.2 Withers",0.2,0.144,28))) {
  q <- p; q$alpha <- nm[[2]]; q$lambda_di <- nm[[3]]; q$Tk <- nm[[4]]
  c3 <- dd_calibration(ref,q$alpha,10,q$lambda_di,q$lambda_dd,q$Tk); pkq <- 1-exp(-q$alpha*K0/47)
  s3 <- suggest_regimens(q,ref); q$lambda_dd <- q$lambda_di/pkq; s4 <- suggest_regimens(q,ref)
  cat(sprintf("%s: DD(0.5) status %s onset %.1f; rows(λdd .5)=%d; matched λdd=%.3f rows=%d (DD-credited %d); HYPNO dDI=%.2f\n", nm[[1]], c3$status, c3$onset_ref, nrow(s3), q$lambda_dd, nrow(s4), sum(s4$d_dd>=0,na.rm=TRUE),
    bed_di(20,2.75,26,q$alpha,10,q$lambda_di,q$Tk)-bed_di(35,2,47,q$alpha,10,q$lambda_di,q$Tk)))
}
