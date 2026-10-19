source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
p <- default_params(); ref <- default_ref()
K0 <- bed_lq(35,2,10); pk0 <- 1-exp(-p$alpha*K0/47)
orig <- bed_dd
FLOOR <- 0
bed_dd <- function(m, d, T, alpha, ab, lambda_dd, c_over_alpha) {
  K <- bed_lq(m, d, ab)
  if (is.na(c_over_alpha)) return(rep(NA_real_, length(K)))
  if (is.infinite(c_over_alpha)) return(K)
  p_kill <- 1 - exp(-alpha * K / T)
  onset <- pmax(FLOOR, c_over_alpha * T / K)
  K - lambda_dd * p_kill * smooth_max0(T - onset) / alpha
}
show <- function(lab, q, ...) { s <- suggest_regimens(q, ref, ...)
  cat(sprintf("%-38s rows=%2d DD-credited=%2d  fx: %s\n", lab, nrow(s), sum(s$d_dd>=0,na.rm=TRUE), paste(sort(s$m),collapse=","))) }
for (fl in c(0,7)) { FLOOR <- fl; cat("\n== onset floor", fl, "days ==\n")
  for (ld in c(0.5, p$lambda_di/pk0, 2)) { q <- p; q$lambda_dd <- ld
    show(sprintf("lambda_dd=%.3f either", ld), q)
    show(sprintf("lambda_dd=%.3f dd-only", ld), q, criterion="dd")
    show(sprintf("lambda_dd=%.3f dd-only, Fowler cap 59", ld), q, criterion="dd", early_cap=59) }
}
cat("\nDD BED of short courses, default lambda_dd 0.5, floor 0 vs 7 (ref DD = 65.0):\n")
cal <- dd_calibration(ref,p$alpha,10,p$lambda_di,0.5,p$Tk)
X <- data.frame(m=c(5,5,10,10,15,20), d=c(6,8,4,5,3.1,2.75)); X$T <- overall_time(X$m)
FLOOR<-0; X$dd0 <- bed_dd(X$m,X$d,X$T,.2,10,.5,cal$c_over_alpha); FLOOR<-7; X$dd7 <- bed_dd(X$m,X$d,X$T,.2,10,.5,cal$c_over_alpha)
X$K <- bed_lq(X$m,X$d,10); X$early <- bed_early(X$m,X$d,X$T,10,.35,log(2)/2.5,7); X$late <- bed_late(X$m,X$d,3)
print(round(X,2)); cat(sprintf("ref early %.2f late %.2f\n", bed_early(35,2,47,10,.35,log(2)/2.5,7), bed_late(35,2,3)))
