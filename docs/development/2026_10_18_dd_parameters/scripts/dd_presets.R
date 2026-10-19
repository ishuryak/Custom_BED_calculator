source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
p <- default_params(); ref <- default_ref(); K0 <- 84
for (nm in list(list("calculator a.2 l.2 Tk28",0.2,0.2,28),list("Withers a.2 l.144 Tk28",0.2,0.144,28),list("Fowler2010 a.35 Tp3 Tk21",0.35,log(2)/3,21),list("Fowler2010 a.35 Tp3 Tk28",0.35,log(2)/3,28),list("Hartley a.3 Tp3 Tk22",0.3,log(2)/3,22))) {
  q <- p; q$alpha <- nm[[2]]; q$lambda_di <- nm[[3]]; q$Tk <- nm[[4]]; pk <- 1-exp(-q$alpha*K0/47)
  for (ld in c(0.5, q$lambda_di/pk)) { q$lambda_dd <- ld
    a <- suggest_regimens(q, ref, min_gain=1, min_days_saved=7); b <- suggest_regimens(q, ref, min_gain=1, min_days_saved=7, early_cap=59)
    hy <- compute_all(20,2.75,26,q,ref); r0 <- compute_all(35,2,47,q,ref)
    cat(sprintf("%-26s λdd=%.3f | cap=std: %2d rows (DD %2d) | cap=59: %2d rows (DD %2d) | HYPNO dDI %+.2f dDD %+.2f\n", nm[[1]], ld, nrow(a), sum(a$d_dd>=1,na.rm=T), nrow(b), sum(b$d_dd>=1,na.rm=T), hy$bed_di-r0$bed_di, hy$bed_dd-r0$bed_dd)) }
}
