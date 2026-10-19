source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
p <- default_params(); ref <- default_ref(); pk0 <- 1-exp(-0.2*84/47)
for (ld in c(0.27,0.3,0.35,0.4,0.5,0.665)) {
  q <- p; q$lambda_dd <- ld; cal <- dd_calibration(ref,.2,10,.2,ld,28)
  s <- suggest_regimens(q, ref, min_gain=1, min_days_saved=7)
  K <- bed_lq(s$m,s$d,10); on <- cal$c_over_alpha*s$T/K; rt <- ld*(1-exp(-.2*K/s$T))
  t <- compute_all(20,2.75,26,q,ref); c <- compute_all(33,2,38,q,ref); r <- compute_all(35,2,47,q,ref)
  cat(sprintf("λdd=%.3f %s | ref: onset %.1f d, rate %.3f/d (Tp %.1f d), loss %.2f Gy/d for %.0f d | suggestions %d, DD-credited %d, onsets %.1f-%.1f d, Tp %.1f-%.1f d, dDD %.1f to %.1f | HYPNO onset %.1f d, test-ctrl dDD %+.2f, test-ref dDD %+.2f\n",
    ld, cal$status, cal$onset_ref, ld*pk0, log(2)/(ld*pk0), ld*pk0/.2, 47-cal$onset_ref, nrow(s), sum(s$d_dd>=1), min(on), max(on), log(2)/max(rt), log(2)/min(rt), min(s$d_dd), max(s$d_dd),
    cal$c_over_alpha*26/70.125, t$bed_dd-c$bed_dd, t$bed_dd-r$bed_dd))
}
