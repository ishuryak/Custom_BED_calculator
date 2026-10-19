source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
p <- default_params(); ref <- default_ref()
for (nm in list(list("calculator",0.2,0.2,28),list("Withers",0.2,0.144,28),list("Fowler2010 Tk21",0.35,log(2)/3,21),list("Hartley Tk22",0.3,log(2)/3,22))) {
  q <- p; q$alpha <- nm[[2]]; q$lambda_di <- nm[[3]]; q$Tk <- nm[[4]]; pk <- 1-exp(-q$alpha*84/47)
  for (ld in c(0.5, q$lambda_di/pk, 2)) { q$lambda_dd <- ld
    t <- compute_all(20,2.75,26,q,ref); c <- compute_all(33,2,38,q,ref)
    cat(sprintf("%-16s λdd=%.3f  HYPNO test-control: dDI %+.2f  dDD %+.2f  (Gy10)\n", nm[[1]], ld, t$bed_di-c$bed_di, t$bed_dd-c$bed_dd)) } }
