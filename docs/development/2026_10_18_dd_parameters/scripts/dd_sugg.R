source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
ref <- default_ref()
for (pt in list(c(0.2,0.5),c(0.30,0.375),c(0.40,0.394))) { p <- default_params(); p$alpha <- pt[1]; p$lambda_dd <- pt[2]
  a <- suggest_regimens(p, ref, min_gain=1, min_days_saved=7); b <- suggest_regimens(p, ref, min_gain=1, min_days_saved=7, early_cap=59)
  cat(sprintf("alpha %.2f ldd %.3f: strict early cap %d rows (DD-credited %d); Fowler cap 59: %d rows (DD-credited %d)\n", pt[1],pt[2], nrow(a), sum(a$d_dd>=1,na.rm=T), nrow(b), sum(b$d_dd>=1,na.rm=T))) }
